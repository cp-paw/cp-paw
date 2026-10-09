// SPDX-License-Identifier: GPL-3.0-or-later

#include <ATen/Context.h>
#include <ATen/Parallel.h>
#include <torch/csrc/autograd/autograd.h>
#include <torch/cuda.h>
#include <torch/script.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int device_cpu = 0;
constexpr int device_cuda = 1;
constexpr int feature_count = 9;

struct SkalaModel {
  torch::jit::script::Module module;
  torch::Device device = torch::Device(torch::kCPU);
  bool features[feature_count] = {};
};

using TensorDict = c10::Dict<std::string, torch::Tensor>;

void set_error(char *buffer, const int capacity, const std::string &message) {
  if (buffer == nullptr || capacity <= 0) {
    return;
  }
  const auto count = std::min(message.size(), static_cast<std::size_t>(capacity - 1));
  std::memcpy(buffer, message.data(), count);
  buffer[count] = '\0';
}

torch::Device select_device(const int device_type, const int device_index) {
  if (device_type == device_cpu) {
    return torch::Device(torch::kCPU);
  }
  if (device_type != device_cuda) {
    throw std::runtime_error("unknown Torch device type");
  }
  if (!torch::cuda::is_available()) {
    throw std::runtime_error("CUDA was requested but LibTorch reports no CUDA device");
  }
  const int index = device_index < 0 ? 0 : device_index;
  if (index >= torch::cuda::device_count()) {
    throw std::runtime_error("requested CUDA device index is out of range");
  }
  return torch::Device(torch::kCUDA, index);
}

void parse_features(const std::string &json, bool features[feature_count]) {
  const std::pair<const char *, int> keys[] = {
      {"density", 0},
      {"grad", 1},
      {"kin", 2},
      {"grid_coords", 3},
      {"grid_weights", 4},
      {"coarse_0_atomic_coords", 5},
      {"atomic_grid_weights", 6},
      {"atomic_grid_sizes", 7},
      {"atomic_grid_size_bound_shape", 8},
  };
  for (const auto &[key, index] : keys) {
    features[index] = json.find(std::string("\"") + key + "\"") != std::string::npos;
  }
}

bool is_differentiable_feature(const std::string &key) {
  return key == "density" || key == "grad" || key == "kin" ||
         key == "grid_coords" || key == "grid_weights" ||
         key == "coarse_0_atomic_coords" || key == "atomic_grid_weights";
}

void require_shape(const TensorDict &features,
                   const char *key,
                   const std::vector<int64_t> &expected) {
  const auto &tensor = features.at(key);
  if (tensor.dim() != static_cast<int64_t>(expected.size())) {
    throw std::runtime_error(std::string("Skala feature '") + key +
                             "' has the wrong rank");
  }
  for (std::size_t i = 0; i < expected.size(); ++i) {
    if (expected[i] >= 0 && tensor.size(static_cast<int64_t>(i)) != expected[i]) {
      throw std::runtime_error(std::string("Skala feature '") + key +
                               "' has an incompatible tensor shape");
    }
  }
}

void validate_feature_shapes(const TensorDict &features) {
  const auto &density = features.at("density");
  const auto &atom_coords = features.at("coarse_0_atomic_coords");
  if (density.dim() != 2 || atom_coords.dim() != 2) {
    throw std::runtime_error("Skala density or atom coordinates have the wrong rank");
  }
  const int64_t npoint = density.size(1);
  const int64_t natom = atom_coords.size(0);
  require_shape(features, "density", {2, npoint});
  require_shape(features, "grad", {2, 3, npoint});
  require_shape(features, "kin", {2, npoint});
  require_shape(features, "grid_coords", {npoint, 3});
  require_shape(features, "grid_weights", {npoint});
  require_shape(features, "coarse_0_atomic_coords", {natom, 3});
  require_shape(features, "atomic_grid_weights", {npoint});
  require_shape(features, "atomic_grid_sizes", {natom});
  require_shape(features, "atomic_grid_size_bound_shape", {-1, 0});
}

torch::Tensor *new_tensor(torch::Tensor tensor) {
  return new torch::Tensor(std::move(tensor));
}

void configure_host_threads() {
#ifdef CPPAW_SKALA_NVHPC_FORTRAN
  static std::once_flag configured;
  std::call_once(configured, [] {
    int thread_count = 1;
    if (const char *value = std::getenv("CPPAW_SKALA_TORCH_THREADS")) {
      char *end = nullptr;
      const long requested = std::strtol(value, &end, 10);
      if (end != value && *end == '\0' && requested > 0) {
        thread_count = static_cast<int>(requested);
      }
    }
    at::set_num_threads(thread_count);
    at::set_num_interop_threads(1);
  });
#endif
}

void configure_reproducibility() {
  static std::once_flag configured;
  std::call_once(configured, [] {
    bool deterministic = false;
    if (const char *value = std::getenv("CPPAW_SKALA_DETERMINISTIC")) {
      const std::string setting(value);
      deterministic = setting == "1" || setting == "true" ||
                      setting == "TRUE" || setting == "on" ||
                      setting == "ON";
    }
    if (!deterministic) {
      return;
    }
#if !defined(_WIN32)
    if (std::getenv("CUBLAS_WORKSPACE_CONFIG") == nullptr) {
      setenv("CUBLAS_WORKSPACE_CONFIG", ":4096:8", 0);
    }
#endif
    at::globalContext().setDeterministicAlgorithms(true, false);
    at::globalContext().setDeterministicCuDNN(true);
    at::globalContext().setAllowTF32CuDNN(false);
    at::globalContext().setAllowTF32CuBLAS(false);
  });
}

} // namespace

extern "C" int cppaw_skala_cuda_available() {
  try {
    return torch::cuda::is_available() ? 1 : 0;
  } catch (...) {
    return 0;
  }
}

extern "C" void *cppaw_skala_model_load(const char *filename,
                                         const int device_type,
                                         const int device_index,
                                         int *features,
                                         char *error,
                                         const int error_capacity) {
  try {
    configure_host_threads();
    configure_reproducibility();
    auto model = std::make_unique<SkalaModel>();
    model->device = select_device(device_type, device_index);
    torch::jit::ExtraFilesMap metadata{{"features", ""}, {"protocol_version", ""}};
    model->module = torch::jit::load(filename, model->device, metadata);
    model->module.eval();
    model->module.to(model->device);
    for (auto parameter : model->module.parameters()) {
      parameter.set_requires_grad(false);
    }
    if (metadata.at("protocol_version") != "2") {
      throw std::runtime_error("unsupported Skala TorchScript protocol version '" +
                               metadata.at("protocol_version") + "'");
    }
    parse_features(metadata.at("features"), model->features);
    for (int i = 0; i < feature_count; ++i) {
      features[i] = model->features[i] ? 1 : 0;
    }
    set_error(error, error_capacity, "");
    return model.release();
  } catch (const std::exception &exception) {
    set_error(error, error_capacity, exception.what());
    return nullptr;
  }
}

extern "C" void cppaw_skala_model_release(void *handle) {
  delete static_cast<SkalaModel *>(handle);
}

extern "C" void *cppaw_skala_dict_create() {
  return new TensorDict();
}

extern "C" void cppaw_skala_dict_release(void *handle) {
  delete static_cast<TensorDict *>(handle);
}

extern "C" int cppaw_skala_dict_insert(void *handle,
                                        const char *key,
                                        void *tensor_handle,
                                        char *error,
                                        const int error_capacity) {
  try {
    auto &dict = *static_cast<TensorDict *>(handle);
    const auto &tensor = *static_cast<torch::Tensor *>(tensor_handle);
    dict.insert(std::string(key), tensor);
    set_error(error, error_capacity, "");
    return 0;
  } catch (const std::exception &exception) {
    set_error(error, error_capacity, exception.what());
    return 1;
  }
}

extern "C" int cppaw_skala_model_evaluate(void *model_handle,
                                            void *input_handle,
                                            double *energy,
                                            void **density_grad,
                                            void **spatial_grad,
                                            void **kin_grad,
                                            void **grid_coord_grad,
                                            void **grid_weight_grad,
                                            void **atom_coord_grad,
                                            void **atomic_weight_grad,
                                            char *error,
                                            const int error_capacity) {
  try {
    auto &model = *static_cast<SkalaModel *>(model_handle);
    const auto &input = *static_cast<TensorDict *>(input_handle);
    validate_feature_shapes(input);

    TensorDict features;
    std::vector<torch::Tensor> leaves;
    std::vector<std::string> leaf_keys;
    for (const auto &entry : input) {
      const std::string key = entry.key();
      torch::Tensor value = entry.value().to(model.device);
      if (is_differentiable_feature(key)) {
        value = value.detach().clone().set_requires_grad(true);
        leaves.push_back(value);
        leaf_keys.push_back(key);
      }
      features.insert(key, value);
    }

    std::vector<c10::IValue> args;
    std::unordered_map<std::string, c10::IValue> kwargs;
    kwargs["mol"] = features;
    torch::Tensor exc_density =
        model.module.get_method("get_exc_density")(args, kwargs).toTensor();
    torch::Tensor exc = (exc_density * features.at("grid_weights")).sum();
    *energy = exc.item<double>();

    // Stationary energy, force, and stress paths need first model derivatives.
    // Keep create_graph=false; response properties require a separate contract.
    const auto gradients = torch::autograd::grad(
        {exc}, leaves, {}, false, false, true);
    std::unordered_map<std::string, torch::Tensor> by_key;
    for (std::size_t i = 0; i < gradients.size(); ++i) {
      if (gradients[i].defined()) {
        by_key.emplace(leaf_keys[i], gradients[i]);
      }
    }

    auto take = [&by_key](const char *key) -> void * {
      const auto found = by_key.find(key);
      if (found == by_key.end()) {
        return nullptr;
      }
      return new_tensor(found->second);
    };
    *density_grad = take("density");
    *spatial_grad = take("grad");
    *kin_grad = take("kin");
    *grid_coord_grad = take("grid_coords");
    *grid_weight_grad = take("grid_weights");
    *atom_coord_grad = take("coarse_0_atomic_coords");
    *atomic_weight_grad = take("atomic_grid_weights");
    set_error(error, error_capacity, "");
    return 0;
  } catch (const std::exception &exception) {
    set_error(error, error_capacity, exception.what());
    return 1;
  }
}

extern "C" int cppaw_skala_tensor_copy_double(void *tensor_handle,
                                                double *destination,
                                                const long long count,
                                                char *error,
                                                const int error_capacity) {
  try {
    if (tensor_handle == nullptr) {
      throw std::runtime_error("Skala did not return a requested feature gradient");
    }
    const auto &tensor = *static_cast<torch::Tensor *>(tensor_handle);
    torch::Tensor host = tensor.detach().to(torch::kCPU).to(torch::kFloat64).contiguous();
    if (host.numel() != count) {
      throw std::runtime_error("Skala feature-gradient shape does not match the CP-PAW buffer");
    }
    std::memcpy(destination, host.data_ptr<double>(),
                static_cast<std::size_t>(count) * sizeof(double));
    set_error(error, error_capacity, "");
    return 0;
  } catch (const std::exception &exception) {
    set_error(error, error_capacity, exception.what());
    return 1;
  }
}

extern "C" void cppaw_skala_tensor_release(void *tensor_handle) {
  delete static_cast<torch::Tensor *>(tensor_handle);
}
