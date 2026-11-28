// SPDX-License-Identifier: GPL-3.0-or-later

#include <ATen/Context.h>
#include <ATen/Parallel.h>
#include <torch/csrc/autograd/autograd.h>
#include <torch/csrc/jit/ir/constants.h>
#include <torch/csrc/jit/passes/dead_code_elimination.h>
#include <torch/csrc/jit/runtime/graph_executor.h>
#include <torch/cuda.h>
#include <torch/script.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int device_cpu = 0;
constexpr int device_cuda = 1;
constexpr int feature_count = 9;

struct SkalaModel {
  struct Constant {
    torch::Tensor host, device;
    uint64_t bytes = 0, age = 0;
  };
  torch::jit::script::Module module;
  torch::Device device = torch::Device(torch::kCPU);
  bool features[feature_count] = {};
  int block = 0;
  uint64_t cache_limit = 64ULL * 1024 * 1024, cache_bytes = 0, age = 0;
  uint64_t cache_hits = 0, cache_misses = 0;
  bool copy_inputs = false;
  std::unordered_map<std::string, Constant> constants;
  std::vector<double> block_energies;
  torch::Tensor spin_directions;
  torch::Tensor spin_response_directions, spin_energy_coefficients;
  std::vector<double> spin_weights;
  int spin_batch_size = 1;
  uint64_t spin_batch_calls = 0, spin_directions_evaluated = 0;
  int gradient_square_fixes = 0;
  bool fixed_executor = false;
};

using TensorDict = c10::Dict<std::string, torch::Tensor>;

bool kind_is(const torch::jit::Node *node, const char *name) {
  return std::strcmp(node->kind().toQualString(), name) == 0;
}

bool real_dtype(const torch::jit::Node *node) {
  const auto &arguments = node->schema().arguments();
  for (std::size_t i = 0; i < arguments.size() && i < node->inputs().size(); ++i) {
    if (arguments[i].name() != "dtype") continue;
    const auto value = torch::jit::toIValue(node->input(i));
    if (!value) return false;
    return value->isNone() || (value->isInt() &&
        (value->toInt() == static_cast<int>(torch::kFloat32) ||
         value->toInt() == static_cast<int>(torch::kFloat64)));
  }
  return !kind_is(node, "aten::to");
}

bool layout_only(const torch::jit::Node *node) {
  const char *names[] = {"aten::permute", "aten::unsqueeze", "aten::reshape",
                        "aten::view", "aten::numpy_T", "aten::slice"};
  return std::any_of(std::begin(names), std::end(names),
                     [&](const char *name) { return kind_is(node, name); });
}

bool real_primitive(torch::jit::Value *value, torch::jit::Value *primitive, int depth) {
  if (value == primitive) return true; // The bridge supplies real Float64 primitives.
  if (depth > 32) return false;
  const auto *node = value->node();
  if (layout_only(node) ||
      ((kind_is(node, "aten::sum") || kind_is(node, "aten::to")) && real_dtype(node))) {
    return real_primitive(node->input(0), primitive, depth+1);
  }
  return false;
}

bool nonnegative_gradient_square(torch::jit::Value *value,
                                 torch::jit::Value *gradient, int depth=0) {
  if (depth > 32) return false;
  const auto *node = value->node();
  if (kind_is(node, "aten::pow")) {
    const auto exponent = torch::jit::toIValue(node->input(1));
    const bool square = exponent && ((exponent->isInt() && exponent->toInt() == 2) ||
                                     (exponent->isDouble() && exponent->toDouble() == 2.0));
    return square && real_primitive(node->input(0), gradient, 0);
  }
  if (layout_only(node) ||
      ((kind_is(node, "aten::sum") || kind_is(node, "aten::to")) && real_dtype(node))) {
    return nonnegative_gradient_square(node->input(0), gradient, depth+1);
  }
  return false;
}

void collect_absolute_nodes(torch::jit::Block *block, std::vector<torch::jit::Node *> &nodes) {
  for (auto *node : block->nodes()) {
    if (kind_is(node, "aten::abs")) nodes.push_back(node);
    for (auto *child : node->blocks()) collect_absolute_nodes(child, nodes);
  }
}

int correct_gradient_square_derivatives(torch::jit::script::Module &model) {
  if (!model.hasattr("semi_local_features")) return 0;
  auto feature = model.attr("semi_local_features").toModule();
  int count = 0;
  for (const auto &method : feature.get_methods()) {
    if (method.name() != "forward" && method.name() != "forward1") continue;
    auto graph = method.graph();
    if (graph->inputs().size() != 4) continue;
    auto *gradient = graph->inputs()[2];
    std::vector<torch::jit::Node *> nodes;
    collect_absolute_nodes(graph->block(), nodes);
    for (auto *absolute : nodes) {
      auto *input = absolute->input(0);
      auto *cast = kind_is(input->node(), "aten::to") ? input->node() : nullptr;
      if (cast && (!real_dtype(cast) || !cast->attributeNames().empty())) continue;
      auto *cat = cast ? cast->input(0)->node() : input->node();
      if (!kind_is(cat, "aten::cat") || cat->inputs().size() != 2 ||
          !kind_is(cat->input(0)->node(), "prim::ListConstruct")) continue;
      const auto pieces = cat->input(0)->node()->inputs();
      if (pieces.size() != 4 || !real_primitive(pieces[0], graph->inputs()[1], 0) ||
          !real_primitive(pieces[2], graph->inputs()[3], 0)) continue;
      std::vector<bool> nonnegative;
      for (auto *piece : pieces) nonnegative.push_back(nonnegative_gradient_square(piece, gradient));
      if (!nonnegative[1] || !nonnegative[3]) continue;

      // abs(sum(g*g)) has an incorrect zero Hessian at g=0 in autograd.
      // Distribute abs over cat, omitting it only for proven real gradient
      // squares. Values and weights are unchanged; density/tau abs stays intact.
      std::vector<torch::jit::Value *> values;
      for (std::size_t i = 0; i < pieces.size(); ++i) {
        auto *piece = pieces[i];
        if (!nonnegative[i]) {
          auto *part = graph->create(c10::Symbol::fromQualString("aten::abs"), {piece});
          part->output()->setType(piece->type());
          part->insertBefore(absolute);
          piece = part->output();
        }
        values.push_back(piece);
      }
      auto *list = graph->create(c10::Symbol::fromQualString("prim::ListConstruct"), values);
      list->output()->setType(cat->input(0)->type());
      list->insertBefore(absolute);
      auto *joined = graph->create(c10::Symbol::fromQualString("aten::cat"),
                                   {list->output(), cat->input(1)});
      joined->output()->setType(cat->output()->type());
      joined->insertBefore(absolute);
      auto *replacement = joined->output();
      if (cast) {
        std::vector<torch::jit::Value *> inputs(cast->inputs().begin(), cast->inputs().end());
        inputs[0] = replacement;
        auto *converted = graph->create(cast->kind(), inputs);
        converted->output()->setType(input->type());
        converted->insertBefore(absolute);
        replacement = converted->output();
      }
      absolute->output()->replaceAllUsesWith(replacement);
      absolute->destroy();
      ++count;
    }
    torch::jit::EliminateDeadCode(graph);
    graph->lint();
  }
  return count;
}

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

bool is_differentiable_feature(const std::string &key, bool geometry) {
  return key == "density" || key == "grad" || key == "kin" ||
         (geometry && (key == "grid_coords" || key == "grid_weights" ||
                       key == "coarse_0_atomic_coords" || key == "atomic_grid_weights"));
}

torch::Tensor constant_value(SkalaModel &model, const std::string &key,
                             const torch::Tensor &input) {
  if (model.block <= 0 || model.cache_limit == 0) {
    return input.to(model.device);
  }
  const std::string id = std::to_string(model.block) + ":" + key;
  auto host = input.to(torch::kCPU);
  auto found = model.constants.find(id);
  if (found != model.constants.end()) {
    auto &saved = found->second;
    if (saved.host.scalar_type() == host.scalar_type() &&
        saved.host.sizes() == host.sizes() && saved.host.strides() == host.strides() &&
        torch::equal(saved.host, host)) {
      saved.age = ++model.age;
      ++model.cache_hits;
      return saved.device;
    }
    model.cache_bytes -= saved.bytes;
    model.constants.erase(found);
  }
  ++model.cache_misses;
  // Count host snapshots and device values conservatively, including CPU mode.
  const uint64_t bytes = 256ULL + 2ULL * host.numel() * host.element_size();
  if (bytes > model.cache_limit) return input.to(model.device);
  while (model.cache_bytes + bytes > model.cache_limit && !model.constants.empty()) {
    auto oldest = std::min_element(model.constants.begin(), model.constants.end(),
        [](const auto &a, const auto &b) { return a.second.age < b.second.age; });
    model.cache_bytes -= oldest->second.bytes;
    model.constants.erase(oldest);
  }
  SkalaModel::Constant saved;
  saved.host = host.clone();
  saved.device = saved.host.to(model.device);
  saved.bytes = bytes;
  saved.age = ++model.age;
  auto value = saved.device;
  model.constants.emplace(id, std::move(saved));
  model.cache_bytes += bytes;
  return value;
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
  const int64_t nspin = density.size(0);
  if (nspin != 2 && nspin != 4) {
    throw std::runtime_error("Skala requires two spin populations or four Pauli components");
  }
  require_shape(features, "density", {nspin, npoint});
  require_shape(features, "grad", {nspin, 3, npoint});
  require_shape(features, "kin", {nspin, npoint});
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

std::pair<torch::Tensor, torch::Tensor> project_spin(const torch::Tensor &pauli,
    const torch::Tensor &axis, const torch::Tensor &response_axis = {}) {
  std::vector<int64_t> shape(pauli.dim(), 1);
  shape[0] = 3;
  const auto spin = (pauli.narrow(0, 1, 3) * axis.reshape(shape)).sum(0);
  const auto half_charge = 0.5 * pauli.select(0, 0);
  const auto half_spin = 0.5 * spin;
  const auto tangent = response_axis.defined()
      ? 0.5 * (pauli.narrow(0, 1, 3) * response_axis.reshape(shape)).sum(0) : half_spin;
  return {torch::stack({half_charge + half_spin, half_charge - half_spin}),
          torch::stack({tangent, -tangent})};
}

torch::Tensor collinear_rows(SkalaModel &model, const TensorDict &features) {
  std::vector<c10::IValue> args;
  std::unordered_map<std::string, c10::IValue> kwargs;
  kwargs["mol"] = features;
  return model.module.get_method("get_exc_density")(args, kwargs).toTensor() *
         features.at("grid_weights");
}

torch::Tensor effective_spin_rows(SkalaModel &model, const TensorDict &projected,
    const std::vector<torch::Tensor> &fields, const std::vector<torch::Tensor> &spin,
    const torch::Tensor &energy_coefficients = {}) {
  const auto weighted = collinear_rows(model, projected);
  if (weighted.dim() != 1) throw std::runtime_error("Skala energy rows must be one-dimensional");
  const auto energy = weighted.sum();
  // Pu et al., Phys. Rev. Research 5, 013036 (2023), Eq. (21).
  // Scaling all spin fields requires a second reverse pass, including geometry.
  const auto first = torch::autograd::grad({energy}, fields, {}, true, true, true);
  auto rows = weighted;
  if (energy_coefficients.defined()) {
    rows = (weighted.reshape({energy_coefficients.numel(), -1}) *
            energy_coefficients.reshape({-1, 1})).reshape({-1});
  }
  for (std::size_t i = 0; i < fields.size(); ++i) {
    if (!first[i].defined()) continue;
    const auto response = first[i] * spin[i];
    rows = rows + (i == 1 ? response.sum(0).sum(0) : response.sum(0));
  }
  return rows;
}

std::pair<torch::Tensor, std::vector<torch::Tensor>> multicollinear_batched(
    SkalaModel &model, const TensorDict &features, const std::vector<torch::Tensor> &leaves) {
  const int64_t npoint = features.at("density").size(1);
  auto rows = torch::zeros_like(features.at("grid_weights"));
  std::vector<torch::Tensor> gradients;
  for (const auto &leaf : leaves) gradients.push_back(torch::zeros_like(leaf));
  const char *keys[] = {"density", "grad", "kin"};
  for (std::size_t q = 0; q < model.spin_weights.size();) {
    const auto count = std::min<std::size_t>(model.spin_batch_size, model.spin_weights.size() - q);
    TensorDict projected;
    for (const auto &entry : features) {
      const std::string key = entry.key();
      if (key == "density" || key == "grad" || key == "kin") continue;
      auto value = entry.value();
      // Each direction keeps complete independent atom blocks. Repetition
      // remains in the graph so geometry adjoints return to the original leaves.
      if (key == "grid_coords" || key == "grid_weights" || key == "coarse_0_atomic_coords" ||
          key == "atomic_grid_weights" || key == "atomic_grid_sizes") {
        std::vector<int64_t> repeats(value.dim(), 1);
        repeats[0] = static_cast<int64_t>(count);
        value = value.repeat(repeats);
      }
      projected.insert(key, value);
    }
    std::vector<torch::Tensor> fields, spin;
    for (const auto *key : keys) {
      std::vector<torch::Tensor> values, tangents;
      for (std::size_t i = 0; i < count; ++i) {
        auto value = project_spin(features.at(key), model.spin_directions[q + i],
            model.spin_response_directions.defined() ? model.spin_response_directions[q + i] : torch::Tensor());
        values.push_back(std::move(value.first));
        tangents.push_back(std::move(value.second));
      }
      fields.push_back(torch::cat(values, -1));
      spin.push_back(torch::cat(tangents, -1));
      projected.insert(key, fields.back());
    }
    const auto coefficients = model.spin_energy_coefficients.defined()
        ? model.spin_energy_coefficients.narrow(0, q, count) : torch::Tensor();
    const auto effective = effective_spin_rows(model, projected, fields, spin, coefficients);
    const auto weights = torch::tensor(std::vector<double>(model.spin_weights.begin() + q,
        model.spin_weights.begin() + q + count), torch::TensorOptions().dtype(torch::kFloat64).device(model.device));
    const auto weighted = (effective.reshape({static_cast<int64_t>(count), npoint}) * weights.unsqueeze(1)).sum(0);
    const auto derivative = torch::autograd::grad({weighted.sum()}, leaves, {}, false, false, true);
    rows.add_(weighted.detach());
    for (std::size_t i = 0; i < gradients.size(); ++i) {
      if (derivative[i].defined()) gradients[i].add_(derivative[i].detach());
    }
    ++model.spin_batch_calls;
    model.spin_directions_evaluated += count;
    q += count;
  }
  return {rows, gradients};
}

std::pair<torch::Tensor, std::vector<torch::Tensor>> multicollinear(
    SkalaModel &model, const TensorDict &features, const std::vector<torch::Tensor> &leaves) {
  if (!model.spin_directions.defined() || model.spin_weights.empty()) {
    throw std::runtime_error("Pauli Skala requires an explicit normalized spin quadrature");
  }
  if (model.spin_batch_size > 1) return multicollinear_batched(model, features, leaves);
  auto rows = torch::zeros_like(features.at("grid_weights"));
  std::vector<torch::Tensor> gradients;
  for (const auto &leaf : leaves) gradients.push_back(torch::zeros_like(leaf));
  const char *keys[] = {"density", "grad", "kin"};
  for (std::size_t q = 0; q < model.spin_weights.size(); ++q) {
    TensorDict projected;
    for (const auto &entry : features) projected.insert(entry.key(), entry.value());
    std::vector<torch::Tensor> fields, spin;
    for (const auto *key : keys) {
      auto projection = project_spin(features.at(key), model.spin_directions[q],
          model.spin_response_directions.defined() ? model.spin_response_directions[q] : torch::Tensor());
      fields.push_back(projection.first);
      spin.push_back(projection.second);
      projected.insert_or_assign(key, fields.back());
    }
    const auto coefficient = model.spin_energy_coefficients.defined()
        ? model.spin_energy_coefficients[q] : torch::Tensor();
    const auto effective_rows = effective_spin_rows(model, projected, fields, spin, coefficient);
    const auto derivative = torch::autograd::grad(
        {effective_rows.sum()}, leaves, {}, false, false, true);
    const double weight = model.spin_weights[q];
    rows.add_(weight * effective_rows.detach());
    for (std::size_t i = 0; i < gradients.size(); ++i) {
      if (derivative[i].defined()) gradients[i].add_(weight * derivative[i].detach());
    }
    ++model.spin_batch_calls;
    ++model.spin_directions_evaluated;
  }
  return {rows, gradients};
}

void set_spin_quadrature(SkalaModel &model, const int count, const double *directions,
    const double *weights, const double *coefficients, const double *responses) {
  if (count < 1 || directions == nullptr || weights == nullptr ||
      (coefficients == nullptr) != (responses == nullptr)) {
    throw std::runtime_error("invalid Skala spin quadrature");
  }
  double total = 0;
  for (int q = 0; q < count; ++q) {
    double norm = 0;
    for (int d = 0; d < 3; ++d) {
      norm += directions[3 * q + d] * directions[3 * q + d];
      if (responses != nullptr && !std::isfinite(responses[3 * q + d])) {
        throw std::runtime_error("Skala spin response directions must be finite");
      }
    }
    if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1e-10 || !std::isfinite(weights[q]) ||
        (coefficients != nullptr && !std::isfinite(coefficients[q]))) {
      throw std::runtime_error("Skala spin quadrature requires finite coefficients, unit directions and weights");
    }
    total += weights[q];
  }
  if (!std::isfinite(total) || std::abs(total - 1.0) > 1e-10) {
    throw std::runtime_error("Skala spin quadrature weights must sum to one");
  }
  auto axes = torch::from_blob(const_cast<double *>(directions), {count, 3}, torch::kFloat64)
                  .clone().to(model.device);
  torch::Tensor response_axes, energy_coefficients;
  if (responses != nullptr) {
    response_axes = torch::from_blob(const_cast<double *>(responses), {count, 3}, torch::kFloat64)
                        .clone().to(model.device);
    energy_coefficients = torch::from_blob(const_cast<double *>(coefficients), {count}, torch::kFloat64)
                              .clone().to(model.device);
  }
  std::vector<double> new_weights(weights, weights + count);
  model.spin_weights = std::move(new_weights);
  model.spin_directions = std::move(axes);
  model.spin_response_directions = std::move(response_axes);
  model.spin_energy_coefficients = std::move(energy_coefficients);
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

bool configure_reproducibility() {
  static std::once_flag configured;
  static bool deterministic = false;
  std::call_once(configured, [] {
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
  return deterministic;
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
    const bool fixed_executor = configure_reproducibility();
    auto model = std::make_unique<SkalaModel>();
    model->fixed_executor = fixed_executor;
    if (const char *value = std::getenv("CPPAW_SKALA_SPIN_BATCH")) {
      char *end = nullptr;
      errno = 0;
      const long count = std::strtol(value, &end, 10);
      if (end == value || *end != '\0' || errno == ERANGE || count < 1 || count > 64) {
        throw std::runtime_error("CPPAW_SKALA_SPIN_BATCH must be an integer from 1 through 64");
      }
      model->spin_batch_size = static_cast<int>(count);
    }
    if (const char *value = std::getenv("CPPAW_SKALA_COPY_INPUTS")) {
      model->copy_inputs = std::string(value) == "1";
    }
    if (const char *value = std::getenv("CPPAW_SKALA_TENSOR_CACHE_MB")) {
      char *end = nullptr;
      errno = 0;
      const auto mb = std::strtoull(value, &end, 10);
      if (value[0] == '-' || end == value || *end != '\0' || errno == ERANGE ||
          mb > UINT64_MAX / (1024ULL * 1024)) {
        throw std::runtime_error("CPPAW_SKALA_TENSOR_CACHE_MB must be a nonnegative integer");
      }
      model->cache_limit = mb * 1024ULL * 1024;
    }
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
    model->gradient_square_fixes = correct_gradient_square_derivatives(model->module);
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

extern "C" int cppaw_skala_model_gradient_fixes(void *handle) {
  return handle == nullptr ? 0 : static_cast<SkalaModel *>(handle)->gradient_square_fixes;
}

extern "C" void cppaw_skala_model_set_block(void *handle, const int block) {
  static_cast<SkalaModel *>(handle)->block = block;
}

extern "C" int cppaw_skala_model_set_spin_quadrature(void *handle, const int count,
    const double *directions, const double *weights, char *error, const int capacity) {
  try {
    if (handle == nullptr) {
      throw std::runtime_error("invalid Skala spin quadrature");
    }
    set_spin_quadrature(*static_cast<SkalaModel *>(handle), count, directions, weights, nullptr, nullptr);
    set_error(error, capacity, "");
    return 0;
  } catch (const std::exception &exception) {
    set_error(error, capacity, exception.what());
    return 1;
  }
}

extern "C" int cppaw_skala_model_set_spin_response_quadrature(void *handle, const int count,
    const double *directions, const double *weights, const double *coefficients,
    const double *responses, char *error, const int capacity) {
  try {
    if (handle == nullptr || coefficients == nullptr || responses == nullptr) {
      throw std::runtime_error("invalid Skala linear spin-response quadrature");
    }
    set_spin_quadrature(*static_cast<SkalaModel *>(handle), count, directions, weights, coefficients, responses);
    set_error(error, capacity, "");
    return 0;
  } catch (const std::exception &exception) {
    set_error(error, capacity, exception.what());
    return 1;
  }
}

extern "C" void cppaw_skala_model_cache_stats(void *handle, int64_t *stats) {
  const auto &model = *static_cast<SkalaModel *>(handle);
  stats[0] = static_cast<int64_t>(model.cache_hits);
  stats[1] = static_cast<int64_t>(model.cache_misses);
  stats[2] = static_cast<int64_t>(model.cache_bytes);
}

extern "C" void cppaw_skala_model_spin_batch_stats(void *handle, int64_t *stats) {
  const auto &model = *static_cast<SkalaModel *>(handle);
  stats[0] = model.spin_batch_size;
  stats[1] = static_cast<int64_t>(model.spin_batch_calls);
  stats[2] = static_cast<int64_t>(model.spin_directions_evaluated);
}

extern "C" int cppaw_skala_model_block_energies(void *handle, const int count,
    double *values, char *error, const int capacity) {
  const auto &energies = static_cast<SkalaModel *>(handle)->block_energies;
  if (count < 1 || static_cast<std::size_t>(count) != energies.size()) {
    set_error(error, capacity, "Skala block-energy extent does not match the last evaluation");
    return 1;
  }
  std::copy(energies.begin(), energies.end(), values);
  return 0;
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

extern "C" int cppaw_skala_model_evaluate_mode(void *model_handle,
                                            void *input_handle,
                                            const int geometry_derivatives,
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
    // Profiling can change FP32 forward/backward arithmetic after warmup even
    // with deterministic kernels. Keep opt-in checks on one execution plan.
    std::optional<torch::jit::GraphOptimizerEnabledGuard> executor_guard;
    if (model.fixed_executor) executor_guard.emplace(false);
    const auto &input = *static_cast<TensorDict *>(input_handle);
    validate_feature_shapes(input);

    TensorDict features;
    std::vector<torch::Tensor> leaves;
    std::vector<std::string> leaf_keys;
    for (const auto &entry : input) {
      const std::string key = entry.key();
      const bool dynamic = key == "density" || key == "grad" || key == "kin";
      torch::Tensor value = dynamic ? entry.value().to(model.device)
                                   : constant_value(model, key, entry.value());
      if (is_differentiable_feature(key, geometry_derivatives != 0)) {
        // Inputs remain alive and read-only until this synchronous call returns.
        value = value.detach();
        if (model.copy_inputs) value = value.clone();
        value.set_requires_grad(true);
        leaves.push_back(value);
        leaf_keys.push_back(key);
      }
      features.insert(key, value);
    }

    torch::Tensor weighted;
    std::vector<torch::Tensor> gradients;
    if (features.at("density").size(0) == 4) {
      auto result = multicollinear(model, features, leaves);
      weighted = std::move(result.first);
      gradients = std::move(result.second);
    } else {
      weighted = collinear_rows(model, features);
      // Collinear stationary quantities need only first model derivatives.
      gradients = torch::autograd::grad({weighted.sum()}, leaves, {}, false, false, true);
    }
    *energy = weighted.sum().item<double>();
    const auto sizes = input.at("atomic_grid_sizes").to(torch::kCPU).contiguous();
    model.block_energies.clear();
    if (sizes.numel() == 1) {
      model.block_energies.push_back(*energy);
    } else {
      std::vector<torch::Tensor> sums;
      int64_t offset = 0;
      const auto *counts = sizes.data_ptr<int64_t>();
      for (int64_t block = 0; block < sizes.numel(); ++block) {
        if (counts[block] < 1 || counts[block] > weighted.numel() - offset) {
          throw std::runtime_error("invalid Skala atom-block extent");
        }
        sums.push_back(weighted.detach().narrow(0, offset, counts[block]).sum());
        offset += counts[block];
      }
      if (offset != weighted.numel()) throw std::runtime_error("incomplete Skala atom-block coverage");
      const auto host = torch::stack(sums).to(torch::kCPU).to(torch::kFloat64).contiguous();
      model.block_energies.assign(host.data_ptr<double>(), host.data_ptr<double>() + host.numel());
    }
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

extern "C" int cppaw_skala_model_evaluate(void *model_handle, void *input_handle,
    double *energy, void **density_grad, void **spatial_grad, void **kin_grad,
    void **grid_coord_grad, void **grid_weight_grad, void **atom_coord_grad,
    void **atomic_weight_grad, char *error, const int error_capacity) {
  return cppaw_skala_model_evaluate_mode(model_handle, input_handle, 1, energy,
      density_grad, spatial_grad, kin_grad, grid_coord_grad, grid_weight_grad,
      atom_coord_grad, atomic_weight_grad, error, error_capacity);
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
