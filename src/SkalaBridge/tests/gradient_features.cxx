// SPDX-License-Identifier: GPL-3.0-or-later

#include <ATen/Parallel.h>
#include <torch/csrc/autograd/autograd.h>
#include <torch/csrc/jit/runtime/graph_executor.h>
#include <torch/script.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" {
void *cppaw_skala_model_load(const char *, int, int, int *, char *, int);
void cppaw_skala_model_release(void *);
int cppaw_skala_model_gradient_fixes(void *);
int cppaw_skala_model_set_spin_quadrature(void *, int, const double *, const double *, char *, int);
int cppaw_skala_model_set_spin_response_quadrature(void *, int, const double *, const double *,
                                                  const double *, const double *, char *, int);
void cppaw_skala_model_spin_batch_stats(void *, int64_t *);
int cppaw_skala_model_evaluate_mode(void *, void *, int, double *, void **, void **, void **,
                                   void **, void **, void **, void **, char *, int);
}

using Inputs = c10::Dict<std::string, torch::Tensor>;
using Model = std::unique_ptr<void, decltype(&cppaw_skala_model_release)>;

bool fixed_executor_requested() {
  const char *raw = std::getenv("CPPAW_SKALA_DETERMINISTIC");
  const std::string value = raw == nullptr ? "" : raw;
  return value == "1" || value == "true" || value == "TRUE" || value == "on" || value == "ON";
}

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

Model load(const std::filesystem::path &file, int expected, int device=0) {
  int features[9];
  char error[4096] = {};
  Model model(cppaw_skala_model_load(file.c_str(), device, 0, features, error, sizeof(error)),
              cppaw_skala_model_release);
  require(model != nullptr, error);
  require(cppaw_skala_model_gradient_fixes(model.get()) == expected, "Unexpected graph rewrite count");
  return model;
}

std::filesystem::path fixture(int index, const std::string &gradient, bool corrected=false,
                              bool with_features=true, bool cast_other=false,
                              const std::string &density="rho") {
  torch::jit::Module model("GradientFixture"+std::to_string(index));
  if (with_features) {
    torch::jit::Module feature("SemiLocalFeatures"+std::to_string(index));
    const std::string body =
        "(self, rho: Tensor, grad: Tensor, tau: Tensor) -> Tensor:\n"
        "    g = "+gradient+"\n"
        "    g2 = torch.sum(g ** 2, dim=[1])\n"
        "    total = torch.sum(torch.sum(g, dim=[0]) ** 2, dim=[0]).unsqueeze(0)\n"
        "    x = torch.cat(["+density+", g2, tau, total], dim=0)\n"+
        (corrected ?
         "    x = torch.cat([torch.abs(rho), g2, torch.abs(tau), total], dim=0)\n"
         "    return torch.log(x + 1.e-5)\n" :
         "    return torch.log(torch.abs(x"+std::string(cast_other ? ".to(rho)" : ".to(dtype=7)")+
         ") + 1.e-5)\n");
    feature.define("def forward"+body+"\ndef forward1"+body);
    model.register_module("semi_local_features", feature);
    model.define(
        "def get_exc_density(self, mol: Dict[str, Tensor]) -> Tensor:\n"
        "    return self.semi_local_features.forward(mol['density'], mol['grad'], mol['kin']).sum(0)\n");
  } else {
    model.define(
        "def get_exc_density(self, mol: Dict[str, Tensor]) -> Tensor:\n"
        "    return mol['density'].sum(0)\n");
  }
  const auto file = std::filesystem::current_path()/("gradient-fixture-"+std::to_string(index)+".fun");
  require(!std::filesystem::exists(file), "Refusing to overwrite a model fixture");
  model.save(file.string(), {{"protocol_version", "2"}, {"features", "[\"density\",\"grad\",\"kin\"]"}});
  return file;
}

Inputs inputs(int nc) {
  const auto options = torch::TensorOptions().dtype(torch::kFloat64);
  Inputs result;
  result.insert("density", torch::zeros({nc, 4}, options));
  result.insert("grad", torch::zeros({nc, 3, 4}, options));
  result.insert("kin", torch::zeros({nc, 4}, options));
  result.insert("grid_coords", torch::zeros({4, 3}, options));
  result.insert("grid_weights", torch::full({4}, .04, options));
  result.insert("coarse_0_atomic_coords", torch::zeros({1, 3}, options));
  result.insert("atomic_grid_weights", torch::full({4}, .05, options));
  result.insert("atomic_grid_sizes", torch::tensor({4}, torch::kInt64));
  result.insert("atomic_grid_size_bound_shape", torch::empty({4, 0}, torch::kInt64));
  if (nc == 2) {
    result.at("density")[0].fill_(.2);
    result.at("density")[1].fill_(-.15);
    result.at("kin")[0].fill_(-.1);
    result.at("kin")[1].fill_(.08);
    result.at("grad")[1][0].fill_(.1);
  } else {
    result.at("density")[0].fill_(.4);
    result.at("density")[3].fill_(.1);
    result.at("kin")[0].fill_(.1);
    result.at("kin")[3].fill_(.02);
    result.at("grad")[0][0].fill_(.1);
    result.at("grad")[3][0].fill_(-.1);
  }
  return result;
}

struct Evaluation {
  double energy;
  std::array<torch::Tensor, 7> gradients;
};

Evaluation evaluate(void *model, Inputs &input) {
  Evaluation result;
  std::array<void *, 7> output = {};
  char error[4096] = {};
  const int status = cppaw_skala_model_evaluate_mode(
      model, &input, 1, &result.energy, &output[0], &output[1], &output[2], &output[3],
      &output[4], &output[5], &output[6], error, sizeof(error));
  for (std::size_t i = 0; i < output.size(); ++i) {
    std::unique_ptr<torch::Tensor> value(static_cast<torch::Tensor *>(output[i]));
    if (value) result.gradients[i] = *value;
  }
  require(status == 0, error);
  require(std::isfinite(result.energy), "Nonfinite fixture energy");
  return result;
}

torch::Tensor original_energy(torch::jit::Module &model, const Inputs &data) {
  std::vector<c10::IValue> args;
  std::unordered_map<std::string, c10::IValue> kwargs{{"mol", data}};
  return (model.get_method("get_exc_density")(args, kwargs).toTensor()*data.at("grid_weights")).sum();
}

void collinear_parity(torch::jit::Module &original, void *corrected, Inputs data,
                      torch::Device device=torch::kCPU) {
  std::vector<torch::Tensor> leaves;
  const char *keys[] = {"density", "grad", "kin", "grid_coords", "grid_weights",
                        "coarse_0_atomic_coords", "atomic_grid_weights"};
  Inputs local;
  for (const auto &item : data) local.insert(item.key(), item.value().to(device).detach());
  for (const char *key : keys) {
    local.at(key).set_requires_grad(true);
    leaves.push_back(local.at(key));
  }
  std::optional<torch::jit::GraphOptimizerEnabledGuard> reference_guard;
  if (fixed_executor_requested()) reference_guard.emplace(false);
  const auto energy = original_energy(original, local);
  const auto gradients = torch::autograd::grad({energy}, leaves, {}, false, false, true);
  const double reference_energy = energy.item<double>();
  reference_guard.reset();
  auto output = evaluate(corrected, data);
  require(std::abs(output.energy-reference_energy) < 1.e-12, "Collinear fixture energy changed");
  double maximum = 0.;
  for (std::size_t i = 0; i < gradients.size(); ++i) {
    const auto expected = gradients[i].defined() ? gradients[i] : torch::zeros_like(leaves[i]);
    const auto actual = output.gradients[i].defined() ? output.gradients[i].to(device) :
                                                      torch::zeros_like(leaves[i]);
    require(torch::isfinite(expected).all().item<bool>() && torch::isfinite(actual).all().item<bool>(),
            "A collinear first derivative is nonfinite");
    const double difference = (actual-expected).abs().max().item<double>();
    maximum = std::max(maximum, difference);
    std::cout << "COLLINEAR_FIELD " << keys[i] << " maximum difference " << difference << '\n';
  }
  std::cout << "COLLINEAR_PARITY maximum first-adjoint difference " << maximum << '\n';
  require(maximum < 1.e-12, "A collinear first derivative changed");
}

void repeatability(void *model, Inputs data) {
  if (!fixed_executor_requested()) return;
  const bool external_mode = torch::jit::getGraphExecutorOptimize();
  const auto first = evaluate(model, data);
  require(torch::jit::getGraphExecutorOptimize() == external_mode, "Model call changed the caller's JIT mode");
  for (int repeat=0; repeat<4; ++repeat) {
    const auto current = evaluate(model, data);
    require(current.energy == first.energy, "Deterministic model energy changed after warmup");
    for (int i=0; i<7; ++i) {
      const auto &a = first.gradients[i], &b = current.gradients[i];
      require(a.defined() == b.defined(), "Repeated model adjoint became unavailable");
      if (a.defined()) {
        require(torch::isfinite(a).all().item<bool>(), "Nonfinite repeated model adjoint");
        require(torch::equal(a, b), "Deterministic model adjoint changed after warmup");
      }
    }
    require(torch::jit::getGraphExecutorOptimize() == external_mode, "Model call leaked its JIT mode");
  }
  Inputs invalid;
  for (const auto &item : data) invalid.insert(item.key(), item.value());
  invalid.insert_or_assign("density", torch::zeros({3, data.at("density").size(1)}, torch::kFloat64));
  bool rejected = false;
  try { evaluate(model, invalid); } catch (const std::exception &) { rejected = true; }
  require(rejected, "Invalid deterministic model input was accepted");
  require(torch::jit::getGraphExecutorOptimize() == external_mode, "Failed model call leaked its JIT mode");
  std::cout << "DETERMINISTIC REPEATED ENERGY AND ALL ADJOINTS UNCHANGED\n";
}

void flat_gradient(void *model) {
  const double axis[] = {0., 0., 1.}, weight[] = {1.};
  char error[4096] = {};
  require(cppaw_skala_model_set_spin_quadrature(model, 1, axis, weight, error, sizeof(error)) == 0,
          error);
  auto data = inputs(4);
  const auto output = evaluate(model, data);
  constexpr double epsilon = 1.e-5, g = .1, spin = -.1;
  const double first = 2*g/(g*g+epsilon);
  const double second = 2*(epsilon-g*g)/std::pow(g*g+epsilon, 2);
  const double expected = 4*.04*(1.5*first+.25*spin*(2/epsilon-second));
  const double actual = output.gradients[1][0][0].sum().item<double>();
  require(std::abs(actual-expected) < 1.e-10, "Flat-gradient Hessian contribution is incorrect");
  std::cout << "FLAT_GRADIENT_ANALYTIC error " << std::abs(actual-expected) << '\n';
}

Inputs complete_inputs() {
  auto data = inputs(4);
  const auto options = torch::TensorOptions().dtype(torch::kFloat64);
  data.insert_or_assign("density", torch::zeros({4, 16}, options));
  data.insert_or_assign("grad", torch::zeros({4, 3, 16}, options));
  data.insert_or_assign("kin", torch::zeros({4, 16}, options));
  data.at("density")[0].fill_(.4);
  data.at("density")[3].fill_(.1);
  data.at("kin")[0].fill_(.1);
  data.at("kin")[3].fill_(.02);
  data.at("grad")[0][0].fill_(.1);
  data.at("grad")[3][0].fill_(-.1);
  auto x = torch::zeros({16, 3}, options), w = torch::zeros({16}, options);
  for (int p = 0; p < 16; ++p) {
    const int j = p%8;
    const double angle = 2*std::acos(-1.)*j/8, radius = .35+.03*((j+1)%3);
    x[p][0] = radius*std::cos(angle);
    x[p][1] = radius*std::sin(angle);
    x[p][2] = 1.4*(p/8)+.02*(j-3);
    w[p] = .04+.001*(j+1);
  }
  auto centers = torch::zeros({2, 3}, options);
  centers[1][2] = 1.4;
  data.insert_or_assign("grid_coords", x);
  data.insert_or_assign("grid_weights", w);
  data.insert_or_assign("coarse_0_atomic_coords", centers);
  data.insert_or_assign("atomic_grid_weights", torch::full({16}, .05, options));
  data.insert_or_assign("atomic_grid_sizes", torch::tensor({8, 8}, torch::kInt64));
  data.insert_or_assign("atomic_grid_size_bound_shape", torch::empty({8, 0}, torch::kInt64));
  return data;
}

std::pair<double, double> scalar_scaling(torch::jit::Module &model, const Inputs &data,
                                       torch::Device device, double shift) {
  std::optional<torch::jit::GraphOptimizerEnabledGuard> reference_guard;
  if (fixed_executor_requested()) reference_guard.emplace(false);
  auto scale = torch::ones({}, torch::TensorOptions().dtype(torch::kFloat64).device(device)
                                                           .requires_grad(true));
  Inputs projected;
  torch::Tensor gradient;
  for (const auto &item : data) {
    auto value = item.value().to(device).detach().clone();
    if (item.key() == "grad") {
      value[0][0].add_(shift);
      value.set_requires_grad(true);
      gradient = value;
    }
    if (item.key() == "density" || item.key() == "grad" || item.key() == "kin") {
      value = torch::stack({.5*(value[0]+scale*value[3]), .5*(value[0]-scale*value[3])});
    }
    projected.insert(item.key(), value);
  }
  const auto energy = original_energy(model, projected);
  const auto response = torch::autograd::grad({energy}, {scale}, {}, true, true)[0];
  const auto effective = energy+response;
  const auto derivative = torch::autograd::grad({effective}, {gradient})[0];
  return {effective.item<double>(), derivative[0][0].sum().item<double>()};
}

void full_model(const std::filesystem::path &file, const std::string &device_name,
                bool diagnostic64) {
  require(device_name == "CPU" || device_name == "CUDA", "Expected CPU or CUDA");
  const auto device = device_name == "CUDA" ? torch::Device(torch::kCUDA, 0) : torch::Device(torch::kCPU);
  auto corrected = load(file, 2, device_name == "CUDA" ? 1 : 0);
  auto original = torch::jit::load(file.string(), device);
  original.eval();
  for (auto parameter : original.parameters()) parameter.set_requires_grad(false);
  auto data = complete_inputs();
  for (bool closed : {true, false}) {
    Inputs populations;
    for (const auto &item : data) {
      auto value = item.value();
      if (item.key() == "density" || item.key() == "grad" || item.key() == "kin") {
        const auto spin = closed ? torch::zeros_like(value[3]) : value[3];
        value = torch::stack({.5*(value[0]+spin), .5*(value[0]-spin)});
      }
      populations.insert(item.key(), value);
    }
    collinear_parity(original, corrected.get(), populations, device);
    repeatability(corrected.get(), populations);
  }
  const double axis[] = {0., 0., 1.}, weights[] = {1.};
  char error[4096] = {};
  require(cppaw_skala_model_set_spin_quadrature(corrected.get(), 1, axis, weights, error, sizeof(error)) == 0,
          error);
  const auto result = evaluate(corrected.get(), data);
  repeatability(corrected.get(), data);
  const double analytic = result.gradients[1][0][0].sum().item<double>();
  const auto baseline = scalar_scaling(original, data, device, 0.);
  require(std::isfinite(analytic) && std::isfinite(baseline.first) && std::isfinite(baseline.second),
          "Nonfinite full-model derivative or scalar-scaling energy");
  // Match the existing scalar-scaling reference's precision-dependent bounds.
  const double reference_limit = diagnostic64 ? 2.e-12 : 1.e-8;
  std::cout << "MODEL_MC_SCALAR_REFERENCE energy difference " << std::abs(result.energy-baseline.first)
            << " limit " << reference_limit << '\n';
  require(std::abs(result.energy-baseline.first) <= reference_limit,
          "Multicollinear energy disagrees with independent scalar scaling");
  std::cout << "MODEL_MC energy " << result.energy << " original analytic " << baseline.second
            << " corrected analytic " << analytic << '\n';
  double last_error = 0.;
  for (double step : {1.e-5, 3.e-6, 1.e-6}) {
    const double numeric = (scalar_scaling(original, data, device, step).first-
                            scalar_scaling(original, data, device, -step).first)/(2*step);
    require(std::isfinite(numeric), "Nonfinite full-model finite difference");
    last_error = std::abs(numeric-analytic);
    std::cout << "MODEL_MC_FD step " << step << " numeric " << numeric << " corrected error " << last_error
              << " original error " << std::abs(numeric-baseline.second) << '\n';
  }
  // The published FP32 network has a different finite-difference noise floor.
  // Only the unchanged-weight FP64 diagnostic is a strict Hessian FD test.
  if (diagnostic64) require(last_error <= 1.e-7, "Full-model flat-gradient derivative failed 1e-7 gate");
  std::cout << "SKALA FULL MODEL GRADIENT CHECKS PASSED " << (diagnostic64 ? "FD64" : "FIRST32") << '\n';
}

void atom_block_spin(const std::string &device) {
  require(device == "CPU" || device == "CUDA", "Expected CPU or CUDA");
  const auto file = std::filesystem::current_path()/"spin-blocks.fun";
  require(!std::filesystem::exists(file), "Refusing to overwrite spin-block fixture");
  torch::jit::Module fixture("AtomicSpinCoupling");
  fixture.define(R"JIT(
def get_exc_density(self, mol: Dict[str, Tensor]) -> Tensor:
    sizes = mol["atomic_grid_sizes"]
    a = mol["atomic_grid_weights"]
    ids = torch.repeat_interleave(torch.arange(sizes.size(0), device=sizes.device), sizes)
    d = mol["grid_coords"] - mol["coarse_0_atomic_coords"].index_select(0, ids)
    rho, grad, tau = mol["density"], mol["grad"], mol["kin"]
    b = rho[0]-rho[1] + a*(tau[0]-tau[1]) + ((grad[0]-grad[1])*d.t()).sum(0)
    moment = torch.zeros_like(sizes).to(b).index_add(0, ids, a*b)
    return a*(rho[0]+rho[1]) + b*moment.index_select(0, ids)
)JIT");
  const char *names[] = {"density","grad","kin","grid_coords","grid_weights",
                         "coarse_0_atomic_coords","atomic_grid_weights"};
  fixture.save(file.string(), {{"protocol_version","2"},
      {"features",R"(["density","grad","kin","grid_coords","grid_weights","coarse_0_atomic_coords","atomic_grid_weights","atomic_grid_sizes","atomic_grid_size_bound_shape"])"}});
  auto model = load(file, 0, device == "CUDA" ? 1 : 0);
  require(std::filesystem::remove(file), "Fixture cleanup failed");
  constexpr int n=5, nq=4;
  const auto opt=torch::TensorOptions().dtype(torch::kFloat64);
  Inputs data;
  std::vector<torch::Tensor> leaves;
  const std::vector<std::vector<int64_t>> shapes={{4,n},{4,3,n},{4,n},{n,3},{n},{2,3},{n}};
  for (size_t j=0; j<shapes.size(); ++j) {
    auto v=torch::empty(shapes[j],opt);
    auto t=(torch::arange(v.numel(),opt)*.73+j).sin().reshape_as(v)*.19;
    if (j == 4 || j == 6) t=t+.28;
    if (j == 0 || j == 2) t[0].add_(.8);
    data.insert(names[j],t.detach().clone());
    leaves.push_back(t.set_requires_grad(true));
  }
  data.insert("atomic_grid_sizes",torch::tensor({2,3},torch::kInt64));
  data.insert("atomic_grid_size_bound_shape",torch::empty({3,0},torch::kInt64));
  const double u=1./std::sqrt(3.);
  double axes[3*nq]={1,0,0, 0,1,0, 0,0,1, u,-u,u};
  double weights[nq]={.12,.25,.35,.28}, coeff[nq]={1.2,-.4,1.1,.7};
  double resp[3*nq]={.1,.4,-.3, .5,-.3,.2, -.2,.1,.8, .3,.2,-.5};
  char error[4096]={};
  require(cppaw_skala_model_set_spin_response_quadrature(model.get(),nq,axes,weights,coeff,resp,error,sizeof(error))==0,error);
  const auto actual=evaluate(model.get(),data);
  // Integrate the quadratic model analytically in the physical Pauli components.
  // Each atom contributes U^T M V with two independent weighted spin moments.
  auto M=torch::zeros({3,3},opt); double charge_weight=0.;
  for (int q=0; q<nq; ++q) {
    auto k=torch::from_blob(axes+3*q,{3},opt), r=torch::from_blob(resp+3*q,{3},opt);
    M.add_(weights[q]*(coeff[q]*torch::outer(k,k)+torch::outer(k,r)+torch::outer(r,k)));
    charge_weight+=weights[q]*coeff[q];
  }
  auto energy=charge_weight*(leaves[0][0]*leaves[4]*leaves[6]).sum();
  for (int atom=0; atom<2; ++atom) {
    const int first=atom == 0 ? 0 : 2, count=atom == 0 ? 2 : 3;
    const auto offset=leaves[3].narrow(0,first,count)-leaves[5][atom];
    const auto a=leaves[6].narrow(0,first,count), w=leaves[4].narrow(0,first,count);
    const auto v=(leaves[0].narrow(0,1,3).narrow(1,first,count)
                 +leaves[2].narrow(0,1,3).narrow(1,first,count)*a
                 +(leaves[1].narrow(0,1,3).narrow(2,first,count)*offset.t()).sum(1));
    energy=energy+torch::dot(torch::matmul(v,w),torch::matmul(M,torch::matmul(v,a)));
  }
  auto expected=torch::autograd::grad({energy},leaves);
  double diff=std::abs(actual.energy-energy.item<double>());
  for (size_t j=0; j<expected.size(); ++j) {
    require(expected[j].abs().max().item<double>()>1.e-6, "Degenerate atom-block reference");
    require(actual.gradients[j].defined() && torch::isfinite(actual.gradients[j]).all().item<bool>(),
            "Missing or nonfinite atom-block derivative");
    diff=std::max(diff,(expected[j]-actual.gradients[j].cpu()).abs().max().item<double>());
  }
  int64_t stats[3]; cppaw_skala_model_spin_batch_stats(model.get(),stats);
  require(stats[0]>=1 && stats[1]==(nq+stats[0]-1)/stats[0] && stats[2]==nq, "Spin batch execution mismatch");
  std::cout << "ATOM BLOCK SPIN: batch " << stats[0] << " calls " << stats[1] << " max error " << diff << '\n';
  require(std::isfinite(diff) && diff<2.e-12, "Nonlocal energy or field/geometry pullback changed");
}

int main(int argc, char **argv) {
  try {
    at::set_num_threads(1);
    std::cout << std::setprecision(17);
    if (argc==3 && std::string(argv[1])=="--atom-block-spin") {
      atom_block_spin(argv[2]);
      return 0;
    }
    if (argc == 4) {
      require(std::string(argv[3]) == "FD64" || std::string(argv[3]) == "FIRST32", "Expected FD64 or FIRST32");
      full_model(argv[1], argv[2], std::string(argv[3]) == "FD64");
      return 0;
    }
    require(argc == 1, "Usage: skala_gradient_features [MODEL CPU|CUDA FD64|FIRST32]");
    require(cppaw_skala_model_gradient_fixes(nullptr) == 0, "Null model count is not zero");
    std::vector<std::filesystem::path> files;
    files.push_back(fixture(1, "grad"));
    {
      auto first = load(files.back(), 2);
      auto original = torch::jit::load(files.back().string());
      collinear_parity(original, first.get(), inputs(2));
      repeatability(first.get(), inputs(2));
      flat_gradient(first.get());
      repeatability(first.get(), inputs(4));
      // A second handle must not inherit an already-executed JIT plan.
      auto second = load(files.back(), 2);
      flat_gradient(second.get());
    }
    files.push_back(fixture(2, "grad", true));
    flat_gradient(load(files.back(), 0).get());
    files.push_back(fixture(3, "grad", false, false));
    load(files.back(), 0);
    files.push_back(fixture(4, "grad.to(dtype=10)"));
    load(files.back(), 0);
    files.push_back(fixture(5, "grad.to(rho)"));
    load(files.back(), 0);
    files.push_back(fixture(6, "grad", false, true, true));
    load(files.back(), 0);
    files.push_back(fixture(7, "grad", false, true, false,
                            "torch.complex(torch.zeros_like(rho), rho)"));
    load(files.back(), 0);
    for (const auto &file : files) require(std::filesystem::remove(file), "Fixture cleanup failed");
    std::cout << "SKALA GRADIENT FEATURE CHECKS PASSED\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
