// Generated by rstantools.  Do not edit by hand.

/*
    bkmrGen is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    bkmrGen is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with bkmrGen.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_fit_poisson_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 24> locations_array__ =
  {" (found before start of program)",
  " (in 'fit_poisson', line 13, column 2 to column 17)",
  " (in 'fit_poisson', line 14, column 2 to column 23)",
  " (in 'fit_poisson', line 15, column 2 to column 20)",
  " (in 'fit_poisson', line 23, column 2 to column 123)",
  " (in 'fit_poisson', line 26, column 2 to column 24)",
  " (in 'fit_poisson', line 27, column 2 to column 24)",
  " (in 'fit_poisson', line 28, column 2 to column 25)",
  " (in 'fit_poisson', line 29, column 2 to column 26)",
  " (in 'fit_poisson', line 2, column 2 to column 17)",
  " (in 'fit_poisson', line 3, column 2 to column 17)",
  " (in 'fit_poisson', line 4, column 2 to column 17)",
  " (in 'fit_poisson', line 5, column 9 to column 10)",
  " (in 'fit_poisson', line 5, column 12 to column 13)",
  " (in 'fit_poisson', line 5, column 2 to column 17)",
  " (in 'fit_poisson', line 6, column 18 to column 19)",
  " (in 'fit_poisson', line 6, column 13 to column 14)",
  " (in 'fit_poisson', line 6, column 2 to column 21)",
  " (in 'fit_poisson', line 7, column 8 to column 9)",
  " (in 'fit_poisson', line 7, column 2 to column 11)",
  " (in 'fit_poisson', line 8, column 2 to column 20)",
  " (in 'fit_poisson', line 13, column 9 to column 10)",
  " (in 'fit_poisson', line 15, column 9 to column 10)",
  " (in 'fit_poisson', line 23, column 9 to column 10)"};
#include <stan_meta_header.hpp>
class model_fit_poisson final : public model_base_crtp<model_fit_poisson> {
private:
  int N;
  int p;
  int d;
  Eigen::Matrix<double,-1,-1> X_data__;
  std::vector<Eigen::Matrix<double,1,-1>> Z;
  std::vector<int> y;
  double rho;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X{nullptr, 0, 0};
public:
  ~model_fit_poisson() {}
  model_fit_poisson(stan::io::var_context& context__, unsigned int
                    random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_fit_poisson_namespace::model_fit_poisson";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 9;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 9;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 9;
      stan::math::check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 10;
      context__.validate_dims("data initialization", "p", "int",
        std::vector<size_t>{});
      p = std::numeric_limits<int>::min();
      current_statement__ = 10;
      p = context__.vals_i("p")[(1 - 1)];
      current_statement__ = 10;
      stan::math::check_greater_or_equal(function__, "p", p, 1);
      current_statement__ = 11;
      context__.validate_dims("data initialization", "d", "int",
        std::vector<size_t>{});
      d = std::numeric_limits<int>::min();
      current_statement__ = 11;
      d = context__.vals_i("d")[(1 - 1)];
      current_statement__ = 11;
      stan::math::check_greater_or_equal(function__, "d", d, 1);
      current_statement__ = 12;
      stan::math::validate_non_negative_index("X", "N", N);
      current_statement__ = 13;
      stan::math::validate_non_negative_index("X", "d", d);
      current_statement__ = 14;
      context__.validate_dims("data initialization", "X", "double",
        std::vector<size_t>{static_cast<size_t>(N), static_cast<size_t>(d)});
      X_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, d,
                   std::numeric_limits<double>::quiet_NaN());
      new (&X) Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_data__.data(), N, d);
      {
        std::vector<local_scalar_t__> X_flat__;
        current_statement__ = 14;
        X_flat__ = context__.vals_r("X");
        current_statement__ = 14;
        pos__ = 1;
        current_statement__ = 14;
        for (int sym1__ = 1; sym1__ <= d; ++sym1__) {
          current_statement__ = 14;
          for (int sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 14;
            stan::model::assign(X, X_flat__[(pos__ - 1)],
              "assigning variable X", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 14;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 15;
      stan::math::validate_non_negative_index("Z", "N", N);
      current_statement__ = 16;
      stan::math::validate_non_negative_index("Z", "p", p);
      current_statement__ = 17;
      context__.validate_dims("data initialization", "Z", "double",
        std::vector<size_t>{static_cast<size_t>(N), static_cast<size_t>(p)});
      Z = std::vector<Eigen::Matrix<double,1,-1>>(N,
            Eigen::Matrix<double,1,-1>::Constant(p,
              std::numeric_limits<double>::quiet_NaN()));
      {
        std::vector<local_scalar_t__> Z_flat__;
        current_statement__ = 17;
        Z_flat__ = context__.vals_r("Z");
        current_statement__ = 17;
        pos__ = 1;
        current_statement__ = 17;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 17;
          for (int sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 17;
            stan::model::assign(Z, Z_flat__[(pos__ - 1)],
              "assigning variable Z", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 17;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 18;
      stan::math::validate_non_negative_index("y", "N", N);
      current_statement__ = 19;
      context__.validate_dims("data initialization", "y", "int",
        std::vector<size_t>{static_cast<size_t>(N)});
      y = std::vector<int>(N, std::numeric_limits<int>::min());
      current_statement__ = 19;
      y = context__.vals_i("y");
      current_statement__ = 20;
      context__.validate_dims("data initialization", "rho", "double",
        std::vector<size_t>{});
      rho = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 20;
      rho = context__.vals_r("rho")[(1 - 1)];
      current_statement__ = 20;
      stan::math::check_greater_or_equal(function__, "rho", rho, 0);
      current_statement__ = 21;
      stan::math::validate_non_negative_index("beta", "d", d);
      current_statement__ = 22;
      stan::math::validate_non_negative_index("h_tilde", "N", N);
      current_statement__ = 23;
      stan::math::validate_non_negative_index("ystar", "N", N);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = d + 1 + N;
  }
  inline std::string model_name() const final {
    return "model_fit_poisson";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_fit_poisson_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,1> beta =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(d, DUMMY_VAR__);
      current_statement__ = 1;
      beta = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(d);
      local_scalar_t__ lambda = DUMMY_VAR__;
      current_statement__ = 2;
      lambda = in__.template read_constrain_lb<local_scalar_t__,
                 jacobian__>(0, lp__);
      Eigen::Matrix<local_scalar_t__,-1,1> h_tilde =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N, DUMMY_VAR__);
      current_statement__ = 3;
      h_tilde = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(N);
      Eigen::Matrix<local_scalar_t__,-1,1> ystar =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N, DUMMY_VAR__);
      current_statement__ = 4;
      stan::model::assign(ystar,
        stan::math::add(
          stan::math::multiply(
            stan::math::cholesky_decompose(
              stan::math::add(stan::math::cov_exp_quad(Z, lambda, rho),
                stan::math::diag_matrix(stan::math::rep_vector(1e-6, N)))),
            h_tilde), stan::math::multiply(X, beta)),
        "assigning variable ystar");
      {
        current_statement__ = 5;
        lp_accum__.add(stan::math::gamma_lpdf<propto__>(lambda, 1, .1));
        current_statement__ = 6;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(beta, 0, 100));
        current_statement__ = 7;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(h_tilde, 0, 1));
        current_statement__ = 8;
        lp_accum__.add(stan::math::poisson_lpmf<propto__>(y,
                         stan::math::exp(ystar)));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_fit_poisson_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,1> beta =
        Eigen::Matrix<double,-1,1>::Constant(d,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      beta = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(d);
      double lambda = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      lambda = in__.template read_constrain_lb<local_scalar_t__,
                 jacobian__>(0, lp__);
      Eigen::Matrix<double,-1,1> h_tilde =
        Eigen::Matrix<double,-1,1>::Constant(N,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 3;
      h_tilde = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(N);
      Eigen::Matrix<double,-1,1> ystar =
        Eigen::Matrix<double,-1,1>::Constant(N,
          std::numeric_limits<double>::quiet_NaN());
      out__.write(beta);
      out__.write(lambda);
      out__.write(h_tilde);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 4;
      stan::model::assign(ystar,
        stan::math::add(
          stan::math::multiply(
            stan::math::cholesky_decompose(
              stan::math::add(stan::math::cov_exp_quad(Z, lambda, rho),
                stan::math::diag_matrix(stan::math::rep_vector(1e-6, N)))),
            h_tilde), stan::math::multiply(X, beta)),
        "assigning variable ystar");
      if (emit_transformed_parameters__) {
        out__.write(ystar);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,1> beta =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(d, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(beta,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,1>>(d),
        "assigning variable beta");
      out__.write(beta);
      local_scalar_t__ lambda = DUMMY_VAR__;
      current_statement__ = 2;
      lambda = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, lambda);
      Eigen::Matrix<local_scalar_t__,-1,1> h_tilde =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(h_tilde,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,1>>(N),
        "assigning variable h_tilde");
      out__.write(h_tilde);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "beta", "double",
        std::vector<size_t>{static_cast<size_t>(d)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "lambda", "double",
        std::vector<size_t>{});
      current_statement__ = 3;
      context__.validate_dims("parameter initialization", "h_tilde",
        "double", std::vector<size_t>{static_cast<size_t>(N)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,1> beta =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(d, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> beta_flat__;
        current_statement__ = 1;
        beta_flat__ = context__.vals_r("beta");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= d; ++sym1__) {
          current_statement__ = 1;
          stan::model::assign(beta, beta_flat__[(pos__ - 1)],
            "assigning variable beta", stan::model::index_uni(sym1__));
          current_statement__ = 1;
          pos__ = (pos__ + 1);
        }
      }
      out__.write(beta);
      local_scalar_t__ lambda = DUMMY_VAR__;
      current_statement__ = 2;
      lambda = context__.vals_r("lambda")[(1 - 1)];
      out__.write_free_lb(0, lambda);
      Eigen::Matrix<local_scalar_t__,-1,1> h_tilde =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> h_tilde_flat__;
        current_statement__ = 3;
        h_tilde_flat__ = context__.vals_r("h_tilde");
        current_statement__ = 3;
        pos__ = 1;
        current_statement__ = 3;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 3;
          stan::model::assign(h_tilde, h_tilde_flat__[(pos__ - 1)],
            "assigning variable h_tilde", stan::model::index_uni(sym1__));
          current_statement__ = 3;
          pos__ = (pos__ + 1);
        }
      }
      out__.write(h_tilde);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"beta", "lambda", "h_tilde"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"ystar"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(d)},
                std::vector<size_t>{},
                std::vector<size_t>{static_cast<size_t>(N)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(N)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= d; ++sym1__) {
      param_names__.emplace_back(std::string() + "beta" + '.' +
        std::to_string(sym1__));
    }
    param_names__.emplace_back(std::string() + "lambda");
    for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
      param_names__.emplace_back(std::string() + "h_tilde" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        param_names__.emplace_back(std::string() + "ystar" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= d; ++sym1__) {
      param_names__.emplace_back(std::string() + "beta" + '.' +
        std::to_string(sym1__));
    }
    param_names__.emplace_back(std::string() + "lambda");
    for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
      param_names__.emplace_back(std::string() + "h_tilde" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        param_names__.emplace_back(std::string() + "ystar" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {}
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"beta\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(d) + "},\"block\":\"parameters\"},{\"name\":\"lambda\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"h_tilde\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "},\"block\":\"parameters\"},{\"name\":\"ystar\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "},\"block\":\"transformed_parameters\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"beta\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(d) + "},\"block\":\"parameters\"},{\"name\":\"lambda\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"h_tilde\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "},\"block\":\"parameters\"},{\"name\":\"ystar\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "},\"block\":\"transformed_parameters\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((d + 1) + N);
    const size_t num_transformed = emit_transformed_parameters * (N);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((d + 1) + N);
    const size_t num_transformed = emit_transformed_parameters * (N);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_fit_poisson_namespace::model_fit_poisson;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_fit_poisson_namespace::profiles__;
}
#endif
#endif
