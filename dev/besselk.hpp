// The formulas and code are based on
// https://github.com/stan-dev/stan/wiki/Stan-Development-Meeting-Agenda/0ca4e1be9f7fc800658bfbd97331e800a4f50011
// Which is in turn based on Equation 26 of Rothwell: Computation of the
// logarithm of Bessel functions of complex argument and fractional order
// https://scholar.google.com/scholar?cluster=2908870453394922596&hl=en&as_sdt=5,33&sciodt=0,33

namespace {

using stan::math::domain_error;
  using stan::math::is_inf;
  using stan::math::is_nan;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;
  using stan::math::var;

  template <typename T_v, typename T_z, typename T_u>
  class inner_integral {
  private:
    T_v v;
    T_z z;

  public:
    typedef typename boost::math::tools::promote_args<T_v, T_z, T_u>::type T_Ret;

    inner_integral(const T_v &v, const T_z &z) : v(v), z(z) {}

    inline T_Ret operator()(const T_u &u, const T_u &) const {
      auto v_ = fabs(v);
      auto v_mhalf = v_ - 0.5;
      auto neg2v_m1 = -2 * v_ - 1;
      auto beta = 16.0 / (2 * v_ + 1);

      T_Ret value;
      T_Ret uB = pow(u, beta);
      T_Ret first
        = beta * exp(-uB) * pow(2 * z + uB, v_mhalf) * boost::math::pow<7>(u);
      T_Ret second = exp(-1.0 / u);
      if (second > 0)
        second = second * pow(u, neg2v_m1) * pow(2 * z * u + 1, v_mhalf);
      value = first + second;

      return value;
    }
  };

  // Uses nested autodiff to get gradient with respect to v
  // Gradient with respect to z can be computed analytically
  // Code modified from gradient_of_f
  template <typename T>
  class inner_integral_grad_v {
  private:
    T v;
    T z;

  public:
    inner_integral_grad_v(const T &v, const T &z) : v(v), z(z) {}

    inline T operator()(const T &u, const T &uc) const {
      double gradient = 0.0;
      start_nested();
      var v_var(stan::math::value_of(v));
      try {
        auto f = inner_integral<var, T, T>(v_var, z);
        var fx = f(u, uc);
        fx.grad();
        gradient = v_var.adj();
        if (is_nan(gradient)) {
          if (fx.val() == 0) {
            gradient = 0;
          } else {
            domain_error("inner_integral_grad_v",
                         "The gradient of inner_integral is nan", 0, "", "");
          }
        }
      } catch (const std::exception &e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
      return gradient;
    }
  };

  double compute_inner_integral_with_gradient(const double &v, const double &z) {
    double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());

    double error;
    double L1;
    size_t levels;

    boost::math::quadrature::tanh_sinh<double> integrator;

    auto f = inner_integral<double, double, double>(v, z);
    double integral = integrator.integrate(f, 0.0, 1.0, relative_tolerance,
                                           &error, &L1, &levels);

    if (error > 1e-6 * L1) {
      domain_error("compute_inner_integral_with_gradient",
                   "error estimate of integral / L1 ", error / L1, "",
                   "is larger than 1e-6");
    }
    return integral;
  }

  var compute_inner_integral_with_gradient(const var &v, const double &z) {
    double integral = compute_inner_integral_with_gradient(
      stan::math::value_of(v), stan::math::value_of(z));

    std::vector<stan::math::var> theta_concat;
    std::vector<double> dintegral_dtheta;

    theta_concat.push_back(v);
    auto f = inner_integral_grad_v<double>(v.val(), z);

    double error;
    double L1;
    size_t levels;
    double condition_number;
    double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());
    boost::math::quadrature::tanh_sinh<double> integrator;

    dintegral_dtheta.push_back(integrator.integrate(
        f, 0.0, 1.0, relative_tolerance, &error, &L1, &levels));

    if (error > 1e-6 * L1) {
      domain_error("compute_inner_integral_with_gradient",
                   "error estimate of integral / L1 ", error / L1, "",
                   "is larger than 1e-6");
    }

    return stan::math::precomputed_gradients(integral, theta_concat,
                                             dintegral_dtheta);
  }

  template <typename T_v, typename T_z>
  typename boost::math::tools::promote_args<T_v, T_z>::type compute_lead(
      const T_v &v, const T_z &z) {
    typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

    using std::exp;
    using std::fabs;
    using std::log;
    using std::pow;

    if (v == 0.5)
      return 0.5 * log(pi() / (2 * z)) - z;
    const T_v v_ = fabs(v);
    const T_Ret lead = 0.5 * log(pi()) - stan::math::lgamma(v_ + 0.5) - v_ * log(2 * z) - z;
    if (is_inf(lead))
      return -z + 0.5 * log(0.5 * pi() / z);

    return lead;
  }

  void check_params(const double& v, const double& z) {
    const char *function = "log_modified_bessel_second_kind_frac";
    if(!std::isfinite(v)) {
      stan::math::domain_error(function,
                               "v must be finite", v, "");
    }
    if(!std::isfinite(z)) {
      stan::math::domain_error(function,
                               "z must be finite", z, "");
    }
    if(z < 0) {
      stan::math::domain_error(function,
                               "z is negative", z, "");
    }
  }

  // Using the first two terms from
  // Sidi & Hoggan 2011 "ASYMPTOTICS OF MODIFIED BESSELFUNCTIONS OF HIGH ORDER"
  // https://pdfs.semanticscholar.org/d659/ecd3e66fefb2c2fafe16bf3086fcf41c573d.pdf
  // Without the second term, the approximation is worse and third term
  // did not help for v in [-1000, 1000].
  // With only the first term the formula reduces to 1.10 of
  // Temme, Journal of Computational Physics, vol 19, 324 (1975)
  // https://doi.org/10.1016/0021-9991(75)90082-0
  template <typename T_v, typename T_z>
  typename boost::math::tools::promote_args<T_v, T_z>::type
    asymptotic_large_v(const T_v &v, const T_z &z) {
      //return 0.5 * (log(stan::math::pi()) - log(2) - log(v)) - v * (log(z) - log(2) - log(v));
      T_v v_ = fabs(v);
      return stan::math::LOG_2 - v_ * (log(z) - stan::math::LOG_2) + stan::math::lgamma(v_)
        + log(1 + (0.25 * boost::math::pow<2>(z) / v) )
        //Third term, currently removed
        //+ (- 0.25 * boost::math::pow<2>(z) +
        //    0.5 * 0.125 * boost::math::pow<4>(z)) / boost::math::pow<2>(v)
        ;
      ;
    }
}  // namespace


const double log_modified_bessel_second_kind_frac_large_v_bound = 80;


template <typename T_v>
T_v log_modified_bessel_second_kind_frac(const T_v &v, const double &z) {
  check_params(value_of(v), value_of(z));
  if(fabs(value_of(v)) > log_modified_bessel_second_kind_frac_large_v_bound) {
    return asymptotic_large_v(v, z);
  }
  T_v lead = compute_lead(v, z);
  T_v Q = compute_inner_integral_with_gradient(v, z);
  return lead + log(Q);
}

template <typename T_v>
var log_modified_bessel_second_kind_frac(const T_v &v, const var &z) {
  check_params(value_of(v), value_of(z));

  if(fabs(value_of(v)) > log_modified_bessel_second_kind_frac_large_v_bound) {
    return asymptotic_large_v(v, z);
  }

  T_v lead = compute_lead(v, z.val());

  T_v Q = compute_inner_integral_with_gradient(v, z.val());

  T_v value = lead + log(Q);

  double value_vm1
    = log_modified_bessel_second_kind_frac(value_of(v) - 1, z.val());
  double gradient_dz
    = -std::exp(value_vm1 - value_of(value)) - value_of(v) / z.val();

  std::vector<var> operands;
  std::vector<double> gradients;
  if (stan::is_var<T_v>::value) {
    // A trick to combine the autodiff gradient with precomputed_gradients
    operands.push_back(value);
    gradients.push_back(1);
  }
  operands.push_back(z);
  gradients.push_back(gradient_dz);

  return precomputed_gradients(value_of(value), operands, gradients);
}

template <typename T0__, typename T1__>
typename boost::math::tools::promote_args<T0__, T1__>::type
  log_besselk_frac(const T0__& v, const T1__& z, std::ostream *msgs) {
    return(log_modified_bessel_second_kind_frac(v,z));


}
