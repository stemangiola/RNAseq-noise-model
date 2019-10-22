template <typename T_a, typename T_b, typename T_c>
typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
  poisson_tweedie_raw(const std::vector<int> &y, const T_a &a, const T_b &b, const T_c &c, std::ostream* pstream__) {

    using T_partials_return = stan::partials_return_type<T_a, T_b, T_c>;

    using stan::is_constant;

    //TODO check input

    double _a = value_of(a);
    double _b = value_of(b);
    double _c = value_of(c);


    double value = poisson_tweedie_raw_lpmf_internal(y, _a, _b, _c, pstream__);

    operands_and_partials<T_a, T_b, T_c> ops_partials(a, b, c);

    if (!is_constant<T_a>::value) {
      auto f = [&](std::complex<double> complex_a) { return poisson_tweedie_raw_lpmf_internal(y, complex_a, _b, _c, pstream__); };
      double grad_a = boost::math::tools::complex_step_derivative(f, _a);
      ops_partials.edge1_.partials_[0] += grad_a;
    }

    if (!is_constant<T_b>::value) {
      auto f = [&](std::complex<double> complex_b) { return poisson_tweedie_raw_lpmf_internal(y, _a, complex_b, _c, pstream__); };
      double grad_b = boost::math::tools::complex_step_derivative(f, _b);
      ops_partials.edge2_.partials_[0] += grad_b;
    }

    if (!is_constant<T_c>::value) {
      auto f = [&](std::complex<double> complex_c) { return poisson_tweedie_raw_lpmf_internal(y, _a, _b, complex_c, pstream__); };
      double grad_c = boost::math::tools::complex_step_derivative(f, _c);
      ops_partials.edge3_.partials_[0] += grad_c;
    }

    return ops_partials.build(value);
  }

double poisson_tweedie_raw_da(const std::vector<int> &y, const double &a, const double &b, const double &c, std::ostream* pstream__) {
  auto f = [&](std::complex<double> complex_a) { return poisson_tweedie_raw_lpmf_internal(y, complex_a, b, c, pstream__); };
  double grad_a = boost::math::tools::complex_step_derivative(f, a);
  return grad_a;
}

/*
 template <typename T_a, typename T_b, typename T_c>
 inline typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
 poisson_tweedie_raw_lpmf(const std::vector<int> &y, const T_a &a, const T_b &b, const T_c &c, std::ostream* pstream__) {
 return poisson_tweedie_raw_lpmf<false>(y, a, b, c);
 }
 */
