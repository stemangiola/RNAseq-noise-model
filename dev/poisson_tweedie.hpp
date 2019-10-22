namespace local_override {
  template <typename _Tp>
  inline int my_is_nan(const std::complex<_Tp> &c) {
    return std::isnan(c.imag()) || std::isnan(c.real());
  }


  inline std::complex<double> my_log1p(const std::complex<double> &x)
  {
    double abs_x = std::abs(x);
    if (!my_is_nan(x)) {
      check_greater("log1p", "x", abs_x, -1);
    }

    if (abs_x > 1e-4)
    {
      return std::log(1.0 + x);
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
    return (-0.5*x + 1.0)*x;
  }

  template<typename T> inline
    T my_expm1(const T &x) {
      return stan::math::expm1(x);
  }

  template<>
  inline std::complex<double> my_expm1<std::complex<double>>(const std::complex<double> &x)
  {
    if (std::abs(x) < 1e-5)
      return x + 0.5*x*x;
    else
      return std::exp(x) - 1.0;
  }

  template<typename T> inline
    T my_log1m(const T &x) {
      return stan::math::log1m(x);
    }

  template<>
  inline std::complex<double> my_log1m<std::complex<double>>(const std::complex<double> &x) {
    if (!my_is_nan(x)) {
      check_less_or_equal("log1m", "x", std::abs(x), 1);
    }
    return my_log1p(-x);
  }

  inline std::complex<double> my_log1m_exp(const std::complex<double> &a) {
    double abs_a = std::abs(a);
    if (abs_a >= 0) {
      return NOT_A_NUMBER;
    } else if (abs_a > -0.693147) {
      return std::log(-my_expm1(a));  // 0.693147 ~= log(2)
    } else {
      return my_log1m(std::exp(a));
    }
  }

  template<typename T1, typename T2>
  inline typename boost::math::tools::promote_args<T1, T2>::type
  my_log_diff_exp(const T1 &x, const T2 &y) {
    stan::math::log_diff_exp(x,y);
  }

  template<>
  inline std::complex<double> my_log_diff_exp<std::complex<double>, std::complex<double>>(const std::complex<double> &x, const std::complex<double> &y) {
    double abs_x = std::abs(x);
    double abs_y = std::abs(y);
    if (std::abs(x - y) < 0) {
      return (std::abs(x) < INFTY && x == y) ? NEGATIVE_INFTY : NOT_A_NUMBER;
    }
    return x + my_log1m_exp(y - x);
  }


}



template <typename T_a, typename T_b, typename T_c>
typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
  poisson_tweedie_raw_lpmf_internal(const std::vector<int> &x, const T_a &a, const T_b &b, const T_c &c, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T_a, T_c>::type T_ac;
    typedef typename boost::math::tools::promote_args<T_a, T_b, T_c>::type T_ret;

    using std::log;
    using std::exp;
    using std::pow;

    //TODO check inputs
    /*
     if(pstream__ != 0) {
     *pstream__ << "Baf:" << a << " " << b << " " << c << std::endl;
     }
     */
    int max_x = *std::max_element(x.begin(), x. end());

    T_b log_b = std::log(b);
    T_c log_c = std::log(c);

    std::vector<T_ac> log_r;
    std::vector<T_ret> log_p;

    log_r.reserve(max_x);
    log_p.reserve(max_x + 1);

    T_ret log_p0;
    if(a == 0.0) {
      log_p0 = b * std::log(1.0 - c);
    } else if (std::abs(a) < 1e-4) {
      //Using Taylor expansion of order 2
      T_c log_1mc = local_override::my_log1m(c);
      log_p0 = log_1mc +
        0.5 * a * boost::math::pow<2>(log_1mc) +
        (boost::math::pow<2>(a) * boost::math::pow<3>(log_1mc))/ 6.0;
    } else {
      log_p0 = b * (std::pow(1.0 - c, a) - 1.0) / a;
    }
    log_p.push_back(log_p0);


    if(max_x > 0) {
      log_r.push_back(std::log(1.0 - a) + log_c);

      log_p.push_back(log_b + log_c + log_p0);

      if(max_x > 1) {
        std::vector<T_ret> summands;
        summands.reserve(max_x);

        for(int k = 1; k < max_x; ++k) {
          T_ac new_log_r = log(k - 1.0 + a) - log((double)k + 1) + log_c + log_r.back();
          log_r.push_back(new_log_r);

          summands.push_back(log_b + log_c + log_p.back());
          for(int j = 1; j<= k; ++j) {
            double log_j = std::log((double)j);
            summands.push_back(log_j + log_r[k - j] + log_p[j]);
          }
          double log_kp1 = std::log((double)k + 1);
          log_p.push_back(-log_kp1 + stan::math::log_sum_exp(summands));
          summands.clear();
        }
      }
    }

    /*
     if(pstream__ != 0) {
     auto print_func = [&](double x) { *pstream__ << x << " "; };
     std::for_each(log_p.begin(), log_p.end(), print_func);
     *pstream__ << std::endl;
     std::for_each(log_r.begin(), log_r.end(), print_func);
     *pstream__ << std::endl;
     }*/

    T_ret lpmf = 0;
    for(int xx : x) {
      lpmf += log_p[xx];
    }

    return lpmf;
  }


template <bool propto = false, typename T_a, typename T_b, typename T_c>
typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
  poisson_tweedie_raw_lpmf(const std::vector<int> &y, const T_a &a, const T_b &b, const T_c &c, std::ostream* pstream__) {

  using stan::partials_return_type;
  using stan::is_constant_all;
  using stan::value_of;

  using T_partials_return = partials_return_type<T_a, T_b, T_c>;


  //TODO check input

  double _a = value_of(a);
  double _b = value_of(b);
  double _c = value_of(c);

  T_partials_return value(poisson_tweedie_raw_lpmf_internal(y, _a, _b, _c, pstream__));

  operands_and_partials<T_a, T_b, T_c> ops_partials(a, b, c);

  if (!is_constant_all<T_a>::value) {
    auto f = [&](std::complex<double> complex_a) { return poisson_tweedie_raw_lpmf_internal(x, complex_a, _b, _c, pstream__); };
    double grad_a = boost::math::tools::complex_step_derivative(f, value_of(a));
    ops_partials.edge1_.partials_[0] += grad_a;
  }

  if (!is_constant_all<T_b>::value) {
    throw std::domain_error("Not supported");
  }

  if (!is_constant_all<T_c>::value) {
    throw std::domain_error("Not supported");
  }

  return ops_partials.build(value);
}


template <typename T_a, typename T_b, typename T_c>
inline typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
  poisson_tweedie_raw_lpmf(const std::vector<int> &y, const T_a &a, const T_b &b, const T_c &c, std::ostream* pstream__) {
  return poisson_tweedie_raw_lpmf<false>(y, a, b, c);
}
