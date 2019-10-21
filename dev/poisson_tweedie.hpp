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

  inline std::complex<double> my_log1m(const std::complex<double> &x) {
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
    /*
     if(pstream__ != 0) {
     *pstream__ << "Baf:" << a << " " << b << " " << c << std::endl;
     }
     */
    int max_x = *std::max_element(x.begin(), x. end());
    T_a _a = a;
    T_b _b = b;
    T_c _c = c;

    T_a log_a = std::log(_a);
    T_b log_b = std::log(_b);
    T_c log_c = std::log(_c);

    std::vector<T_ac> log_r;
    std::vector<T_ret> log_p;

    log_r.reserve(max_x);
    log_p.reserve(max_x + 1);

    T_ret log_p0;
    if(_a == 0.0) {
      log_p0 = _b * std::log(1.0 - _c);
    } else {
      log_p0 = _b * (std::pow(1.0 - _c, _a) - 1.0) / _a;
      //log_p0 = b * (local_override::my_expm1(_a * log(1.0 - _c) / _a));
    }
    log_p.push_back(log_p0);


    if(max_x > 0) {
      log_r.push_back(std::log(1.0 - _a) + log_c);

      log_p.push_back(log_b + log_c + log_p0);

      if(max_x > 1) {
        std::vector<T_ret> summands;
        summands.reserve(max_x);

        for(int k = 1; k < max_x; ++k) {
          T_ac new_log_r = log(k - 1.0 + _a) - log((double)k + 1) + log_c + log_r.back();
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

    //return {lpmf, dlog_p_da };
    return lpmf;
  }


template <typename T_a, typename T_b, typename T_c>
typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
  poisson_tweedie_raw(const std::vector<int> &x, const T_a &a, const T_b &b, const T_c &c, std::ostream* pstream__);

template <>
double
  poisson_tweedie_raw<double, double, double>(const std::vector<int> &x, const double &a, const double &b, const double &c, std::ostream* pstream__) {

    return poisson_tweedie_raw_lpmf_internal(x, a, b, c, pstream__);
  }

/*
double poisson_tweedie_raw_da(const std::vector<int> &x, const double &a, const double &b, const double &c, std::ostream* &pstream__) {
  stan::math::var _a(a);
  return poisson_tweedie_raw_lpmf_internal(x, _a, b, c, pstream__).dlpmf_da;
}
 */

double poisson_tweedie_raw_da(const std::vector<int> &x, const double &a, const double &b, const double &c, std::ostream* pstream__) {
  auto f = [&](std::complex<double> _a) { return poisson_tweedie_raw_lpmf_internal(x, _a, b, c, pstream__); };
  return boost::math::tools::complex_step_derivative(f, a);
}

