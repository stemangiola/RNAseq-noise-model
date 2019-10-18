#include<algorithm>

struct TweedieData {
  double lpmf;
  double dlpmf_da;
};

template <bool propto = true, typename T_a, typename T_b, typename T_c>
TweedieData
  poisson_tweedie_raw_lpmf_internal(const std::vector<int> &x, const T_a &a, const T_b &b, const T_c &c, std::ostream* &pstream__) {
    typedef typename boost::math::tools::promote_args<T_a, T_b, T_c>::type T_ret;

    using stan::is_var;
    /*
     if(pstream__ != 0) {
     *pstream__ << "Baf:" << a << " " << b << " " << c << std::endl;
     }
     */
    int max_x = *std::max_element(x.begin(), x. end());
    double _a = value_of(a);
    double _b = value_of(b);
    double _c = value_of(c);

    double log_b = std::log(_b);
    double log_c = std::log(_c);

    std::vector<double> log_r;
    std::vector<double> log_p;
    std::vector<double> dlog_r_da;
    std::vector<double> log_dlog_p_da;

    log_r.reserve(max_x);
    log_p.reserve(max_x + 1);
    dlog_r_da.reserve(max_x);

    if(stan::is_var<T_a>::value) {
      log_dlog_p_da.reserve(max_x + 1);
    }

    double log_p0;
    if(_a == 0) {
      log_p0 = _b * std::log(1 - _c);
    } else {
      log_p0 = _b * (std::pow(1 - _c, _a) - 1) / _a;
    }
    log_p.push_back(log_p0);

    if(stan::is_var<T_a>::value) {
      double log_dlog_p0_da;
      if(_a == 0) {
        log_dlog_p0_da = std::log(0.5) + log_b + 2 * std::log(std::abs(std::log(1 - _c)));
      } else {
        double one_m_c_to_a = std::pow((1 - _c),_a);
        log_dlog_p0_da = log_b + std::log(-one_m_c_to_a + _a * one_m_c_to_a * std::log(1 - _c) + 1)
          - 2* std::log(std::abs(_a));
      }
      log_dlog_p_da.push_back(log_dlog_p0_da);
    }

    if(max_x > 0) {
      log_r.push_back(std::log(1 - _a) + log_c);

      log_p.push_back(log_b + log_c + log_p0);

      if(stan::is_var<T_a>::value) {
        dlog_r_da.push_back(- 1 / (1- _a));
        log_dlog_p_da.push_back(log_dlog_p_da.back());
      }

      if(max_x > 1) {
        std::vector<double> summands;
        std::vector<double> summands_da;
        summands.reserve(max_x);

        if(stan::is_var<T_a>::value) {
          summands_da.reserve(max_x * 2 - 1);
        }
        for(int k = 1; k < max_x; ++k) {
          log_r.push_back(log(k - 1 + _a) - log((double)k + 1) + log_c + log_r.back());
          dlog_r_da.push_back(1 / (k - 1 + _a) + dlog_r_da.back());

          summands.push_back(log_b + log_c + log_p.back());
          if(stan::is_var<T_a>::value) {
            dlog_r_da.push_back(dlog_r_da.back());
            summands_da.push_back(log_b + log_c + log_p.back() + std::exp(log_dlog_p_da.back()));
          }
          for(int j = 1; j<= k; ++j) {
            double log_j = std::log((double)j);
            summands.push_back(log_j + log_r[k - j] + log_p[j]);
            if(stan::is_var<T_a>::value) {
              double base = log_j + log_p[j] + log_r[k - j];
              summands_da.push_back(base + dlog_r_da[k - j]);
              summands_da.push_back(base + std::exp(log_dlog_p_da[j]));
            }
          }
          double log_kp1 = std::log((double)k + 1);
          log_p.push_back(-log_kp1 + stan::math::log_sum_exp(summands));
          if(stan::is_var<T_a>::value) {
            log_dlog_p_da.push_back(-log_kp1 - log_p.back() + stan::math::log_sum_exp(summands_da));
          }
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

    double lpmf = 0;
    double dlog_p_da = 0;
    std::for_each(x.begin(), x.end(), [&](int x) {
      lpmf += log_p[x];
      if(stan::is_var<T_a>::value) {
        dlog_p_da += std::exp(log_dlog_p_da[x]);
      }
    });

    return {lpmf, dlog_p_da };
  }


template <bool propto = true, typename T_a, typename T_b, typename T_c>
typename boost::math::tools::promote_args<T_a, T_b, T_c>::type
  poisson_tweedie_raw_lpmf(const std::vector<int> &x, const T_a &a, const T_b &b, const T_c &c, std::ostream* &pstream__);

template <bool propto = true>
double
  poisson_tweedie_raw_lpmf(const std::vector<int> &x, const double &a, const double &b, const double &c, std::ostream* &pstream__) {

    return poisson_tweedie_raw_lpmf_internal(x, a, b, c, pstream__).lpmf;
  }


double poisson_tweedie_raw_da(const std::vector<int> &x, const double &a, const double &b, const double &c, std::ostream* &pstream__) {
  stan::math::var _a(a);
  return poisson_tweedie_raw_lpmf_internal(x, _a, b, c, pstream__).dlpmf_da;
}
