#include <stan/math/rev/scal/meta/operands_and_partials.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>

#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

#include <stan/math/prim/scal/err/check_greater.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>

#include <stan/math/prim/arr/fun/log_sum_exp.hpp>

#include <boost/math/tools/numerical_differentiation.hpp>


#include<algorithm>

namespace stan { namespace math {


namespace local_override {

template <typename T>
inline T my_log_sum_exp(const std::vector<T> &x) {
  return stan::math::log_sum_exp(x);
}


template<>
inline std::complex<double> my_log_sum_exp<std::complex<double>>(const std::vector<std::complex<double>>& x) {
  using std::exp;
  using std::log;
  using std::numeric_limits;
  double r_max = -numeric_limits<double>::infinity();
  for (std::complex<double> xx : x) {
    double xx_abs = std::abs(xx);
    if (xx_abs > r_max) {
      r_max = xx_abs;
    }
  }

  std::complex<double> sum = 0.0;
  for (size_t ii = 0; ii < x.size(); ii++) {
    if (std::abs(x[ii]) != -numeric_limits<double>::infinity()) {
      sum += exp(x[ii] - r_max);
    }
  }

  return r_max + log(sum);
}


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
      log_p0 = b * local_override::my_log1m(c);
    } else if (std::abs(a) < 1e-4) {
      //Using Taylor expansion of order 2
      T_c log_1mc = local_override::my_log1m(c);
      log_p0 = b * (
        log_1mc +
        0.5 * a * boost::math::pow<2>(log_1mc) +
        (boost::math::pow<2>(a) * boost::math::pow<3>(log_1mc))/ 6.0
      );
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
          //TODO - first term can be log(0) if a <= 0
          T_ac new_log_r = log(k - 1.0 + a) - log((double)k + 1) + log_c + log_r.back();
          log_r.push_back(new_log_r);

          summands.push_back(log_b + log_c + log_p.back());
          for(int j = 1; j<= k; ++j) {
            double log_j = std::log((double)j);
            summands.push_back(log_j + log_r[k - j] + log_p[j]);
          }
          double log_kp1 = std::log((double)k + 1);
          log_p.push_back(-log_kp1 + local_override::my_log_sum_exp(summands));
          summands.clear();
        }
      }
    }
/*

    if(pstream__ != 0) {
       auto print_func = [&](T_ret x) { *pstream__ << x << " "; };
       std::for_each(log_p.begin(), log_p.end(), print_func);
       *pstream__ << std::endl;
       std::for_each(log_r.begin(), log_r.end(), print_func);
       *pstream__ << std::endl;
    }
*/
    T_ret lpmf = 0;
    for(int xx : x) {
      lpmf += log_p[xx];
    }

    return lpmf;
  }



}} //namespace

