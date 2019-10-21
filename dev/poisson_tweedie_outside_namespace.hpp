#include <stan/math/prim/scal/fun/constants.hpp>

#include<algorithm>

namespace stan { namespace math {
inline std::complex<double> log_sum_exp(const std::vector<std::complex<double>>& x) {
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


// template <typename T1, typename T2>
// inline return_type_t<std::complex<T1>, std::complex<T2>> log_diff_exp(const std::complex<T1> &x, const std::complex<T2> &y) {
//   T1 abs_x = std::abs(x);
//   T2 abs_y = std::abs(y);
//   if (std::abs(x - y) < 0) {
//     return (std::abs(x) < INFTY && x == y) ? NEGATIVE_INFTY : NOT_A_NUMBER;
//   }
//   return x + log1m_exp(y - x);
// }

}} //namespace

#include <boost/math/tools/numerical_differentiation.hpp>
