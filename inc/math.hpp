#ifndef H_MATH
#define H_MATH

#include <cmath>

namespace Math{
  const double pi(4.0*atan(1.0));

  double GetDeltaPhi(const double, const double);
  double GetAbsDeltaPhi(const double, const double);
  double GetDeltaR(const double, const double, const double, const double);
}

#endif
