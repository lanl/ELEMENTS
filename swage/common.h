#ifndef COMMON_H
#define COMMON_H

#include <cmath>
#include <cstdio>
#include <complex>
#include <limits>
#include <algorithm>

#define NUM_EPS std::numeric_limits<RealNumber>::epsilon()
#define NUM_MIN std::numeric_limits<RealNumber>::min()
#define NUM_MAX std::numeric_limits<RealNumber>::max()

typedef double RealNumber;
typedef std::complex<double> ComplexNumber;
typedef size_t SizeType;
#ifdef USE_COMPLEX_NUMBERS
typedef ComplexNumber NumType;
#else
typedef RealNumber NumType;
#endif

namespace common {
  // Definitions used to ensure compatibility when switching from real number
  // type to complex number type to test derivative implementations via the
  // complex step method
  inline RealNumber real(RealNumber number) { return number; }
  inline RealNumber real(ComplexNumber number) { return number.real(); }

  inline RealNumber imag(RealNumber number) { return 0.0; }
  inline RealNumber imag(ComplexNumber number) { return number.imag(); }

  inline RealNumber abs(RealNumber number) { return std::abs(number); }
  inline RealNumber abs(ComplexNumber number) { 
      return std::abs(number.real()); }
}
#endif // COMMON_H
