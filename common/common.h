#pragma once

#include "matar.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

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

  // Converting between representations of array indices
  void base_10_to_mixed_radix(const SizeType &Nb, const SizeType *b, 
      SizeType x, SizeType *y);
  SizeType mixed_radix_to_base_10(const SizeType &Nb, const SizeType *b, 
      SizeType *x);

  // Encoding and decoding a requested partial derivative to and from an
  // unsigned integer
  SizeType encode_partial_derivative(const SizeType &nx, const SizeType &ny, 
      const SizeType &nz);
  void decode_partial_derivative(SizeType e, SizeType &nx, SizeType &ny, 
      SizeType &nz);
}
