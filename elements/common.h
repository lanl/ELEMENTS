#pragma once
#include <cmath>
#include <cstdio>
#include <complex>
#include <limits>
#include <algorithm>

typedef double RealScalar;
typedef std::complex<double> ComplexScalar;
typedef size_t SizeType;

#define SCALAR_EPSILON std::numeric_limits<RealScalar>::epsilon()
#define SCALAR_MIN std::numeric_limits<RealScalar>::min()
#define SCALAR_MAX std::numeric_limits<RealScalar>::max()

#ifdef USE_COMPLEX_SCALAR
typedef ComplexScalar Scalar;
#else
typedef RealScalar Scalar;
#endif

namespace common {
  inline RealScalar real(RealScalar scalar) { return scalar; }
  inline RealScalar real(ComplexScalar scalar) { return scalar.real(); }

  inline RealScalar imag(RealScalar scalar) { return 0.0; }
  inline RealScalar imag(ComplexScalar scalar) { return scalar.imag(); }

  inline RealScalar abs(RealScalar scalar) { return std::abs(scalar); }
  inline RealScalar abs(ComplexScalar scalar) { 
      return std::abs(scalar.real()); }
}
