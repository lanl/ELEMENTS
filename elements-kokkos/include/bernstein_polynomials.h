#pragma once

#include "common.h"

// Bernstein polynomial of order n:
// B(n,v) = (0.5)^n * (C^n_v) * (1-x)^(n-v) * (1+x)^v
// for x \in [-1, 1]

namespace bernstein{
  // Bernstein polynomials
  // n is the order, v is the "index", X \in [-1, 1]
  template <typename NumType> NumType eval(const int n, const int v, const NumType X);
    
  template <typename NumType> NumType eval_der(const int n, const int v, const NumType X);

  //Bernstein Approximations
  template <typename NumType>
  NumType eval_approx(const SizeType N, const NumType *c, const NumType X);

  template <typename NumType>
  NumType eval_der_approx(const SizeType N, const NumType *c, const NumType X);
}
