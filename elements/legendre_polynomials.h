#pragma once

#include "jacobi_polynomials.h"
#include "common.h"

namespace legendre {
  /* Polynomial and derivative evaluation */
  template <typename NumType>
  NumType evaluate(int n, NumType ksi);

  template <typename NumType>
  NumType evaluate_derivative(int n, NumType ksi);
}
