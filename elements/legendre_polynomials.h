#pragma once

#include "jacobi_polynomials.h"
#include "common.h"

namespace legendre {
  /* Polynomial and derivative evaluation */
  NumType evaluate(int n, NumType ksi);
  NumType evaluate_derivative(int n, NumType ksi);
}
