#pragma once

#include "common.h"

namespace jacobi {
  /* Recurrence relation parameters */
  template <typename NumType>
  inline NumType a(NumType alpha, NumType beta, int n) {
    if (n == 1) return 0.5*(alpha + beta) + 1.0;
    return (2.0*double(n) + alpha + beta - 1.0)*(2.0*double(n) + alpha + beta)
        /(2.0*double(n)*(double(n) + alpha + beta));
  };

  template <typename NumType>
  inline NumType b(NumType alpha, NumType beta, int n) {
    if (n == 1) return 0.5*(alpha - beta); 
    return (alpha*alpha - beta*beta)*(2*double(n) + alpha + beta - 1.0)
        /(2.0*double(n)*(double(n) + alpha + beta)
            *(2.0*double(n) + alpha + beta - 2.0));
  };

  template <typename NumType>
  inline NumType c(NumType alpha, NumType beta, int n) {
    if (n == 1) return 0.0;
    return (double(n) + alpha - 1.0)*(double(n) + beta - 1.0)
       *(2.0*double(n) + alpha + beta)
       /(double(n)*(double(n) + alpha + beta)
           *(2.0*double(n) + alpha + beta - 2.0));
  };

  /* Polynomial and polynomial derivative evaluation */
  template <typename NumType>
  NumType eval(int n, NumType alpha, NumType beta, NumType X);

  template <typename NumType>
  NumType eval_der(int n, NumType alpha, NumType beta, NumType X);
}
