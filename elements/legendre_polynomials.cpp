#include "legendre_polynomials.h"

namespace legendre {
  NumType evaluate(int n, NumType ksi) {
    if (n == -1) return 0.0;
    if (n == 0) return 1.0;
    return (1.0/double(n))*((2.0*double(n) - 1.0)*ksi
        *evaluate(n - 1, ksi)
        - (double(n) - 1.0)*evaluate(n - 2, ksi));
  };

  NumType evaluate_derivative(int n, NumType ksi) {
    return jacobi::evaluate_derivative(n, 0.0, 0.0, ksi);
  };
}
