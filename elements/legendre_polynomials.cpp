#include "legendre_polynomials.h"

namespace legendre {
  template <typename NumType>
  NumType evaluate(int n, NumType ksi) {
    if (n == -1) return 0.0;
    if (n == 0) return 1.0;
    return (1.0/double(n))*((2.0*double(n) - 1.0)*ksi
        *evaluate(n - 1, ksi)
        - (double(n) - 1.0)*evaluate(n - 2, ksi));
  };

  template <typename NumType>
  NumType evaluate_derivative(int n, NumType ksi) {
    return jacobi::evaluate_derivative(n, NumType(0.0), NumType(0.0), ksi);
  };

  // Explicit instatiations of template functions
  template Real evaluate(int n,  Real ksi);
  template Complex evaluate(int n,  Complex ksi);

  template Real evaluate_derivative(int n, Real ksi);
  template Complex evaluate_derivative(int n, Complex ksi);
}
