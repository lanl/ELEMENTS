#include "jacobi_polynomials.h"

namespace jacobi {
  template <typename NumType>
  NumType evaluate(int n, NumType alpha, NumType beta, NumType ksi) {
    if (n == -1) return 0.0;
    if (n == 0) return 1.0;
    return (a(alpha, beta, n)*ksi + b(alpha, beta, n))
        *evaluate(n - 1, alpha, beta, ksi) 
        - c(alpha, beta, n)*evaluate(n - 2, alpha, beta, ksi);
  }

  template <typename NumType>
  NumType evaluate_derivative(int n, NumType alpha, NumType beta, NumType ksi) {
    if (n == 0) return 0.0;
    return 0.5*(double(n) + alpha + beta + 1.0)
        *evaluate(n-1, alpha + 1.0, beta + 1.0, ksi);
  }

  // Explicit instantiations of template functions
  template Real evaluate(int n, Real alpha, Real beta, Real ksi);
  template Complex evaluate(int n, Complex alpha, Complex beta, Complex ksi);

  template Real evaluate_derivative(int n, Real alpha, Real beta, Real ksi);
  template Complex evaluate_derivative(int n, Complex alpha, Complex beta, 
      Complex ksi);
}
