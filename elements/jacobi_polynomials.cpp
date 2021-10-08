#include "jacobi_polynomials.h"

namespace jacobi {
  template <typename NumType>
  NumType eval(int n, NumType alpha, NumType beta, NumType X) {
    if (n == -1) return 0.0;
    if (n == 0) return 1.0;
    return (a(alpha, beta, n)*X + b(alpha, beta, n))
        *eval(n - 1, alpha, beta, X) 
        - c(alpha, beta, n)*eval(n - 2, alpha, beta, X);
  }

  template <typename NumType>
  NumType eval_der(int n, NumType alpha, NumType beta, NumType X) {
    if (n == 0) return 0.0;
    return 0.5*(double(n) + alpha + beta + 1.0)
        *eval(n-1, alpha + 1.0, beta + 1.0, X);
  }

  // Explicit instantiations of template functions
  template Real eval(int n, Real alpha, Real beta, Real X);
  template Complex eval(int n, Complex alpha, Complex beta, Complex X);

  template Real eval_der(int n, Real alpha, Real beta, Real X);
  template Complex eval_der(int n, Complex alpha, Complex beta, Complex X);
}
