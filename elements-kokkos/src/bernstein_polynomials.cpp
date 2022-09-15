#include "bernstein_polynomials.h"


namespace bernstein {
  template <typename NumType> 
  NumType eval( const int n, const int v, const NumType X) {
    if ( n == 0 && v != 0 ) return 0;
    if ( n == 0 && v == 0 ) return 1;
    if ( n < v ) return 0;
    return 0.5*((1.0 - X)*eval(n-1, v, X) + (1.0 + X)*eval(n-1, v-1, X)); 
  };
  template <typename NumType>
  NumType eval_der( const int n, const int v, const NumType X) {
    if ( n == 0) return 0;
    if ( n < v) return 0;
    if ( n == 1 && v == 0) return -0.5;
    if ( n == 1 && v == 1) return 0.5;
    return 0.5*((1.0 - X)*eval_der(n-1,v, X) + (1.0 + X)*eval_der(n-1, v-1, X) 
        + eval(n-1, v-1, X) - eval(n-1, v, X));
  };

  template <typename NumType>
  NumType eval_approx(const SizeType N, const NumType *c, const NumType X) {
    NumType sum = 0.0;
    for (SizeType j = 0; j < N; j++) {
      sum += c[j]*eval(N, j, X);
    }
    return sum;
  }

  template <typename NumType>
  NumType eval_der_approx(const SizeType N, const NumType *c, const NumType X) { 
    NumType sum = 0.0;
    for (SizeType j = 0; j < N; j ++) { 
      sum += c[j]*eval_der( N, j, X);
    }
    return sum;
  };


  // Explicit instantiations of template functions
  template Real eval(const int n, const int v, const Real X);
  template Complex eval(const int n, const int v, const Complex X);
  
  template Real eval_der(const int n, const int v, const Real X);
  template Complex eval_der(const int n, const int v, const Complex X);
  
  template Real eval_approx(const SizeType N, const Real *c, const Real X);
  template Complex eval_approx(const SizeType N, const Complex *c, 
      const Complex X);
  
  template Real eval_der_approx(const SizeType N, const Real *c, const Real X);
  template Complex eval_der_approx(const SizeType N, const Complex *c, 
      const Complex X);
}
