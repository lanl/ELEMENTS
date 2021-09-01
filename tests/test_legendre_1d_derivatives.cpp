#include "legendre_polynomials.h"
#include "jacobi_polynomials.h"

#include <iomanip>

using namespace std;

/** Test Jacobi polynomial implementation */
int main() {
  cout.precision(15);

  // Test parameters (alpha = beta = 0 means Legendre)
  const size_t 
  num_coords = 1001, 
  degree     = 7;

  const NumType
  alpha = -0.5, 
  beta  = -0.5;

  const RealNumber h = 1e-30;

  // Define grid
  NumType ksi = -1.0 + NumType(0.0, h);
  NumType delta_ksi = 2.0/(num_coords - 1);

  // Define Jacobi polynomial and its derivative
  cout << "Jacobi polynomial" << "\n"
       << setw(10) << "alpha:"
       << setw(25) << common::real(alpha) << "\n"
       << setw(10) << "beta:"
       << setw(25) << common::real(beta) << "\n"
       << setw(10) << "degree:" 
       << setw(25) << degree << endl;
  auto P = [&degree, &alpha, &beta](NumType ksi) {
    return legendre::evaluate(degree, ksi);
    //return jacobi::evaluate(degree, alpha, beta, ksi);
  };

  auto dP = [&degree, &alpha, &beta](NumType ksi) {
    return legendre::evaluate_derivative(degree, ksi);
    //return jacobi::evaluate_derivative(degree, alpha, beta, ksi);
  };
  
  // Evaluate Jacobi polynomial and its derivative over grid
  cout << setw(25) << "ksi"
       << setw(25) << "eval"
       << setw(25) << "deriv"
       << setw(25) << "deriv error"
       << endl;

  for (int i = 0; i < num_coords; i++) {
    cout << setw(25) << common::real(ksi)
         << setw(25) << common::real(P(ksi))
         << setw(25) << common::real(dP(ksi))
         << setw(25) << fabs(common::real(dP(ksi)) - common::imag(P(ksi))/h)
       << endl;
    ksi += delta_ksi;
  }

  return 0;
}
