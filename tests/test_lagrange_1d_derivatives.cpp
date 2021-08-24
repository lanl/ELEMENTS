#include "lagrange_polynomials.h"

#include <iostream>
#include <iomanip>

// Test function
NumType f(NumType x) {
  return sin(x);
};

int main() {
  std::cout.precision(15);

  SizeType N = 40, Ne = 20;
  NumType zl = -10.0, zr = 10.0;

  // Construct a Lagrange interpolant of a function using Chebyshev points
  NumType z[N];
  lagrange::chebyshev_points(N, zl, zr, z);

  NumType F[N];
  for (int i = 0; i < N; i++)
    F[i] = f(z[i]);

  NumType w[N];
  lagrange::compute_barycentric_weights(N, z, w);

  // Create an array of equispaced points
  NumType xe[Ne];
  lagrange::equispaced_points(Ne, zl, zr, xe);

  // Work array for intermediate coefficients
  NumType G[N];

  // At each equispaced point...
  std::cout << std::setw(20) << "x"
            << std::setw(20) << "p(x)"
            << std::setw(20) << "dp(x)"
            << std::setw(20) << "dp(x), CS"
            << std::setw(20) << "dp(x) error"
            << std::setw(20) << "d2p(x)"
            << std::setw(20) << "d2p(x), CS"
            << std::setw(20) << "d2p(x) error"
            << std::endl;
  std::cout << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::setw(20) << "------------"
            << std::endl;
  for (int i = 0; i < Ne; i++) {
    // check if the point is coincident with a vertex
    SizeType ic;
    lagrange::find_coincident_vertex(N, z, xe[i], ic);

    // apply a small perturbation to its imaginary part
    RealNumber h = 1e-30;
    xe[i] = ComplexNumber(common::real(xe[i]), h);

    // evaluate the interpolant
    NumType p;
    lagrange::evaluate_1d(0, ic, N, z, w, xe[i], F, G, &p);

    // evaluate the first derivative of the interpolant
    NumType dp;
    lagrange::evaluate_1d(1, ic, N, z, w, xe[i], F, G, &dp);

    // approximate the first derivative (using the complex step approximation)
    RealNumber Dp = common::imag(p)/h; 

    // measure absolute error in first derivative
    RealNumber E1 = abs(common::real(dp) - Dp);

    // evaluate the second derivative of the interpolant
    NumType d2p;
    lagrange::evaluate_1d(2, ic, N, z, w, xe[i], F, G, &d2p);

    // approximate the second derivative
    RealNumber D2p = common::imag(dp)/h;

    // measure absolute error in second derivative
    RealNumber E2 = abs(common::real(d2p) - D2p);

    // print evaluations and errors
    std::cout << std::setw(20) << common::real(xe[i])
              << std::setw(20) << common::real(p)
              << std::setw(20) << common::real(dp)
              << std::setw(20) << common::real(Dp)
              << std::setw(20) << E1
              << std::setw(20) << common::real(d2p)
              << std::setw(20) << common::real(D2p)
              << std::setw(20) << E2
              << std::endl;
  }

  return 0;
}
