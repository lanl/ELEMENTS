#include "lagrange_polynomials.h"

#include <iostream>
#include <iomanip>

/* 
 * Test derivatives of 3D Lagrange interpolation by computing derivatives of
 * the spatial coordinates with respect to the reference coordinates
 * (ultimately what is needed to compute the Jacobian)
 */
int main() {
  std::cout.precision(15);

  const SizeType Np = 3;  // polynomial order

  // Define vertex reference coordinates (only in 1D, assume equispaced)
  const SizeType Nv = Np + 1;
  NumType zl = -1.0, zr = 1.0;
  NumType z[Nv];
  lagrange::equispaced_points(Nv, zl, zr, z);

  // Evaluate the barycentric node weights
  NumType w[Nv];
  lagrange::compute_barycentric_weights(Nv, z, w);

  // Define vertex spatial coordinates equal to reference coordinates...
  NumType cx[Nv*Nv*Nv];
  NumType fy[Nv*Nv*Nv];
  NumType fz[Nv*Nv*Nv];

  for (SizeType k = 0; k < Nv; k++) {
    for (SizeType j = 0; j < Nv; j++) {
      for (SizeType i = 0; i < Nv; i++) {
        NumType x[3];
        cx[i + j*Nv + k*Nv*Nv] = z[i];
        fy[i + j*Nv + k*Nv*Nv] = z[j];
        fz[i + j*Nv + k*Nv*Nv] = z[k];
      }
    }
  }

  // ...except for a perturbation to the first vertex
  cx[0] = -2.0;
  fy[0] = -2.0;
  fz[0] = -2.0;


  // Create equispaced "quadrature" points
  const SizeType Nq = 3;  // 1D quadrature points

  NumType ksi_c[Nq];
  NumType eta_c[Nq];
  NumType zet_c[Nq];
  lagrange::equispaced_points(Nq, zl, zr, ksi_c);
  lagrange::equispaced_points(Nq, zl, zr, eta_c);
  lagrange::equispaced_points(Nq, zl, zr, zet_c);

  NumType co[3*Nv+1]; // intermediate coefficients array

  // At each quadrature point...
  std::cout << std::setw(25) << "ksi"
            << std::setw(25) << "x"
            << std::setw(25) << "dx/dksi"
            << std::setw(25) << "dx/dksi, CS"
            << std::setw(25) << "dx/dksi, error"
            << std::endl;
  std::cout << std::setw(25) << "---------------"
            << std::setw(25) << "---------------"
            << std::setw(25) << "---------------"
            << std::setw(25) << "---------------"
            << std::setw(25) << "---------------"
            << std::endl;
  for (SizeType k = 0; k < Nq; k++) {
    for (SizeType j = 0; j < Nq; j++) {
      for (SizeType i = 0; i < Nq; i++) {
          // perturb the imaginary part of the first component of ksi
          NumType ksi[3];

          ksi[0] = ksi_c[i];
          ksi[1] = eta_c[j];
          ksi[2] = zet_c[k];

          RealNumber h = 1e-30;
          ksi[0] = ComplexNumber(common::real(ksi[0]), h);

          // evaluate x, y, z for the given ksi
          lagrange::evaluate_3d(0, Nv, z, w, ksi, cx, co);
          NumType x = co[3*Nv];

          // evaluate the derivative of x with respect to ksi
          lagrange::evaluate_3d(1, Nv, z, w, ksi, cx, co);
          NumType dxdksi = co[3*Nv];

          // approximate the derivative (using the complex step approximation)
          RealNumber DxDksi = common::imag(x)/h; 

          // measure absolute error in first derivative
          RealNumber E11 = std::fabs(common::real(dxdksi) - DxDksi);

          // print evaluations and errors
          std::cout << std::setw(25) << common::real(ksi[0])
                    << std::setw(25) << common::real(x)
                    << std::setw(25) << common::real(dxdksi)
                    << std::setw(25) << DxDksi
                    << std::setw(25) << E11
                    << std::endl;
      }
    }
  }
  
  return 0;
}
