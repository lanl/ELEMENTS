#pragma once

#include "jacobi_polynomials.h"
#include "blas_lapack_c_interface.h"
#include "error.h"

/* 
 * For a quadrature rule based on the interpolation of the integrand by a
 * n-degree Jacobi polynomial, compute the quadrature points and weights using
 * the algorithm presented in "Calculation of Guass Quadrature Rules" (Golub
 * and Welsh, 1969)
 */
ErrorCode compute_gauss_jacobi_quadrature_rule(
    size_t n, NumType alpha, NumType beta, 
    CArray<NumType> &points, CArray<NumType> &weights) {
  if (n > 1) {
    // Populate arrays corresponding to the diagonal and subdiagonal (or
    // superdiagonal) of the symmetric tridiagonal matrix
    RealNumber diag[n], subdiag[n-1];

    for (int k = 0; k < n; k++) {
      RealNumber 
      a_k = bojador::real(jacobi::a(alpha, beta, k+1)),
      b_k = bojador::real(jacobi::b(alpha, beta, k+1));

      diag[k] = -b_k/a_k;
    }

    for (int k = 0; k < n-1; k++) {
      RealNumber 
      a_k   = bojador::real(jacobi::a(alpha, beta, k+1)),
      a_kp1 = bojador::real(jacobi::a(alpha, beta, k+2)),
      c_kp1 = bojador::real(jacobi::c(alpha, beta, k+2));

      subdiag[k] = sqrt(c_kp1/(a_k*a_kp1));
    }

    // Use dstev from LAPACK to compute eigenvectors and eigenvalues of the
    // tridiagonal matrix
    const char compute_eigenvectors = 'V';
    int N = n;
    RealNumber eigenvectors[n*n], work_array[2*N-2];
    int info = 0;
    lapack_stev(
        &compute_eigenvectors, 
        &N, 
        diag, subdiag, 
        eigenvectors, &N, 
        work_array, 
        &info);

    // TODO Check info for error

    // Compute zeroth moment (integral of weight function times 1); this integral
    // can be calculated analytically using the beta function (discovered this in
    // the SciPy source code); I implemented the beta function here in terms of
    // the gamma function
    RealNumber mu0 = pow(2.0, bojador::real(alpha + beta + 1.0))
        *std::tgamma(bojador::real(alpha + 1.0))
        *std::tgamma(bojador::real(beta + 1.0))
        /std::tgamma(bojador::real(alpha + beta + 2.0));

    // Extract quadrature points from the eigenvalues and weights from the
    // magnitude of the eigenvector
    for (int j = 0; j < n; j++) {
      points(j) = diag[j]; 

      RealNumber q_0j = eigenvectors[0+n*j];
      weights(j) = q_0j*q_0j*mu0;

    }
  } else if (n == 1) {
    points(0) = 0;
    weights(0) = 2;
  } else if (n == 0) {
    return QUADRATURE_ERROR;
  } 

  return 0;
}
