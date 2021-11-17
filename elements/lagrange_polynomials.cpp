#include "lagrange_polynomials.h"

namespace lagrange {
  /*
   * Assess equality of two double precision floating point numbers of the same
   * order of magnitude
   *
   * Paramters
   * ---------
   * a : a double precision number
   * b : another double precision number
   *
   * Returns
   * -------
   *   a boolean 
   */
  inline bool almost_equal(NumType a, NumType b) {
    return common::abs(a - b) < 2.0*NUM_EPS;
  }


  /*
   * Branchless choice between two SizeType variables based on a logical
   * condition
   *
   * Parameters
   * ----------
   * c : condition 
   * a : to be returned if condition is true
   * b : to be returned if condition is false
   *
   * Returns
   * -------
   *   decision between the two variables based on the condition
   */
  inline SizeType branchless_choice(bool c, SizeType a, SizeType b) {
    return a ^ ((b ^ a) & -(!c));
  }

  /*
   * Fill an array with points equally spaced between end points (inclusive)
   *
   * Parameters
   * ----------
   * N  : number of points
   * zl : left end point 
   * zr : right end point 
   *
   * Returns
   * -------
   * z : array of equispaced points
   */
  void equispaced_points(SizeType N, NumType &zl, NumType &zr, 
      NumType *z) {
    for (SizeType i = 0; i < N; i++)
      z[i] = zl + double(i)/double(N - 1)*(zr - zl);
  }

  /*
   * Fill an array with Chebyshev points of the second kind in the interval
   * defined by the specified end points. The symmetry-preserving technique from
   * Chebfun (chebtech1/chebpts.m) is used.
   *
   * Parameters
   * ----------
   * N  : number of points
   * zl : left end point 
   * zr : right end point 
   *
   * Returns
   * -------
   * z : array of Chebyshev points
   */
  void chebyshev_points(SizeType N, NumType &zl, NumType &zr, 
      NumType *z) {
    // Evaluate the points using sine function to preserve symmetry
    NumType f = 0.5*M_PI/double(N);
    for (SizeType i = 0; i < N; i++) {
      int j = (-1*int(N) + 1) + 2*int(i);
      z[i] = sin(f*double(j));
    }

    // Scale the points to fit the domain
    for (int i = 0; i < N; i++)
      z[i] = 0.5*(1.0 - z[i])*zl + 0.5*(1.0 + z[i])*zr;
  }

  /*
   * If the input coordinate is coincident with any of the input vertices,
   * return the index of the vertex with which it coincides
   *
   * It is typically bad practice to check if two floating point numbers, say x
   * and y, are equal by evaluating x - y == 0. But that is precisely what the
   * authors suggest to do in "Barycentric Lagrange Interpolation" (Berrut and
   * Trefethen, 2004). It is also what the authors of Nektar++ do in their
   * implementation of barycentric Lagrange interpolation. I'm not sure what to
   * do
   *
   * Parameters
   * ----------
   * N : number of vertices
   * z : vertex coordinates
   * x : coordinate to test
   * 
   * Returns
   * -------
   * k : index of coincident vertex
   */
  void find_coincident_vertex(const SizeType &N, const NumType *z, 
      const NumType &x, SizeType &k) {
    // The largest number an unsigned integer can be
    k = -1;

    // Adding to -1 in unsigned rolls over zero
    for (SizeType j = 0; j < N; j++) {
      //k += (j + 1)*SizeType(z[j] - x == 0.0);
      k += (j + 1)*SizeType(lagrange::almost_equal(z[j], x));
    }
  }

  /*
   * Calculation of barycentric weights of Lagrange interpolant vertices
   *
   * Parameters
   * ----------
   * N : number of vertices 
   * z : vertex coordinates
   *
   * Returns
   * -------
   * w : barycentric vertex weights
   *
   * Notes
   * -----
   * If one evaluates the interpolant using the barycentric formula of the
   * second kind, the barycentric weights may be scaled by an arbitrary scalar
   * (the maximum weight, for example). Since I am using the barycentric
   * formula of the first kind (also called the modified Lagrange formula), I
   * cannot and do not scale the weights 
   */
  void compute_barycentric_weights(const SizeType &N, const NumType *z, 
      NumType *w) {
    // Compute weights
    for (SizeType j = 0; j < N; j++) {
      w[j] = 1.0;
      for (SizeType k = 0; k < N; k++) {
        // TODO optimization: make branchless? Trefethen uses log o exp trick
        if (j != k) w[j] *= 1.0/(z[j] - z[k]);
      }
    }
  };

  /*
   * Evaluation of (derivative of) 1D Lagrange interpolant 
   *
   * The implementation of the evaluation of the interpolant uses the modified
   * Lagrange formula (barycentric formula of the first kind). The evaluation
   * of derivatives (up to order 2) is based on the formula for the derivative
   * of a rational interpolant given in "Some New Aspects of Rational
   * Interpolation" (Schneider and Werner, 1986).
   *
   * Parameters
   * ----------
   * ND  : order of derivative
   * ic  : index of vertex coincident with coordinate above (-1 if not coincident)
   * Nv  : number of vertices 
   * z   : vertex coordinates
   * w   : barycentric vertex weights
   * x   : coordinate at which to evaluate interpolant
   * c0  : input coefficients (function values at vertices)
   * ci  : work array for intermediate coefficients
   *
   * Returns
   * -------
   * co : output coefficient (evaluation of interpolant or its derivative)
   */
  void evaluate_1d(
      const SizeType &ND, 
      const SizeType &ic, 
      const SizeType &Nv, 
      const NumType *z, 
      const NumType *w, 
      const NumType &x, 
      const NumType *c0, 
      NumType *ci, 
      NumType *co) {
    // Copy input coefficients into the intermediate coefficients
    std::copy(c0, c0+Nv, ci);

    NumType M = 1.0;  // factorial(n)
    NumType L = 1.0;  // nodal polynomial

    // Evaluate the interpolant
    if (ic < Nv) {  // coincident
      *co = ci[ic];
    } else {  // non-coincident
      // Loop over vertices
      *co = 0.0;
      for (SizeType j = 0; j < Nv; j++) {
        // Contribution to interpolant evaluation
        *co += w[j]*ci[j]/(x - z[j]);

        // Contribution to scalar (nodal polynomial evaluation)
        L *= (x - z[j]);
      }

      // Apply scaling
      *co *= L;
    }

    // Evaluate derivatives, building up to specified order
    NumType dnp = *co;  // initialize with n = 0, the evaluation
    if (ic < Nv) {  // coincident
      for (SizeType n = 1; n <= ND; n++) {
        // Zero the output
        *co = 0.0;

        for (SizeType j = 0; j < Nv; j++) {
          // Calculate divided difference and store in intermediate coefficients
          NumType sx = 1.0/(z[ic] - z[j]);
          ci[j] = sx*(dnp - ci[j]);

          // Update output coefficient with contribution from vertex
          // TODO optimization: make branchless? Trefethen uses log o exp trick
          if (j != ic) *co += w[j]*ci[j];
        }

        // Scale the output and copy for use in calculating next order
        M *= n;
        *co *= -M/w[ic];
        dnp = *co;
      }
    } else {  // non-coincident
      for (SizeType n = 1; n <= ND; n++) {
        // Zero the output
        *co = 0.0;

        for (SizeType j = 0; j < Nv; j++) {
          // Calculate divided difference and store in intermediate coefficients
          NumType sx = 1.0/(x - z[j]);
          ci[j] = sx*(dnp - ci[j]);

          // Update output coefficient with contribution from vertex
          *co += w[j]*ci[j]/(x - z[j]);
        }

        // Scale the output and copy for use in calculating next order
        M *= n;
        *co *= L*M;
        dnp = *co;
      }
    }
  };

  /*
   * Evaluation of 3D tensor product Lagrange interpolant and its partial
   * derivatives
   *
   * This implementation uses a dimension-by-dimension approach that reuses the
   * 1D interpolant evaluation kernel. It requires a work array of size 3N to
   * hold intermediate coefficients.
   *
   * Parameters
   * ----------
   * pde : partial derivative encoding
   * N   : number of vertices 
   * z   : vertex coordinates
   * w   : barycentric vertex weights
   * f   : values of function at vertices
   * x   : coordinates at which to evaluate interpolants
   *
   * Returns
   * -------
   * g  : work array of intermediate coefficients
   * p : evaluation of polynomial interpolants
   */
  void evaluate_3d(
      const SizeType &pde,
      const SizeType &N,
      const NumType *z, 
      const NumType *w, 
      const NumType *x, 
      const NumType *f, 
      NumType *g) {  
    // Check the coincidence of the coordinates with the nodes
    SizeType ix, iy, iz;
    find_coincident_vertex(N, z, x[0], ix);
    find_coincident_vertex(N, z, x[1], iy);
    find_coincident_vertex(N, z, x[2], iz);

    // Decode partial derivative information
    SizeType NDx, NDy, NDz;
    common::decode_partial_derivative(pde, NDx, NDy, NDz);

    const NumType *c0;
    NumType *ci, *co;

    // Collapse second dimension into coefficients for third dimension
    for (int k = 0; k < N; k++) {
      // Collapse first dimension into coefficients for second dimension
      for (int j = 0; j < N; j++) {
        c0 = f+j*N+k*N*N;
        ci = g;
        co = g+N+j;
        evaluate_1d(NDx, ix, N, z, w, x[0], c0, ci, co);
      }

      c0 = g+N;
      ci = g;
      co = g+2*N+k;
      evaluate_1d(NDy, iy, N, z, w, x[1], c0, ci, co); 
    }

    // Evaluate 3D interpolant as 1D interpolant
    c0 = g+2*N;
    ci = g;
    co = g+3*N;
    evaluate_1d(NDz, iz, N, z, w, x[2], c0, ci, co);
  };

  /*
   * Evaluate determinant of the Jacobian of the mapping from reference
   * coordinates to spatial coordinates
   *
   * Parameters
   * ----------
   * N : number of vertices 
   * z : vertex coordinates
   * w : barycentric vertex weights
   * x : x spatial coordinate at vertices
   * g : work array of intermediate coefficients
   * X : reference coordinate at which Jacobian is to be evaluated
   *
   * Returns
   * -------
   * J : determinant of Jacobian of mapping
   */
  void compute_jacobian_determinant(
      SizeType &N, 
      NumType *z, 
      NumType *w, 
      NumType *cx, 
      NumType *cy, 
      NumType *cz, 
      NumType *g, 
      NumType X[3], 
      NumType &J) {
    // Evaluate Jacobian (denoted F)
    NumType F[9];
    for (SizeType j = 0; j < 3; j++) {
      SizeType pde = 1 << j;  // partial derivative code (dx, dy, or dz)

      evaluate_3d(pde, N, z, w, X, cx, g);
      F[0+3*j] = g[3*N];
      evaluate_3d(pde, N, z, w, X, cy, g);
      F[1+3*j] = g[3*N];
      evaluate_3d(pde, N, z, w, X, cz, g);
      F[2+3*j] = g[3*N];
    }

    // Calculate Jacobian determinant
    J = F[0]*(F[4]*F[8] - F[5]*F[7]) 
        - F[3]*(F[1]*F[8] - F[2]*F[7]) 
        + F[6]*(F[1]*F[5] - F[2]*F[4]);
  }
}
