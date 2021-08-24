#ifndef LAGRANGE_POLYNOMIALS_H
#define LAGRANGE_POLYNOMIALS_H

#include "common.h"

namespace lagrange {
  /*** Utilities ****/

  // Possible optimization
  inline SizeType branchless_choice(bool c, SizeType a, SizeType b);

  // Identifying singularities
  inline bool almost_equal(NumType a, NumType b);
  void find_coincident_vertex(const SizeType &N, const NumType *z, 
      const NumType &x, SizeType &k);

  // Converting between representations of array indices
  void base_10_to_mixed_radix(const SizeType &Nb, const SizeType *b, 
      SizeType x, SizeType *y);
  void mixed_radix_to_base_10(const SizeType &Nb, const SizeType *b, 
      SizeType *x, SizeType &y);

  // Point distributions
  void equispaced_points(SizeType N, NumType &zl, NumType &zr, NumType *z);
  void chebyshev_points(SizeType N, NumType &zl, NumType &zr, NumType *z);

  // Encoding and decoding a requested partial derivative to and from an
  // unsigned integer
  void encode_partial_derivative(const SizeType &nx, const SizeType &ny, 
      const SizeType &nz, SizeType &e);
  void decode_partial_derivative(SizeType e, SizeType &nx, SizeType &ny, 
      SizeType &nz);

  
  /*** Barycentric interpolation ***/
  void compute_barycentric_weights(const SizeType &N, const NumType *z, 
      NumType *w);
  void evaluate_1d(const SizeType &ND, const SizeType &ic, const SizeType &Nv, 
      const NumType *z, const NumType *w, const NumType &x, const NumType *c0, 
      NumType *ci, NumType *co);
  void evaluate_3d(const SizeType &pde, const SizeType &N, const NumType *z, 
      const NumType *w, const NumType *x, const NumType *f, NumType *g);


  /*** Geometric quantities from interpolants ***/
  void compute_jacobian_determinant(SizeType &N, NumType *z, NumType *w, 
      NumType *cx, NumType *cy, NumType *cz, NumType *g, NumType X[3], 
      NumType &J);
}
#endif // LAGRANGE_POLYNOMIALS_H
