#pragma once

#include "common.h"

namespace lagrange {
  // Comparing floating point numbers
  inline bool almost_equal(NumType a, NumType b);

  // Possible optimization
  inline SizeType branchless_choice(bool c, SizeType a, SizeType b);

  // Point distributions
  void equispaced_points(SizeType N, NumType &zl, NumType &zr, NumType *z);
  void chebyshev_points(SizeType N, NumType &zl, NumType &zr, NumType *z);

  // Identifying singularities
  void find_coincident_vertex(const SizeType &N, const NumType *z, 
      const NumType &x, SizeType &k);

  // Barycentric interpolation
  void compute_barycentric_weights(const SizeType &N, const NumType *z, 
      NumType *w);
  void evaluate_1d(const SizeType &ND, const SizeType &ic, const SizeType &Nv, 
      const NumType *z, const NumType *w, const NumType &x, const NumType *c0, 
      NumType *ci, NumType *co);
  void evaluate_3d(const SizeType &pde, const SizeType &N, const NumType *z, 
      const NumType *w, const NumType *x, const NumType *f, NumType *g);

  // Geometric quantities from interpolants
  void compute_jacobian_determinant(SizeType &N, NumType *z, NumType *w, 
      NumType *cx, NumType *cy, NumType *cz, NumType *g, NumType X[3], 
      NumType &J);
}
