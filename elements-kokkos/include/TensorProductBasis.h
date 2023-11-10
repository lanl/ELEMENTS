#pragma once


#include "common.h"


/*
 * Interface class for 3D tensor product basis, to be called for evaluating: 
 *  - basis functions and their gradients;
 *  - approximations in the basis and their derivatives; and,
 *  - for choices of basis that are typically used to represent the geometry of
 *    the mesh (Lagrange and Bernstein), the Jacobian of the spatial mapping.
 */
template <typename NumType>
struct TensorProductBasis {
  static const SizeType Nd = 3;
  SizeType Np;

  // Arrays for converting from flat to multidimensional indices
  SizeType ijk[Nd];
  SizeType rad[Nd];

  TensorProductBasis(const SizeType);
  ~TensorProductBasis() {};

  // Basis functions and basis function gradients
  virtual NumType eval_basis(const SizeType, const NumType *) = 0;
  virtual void eval_grad_basis(const SizeType, const NumType *, NumType *) = 0;

  // Function approximation over element
  virtual NumType eval_approx(const NumType *, const NumType *) = 0;
  virtual void eval_grad_approx(const NumType *, const NumType *, NumType *) = 0;

  // Jacobian of spatial mapping
  virtual void eval_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *, NumType *);
  virtual NumType eval_det_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *);
  virtual void eval_inv_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *, NumType *);
};
