#pragma once

#include "TensorProductBasis.h"
#include "common.h"


template <typename NumType>
struct BernsteinProductBasis : public TensorProductBasis { 
  SizeType N;
  SizeType Ne;
  
  // Work array for intermediate coefficients
  NumType *C;

  BernsteinProductBasis(const SizeType);
  ~BernsteinProductBasis();
  
  // Basis functions and basis function gradients
  NumType eval_basis(const SizeType, const NumType *);
  void eval_grad_basis(const SizeType, const NumType *, NumType *);

  // Function approximation over element
  NumType eval_approx(const NumType *, const NumType *); 
  void eval_grad_approx(const NumType *, const NumType *, NumType *);

  // Jacobian of spatial mapping
  void eval_jac(const NumType *, const NumType *, const NumType *,
                const NumType *, NumType *);
  NumType eval_det_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *);
  void eval_inv_jac(const NumType *, const NumType *, const NumType *,
       const NumType *, NumType *);
};
