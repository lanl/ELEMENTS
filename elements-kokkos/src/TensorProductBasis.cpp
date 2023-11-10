#include "TensorProductBasis.h"

#include "error.h"


template <typename NumType>
TensorProductBasis<NumType>::TensorProductBasis(const SizeType order) : Np(order) {}

template <typename NumType>
void TensorProductBasis<NumType>::eval_jac(const NumType *, const NumType *, 
    const NumType *, const NumType *, NumType *) {
  std::string err_msg(
      "Error: Jacobian routines not implemented for this basis type.");
  throw NotImplementedError(err_msg);
};

template <typename NumType>
NumType TensorProductBasis<NumType>::eval_det_jac(const NumType *, 
    const NumType *, const NumType *, const NumType *) {
  std::string err_msg(
      "Error: Jacobian routines not implemented for this basis type.");
  throw NotImplementedError(err_msg);
  return 0.0;
}

template <typename NumType>
void TensorProductBasis<NumType>::eval_inv_jac(const NumType *, 
    const NumType *, const NumType *, const NumType *, NumType *) {
  std::string err_msg(
      "Error: Jacobian routines not implemented for this basis type.");
  throw NotImplementedError(err_msg);
}

// Explicit instantiation of template class
template class TensorProductBasis<Real>;
template class TensorProductBasis<Complex>;
