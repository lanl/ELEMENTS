#include "TensorProductBasis.h"

#include "error.h"


void TensorProductBasis::eval_jac(const NumType *, const NumType *, 
    const NumType *, const NumType *, NumType *) {
  std::string err_msg(
      "Error: Jacobian routines not implemented for this basis type.")
  throw NotImplementedError(err_msg);
};

NumType TensorProductBasis::eval_det_jac(const NumType *, const NumType *, 
    const NumType *, const NumType *) {
  std::string err_msg(
      "Error: Jacobian routines not implemented for this basis type.")
  throw NotImplementedError(err_msg);
  return 0.0;
}

void TensorProductBasis::eval_inv_jac(const NumType *, const NumType *, 
    const NumType *, const NumType *, NumType *) {
  std::string err_msg(
      "Error: Jacobian routines not implemented for this basis type.")
  throw NotImplementedError(err_msg);
}

// Explicit instantiation of template class
template class TensorProductBasis<Real>;
template class LegendreProductBasis<Complex>;
