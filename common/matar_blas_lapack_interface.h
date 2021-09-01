#pragma once

#include "c_blas_lapack_interface.h"

#include "error.h"
#include "common.h"

// Interface to selected BLAS routines for MATAR CArrays
namespace matar2blas {
  void transpose(const CArray<NumType> &, CArray<NumType> &);

  void matvec(const CArray<NumType> &, 
      const CArray<NumType> &, CArray<NumType> &);

  void matmul(const CArray<NumType> &, 
      const CArray<NumType> &, CArray<NumType> &);
}

// Interface to selected LAPACK routines for MATAR CArrays
namespace matar2lapack {
  void invert(const CArray<NumType> &, CArray<NumType> &);

  void eig_sym_tri(const CArray<NumType> &diag, 
    const CArray<NumType> &subdiag, CArray<NumType> &eigvals, 
    CArray<NumType> &eigvecs);
}
