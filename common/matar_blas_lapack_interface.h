#pragma once

#include "c_blas_lapack_interface.h"
#include "error.h"

#include "matar.h"

ErrorCode transpose_2d_matar_carray(const CArray<NumType> &, CArray<NumType> &);

ErrorCode matvec_matar_carray_blas(const CArray<NumType> &, 
    const CArray<NumType> &, CArray<NumType> &);

ErrorCode multiply_2d_matar_carray_blas(const CArray<NumType> &, 
    const CArray<NumType> &, CArray<NumType> &);

ErrorCode invert_2d_matar_carray_lapack(const CArray<NumType> &, 
    CArray<NumType> &);
