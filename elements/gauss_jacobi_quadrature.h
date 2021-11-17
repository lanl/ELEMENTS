#pragma once

#include "jacobi_polynomials.h"
#include "matar_blas_lapack_interface.h"
#include "error.h"
#include "common.h"

void compute_gauss_jacobi_quadrature_rule(SizeType n, NumType alpha, 
    NumType beta, CArray<NumType> &points, CArray<NumType> &weights);
