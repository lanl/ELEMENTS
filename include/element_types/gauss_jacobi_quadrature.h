#pragma once

#include "element_types/jacobi_polynomials.h"

#include "common/matar_blas_lapack_interface.h"
#include "common/error.h"
#include "common/common.h"

void compute_gauss_jacobi_quadrature_rule(SizeType n, Real alpha, 
    Real beta, CArray<Real> &points, CArray<Real> &weights);
