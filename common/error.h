#pragma once

#include "common.h"

typedef int ErrorCode;

enum ErrorCodeType {
  NOT_A_VECTOR_ERROR      = 100,
  NOT_A_MATRIX_ERROR      = 200,
  NON_MATCHING_DIMS_ERROR = 300,
  NON_SQUARE_MATRIX_ERROR = 400,
  INCOMPATIBLE_DIMS_ERROR = 500,
  LINEAR_ALGEBRA_ERROR    = 600,
  VANDERMONDE_ERROR       = 700,
  SOLUTION_ERROR          = 800,
  VTK_IO_ERROR            = 900,
  QUADRATURE_ERROR        = 1000
};

void error_handler(ErrorCode);

void check_error(ErrorCode, const char *, const char *, int);

#define check(error) (check_error(error, __FILE__, __FUNCTION__, __LINE__))
