#include "error.h"

void error_decoder(ErrorCode error) {  
  switch (error) {
    case (NOT_A_MATRIX_ERROR):
      std::cerr << "Error: input/ouput arrays are not all 2D." << std::endl;
      break;
    case (NON_MATCHING_DIMS_ERROR):
      std::cerr << "Error: matrix dimensions don't match." << std::endl;
      break;
    case (NON_SQUARE_MATRIX_ERROR):
      std::cerr << "Error: matrices are not square." << std::endl;
      break;
    case (INCOMPATIBLE_DIMS_ERROR):
      std::cerr << "Error: matrix dimensions are incompatible" << std::endl;
      break;
    case (LINEAR_ALGEBRA_ERROR + 1):
      std::cerr << "Error: unsuccesful LU factorization (LAPACK)" << std::endl;
      break;
    case (LINEAR_ALGEBRA_ERROR + 2):
      std::cerr << "Error: unsuccesful LU solution (LAPACK)" << std::endl;
      break;
    case (VANDERMONDE_ERROR + 1):
      std::cerr << "Error: the input order must be >= 0." << std::endl;
      break;
    case (VANDERMONDE_ERROR + 2):
      std::cerr << "Error: incorrect numbers of dimensions "
                << "in input arrays" 
                << std::endl;
      break;
    case (SOLUTION_ERROR + 1):
      std::cerr << "Error: the leading dimension of the output array " 
                << "must be equal to (order + 1)." 
                << std::endl;
      break;
    case (VTK_IO_ERROR + 1):
      std::cerr << "Error: the input mesh file name is missing the supported " 
                << "extensions (.vtk, .vtu)" 
                << std::endl;
      break;
    case (VTK_IO_ERROR + 2):
      std::cerr << "Error: the number of nodes for the current VTK cell does " 
                << "not match the number of nodes expected for this element "
                << "in SWAGE" 
                << std::endl;
      break;
    default:
      std::cerr << "Error: unknown error." << std::endl;
  }
}

void check_error(const char *file, const char *func, int line, ErrorCode error) {
  if (error) {
    std::cerr << std::endl;
    error_decoder(error);
    std::cerr << "In file:     " << file << std::endl;
    std::cerr << "In function: " << func << std::endl;
    std::cerr << "On line:     " << line << std::endl;
    std::cerr << std::endl;
    exit(error);
  }
}
