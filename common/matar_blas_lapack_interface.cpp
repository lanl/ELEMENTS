#include "matar_blas_lapack_interface.h"
#include "c_blas_lapack_interface.h"

/*
 * Transpose a 2D MATAR CArray and return the result in B
 *
 * B = A^{T}
 *
 */
ErrorCode transpose_2d_matar_carray(const CArray<NumType> &A, 
    CArray<NumType> &B) {
  // Throw error if MATAR CArrays are not both 2D, i.e. they're not matrices
  int 
  num_dim_a = A.order(),
  num_dim_b = B.order();

  bool both_are_2d = num_dim_a == 2 and num_dim_b == 2;

  if (not both_are_2d) return(NOT_A_MATRIX_ERROR);

  // Throw error if the dimensions of A and B do not match
  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  m_b = B.dims(0),
  n_b = B.dims(1);

  bool matching_dims = m_a == n_b and n_a == m_b;
  if (not matching_dims) return (INCOMPATIBLE_DIMS_ERROR);

  // Put the transpose of A into B
  for (int i = 0; i < m_a; i++)
    for (int j = 0; j < n_a; j++)
      B(j,i) = A(i,j);

  return(0);
}

/*
 * Compute the matrix-vector multiplication of two MATAR CArrays, A (2D) and x
 * (1D), using BLAS's matrix-vector multiplication routine (dgemv) and return
 * the result in y
 *
 * y = A \times x
 *
 */
ErrorCode matvec_matar_carray_blas(const CArray<NumType> &A, 
    const CArray<NumType> &x, CArray<NumType> &y) {
  // Check if A is 2D, i.e. that it's a matrix
  int num_dim_a = A.order();
  bool a_is_2d = num_dim_a == 2;
  if (not a_is_2d) return(NOT_A_MATRIX_ERROR);

  // Check that x and y are both 1D, i.e. that they're vectors
  int num_dim_x = x.order(), num_dim_y = y.order();
  bool x_and_y_are_1d = num_dim_x == 1 and num_dim_y == 1;
  if (not x_and_y_are_1d) return(NOT_A_VECTOR_ERROR);

  // Check that the dimensions of A, x, and y are compatible for matvec
  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  n_x = x.dims(0),
  m_y = y.dims(0);
  bool compatible_dims = m_y == m_a and n_a == n_x;
  if (not compatible_dims) return(INCOMPATIBLE_DIMS_ERROR);
  
  int k = n_a;

  /*
   * Copy the contents of A into regular, unwrapped arrays in column major
   * order, which what the BLAS routine requires
   */
  NumType *a = new NumType[m_a*n_a];
  for (int i = 0; i < m_a; i++)
    for (int j = 0; j < n_a; j++)
      a[i+n_a*j] = A(i,j);

  /*
   * Use dgemv from BLAS to compute the matrix-vector multiplication. Note that
   * dgemv computes 
   *
   *   y = \alpha A \times x + \beta y. 
   *
   * In this case, we set \alpha = 1 and \beta = 0.
   */
  const char if_transpose = 'N';  // whether to tranpose matrix A

  NumType 
  alpha = 1.0,  // scaling of matrix-vector multiplication
  beta  = 0.0;  // Scaling of input/output vector y

  int 
  lda  = m_a,  // leading dimension of matrix A
  incx = 1,    // increment of vector x
  incy = 1;    // increment of vector y

  blas_gemv( 	
    &if_transpose,
		&m_a, &n_a,
		&alpha,
		a, &lda,
		x.get_pointer(), &incx,
		&beta,
		y.get_pointer(), &incy 	
	);

  delete[] a;

  return(0);
}

/*
 * Compute the multiplication of two 2D MATAR CArrays, A and B, using BLAS's
 * matrix-matrix multiplication routine (dgemm) and return the result in C
 *
 * C = A \times B
 *
 */
ErrorCode multiply_2d_matar_carray_blas(const CArray<NumType> &A, 
    const CArray<NumType> &B, CArray<NumType> &C) {
  // Throw error if MATAR CArrays are not all 2D, i.e. they're not matrices
  int 
  num_dim_a = A.order(),
  num_dim_b = B.order(),
  num_dim_c = C.order();
  bool all_are_2d = num_dim_a == 2 and num_dim_b == 2 and num_dim_c == 2;
  if (not all_are_2d) return(NOT_A_MATRIX_ERROR);

  // Throw error if dimensions of A and B aren't compatible for multiplication
  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  m_b = B.dims(0),
  n_b = B.dims(1);
  bool compatible_dims = n_a == m_b;
  if (not compatible_dims) return(INCOMPATIBLE_DIMS_ERROR);
  
  int k = n_a;

  /*
   * Copy the contents of the input MATAR CArrays into regular, unwrapped
   * arrays in column major order, which what the BLAS routine requires
   */
  NumType *a = new NumType[m_a*n_a];
  for (int i = 0; i < m_a; i++)
    for (int j = 0; j < n_a; j++)
      a[i+n_a*j] = A(i,j);

  NumType *b = new NumType[m_b*n_b];
  for (int i = 0; i < m_b; i++)
    for (int j = 0; j < n_b; j++)
      b[i+n_b*j] = B(i,j);

  /*
   * Use dgemm from BLAS to compute the matrix-matrix multiplication. Note that
   * dgemm computes 
   *
   *   C = \alpha A \times B + \beta C. 
   *
   * In this case, we set \alpha = 1 and \beta = 0.
   */
  const char if_transpose_a = 'N';  // whether to tranpose matrix A
  const char if_transpose_b = 'N';  // whether to tranpose matrix B

  NumType 
  alpha = 1.0,  // scaling of matrix-matrix multiplication
  beta  = 0.0;  // Scaling of input/output matrix C

  int 
  lda = m_a,  // leading dimension of matrix A
  ldb = m_b,  // leading dimension of matrix B
  ldc = m_a;  // leading dimension of matrix C

  NumType *c = new NumType[m_a*n_b];

  blas_gemm(
      &if_transpose_a, &if_transpose_b, 
      &m_a, &n_b, &k, 
      &alpha, 
      a, &lda, 
      b, &ldb, 
      &beta, 
      c, &lda
    );

  // Copy the result into the output MATAR CArray
  for (int j = 0; j < n_b; j++) {
    for (int i = 0; i < m_a; i++) {
      C(i,j) = c[i+n_b*j];
    }
  }

  delete[] a, b, c;

  return(0);
}

/*
 * Compute the inverse of a 2D MATAR CArray A using LAPACK's LU factorization
 * and solution routines (dgetrf/dgetrf) and return the result in B
 *
 * B = A^{-1}
 *
 */
ErrorCode invert_2d_matar_carray_lapack(const CArray<NumType> &A, CArray<NumType> &B) {
  // Throw error if MATAR CArrays are not both 2D, i.e. they're not matrices
  int num_dim_a = A.order();
  int num_dim_b = B.order();
  bool both_are_2d = num_dim_a == 2 and num_dim_b == 2;
  if (not both_are_2d) return(NOT_A_MATRIX_ERROR);
  
  // Throw error if matrices do not have same dimensions
  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  m_b = B.dims(0),
  n_b = B.dims(1);
  bool matrices_have_same_dims = m_a == m_b and n_a == n_b;
  if (not matrices_have_same_dims) return(NON_MATCHING_DIMS_ERROR);
  
  // Throw error if matrices are not square
  int m = m_a;
  int n = n_a;
  bool matrices_are_square = m == n;
  if (not matrices_are_square) return(NON_SQUARE_MATRIX_ERROR);

  /*
   * Copy the contents of the input MATAR CArray into a regular, unwrapped
   * array in column major order, which what the LAPACK routine requires
   */
  NumType *a = new NumType[m*n];
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      a[i+n*j] = A(i,j);

  // Compute LU factorization of matrix A
  int 
  lda = m,      // leading dimension of matrix A
  info_factor;  // error code

  int *ipiv = new int[n];  // pivot indices in LU factorization

  lapack_getrf(&m, &n, a, &lda, ipiv, &info_factor);

  // Throw error if unsuccessful factorization
  if (info_factor != 0) return(LINEAR_ALGEBRA_ERROR + 1);

  // Populate right-hand side array (identity since A = L*U, L*U*A^{-1} = I)
  NumType *b = new NumType[m*n];
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      b[i+n*j] = (i == j) ? 1 : 0;

  // Use the LU factorization to compute the inverse
  const char if_transpose = 'N';  // whether A should be transposed

  int 
  nrhs  = n,   // number of right-hand sides
  ldb   = m,   // leading dimension of matrix B
  info_solve;  // error code 

  lapack_getrs(&if_transpose, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info_solve);

  // Throw error if unsuccessful solution
  if (info_solve < 0) return(LINEAR_ALGEBRA_ERROR + 2);

  // Copy the result into the output MATAR CArray
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      B(i,j) = b[i+n*j];
    }
  }
  
  delete[] a, ipiv, b;

  return(0);
}
