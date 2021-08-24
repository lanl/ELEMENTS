#pragma once
#include "common.h"

/* 
 * BLAS/LAPACK naming convention help
 *
 * GEMV  - general matrix-vector multiplication
 * GEMM  - general matrix-matrix multiplication
 * GETRF - general triangular (LU) factorization
 * GETRS - general triangular (LU) solution
 * STEV  - symmetric tridiagonal eigenvalues/vectors
 */

// Change definitions of routines
#ifdef USE_COMPLEX_NUMBERS
#define blas_gemv zgemv_
#define blas_gemm zgemm_
#define lapack_getrf zgetrf_
#define lapack_getrs zgetrs_
#define lapack_stev dstev_
#else
#define blas_gemv dgemv_
#define blas_gemm dgemm_       
#define lapack_getrf dgetrf_
#define lapack_getrs dgetrs_
#define lapack_stev dstev_
#endif

extern "C" {
  /* Real number routines */
  void dgemv_( 	
    const char *,
		int *, int *,
		NumType *,
		NumType *, int *,
		NumType *, int *,
		NumType *,
		NumType *, int * 	
	);

  void dgemm_(
    const char *, const char *, 
    int	*, int	*, int	*, 
    NumType *, NumType *, int *, 
    NumType *, int *, 
    NumType *, NumType *, int *
  );

  void dgetrf_(int*, int *, NumType *, int	*, int *, int	*);

  void dgetrs_(
    const char *, 
    int *, int *, 
    NumType *, int *, 
    int *, 
    NumType *, int *, 
    int *);

  void dstev_(
    const char *,
		int *,
		RealNumber *,
		RealNumber *,
		RealNumber *, int *,
		RealNumber *,
		int *
	);

  // Complex number routines
  void zgemv_( 	
    const char *,
		int *, int *,
		NumType *,
		NumType *, int *,
		NumType *, int *,
		NumType *,
		NumType *, int * 	
	);

  void zgemm_(
    const char *, const char *, 
    int	*, int	*, int	*, 
    NumType *, NumType *, int *, 
    NumType *, int *, 
    NumType *, NumType *, int *
  );

  void zgetrf_(int*, int *, NumType *, int	*, int *, int	*);

  void zgetrs_(
    const char *, 
    int *, int *, 
    NumType *, int *, 
    int *, 
    NumType *, int *, 
    int *);
}
