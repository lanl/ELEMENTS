# Testing Lagrange polynomials
add_executable(
  test_lagrange_1d_derivatives 
  test_lagrange_1d_derivatives.cpp
)

add_executable(
  test_slam_math
  test_slam_math.cpp
)
target_link_libraries(
  slam
)

target_link_libraries(
  test_lagrange_1d_derivatives 
  common
  polynomials_complex
)

target_compile_options(
  test_lagrange_1d_derivatives 
  PRIVATE
  -DUSE_COMPLEX_NUMBERS
)

add_executable(
  test_lagrange_3d_derivatives 
  test_lagrange_3d_derivatives.cpp
)
target_link_libraries(
  test_lagrange_3d_derivatives 
  common
  polynomials_complex
)
target_compile_options(
  test_lagrange_3d_derivatives 
  PRIVATE
  -DUSE_COMPLEX_NUMBERS
)

add_executable(
  test_lagrange_jacobian_determinant 
  test_lagrange_jacobian_determinant.cpp
)
target_link_libraries(
  test_lagrange_jacobian_determinant 
  common
  elements
)


# Test Legendre polynomials
add_executable(
  test_legendre_1d_derivatives 
  test_legendre_1d_derivatives.cpp
)
target_link_libraries(
  test_legendre_1d_derivatives 
  common
  polynomials_complex
)
target_compile_options(
  test_legendre_1d_derivatives 
  PRIVATE
  -DUSE_COMPLEX_NUMBERS
)


# Test Gauss-Jacobi quadrature
# (Only if BLAS/LAPACK interface enabled)
if (ENABLE_BLAS_LAPACK)
  add_executable(
    test_gauss_jacobi_quadrature
    test_gauss_jacobi_quadrature.cpp
  )
  target_link_libraries(
    test_gauss_jacobi_quadrature
    common
    elements
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  )
endif()


# Test capability to read mesh and write solution using VTK
