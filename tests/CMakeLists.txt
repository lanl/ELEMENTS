# Test SLAM
add_executable(
  test_slam_math
  test_slam_math.cpp
)

target_link_libraries(
  slam
)

# Test Lagrange polynomials
add_executable(
  test_lagrange_polynomials
  test_lagrange_polynomials.cpp
)

target_link_libraries(
  test_lagrange_polynomials
  elements
)

# Test Lagrange element
add_executable(
  test_lagrange_element
  test_lagrange_element.cpp
)
target_link_libraries(
  test_lagrange_element
  elements
)

# Test Legendre polynomials
add_executable(
  test_legendre_polynomials
  test_legendre_polynomials.cpp
)
target_link_libraries(
  test_legendre_polynomials
  elements
)

# Test Legendre element
add_executable(
  test_legendre_element
  test_legendre_element.cpp
)
target_link_libraries(
  test_legendre_element
  elements
)

# Test Gauss-Jacobi quadrature
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
