if (ENABLE_BLAS_LAPACK)
  add_library(
    elements
    elements.cpp
    lagrange_polynomials.cpp
    lagrange_element.cpp
    legendre_polynomials.cpp
    legendre_element.cpp
    jacobi_polynomials.cpp
    gauss_jacobi_quadrature.cpp
    point_distributions.cpp
  )
else()
  add_library(
    elements
    elements.cpp
    lagrange_polynomials.cpp
    lagrange_element.cpp
    legendre_polynomials.cpp
    legendre_element.cpp
    jacobi_polynomials.cpp
    point_distributions.cpp
  )
endif()
target_link_libraries(
  elements
  common
)

target_include_directories(elements PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}
)
