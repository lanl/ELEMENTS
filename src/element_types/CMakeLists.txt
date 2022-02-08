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

target_include_directories(
  elements
    PUBLIC
      ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
      ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)

target_link_libraries(
  elements
    common
)

