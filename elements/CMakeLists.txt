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

if (MATAR_WITH_KOKKOS)
  target_compile_definitions(
    elements
    PUBLIC
      MATAR_WITH_KOKKOS
      HAVE_KOKKOS=1
      HAVE_OPENMP=1
  )

  if (CUDA)
    target_compile_definitions(elements PUBLIC -DHAVE_CUDA=1)
  elseif (HIP)
    target_compile_definitions(elements PUBLIC -DHAVE_HIP=1)
  elseif (OPENMP)
    target_compile_definitions(elements PUBLIC -DHAVE_OPENMP=1)
  elseif (THREADS)
    target_compile_definitions(elements PUBLIC -DHAVE_THREADS=1)
  else()
    message(FATAL_ERROR "Please specify whether you are using CUDA, HIP, OPENMP, etc.")
  endif()

  target_link_libraries(
    elements
      common
      matar
      Kokkos::kokkos
  )
else()
  target_link_libraries(
    elements
      common
      matar
  )
endif()
