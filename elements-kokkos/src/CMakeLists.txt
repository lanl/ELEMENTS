add_library(
  elements_kokkos
    HexRef.cpp
    LagrangeProductBasis.cpp  
    LegendreProductBasis.cpp  
    BernsteinProductBasis.cpp  
    TensorProductBasis.cpp
    lagrange_polynomials.cpp  
    legendre_polynomials.cpp  
    jacobi_polynomials.cpp  
    bernstein_polynomials.cpp  
    points_and_weights.cpp
)

target_compile_definitions(
  elements_kokkos
  PUBLIC
    MATAR_WITH_KOKKOS
    HAVE_KOKKOS=1
    HAVE_OPENMP=1
)

target_include_directories(
  elements_kokkos
  PUBLIC
    ${ELEMENTS_KOKKOS_INCLUDES}
)


if (CUDA)
  target_compile_definitions(elements_kokkos PUBLIC -DHAVE_CUDA=1)
elseif (HIP)                        
  target_compile_definitions(elements_kokkos PUBLIC -DHAVE_HIP=1)
elseif (OPENMP)                    
  target_compile_definitions(elements_kokkos PUBLIC -DHAVE_OPENMP=1)
elseif (THREADS)                  
  target_compile_definitions(elements_kokkos PUBLIC -DHAVE_THREADS=1)
else()
  message(FATAL_ERROR "Please specify whether you are using CUDA, HIP, OPENMP, etc.")
endif()

target_link_libraries(
  elements_kokkos
  PUBLIC
    common
    matar
    Kokkos::kokkos
)
