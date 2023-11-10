add_executable(hexref hexref.cpp)
target_compile_definitions(
  hexref 
  PUBLIC 
    MATAR_WITH_KOKKOS 
    HAVE_KOKKOS=1 
    HAVE_OPENMP=1
)
target_include_directories(hexref PUBLIC ${ELEMENTS_KOKKOS_INCLUDES})
target_link_libraries(hexref PUBLIC elements_kokkos)
