
add_executable(mtest main_kokkos.cpp)

target_link_libraries (mtest matar ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

include_directories(virtualFcnKokkos)
add_subdirectory(virtualFcnKokkos)
add_definitions(-DHAVE_KOKKOS=1 -DHAVE_OPENMP=1)

