
add_executable(virttestkokkos child.cpp  child.hpp  inherited_inits.cpp  inherited_inits.hpp  kokkos_alias.h  main_kokkos_vfcn.cpp  parents.h)
target_link_libraries (virttestkokkos matar ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
add_definitions(-DHAVE_KOKKOS=1 -DHAVE_OPENMP=1)