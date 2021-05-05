set(Solver_SRC_H ../Include/header.h ../Include/state.h ../Include/utilities.h ../Include/math_util.h)
set(Class_SRC Solver.h ../Include/Simulation_Parameters.h ../include/Simulation_Parameters.cpp)
set(Laplace_SRC Pseudo_Laplace/Pseudo_Laplacian.h Pseudo_Laplace/Pseudo_Laplacian.cpp)
set(Static_SRC Static_Solver/Static_Solver.h Static_Solver/Static_Solver.cpp)
set(Static_Parallel_SRC Static_Solver_Parallel/Static_Solver_Parallel.h Static_Solver_Parallel/Static_Solver_Parallel.cpp)
set(Solver_SRC_CXX main.cpp Solver.cpp ../Include/math_util.cpp)

find_package(MPI REQUIRED)
SET(CMAKE_PREFIX_PATH ../Trilinos/build/lib/cmake/Trilinos ${CMAKE_PREFIX_PATH})
#SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos 13.1 REQUIRED)

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

# Make sure to use same compilers and flags as Trilinos
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
SET(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

include_directories(../Include)
include_directories(Static_Solver)
include_directories(Static_Solver_Parallel)
include_directories(Pseudo_Laplace)
include_directories(../MATAR/src)

INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_executable(main ${Solver_SRC_CXX} ${Solver_SRC_H} ${Class_SRC} ${Laplace_SRC} ${Static_SRC} ${Static_Parallel_SRC})

target_link_libraries (main swage)
target_link_libraries (main elements)
target_link_libraries (main matar)
target_include_directories(main PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_compile_options (main PRIVATE ${CMAKE_CXX_COMPILE_FLAGS})
target_link_libraries (main ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

add_definitions(-DHAVE_KOKKOS=1 -DTRILINOS_INTERFACE=1)
