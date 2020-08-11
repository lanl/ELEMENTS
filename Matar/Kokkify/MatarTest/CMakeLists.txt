cmake_minimum_required(VERSION 3.1)
project(mtest CXX C)
#project(unittest_esn CXX C)
#set(CMAKE_C_FLAGS "-g -G")
#set(CMAKE_CXX_FLAGS "-g -G")
# Enable verbose
# set(CMAKE_C_FLAGS "-g -G") # VERBOSE=1")
# set(CMAKE_CXX_FLAGS "-g -G") # VERBOSE=1")

set (CMAKE_CXX_STANDARD 17)

if (KOKKOS)
  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "")
  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
  endif ()
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
endif ()


# adds -g and -02
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

find_package(Vector)
if (CMAKE_VECTOR_NOVEC)
  set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
  set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
endif (CMAKE_VECTOR_NOVEC)
if (CMAKE_VECTOR_VERBOSE)
  set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
  set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
  set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")

#add Likwid module path here

set(mtest_CXX_SRCS main.cpp pseudo_mesh.cpp child.cpp inherited_inits.cpp)
set(mtest_H_SRCS pseudo_mesh.hpp matar.h kokkos_alias.h child.hpp parents.h inherited_inits.hpp)
set(mtest_SRCS ${mtest_CXX_SRCS} ${mtest_H_SRCS})
add_executable(mtest ${mtest_SRCS})

# Executable that runs CPU stream benchmark on MATAR's Kokkos-specific classes
#add_executable(mbench main_benchmark.cpp matar.h kokkos_alias.h)

include_directories(./)

add_subdirectory(Kokkos/kokkos)

target_link_libraries(mtest Kokkos::kokkos)
#target_link_libraries(main_pmcb Kokkos::kokkos)
