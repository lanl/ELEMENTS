cmake_minimum_required(VERSION 3.1)
project(mtest CXX C)
#project(unittest_esn CXX C)
#set(CMAKE_C_FLAGS "-g -G")
#set(CMAKE_CXX_FLAGS "-g -G")
# Enable verbose
# set(CMAKE_C_FLAGS "-g -G") # VERBOSE=1")
# set(CMAKE_CXX_FLAGS "-g -G") # VERBOSE=1")

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


add_executable(mtest main.cpp pseudo_mesh.hpp pseudo_mesh.cpp matar.h kokkos_alias.h)
#add_executable(dummy dummy.cpp matar_test.h kokkos_alias.h)
#add_executable(mtest main.cpp pseudo_mesh.hpp pseudo_mesh.cpp pseudo_var.hpp pseudo_var.cpp matar.h)
#add_executable(unittest_esn main_unittest_esn.cpp pseudo_var.hpp pseudo_var.cpp matar.h)

include_directories(./)

add_subdirectory(Kokkos/kokkos)

target_link_libraries(mtest Kokkos::kokkos)
#target_link_libraries(dummy Kokkos::kokkos)
#target_link_libraries(unittest_esn Kokkos::kokkos)
