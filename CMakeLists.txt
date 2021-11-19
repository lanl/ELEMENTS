################################################################################
# Preliminaries
################################################################################
project(ELEMENTS)
cmake_minimum_required(VERSION 3.2)

set (CMAKE_CXX_STANDARD 14)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# Compiler flags
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")

################################################################################
# Modules, macros, and packages
################################################################################
# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

# Vector module
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

# TargetDistclean macro
INCLUDE(cmake/Modules/TargetDistclean.cmake OPTIONAL)

# Optional BLAS/LAPACK dependency (disabled by default)
option(ENABLE_BLAS_LAPACK "Build ELEMENTS with BLAS/LAPACK" OFF)
if (ENABLE_BLAS_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# Optional VTK dependency (disabled by default)
option(ENABLE_VTK "Build ELEMENTS with VTK" OFF)
if (ENABLE_VTK)
  # Have user provide VTK_DIR flag to identify VTK installation location
  find_package(VTK REQUIRED NO_MODULE)
endif()

################################################################################
# Build organization
################################################################################
# MATAR submodule
# (Note: the way the MATAR submodule is handled may need to change in the
# future. Currently, it is deployed as a header-only library, which is why it
# seemed sufficient to simply include the MATAR header. However, the
# requirements of linking Kokkos and other parts of Trilinos may change this.
# Please update this as needed.)
add_subdirectory(matar)
include_directories(matar/src)

# Includes
include_directories(common)
include_directories(elements)
include_directories(geometry)
include_directories(slam)
include_directories(swage)
include_directories(io)

# Libraries
add_subdirectory(common)
add_subdirectory(elements)
add_subdirectory(geometry)
add_subdirectory(slam)
add_subdirectory(swage)
add_subdirectory(io)

# Tests
add_subdirectory(tests)

# Examples
add_subdirectory(examples)

# mesh tools
add_subdirectory(mesh_tools)

# Install headers in include/ subdirectory of the specified install directory
# (Note: This is an inelegant way of installing the headers. Resursing through
# subdirectories looking for anything that matches the pattern "*.h" is
# indiscriminate, and the resulting include/ directory isn't affected by
# distclean. If there's a better way to do this, please introduce it here.)
file(GLOB_RECURSE ALL_HEADERS *.h)
install(FILES ${ALL_HEADERS} DESTINATION include)

# Install libraries in lib/ subdirectory of the specified install directory
install(TARGETS elements DESTINATION lib)
install(TARGETS geometry DESTINATION lib)
install(TARGETS slam DESTINATION lib)
install(TARGETS swage DESTINATION lib)
install(TARGETS common DESTINATION lib)
