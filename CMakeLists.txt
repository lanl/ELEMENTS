################################################################################
# Preliminaries
################################################################################
project(ELEMENTS)
cmake_minimum_required(VERSION 3.2)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 14)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# Compiler flags
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} 
    ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")


################################################################################
# Modules, macros, and packages
################################################################################
# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

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

# Optional Doxygen and Sphinx dependencies
option(WITH_DOCS "Build the documentation locally" OFF)
if (WITH_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

# MATAR submodule
# (Note: the way the MATAR submodule is handled may need to change in the
# future. Currently, it is deployed as a header-only library, which is why it
# is sufficient to simply include the MATAR header. However, the requirements
# of linking Kokkos and other parts of Trilinos may change this.  Please update
# this as needed.)
add_subdirectory(submodules/matar)
include_directories(submodules/matar/src)

# GoogleTest submodule
add_subdirectory(submodules/googletest)


################################################################################
# Build
################################################################################
# Includes
include_directories(include)

# Source
add_subdirectory(src)

# Tests and examples
add_subdirectory(tests)
add_subdirectory(examples)

# Documentation
if (WITH_DOCS)
  add_subdirectory(docs)
endif()

# Install headers in include/ subdirectory of the specified install directory
install(
  DIRECTORY include/ 
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
)

# Install libraries in lib/ subdirectory of the specified install directory
install(
  TARGETS 
    common
    elements 
    geometry 
    slam 
    swage
  DESTINATION 
    lib
)
if (ENABLE_VTK AND VTK_FOUND)
  install(TARGETS swage2vtk DESTINATION lib)
endif()
