################################################################################
# Preliminaries
################################################################################
project(ELEMENTS)

# CMake version > 3.8 required for C++17 standard
cmake_minimum_required(VERSION 3.8)

# C++17 standard required for MATAR
set(CMAKE_CXX_STANDARD 17)

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

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Public header files
set(ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include)

################################################################################
# Dependencies
################################################################################
# Optional BLAS/LAPACK dependency (disabled by default)
option(WITH_BLAS_LAPACK "Build ELEMENTS with BLAS/LAPACK" OFF)
if (WITH_BLAS_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# Optional VTK dependency (disabled by default)
option(WITH_VTK "Build ELEMENTS with VTK" OFF)
if (WITH_VTK)
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
set(MATAR_PUBLIC_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/submodules/matar/src)

# GoogleTest submodule
option(WITH_GTEST "Use GoogleTest as testing framework" OFF)
if (WITH_GTEST)
  add_subdirectory(submodules/googletest)
endif()


################################################################################
# Build & install
################################################################################
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
# Install non-example, non-test executables in bin/
install(
  TARGETS 
    common
    elements 
    geometry 
    slam 
    swage
    Meshbox
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

if (WITH_VTK AND VTK_FOUND)
  install(TARGETS swage2vtk DESTINATION lib)
endif()
