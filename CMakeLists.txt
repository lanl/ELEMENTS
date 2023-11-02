################################################################################
# Preliminaries
################################################################################
project(ELEMENTS)

# CMake version > 3.8 required for C++17 standard
cmake_minimum_required(VERSION 3.8)

# C++17 standard required for MATAR
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)


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
  # Have user provide VTK_DIR flag to identify VTK installation location and
  # VTK_SOURCE_DIR to identify VTK source location
  find_package(VTK REQUIRED NO_MODULE)
endif()

# Optional Doxygen and Sphinx dependencies
option(WITH_DOCS "Build the documentation locally" OFF)
if (WITH_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

add_subdirectory(matar)
include_directories(matar/src/include)


################################################################################
# Build
################################################################################
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
if (BUILD_TESTS)
  add_subdirectory(tests)
endif()

# Examples
add_subdirectory(examples)

# Documentation
if (WITH_DOCS)
  add_subdirectory(docs)
endif()

# mesh tools
add_subdirectory(mesh_tools)


################################################################################
# Install
################################################################################
# Install headers in include/ subdirectory of the specified install directory
# (Note: This is an inelegant way of installing the headers. Resursing through
# subdirectories looking for anything that matches the pattern "*.h" is
# indiscriminate, and the resulting include/ directory isn't affected by
# distclean. If there's a better way to do this, please introduce it here.)
add_library(Elements SHARED src/elements_lib_dummy.cpp)

file(GLOB_RECURSE ALL_HEADERS elements/*.h geometry/*.h slam/*.h swage/*.h common/*.h swage2vtk/*.h)
install(FILES ${ALL_HEADERS} DESTINATION include)

target_include_directories(Elements PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Elements>  
    $<INSTALL_INTERFACE:include/Elements>
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/ElementsConfigVersion.cmake"
    VERSION 1.0
    COMPATIBILITY AnyNewerVersion
)

install(TARGETS Elements
    EXPORT ElementsTargets
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development
    RUNTIME DESTINATION bin COMPONENT Runtime
    PUBLIC_HEADER DESTINATION include COMPONENT Development
    BUNDLE DESTINATION bin COMPONENT Runtime
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/ElementsConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/ElementsConfig.cmake"
    INSTALL_DESTINATION lib/cmake/Elements
)

set_target_properties(Elements
    PROPERTIES
    PUBLIC_HEADER "${ALL_HEADERS}"
)

install(EXPORT ElementsTargets DESTINATION lib/cmake/Elements)
install(FILES "${PROJECT_BINARY_DIR}/ElementsConfigVersion.cmake"
    "${PROJECT_BINARY_DIR}/ElementsConfig.cmake"
        DESTINATION lib/cmake/elements)
    #install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/include/ DESTINATION include)

# Install libraries in lib/ subdirectory of the specified install directory
install(TARGETS elements EXPORT ElementsTargets DESTINATION lib)
install(TARGETS geometry EXPORT ElementsTargets DESTINATION lib)
install(TARGETS slam EXPORT ElementsTargets DESTINATION lib)
install(TARGETS swage EXPORT ElementsTargets DESTINATION lib)
if (WITH_VTK AND VTK_FOUND)
    install(TARGETS swage2vtk EXPORT ElementsTargets DESTINATION lib)
endif()
install(TARGETS common EXPORT ElementsTargets DESTINATION lib)

target_link_libraries(Elements elements geometry slam swage common)
#target_link_libraries(Elements geometry)
if (WITH_VTK AND VTK_FOUND)
    target_link_libraries(Elements PUBLIC swage2vtk)
endif()
target_link_libraries(Elements matar)
if(Matar_ENABLE_KOKKOS)
    target_link_libraries(Elements Kokkos::kokkos)
endif()

