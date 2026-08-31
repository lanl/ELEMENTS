cmake_minimum_required(VERSION 3.22)
project(ELEMENTS VERSION 1.0.0 LANGUAGES CXX C)

# --- Global Compiler Fixes (Scoped to this project) ---
# These settings apply to ELEMENTS and all its dependencies (Scotch, etc.)
# but DO NOT leak up to the parent project (the solver).
# Force C11 Standard (Required for PT-Scotch to build with GCC 15)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

# --- Check Scope for Options ---
# We still use this check for things like "Should we build examples?"
# because a solver user likely doesn't want to compile your examples.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(ELEMENTS_IS_TOP_LEVEL ON)
else()
    set(ELEMENTS_IS_TOP_LEVEL OFF)
endif()

# --- Check if building documentation ---
# The default is OFF, but this option can be set to ON from the terminal/cmake command line.
option(ELEMENTS_BUILD_DOCS "Build documentation" OFF)
# Setting -DELEMENTS_BUILD_DOCS=ON from the terminal or cmake-gui will overwrite the default.


# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)


# --- Optional Doxygen and Sphinx dependencies ---
if (ELEMENTS_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()


# --- Build documentation ---
if (ELEMENTS_BUILD_DOCS)
  add_subdirectory(docs)
endif()


# --- Set build type if not already set ---
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# --- Find MPI (Required) ---
find_package(MPI REQUIRED)
if(NOT MSVC)
    if(MPI_CXX_VERSION LESS 3)
        message(WARNING "ELEMENTS requires MPI version 3+. Found ${MPI_CXX_VERSION}.")
    endif()
endif()

# --- Options ---
#
# Backend selection: ELEMENTS is built on Kokkos (via MATAR), and the backend
# is expressed as a host/device pair, matching MATAR's MATAR_HOST_BACKEND /
# MATAR_DEVICE_BACKEND front door. Leaving both empty keeps the historical
# default of a serial build.
set(ELEMENTS_HOST_BACKEND   "" CACHE STRING "Host-side backend: serial|openmp|pthreads")
set(ELEMENTS_DEVICE_BACKEND "" CACHE STRING "Device-side backend: serial|openmp|pthreads|cuda|hip|sycl")
set_property(CACHE ELEMENTS_HOST_BACKEND   PROPERTY STRINGS serial openmp pthreads)
set_property(CACHE ELEMENTS_DEVICE_BACKEND PROPERTY STRINGS serial openmp pthreads cuda hip sycl)

option(ELEMENTS_ENABLE_GPU_AWARE_MPI "Assume the MPI implementation is GPU-aware" OFF)

# Precision tiers, forwarded straight through to MATAR's MATAR_REAL /
# MATAR_HIGH_REAL / MATAR_LOW_REAL (see matar/CMakeLists.txt for the full
# semantics and validation of these values).
set(ELEMENTS_REAL      "double" CACHE STRING "Precision of real_t: double|float|half|bfloat16|quad")
set(ELEMENTS_HIGH_REAL "double" CACHE STRING "Precision of high_real_t: double|float|quad")
set(ELEMENTS_LOW_REAL  "double" CACHE STRING "Precision of low_real_t: double|float|half|bfloat16")

option(ELEMENTS_INSTALL "Generate install/export rules" ${ELEMENTS_IS_TOP_LEVEL})

option(ELEMENTS_BUILD_EXAMPLES "Build examples" ${ELEMENTS_IS_TOP_LEVEL})
option(ELEMENTS_BUILD_TESTS "Build unit tests" OFF)

# --- Deprecated backend options ---
# Map the old per-backend booleans onto the new host/device strings. These are
# intentionally NOT declared with option() -- that would define every one of
# them unconditionally and fire every deprecation warning on every build.
# Only variables the user actually set on the command line are DEFINED here.
if(DEFINED ELEMENTS_ENABLE_SERIAL)
    message(DEPRECATION "ELEMENTS_ENABLE_SERIAL is deprecated; use ELEMENTS_DEVICE_BACKEND=serial (and/or ELEMENTS_HOST_BACKEND=serial)")
    if(ELEMENTS_ENABLE_SERIAL AND NOT ELEMENTS_DEVICE_BACKEND)
        set(ELEMENTS_DEVICE_BACKEND "serial")
    endif()
endif()
if(DEFINED ELEMENTS_ENABLE_OPENMP)
    message(DEPRECATION "ELEMENTS_ENABLE_OPENMP is deprecated; use ELEMENTS_DEVICE_BACKEND=openmp (and/or ELEMENTS_HOST_BACKEND=openmp)")
    if(ELEMENTS_ENABLE_OPENMP)
        set(ELEMENTS_DEVICE_BACKEND "openmp")
    endif()
endif()
if(DEFINED ELEMENTS_ENABLE_PTHREADS)
    message(DEPRECATION "ELEMENTS_ENABLE_PTHREADS is deprecated; use ELEMENTS_DEVICE_BACKEND=pthreads (and/or ELEMENTS_HOST_BACKEND=pthreads)")
    if(ELEMENTS_ENABLE_PTHREADS)
        set(ELEMENTS_DEVICE_BACKEND "pthreads")
    endif()
endif()
if(DEFINED ELEMENTS_ENABLE_CUDA)
    message(DEPRECATION "ELEMENTS_ENABLE_CUDA is deprecated; use ELEMENTS_DEVICE_BACKEND=cuda")
    if(ELEMENTS_ENABLE_CUDA)
        set(ELEMENTS_DEVICE_BACKEND "cuda")
    endif()
endif()
if(DEFINED ELEMENTS_ENABLE_HIP)
    message(DEPRECATION "ELEMENTS_ENABLE_HIP is deprecated; use ELEMENTS_DEVICE_BACKEND=hip")
    if(ELEMENTS_ENABLE_HIP)
        set(ELEMENTS_DEVICE_BACKEND "hip")
    endif()
endif()

# Historical default: a bare `cmake ..` with nothing specified builds serial.
if(NOT ELEMENTS_HOST_BACKEND AND NOT ELEMENTS_DEVICE_BACKEND
   AND NOT DEFINED Kokkos_ENABLE_SERIAL AND NOT DEFINED Kokkos_ENABLE_OPENMP
   AND NOT DEFINED Kokkos_ENABLE_THREADS AND NOT DEFINED Kokkos_ENABLE_CUDA
   AND NOT DEFINED Kokkos_ENABLE_HIP AND NOT DEFINED Kokkos_ENABLE_SYCL)
    set(ELEMENTS_DEVICE_BACKEND "serial")
endif()

include(FetchContent)

# Kokkos 5's CUDA backend needs CMake >= 3.25.2 for C++20 support in the CUDA
# language; catch this before MATAR's own configure gets further into it.
if((ELEMENTS_DEVICE_BACKEND STREQUAL "cuda" OR Kokkos_ENABLE_CUDA) AND CMAKE_VERSION VERSION_LESS 3.25.2)
    message(FATAL_ERROR
        "A CUDA backend requires CMake >= 3.25.2 (found ${CMAKE_VERSION}) for "
        "C++20 support in the CUDA language.")
endif()

# --- Dependency: MATAR (git submodule; bundles Kokkos as a nested submodule) ---
if(NOT TARGET matar::matar)
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/matar/CMakeLists.txt)
        message(FATAL_ERROR
            "The MATAR submodule is missing. Run:\n"
            "    git submodule update --init --recursive")
    endif()

    set(MATAR_ENABLE_KOKKOS        ON  CACHE BOOL "" FORCE)   # ELEMENTS is Kokkos-based
    set(MATAR_ENABLE_MPI           ON  CACHE BOOL "" FORCE)   # ELEMENTS requires MPI
    set(MATAR_ENABLE_GPU_AWARE_MPI ${ELEMENTS_ENABLE_GPU_AWARE_MPI} CACHE BOOL "" FORCE)
    set(MATAR_HOST_BACKEND   "${ELEMENTS_HOST_BACKEND}"   CACHE STRING "" FORCE)
    set(MATAR_DEVICE_BACKEND "${ELEMENTS_DEVICE_BACKEND}" CACHE STRING "" FORCE)
    set(MATAR_REAL      "${ELEMENTS_REAL}"      CACHE STRING "" FORCE)
    set(MATAR_HIGH_REAL "${ELEMENTS_HIGH_REAL}" CACHE STRING "" FORCE)
    set(MATAR_LOW_REAL  "${ELEMENTS_LOW_REAL}"  CACHE STRING "" FORCE)
    set(MATAR_BUILD_EXAMPLES   OFF CACHE BOOL "" FORCE)
    set(MATAR_BUILD_TESTS      OFF CACHE BOOL "" FORCE)
    set(MATAR_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
    set(MATAR_INSTALL ${ELEMENTS_INSTALL} CACHE BOOL "" FORCE)

    add_subdirectory(matar)
endif()

# --- Dependency: Scotch / PT-Scotch ---
if(NOT TARGET ptscotch)
    set(SCOTCH_MPI ON CACHE BOOL "Enable MPI for Scotch" FORCE)
    set(SCOTCH_BUILD_TESTS OFF CACHE BOOL "Disable Scotch Tests" FORCE)
    set(BUILD_FORTRAN OFF CACHE BOOL "Disable Scotch Fortran interface" FORCE)
    if(WIN32)
        # Scotch does not ship pthreads on Windows; ensure all pthread toggles are off
        set(SCOTCH_PTHREAD OFF CACHE BOOL "Disable Scotch pthread on Windows" FORCE)
        set(SCOTCH_PTHREAD_MPI OFF CACHE BOOL "Disable Scotch pthread MPI on Windows" FORCE)
        set(COMMON_PTHREAD OFF CACHE BOOL "Disable common pthread on Windows" FORCE)
    endif()

    FetchContent_Declare(
        Scotch
        GIT_REPOSITORY https://gitlab.inria.fr/scotch/scotch.git
        GIT_TAG        v7.0.4
        GIT_SHALLOW    TRUE
        PATCH_COMMAND  ${CMAKE_COMMAND} -DFileToPatch=CMakeLists.txt -P "${PROJECT_SOURCE_DIR}/cmake/PatchScotch.cmake"
    )
    FetchContent_MakeAvailable(Scotch)

    # FIX: GCC 15 Warning Suppressions (Target-Specific)
    # Even though C11 is set globally above, we still need to silence the
    # specific GCC 15 warnings for Scotch that are not fixed by C11 alone.
    if(NOT MSVC)
        set(_scotch_targets 
            scotch ptscotch scotcherr scotcherrexit
            amk_ccc amk_fft2 amk_grf amk_hy amk_m2 amk_p2 
            gbase gmap gord gout gpart gscat gtst 
            mcv mmk_m2 mmk_m3 mord mtst
            dggath dgmap dgord dgpart dgscat dgtst
        )

        foreach(_target ${_scotch_targets})
            if(TARGET ${_target})
                # Silence annoying warnings for legacy C code
                target_compile_options(${_target} PRIVATE 
                    $<$<COMPILE_LANGUAGE:C>:-Wno-error=incompatible-pointer-types>
                    $<$<COMPILE_LANGUAGE:C>:-Wno-error=int-conversion>
                    $<$<COMPILE_LANGUAGE:C>:-Wno-error=implicit-int>
                    $<$<COMPILE_LANGUAGE:C>:-Wno-incompatible-pointer-types>
                    $<$<COMPILE_LANGUAGE:C>:-Wno-int-conversion>
                    $<$<COMPILE_LANGUAGE:C>:-Wno-implicit-int>
                )
            endif()
        endforeach()
    endif()
endif()

# --- Add ELEMENTS Source ---
add_subdirectory(src)

# --- Install / export ---
if(ELEMENTS_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(EXPORT ElementsTargets
            NAMESPACE ELEMENTS::
            FILE ElementsTargets.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Elements)

    configure_package_config_file(
        cmake/ElementsConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/ElementsConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Elements)

    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/ElementsConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)

    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/ElementsConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/ElementsConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Elements)
endif()

# --- Optional Parts ---
if(ELEMENTS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(ELEMENTS_BUILD_TESTS)
    if(EXISTS ${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt)
        enable_testing()
        add_subdirectory(tests)
    else()
        message(WARNING "ELEMENTS_BUILD_TESTS is ON but no tests/ directory exists; skipping.")
    endif()
endif()
