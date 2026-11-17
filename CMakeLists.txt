cmake_minimum_required(VERSION 3.22)  # Kokkos 5 (bundled via MATAR) minimum; CUDA needs 3.25.2+
project(ELEMENTS VERSION 1.0.0 LANGUAGES CXX C)

# --- Global Compiler Fixes (Scoped to this project) ---
# These settings apply to ELEMENTS and all its dependencies (Scotch, etc.)
# but DO NOT leak up to the parent project (the solver).
# Force C11 Standard (Required for PT-Scotch to build with GCC 15)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

# C++20 is required by Kokkos 5 (bundled via MATAR)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
option(ELEMENTS_ENABLE_SERIAL "Enable Serial support in Kokkos" ON)
option(ELEMENTS_ENABLE_OPENMP "Enable OpenMP support in Kokkos" OFF)
option(ELEMENTS_ENABLE_PTHREADS "Enable PThreads support in Kokkos" OFF)
option(ELEMENTS_ENABLE_CUDA "Enable CUDA support in Kokkos" OFF)
option(ELEMENTS_ENABLE_HIP "Enable HIP support in Kokkos" OFF)

option(ELEMENTS_BUILD_EXAMPLES "Build examples" ${ELEMENTS_IS_TOP_LEVEL})
option(ELEMENTS_BUILD_TESTS "Build unit tests" OFF)

# --- Host / device backend selection (mirrors MATAR_HOST_BACKEND / MATAR_DEVICE_BACKEND) ---
# The device backend runs MATAR's standard macros (FOR_ALL, DO_ALL, ...); the
# host backend runs the _HOST twins (FOR_ALL_HOST, ...). Choosing different
# backends for the two lets CPU work overlap in-flight GPU kernels, e.g.
#
#   cmake -B build -DELEMENTS_DEVICE_BACKEND=cuda -DELEMENTS_HOST_BACKEND=openmp
#
# These are forwarded verbatim to MATAR, which validates the pair and maps it
# onto the Kokkos_ENABLE_* variables. Leaving them empty keeps the
# ELEMENTS_ENABLE_* behavior below.
set(ELEMENTS_HOST_BACKEND   "" CACHE STRING "Host backend for the _HOST macros: serial|openmp|pthreads")
set(ELEMENTS_DEVICE_BACKEND "" CACHE STRING "Device backend for FOR_ALL and friends: serial|openmp|pthreads|cuda|hip|sycl")
set_property(CACHE ELEMENTS_HOST_BACKEND   PROPERTY STRINGS "" serial openmp pthreads)
set_property(CACHE ELEMENTS_DEVICE_BACKEND PROPERTY STRINGS "" serial openmp pthreads cuda hip sycl)

# --- Floating-point precision selection (mirrors MATAR_REAL / HIGH_REAL / LOW_REAL) ---
# These choose what MATAR's three tier types mean for this build:
#   real_t       the default working precision for field data
#   high_real_t  fields that must stay accurate (coordinates, conserved sums)
#   low_real_t   tolerant or bulk-storage fields
#
#   cmake -B build -DELEMENTS_REAL=float -DELEMENTS_HIGH_REAL=double
#
# Forwarded verbatim to MATAR, which validates each value. Leaving one empty
# keeps MATAR's default of double for that tier. high_real_t may not drop
# below float. Only double and float are built and tested by ELEMENTS CI.
set(ELEMENTS_REAL      "" CACHE STRING "Precision of real_t: double|float|half|bfloat16|quad")
set(ELEMENTS_HIGH_REAL "" CACHE STRING "Precision of high_real_t: double|float|quad")
set(ELEMENTS_LOW_REAL  "" CACHE STRING "Precision of low_real_t: double|float|half|bfloat16")
set_property(CACHE ELEMENTS_REAL      PROPERTY STRINGS "" double float half bfloat16 quad)
set_property(CACHE ELEMENTS_HIGH_REAL PROPERTY STRINGS "" double float quad)
set_property(CACHE ELEMENTS_LOW_REAL  PROPERTY STRINGS "" double float half bfloat16)

include(FetchContent)

# --- Dependency: Kokkos ---
# Kokkos is NOT fetched here: the MATAR submodule bundles it and builds it from
# the standard Kokkos cache variables. A parent project that already provides
# Kokkos::kokkos wins -- MATAR detects that target and reuses it.
if(ELEMENTS_HOST_BACKEND OR ELEMENTS_DEVICE_BACKEND)
    # MATAR owns the Kokkos_ENABLE_* mapping and the validation (a host backend
    # cannot be a GPU space, openmp and pthreads cannot coexist, the device
    # backend must outrank the host one, ...), so do not also set them from the
    # ELEMENTS_ENABLE_* toggles.
    set(MATAR_HOST_BACKEND   "${ELEMENTS_HOST_BACKEND}"   CACHE STRING "Set by ELEMENTS_HOST_BACKEND"   FORCE)
    set(MATAR_DEVICE_BACKEND "${ELEMENTS_DEVICE_BACKEND}" CACHE STRING "Set by ELEMENTS_DEVICE_BACKEND" FORCE)
    # MATAR reports any problem with this pair against the MATAR_* names, so
    # record the mapping here.
    message(STATUS
        "ELEMENTS backends: device='${ELEMENTS_DEVICE_BACKEND}' host='${ELEMENTS_HOST_BACKEND}'"
        " -> MATAR_DEVICE_BACKEND / MATAR_HOST_BACKEND")

    foreach(_legacy OPENMP PTHREADS CUDA HIP)
        if(ELEMENTS_ENABLE_${_legacy})
            message(WARNING
                "ELEMENTS_ENABLE_${_legacy}=ON is ignored because ELEMENTS_HOST_BACKEND / "
                "ELEMENTS_DEVICE_BACKEND are set; the backend pair is taken from those. "
                "Use one scheme or the other.")
        endif()
    endforeach()
elseif(NOT TARGET Kokkos::kokkos)
    set(Kokkos_ENABLE_SERIAL   ${ELEMENTS_ENABLE_SERIAL}   CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_OPENMP   ${ELEMENTS_ENABLE_OPENMP}   CACHE BOOL "" FORCE)
    # Kokkos 5 spells the std::threads backend THREADS; PTHREADS no longer exists.
    set(Kokkos_ENABLE_THREADS  ${ELEMENTS_ENABLE_PTHREADS} CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA     ${ELEMENTS_ENABLE_CUDA}     CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_HIP      ${ELEMENTS_ENABLE_HIP}      CACHE BOOL "" FORCE)
endif()

# --- Dependency: MATAR (and Kokkos 5, bundled as a MATAR submodule) ---
# A parent project (e.g. a solver) that already provides matar::matar is
# respected; otherwise it comes from the matar submodule.
if(NOT TARGET matar::matar)
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/matar/CMakeLists.txt OR
       (NOT TARGET Kokkos::kokkos AND NOT EXISTS ${PROJECT_SOURCE_DIR}/matar/src/Kokkos/kokkos/CMakeLists.txt))
        message(FATAL_ERROR
            "The MATAR/Kokkos submodules are missing. Run:\n"
            "    git submodule update --init --recursive")
    endif()

    set(MATAR_ENABLE_KOKKOS  ON  CACHE BOOL "" FORCE)
    # ELEMENTS requires MPI, and its examples use MATAR's MPICArrayKokkos and
    # CommunicationPlan, which only exist in an MPI-enabled MATAR.
    set(MATAR_ENABLE_MPI     ON  CACHE BOOL "" FORCE)
    set(MATAR_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
    set(MATAR_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(MATAR_INSTALL        ON  CACHE BOOL "" FORCE)

    # Forward only the tiers the caller actually set, so an unset tier keeps
    # MATAR's own default of double rather than being pinned here.
    foreach(_tier REAL HIGH_REAL LOW_REAL)
        if(ELEMENTS_${_tier})
            set(MATAR_${_tier} "${ELEMENTS_${_tier}}" CACHE STRING "Set by ELEMENTS_${_tier}" FORCE)
        endif()
    endforeach()
    if(ELEMENTS_REAL OR ELEMENTS_HIGH_REAL OR ELEMENTS_LOW_REAL)
        message(STATUS
            "ELEMENTS precision: real='${ELEMENTS_REAL}' high_real='${ELEMENTS_HIGH_REAL}'"
            " low_real='${ELEMENTS_LOW_REAL}' -> MATAR_REAL / MATAR_HIGH_REAL / MATAR_LOW_REAL")
    endif()

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

    # Silence warnings from Scotch's own (third-party) C sources and CMake
    # files. The Scotch subdirectory captures these values when it is added,
    # so restoring them afterwards keeps ELEMENTS' own flags untouched.
    set(_elements_saved_c_flags "${CMAKE_C_FLAGS}")
    if(MSVC)
        string(APPEND CMAKE_C_FLAGS " /w")
    else()
        string(APPEND CMAKE_C_FLAGS " -w")
    endif()
    set(CMAKE_POLICY_DEFAULT_CMP0177 NEW)  # install() DESTINATION normalization
    FetchContent_MakeAvailable(Scotch)
    set(CMAKE_C_FLAGS "${_elements_saved_c_flags}")
    unset(CMAKE_POLICY_DEFAULT_CMP0177)

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

# --- 9. Optional Parts ---
if(ELEMENTS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()



if(ELEMENTS_BUILD_TESTS AND EXISTS ${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt)
    enable_testing()
    add_subdirectory(tests)
endif()