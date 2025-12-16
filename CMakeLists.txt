cmake_minimum_required(VERSION 3.20)
project(ELEMENTS VERSION 1.0.0 LANGUAGES CXX C)


# --- GLOBAL COMPILER FIXES (GCC 15 Compatibility) ---

# 1. Force C11 Standard Globally
#    This fixes the issue where "int func()" meant "0 arguments" in C23 
#    but "unspecified arguments" in older C.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON) # Enable GNU extensions (gnu11)

# 2. Globally Downgrade GCC 15 Strictness Errors to Warnings
#    Applies only to C files. Fixes "incompatible pointer type" errors in Scotch.
add_compile_options(
    $<$<COMPILE_LANGUAGE:C>:-Wno-error=incompatible-pointer-types>
    $<$<COMPILE_LANGUAGE:C>:-Wno-error=int-conversion>
    $<$<COMPILE_LANGUAGE:C>:-Wno-error=implicit-int>
    # Optional: Silence the warnings completely so they don't spam the log
    $<$<COMPILE_LANGUAGE:C>:-Wno-incompatible-pointer-types>
    $<$<COMPILE_LANGUAGE:C>:-Wno-int-conversion>
    $<$<COMPILE_LANGUAGE:C>:-Wno-implicit-int>
)


# --- 1. Find MPI  ---
find_package(MPI)
if(NOT MPI_FOUND)
    message(WARNING "MPI not found! You may not be able to run parallel examples.")
else()
    if(MPI_CXX_VERSION LESS 3)
        message(WARNING "MPI version ${MPI_CXX_VERSION} found, but version 3 or higher is required!")
    endif()
endif()


# --- Options for Users ---
option(ELEMENTS_BUILD_TESTS "Build unit tests" OFF) ## WARNING< NO TEST IMPLEMENTED YET
option(ELEMENTS_BUILD_EXAMPLES "Build examples" ON) 

option(ELEMENTS_ENABLE_SERIAL "Enable Serial support in Kokkos" ON)
option(ELEMENTS_ENABLE_OPENMP "Enable OpenMP support in Kokkos" ON)
option(ELEMENTS_ENABLE_PTHREADS "Enable PThreads support in Kokkos" OFF)
option(ELEMENTS_ENABLE_CUDA "Enable CUDA support in Kokkos" OFF)
option(ELEMENTS_ENABLE_HIP "Enable HIP support in Kokkos" OFF)

# --- Dependency: Kokkos ---
include(FetchContent)

# Forward options to Kokkos BEFORE FetchContent
# We use ordinary variables first; if Kokkos checks CACHE, we force it.
set(Kokkos_ENABLE_SERIAL   ${ELEMENTS_ENABLE_SERIAL}   CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP   ${ELEMENTS_ENABLE_OPENMP}   CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_PTHREADS ${ELEMENTS_ENABLE_PTHREADS} CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA     ${ELEMENTS_ENABLE_CUDA}     CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_HIP      ${ELEMENTS_ENABLE_HIP}      CACHE BOOL "" FORCE)
FetchContent_Declare(
    Kokkos
    GIT_REPOSITORY https://github.com/kokkos/kokkos
    GIT_TAG        4.3.00 # Pin a specific version for stability
)
FetchContent_MakeAvailable(Kokkos)

# --- Dependency: MATAR ---
# Assuming MATAR is a CMake project. If it's just headers, this works fine if MATAR has a CMakeLists.txt
FetchContent_Declare(
    MATAR
    GIT_REPOSITORY https://github.com/lanl/MATAR # Replace with actual URL
    GIT_TAG        main   # Pin a specific commit/tag
)
FetchContent_MakeAvailable(MATAR)


# Get the directory where MATAR was downloaded
FetchContent_GetProperties(MATAR SOURCE_DIR MATAR_SOURCE_DIR)

# --- MATAR Patching ---
if(TARGET matar)
    if(NOT TARGET matar::matar)
        add_library(matar::matar ALIAS matar)
    endif()

    FetchContent_GetProperties(MATAR SOURCE_DIR MATAR_SOURCE_DIR)
    
    get_target_property(_matar_includes matar INTERFACE_INCLUDE_DIRECTORIES)

    # 1. Fix Include Path (Existing Fix)
    if(NOT "${MATAR_SOURCE_DIR}/src/include" IN_LIST _matar_includes)
        message(STATUS "Patching MATAR target with missing include directory")
        target_include_directories(matar INTERFACE 
            $<BUILD_INTERFACE:${MATAR_SOURCE_DIR}/src/include>
        )
    endif()

    # 2. Fix Definitions (NEW FIX)
    # This adds -DHAVE_KOKKOS to the compile line for anything using matar
    target_compile_definitions(matar INTERFACE HAVE_KOKKOS)

    # 3. Fix MPI Definitions (Highly Recommended)
    # Since you are using MPI, MATAR likely needs HAVE_MPI to enable its MPI types
    if(MPI_FOUND)
        target_compile_definitions(matar INTERFACE HAVE_MPI)
    endif()
endif()


# --- Dependency: Scotch / PT-Scotch ---
# We force SCOTCH_MPI to ON to build PT-Scotch.
set(SCOTCH_MPI ON CACHE BOOL "Enable MPI for Scotch" FORCE)
set(SCOTCH_BUILD_TESTS OFF CACHE BOOL "Disable Scotch Tests" FORCE)

FetchContent_Declare(
    Scotch
    GIT_REPOSITORY https://gitlab.inria.fr/scotch/scotch.git
    GIT_TAG        v7.0.4
    # FIX: Robustly replace CMAKE_SOURCE_DIR with CMAKE_CURRENT_SOURCE_DIR.
    # This prevents shell syntax errors and fixes the "include not found" bug.
    PATCH_COMMAND  sed -i "s/CMAKE_SOURCE_DIR/CMAKE_CURRENT_SOURCE_DIR/g" CMakeLists.txt
)
FetchContent_MakeAvailable(Scotch)

# --- FIX: Compatibility for Scotch with GCC 15 ---
# We apply C11 standard and warning suppressions to ALL Scotch libraries and executables.

# # 1. Libraries
# set(_scotch_targets 
#     scotch ptscotch scotcherr scotcherrexit
# )

# # 2. Sequential Utilities (amk_m2, gmap, etc.)
# list(APPEND _scotch_targets 
#     amk_ccc amk_fft2 amk_grf amk_hy amk_m2 amk_p2 
#     gbase gmap gord gout gpart gscat gtst 
#     mcv mmk_m2 mmk_m3 mord mtst
# )

# # 3. Parallel Utilities (dggath, dgmap, etc.)
# list(APPEND _scotch_targets 
#     dggath dgmap dgord dgpart dgscat dgtst
# )

# foreach(_target ${_scotch_targets})
#     if(TARGET ${_target})
#         # Force C11 Standard (Fixes "too many arguments" errors)
#         set_target_properties(${_target} PROPERTIES C_STANDARD 11)
        
#         # Silence strict GCC 15 warnings
#         target_compile_options(${_target} PRIVATE 
#             $<$<COMPILE_LANGUAGE:C>:-Wno-incompatible-pointer-types>
#             $<$<COMPILE_LANGUAGE:C>:-Wno-int-conversion>
#             $<$<COMPILE_LANGUAGE:C>:-Wno-implicit-int>
#         )
#     endif()
# endforeach()

# --- Add ELEMENTS Source ---
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/swage
    ${CMAKE_CURRENT_SOURCE_DIR}/elements
)
add_subdirectory(swage)

# --- Optional Tests --- ## WARNING: NOTE YET IMPLEMENTED
if(ELEMENTS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


# --- Optional Examples ---
if(ELEMENTS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()