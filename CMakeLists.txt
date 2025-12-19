cmake_minimum_required(VERSION 3.20)
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
if(MPI_CXX_VERSION LESS 3)
    message(WARNING "ELEMENTS requires MPI version 3+. Found ${MPI_CXX_VERSION}.")
endif()

# --- Options ---
option(ELEMENTS_ENABLE_SERIAL "Enable Serial support in Kokkos" ON)
option(ELEMENTS_ENABLE_OPENMP "Enable OpenMP support in Kokkos" ON)
option(ELEMENTS_ENABLE_PTHREADS "Enable PThreads support in Kokkos" OFF)
option(ELEMENTS_ENABLE_CUDA "Enable CUDA support in Kokkos" OFF)
option(ELEMENTS_ENABLE_HIP "Enable HIP support in Kokkos" OFF)

option(ELEMENTS_BUILD_EXAMPLES "Build examples" ${ELEMENTS_IS_TOP_LEVEL}) 
option(ELEMENTS_BUILD_TESTS "Build unit tests" OFF)

include(FetchContent)

# --- Dependency: Kokkos ---
if(NOT TARGET Kokkos::kokkos)
    set(Kokkos_ENABLE_SERIAL   ${ELEMENTS_ENABLE_SERIAL}   CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_OPENMP   ${ELEMENTS_ENABLE_OPENMP}   CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_PTHREADS ${ELEMENTS_ENABLE_PTHREADS} CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA     ${ELEMENTS_ENABLE_CUDA}     CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_HIP      ${ELEMENTS_ENABLE_HIP}      CACHE BOOL "" FORCE)
    
    FetchContent_Declare(
        Kokkos
        GIT_REPOSITORY https://github.com/kokkos/kokkos
        GIT_TAG        4.3.00
    )
    FetchContent_MakeAvailable(Kokkos)
endif()

# --- Dependency: MATAR ---
if(NOT TARGET matar::matar)
    FetchContent_Declare(
        MATAR
        GIT_REPOSITORY https://github.com/lanl/MATAR
        GIT_TAG        main
    )
    FetchContent_MakeAvailable(MATAR)

    if(TARGET matar)
        if(NOT TARGET matar::matar)
            add_library(matar::matar ALIAS matar)
        endif()
        
        FetchContent_GetProperties(MATAR SOURCE_DIR MATAR_SOURCE_DIR)
        get_target_property(_matar_includes matar INTERFACE_INCLUDE_DIRECTORIES)
        
        if(NOT "${MATAR_SOURCE_DIR}/src/include" IN_LIST _matar_includes)
            target_include_directories(matar INTERFACE 
                $<BUILD_INTERFACE:${MATAR_SOURCE_DIR}/src/include>
            )
        endif()
        
        target_compile_definitions(matar INTERFACE HAVE_KOKKOS)
        if(MPI_FOUND)
            target_compile_definitions(matar INTERFACE HAVE_MPI)
        endif()
    endif()
endif()

# --- Dependency: Scotch / PT-Scotch ---
if(NOT TARGET ptscotch)
    set(SCOTCH_MPI ON CACHE BOOL "Enable MPI for Scotch" FORCE)
    set(SCOTCH_BUILD_TESTS OFF CACHE BOOL "Disable Scotch Tests" FORCE)

    FetchContent_Declare(
        Scotch
        GIT_REPOSITORY https://gitlab.inria.fr/scotch/scotch.git
        GIT_TAG        v7.0.4
        PATCH_COMMAND  ${CMAKE_COMMAND} -DFileToPatch=CMakeLists.txt -P "${PROJECT_SOURCE_DIR}/cmake/PatchScotch.cmake"
    )
    FetchContent_MakeAvailable(Scotch)

    # FIX: GCC 15 Warning Suppressions (Target-Specific)
    # Even though C11 is set globally above, we still need to silence the
    # specific GCC 15 warnings for Scotch that are not fixed by C11 alone.
    set(_scotch_targets 
        scotch ptscotch scotcherr scotcherrexit
        amk_ccc amk_fft2 amk_grf amk_hy amk_m2 amk_p2 
        gbase gmap gord gout gpart gscat gtst 
        mcv mmk_m2 mmk_m3 mord mtst
        dggath dgmap dgord dgpart dgscat dgtst
    )

    if(NOT MSVC)
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



if(ELEMENTS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()