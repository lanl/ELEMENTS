# ELEMENTS

[![Linux](https://github.com/lanl/ELEMENTS/actions/workflows/Linux.yaml/badge.svg)](https://github.com/lanl/ELEMENTS/actions/workflows/Linux.yaml)
[![MacOS](https://github.com/lanl/ELEMENTS/actions/workflows/MacOS.yaml/badge.svg)](https://github.com/lanl/ELEMENTS/actions/workflows/MacOS.yaml)

## What is ELEMENTS?

The C++ **ELEMENTS** library (LANL open-source code number C20058) is a collection of sub-libraries to support implementing a diverse range of numerical methods on low and high-order meshes.  The **ELEMENTS** library can be used for research and development of both continuous and discontinuous finite element methods, as well as, finite volume methods to solve a diverse range of partial differential equations. The **ELEMENTS** library includes the following sub-libraries:  **MATAR** contains the routines to support dense and sparse **mat**rices and **ar**rays, **SLAM** contains the interfaces to **s**olvers, **l**inear **a**lgebra, and **m**athematical routines or external packages (e.g., Trilinos),  **elements** contains the mathematical functions to support a large range of elements types including serendipity elements, **SWAGE** contains the routines and data-structures to support unstructured arbitrary-order 3D meshes that move or remain stationary, and **geometry** combines together **SWAGE** and **elements**.  The **ELEMENTS** libary is designed to support Lagrangian (mesh moves) solid dynamics and mechanics codes, Eulerian (mesh is stationary) fluid dynamics codes, and many other code applications.  

<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/codeStructureELEMENTS.png" width="400">
<p align="center">Fig. Code structure layout
  
<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/TaylorGreenVortex-t0.png" width="400"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/TaylorGreenVortex-tEnd.png" width="400">
<p align="center">Fig. A high-order 3D mesh deforming in the Taylor-Green vortex


# Cloning the code
ELEMENTS carries MATAR as a git submodule, and MATAR in turn carries Kokkos 5.2.1 as its own
nested submodule, so a recursive clone is required.

If the user has set up ssh keys with GitHub, type
```
git clone --recursive ssh://git@github.com/lanl/ELEMENTS.git
```
The code can also be cloned using
```
git clone --recursive https://github.com/lanl/ELEMENTS.git
```
If you already have a non-recursive clone, populate the submodules with:
```
git submodule update --init --recursive
```

## Requirements
- CMake >= 3.22 (>= 3.25.2 if building a CUDA backend, for C++20 support in the CUDA language)
- A C++20 compiler (required transitively by Kokkos 5, via MATAR)
- MPI (required)

## Getting started

To build the examples locally:

1. Configure into a separate build tree:
```
mkdir build
cd build
cmake ..

2. Build in parallel using <num_cores>
make -j<num_cores>
```

Or, using the provided presets, which also select a Kokkos backend:
```
cmake --preset serial
cmake --build --preset serial
```
Other presets: `openmp`, `pthreads`, `cuda`, `hip`, and their `-debug` variants, plus
`cuda-hostomp` (CUDA device backend with a concurrent OpenMP host backend). See
`CMakePresets.json`.

Key CMake options (defaults in parentheses):
- `CMAKE_BUILD_TYPE` (`RelWithDebInfo` if unset)
- `ELEMENTS_BUILD_DOCS` (`OFF`)
- `ELEMENTS_BUILD_EXAMPLES` (`ON` when ELEMENTS is the top-level project, otherwise `OFF`)
- `ELEMENTS_BUILD_TESTS` (`OFF`)
- `ELEMENTS_HOST_BACKEND` (unset): host-side Kokkos backend for the `_HOST` macros — `serial`, `openmp`, or `pthreads`
- `ELEMENTS_DEVICE_BACKEND` (unset, defaults to `serial` if nothing else is specified): device-side backend — `serial`, `openmp`, `pthreads`, `cuda`, `hip`, or `sycl`
- `ELEMENTS_ENABLE_GPU_AWARE_MPI` (`OFF`)
- `ELEMENTS_REAL` / `ELEMENTS_HIGH_REAL` / `ELEMENTS_LOW_REAL` (`double`): precision tiers, forwarded to MATAR's `real_t` / `high_real_t` / `low_real_t`
- `ELEMENTS_INSTALL` (`ON` when ELEMENTS is the top-level project, otherwise `OFF`)

The older `ELEMENTS_ENABLE_SERIAL` / `ELEMENTS_ENABLE_OPENMP` / `ELEMENTS_ENABLE_PTHREADS` /
`ELEMENTS_ENABLE_CUDA` / `ELEMENTS_ENABLE_HIP` booleans still work and map onto
`ELEMENTS_DEVICE_BACKEND`, but emit a deprecation warning; prefer the backend options above.

## Using ELEMENTS from another CMake project

### Option 1: git submodule (recommended)

Vendor ELEMENTS into your repository as a submodule (which pulls in MATAR and Kokkos as nested
submodules):
```
git submodule add https://github.com/lanl/ELEMENTS.git external/ELEMENTS
git submodule update --init --recursive
```

```cmake
cmake_minimum_required(VERSION 3.22)
project(myapp LANGUAGES CXX)

# Set backend/build options BEFORE ELEMENTS is added.
set(ELEMENTS_DEVICE_BACKEND "openmp" CACHE STRING "" FORCE)
set(ELEMENTS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ELEMENTS_BUILD_TESTS    OFF CACHE BOOL "" FORCE)

add_subdirectory(external/ELEMENTS)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE ELEMENTS)
```

### Option 2: FetchContent

```cmake
include(FetchContent)

# Configure ELEMENTS before it is fetched (optional overrides)
set(ELEMENTS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ELEMENTS_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(ELEMENTS_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
# Choose a backend as needed; see the option table above
set(ELEMENTS_DEVICE_BACKEND "openmp" CACHE STRING "" FORCE)

FetchContent_Declare(ELEMENTS
  GIT_REPOSITORY https://github.com/lanl/ELEMENTS.git
  GIT_TAG        main            # or a release tag/commit
)
FetchContent_MakeAvailable(ELEMENTS)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE ELEMENTS)
```
`FetchContent` updates submodules recursively by default, so MATAR and its bundled Kokkos are
fetched and built automatically — nothing extra is needed.

`ELEMENTS` is an INTERFACE target that propagates its include directories and required dependencies
(Kokkos and MATAR via `matar::matar`, MPI, and PT-Scotch in the build tree). Set the options before
`add_subdirectory` / `FetchContent_MakeAvailable` to control which backend builds.

Note: PT-Scotch is fetched at configure time and is a build-tree-only dependency — it is not
re-exported by ELEMENTS' own install/export rules (`ELEMENTS_INSTALL`), so a consumer of an
*installed* ELEMENTS must provide PT-Scotch itself if it needs `decomp_utilities`.

To learn more about ELEMENTS and how to get started using it, please see the [ELEMENTS documentation](https://lanl.github.io/ELEMENTS/).


## Examples
ELEMENTS has some small examples. Enable them at configure time (default when ELEMENTS is the top-level project) with `-DELEMENTS_BUILD_EXAMPLES=ON`, then build as shown above. Executables live in your build tree under `examples/<name>/`.

### Mesh decomposition (`examples/decomp_example`)
- Demonstrates building an arbitrary order 3D box (or 2D polar) mesh on rank 0, partitioning with PT-Scotch, and exchanging element/node data across ranks using Swage communication plans.
- Build target: `mesh_decomp`.
- Run (from the build directory): `mpirun -n <num_ranks> examples/decomp_example/mesh_decomp`.
- Dependencies: MPI (required), Kokkos and MATAR (transitively provided), PT-Scotch for parallel partitioning.

## Code formatting

ELEMENTS uses MATAR's own clang-format style and macro-aware post-processor directly from the
`matar/` submodule (no duplicated config). Format ELEMENTS' own source with:
```
formatting/format.sh            # format src/, examples/, tests/ in place
formatting/format.sh --check    # report non-conforming files without modifying them
```
See [`formatting/README.md`](formatting/README.md) for details.

## How to cite

If you use the ELEMENTS library in your work, please cite the following in any pursuant research papers.

```
@article{MOORE2019100257,
  title = "{ELEMENTS: A high-order finite element library in C++}",
  journal = {SoftwareX},
  volume = {10},
  pages = {100257},
  year = {2019},
  issn = {2352-7110},
  doi = {https://doi.org/10.1016/j.softx.2019.100257},
  url = {https://www.sciencedirect.com/science/article/pii/S235271101930113X},
  author = {Jacob L. Moore and Nathaniel R. Morgan and Mark F. Horstemeyer},
  keywords = {Element Library, C++, High-order elements, Spectral elements, Serendipity elements}
}
```

```
@proceedings{10.1115/DETC2022-89562,
    title = {SWAGE: A 3D Arbitrary-Order Element Mesh Library to Support Diverse Numerical Methods},
    volume = {Volume 2: 42nd Computers and Information in Engineering Conference (CIE)},
    series = {International Design Engineering Technical Conferences and Computers and Information in Engineering Conference},
    pages = {V002T02A029},
    year = {2022},
    month = {08},
    doi = {10.1115/DETC2022-89562},
    url = {https://doi.org/10.1115/DETC2022-89562},
    author = {Morgan, Nathaniel R. and Moore, Jacob and Kiviaho, Jan and Diaz, Adrian},
    eprint = {https://asmedigitalcollection.asme.org/IDETC-CIE/proceedings-pdf/IDETC-CIE2022/86212/V002T02A029/6942954/v002t02a029-detc2022-89562.pdf},
}
```
