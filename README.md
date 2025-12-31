# ELEMENTS

[![Linux](https://github.com/lanl/ELEMENTS/actions/workflows/Linux.yaml/badge.svg)](https://github.com/lanl/ELEMENTS/actions/workflows/Linux.yaml)
[![MacOS](https://github.com/lanl/ELEMENTS/actions/workflows/MacOS.yaml/badge.svg)](https://github.com/lanl/ELEMENTS/actions/workflows/MacOS.yaml)

## What is ELEMENTS?

The C++ **ELEMENTS** library is a collection of sub-libraries to support implementing a diverse range of numerical methods on low and high-order meshes.  The **ELEMENTS** library can be used for research and development of both continuous and discontinuous finite element methods, as well as, finite volume methods to solve a diverse range of partial differential equations. The **ELEMENTS** library includes the following sub-libraries:  **MATAR** contains the routines to support dense and sparse **mat**rices and **ar**rays, **SLAM** contains the interfaces to **s**olvers, **l**inear **a**lgebra, and **m**athematical routines or external packages (e.g., Trilinos),  **elements** contains the mathematical functions to support a large range of elements types including serendipity elements, **SWAGE** contains the routines and data-structures to support unstructured arbitrary-order 3D meshes that move or remain stationary, and **geometry** combines together **SWAGE** and **elements**.  The **ELEMENTS** libary is designed to support Lagrangian (mesh moves) solid dynamics and mechanics codes, Eulerian (mesh is stationary) fluid dynamics codes, and many other code applications.  

<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/codeStructureELEMENTS.png" width="400">
<p align="center">Fig. Code structure layout
  
<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/TaylorGreenVortex-t0.png" width="400"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/TaylorGreenVortex-tEnd.png" width="400">
<p align="center">Fig. A high-order 3D mesh deforming in the Taylor-Green vortex

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

Key CMake options (defaults in parentheses):
- `CMAKE_BUILD_TYPE` (`RelWithDebInfo` if unset)
- `ELEMENTS_BUILD_DOCS` (`OFF`)
- `ELEMENTS_BUILD_EXAMPLES` (`ON` when ELEMENTS is the top-level project, otherwise `OFF`)
- `ELEMENTS_BUILD_TESTS` (`OFF`)
- `ELEMENTS_ENABLE_SERIAL` (`ON`)
- `ELEMENTS_ENABLE_OPENMP` (`ON`)
- `ELEMENTS_ENABLE_PTHREADS` (`OFF`)
- `ELEMENTS_ENABLE_CUDA` (`OFF`)
- `ELEMENTS_ENABLE_HIP` (`OFF`)
 
## Using ELEMENTS from another CMake project

You can add ELEMENTS to your own CMake build with `FetchContent`:

```cmake
include(FetchContent)

# Configure ELEMENTS before it is fetched (optional overrides)
set(ELEMENTS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ELEMENTS_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(ELEMENTS_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
# Choose backends as needed; defaults are listed above
# set(ELEMENTS_ENABLE_OPENMP ON  CACHE BOOL "" FORCE)
# set(ELEMENTS_ENABLE_CUDA   ON  CACHE BOOL "" FORCE)

FetchContent_Declare(ELEMENTS
  GIT_REPOSITORY https://github.com/lanl/ELEMENTS.git
  GIT_TAG        main            # or a release tag/commit
)
FetchContent_MakeAvailable(ELEMENTS)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE ELEMENTS)
```

`ELEMENTS` is an INTERFACE target that propagates its include directories and required dependencies (Kokkos, MATAR, MPI, Scotch). Set the options before `FetchContent_MakeAvailable` to control which backends build.

To learn more about ELEMENTS and how to get started using it, please see the [ELEMENTS documentation](https://lanl.github.io/ELEMENTS/).


## Examples
ELEMENTS has some small examples. Enable them at configure time (default when ELEMENTS is the top-level project) with `-DELEMENTS_BUILD_EXAMPLES=ON`, then build as shown above. Executables live in your build tree under `examples/<name>/`.

### Mesh decomposition (`examples/decomp_example`)
- Demonstrates building an arbitrary order 3D box (or 2D polar) mesh on rank 0, partitioning with PT-Scotch, and exchanging element/node data across ranks using Swage communication plans.
- Build target: `mesh_decomp`.
- Run (from the build directory): `mpirun -n <num_ranks> examples/decomp_example/mesh_decomp`.
- Dependencies: MPI (required), Kokkos and MATAR (transitively provided), PT-Scotch for parallel partitioning.

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
