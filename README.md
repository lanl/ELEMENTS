# ELEMENTS

The C++ **ELEMENTS** library is a collection of libraries to support implementing a diverse range of numerical methods on low and high-order meshes.  The **ELEMENTS** library can be used for research and development of both continuous and discontinuous finite element methods, as well as, finite volume methods to solve a diverse range of partial differential equations. The **ELEMENTS** library includes the following libraries:  **MATAR** contains the routines to support dense and sparse **mat**rices and **ar**rays, **SLAM** contains the interfaces to **s**olvers, **l**inear **a**lgebra, and **m**athematical routines or external packages (e.g., Trilinos),  **elements** contains the mathematical functions to support a large range of elements types including serendiptiy elements, **SWAGE** contains the routines and data-structures to support unstructured arbitrary-order 3D meshes that move or remain stationary, and **geometry** combines **SWAGE** with **elements**.  The **ELEMENTS** libary is designed to support Lagrangian (mesh moves) solid dynamics and mechanics codes, Eulerian (mesh is stationary) fluid dynamics codes, and many other code applications.  

<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/master/examples/figures/codeStructureELEMENTS.png" width="400">
  
## MATAR
The **MATAR** library is designed for dense and sparse data representations, and follows the data-oriented programming approach for highly efficient calculations.  **MATAR** leverages **Kokkos** for performance portability over diverse architectures. The data representations developed in **MATAR** and the numerical tools in **ELEMENTS** are designed for performance, portability, and productivity (i.e., ease of use).  **MATAR** is stored in a separate repository as it can aid many applications, as such, it is included as a GitHub submodule in the **ELEMENTS** repository.


## SLAM
The **SLAM** libary is a soon to be released capability to seemlessly allow users to solve linear systems (dense or sparse), to perform linear algebra opperations, or to apply mathematical opperators to data stored in **MATAR**.  A prototype of **SLAM** exists and is being developed.


## elements
A mesh is composed of non-overlapping elements, where the **elements** library contains the mathematical functions to support a very broad range of element types including: 

* linear, quadratic, and cubic serendipity elements in 2D and 3D; 
* arbitrary-order tensor product elements in 2D and 3D;
* arbitrary-order spectral elements; and 
* a linear 4D element. 

The **elements** library has functions to calculate quantities that are commonly used in finite element methods (both continuous and discontinous) such as a basis function, gradient of a basis function, the Jacobian matrix, the inverse Jacobian matrix, the determinant of the Jacobiam matrix, and a physical position inside the element, to name a few examples. The **elements** library also supports both Gauss-Legendre and Gauss-Lobatto quadrature rules up to 8 quadrature points in each coordinate direction. 


## SWAGE
The **SWAGE** library contains a large suite of mesh data structures, a rich set of index spaces for different geometric entities, and many connectivity arrays between the various index spaces.  This library supports both unstructured linear and arbitrary-order meshes.  **SWAGE** is designed to support diverse range of methods that arise in computational physics and engineering.

## geometry
The **geometry** library combines **SWAGE** with **elements** to deliver the required capabilities to implement a large range of numerical methods on unstructured linear or high-order meshes.

## examples
The examples folder contains a simple code to calculate an average from the cells to nodes and then back.  A figures folder within the examples folder contains diagrams to illustrate, for instance, the code project layout, geometric index spaces, and various capabilties supported by the **ELEMENTS** library.  


## Building & Installation

To build **ELEMENTS** in place, simply enter
```
cmake ..; make
```
at the command line.

To build **ELEMENTS** in a separate build directory, create the build directory somewhere in your filesystem.
```
mkdir build
```
Then create a configuration script like, for example,
```
#!/bin/bash
ELEMENTS_DIR=..  # relative path of ELEMENTS repository
cmake \
  -DCMAKE_INSTALL_PREFIX=`pwd` \
  -DCMAKE_CXX_FLAGS=-g \
  ${ELEMENTS_DIR}
```
where the installation directory is configured to be the directory in which the configuration script is run.
Then enter
```
./my_config.sh; make
```
at the command line (assuming you named your configuration script `my_config.sh`).

To install the **ELEMENTS** libraries and header files (with an in-place build or otherwise), enter
```
make install
```
at the command line.
This will copy the **ELEMENTS** libraries to the directory specified by the `CMAKE_INSTALL_PREFIX` variable.
In particular, the libraries will be copied to a `lib/` subdirectory and the header files will be copied to a `include/` subdirectory.
If the install prefix is not specified in a configuration file, CMake will use the default, which on Linux is `/usr/local/` and the copy will require administrative privileges.

Warning: if linking against the `libelements` library, you must define a global variable `elem` that specifies the class of element.

## Citation
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
```
