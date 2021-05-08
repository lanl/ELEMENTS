# ELEMENTS

ELEMENTS is a library of mathematical functions for supporting a very broad range of element types including: 

* linear, quadratic, and cubic serendipity elements in 2D and 3D; 
* high-order spectral elements; and 
* a linear 4D element. 

The ELEMENTS library can be used for research and development of both continuous and discontinuous finite element methods for solving a diverse range of partial differential equations. 
The library has functions for calculating quantities that are commonly used in finite element methods such as the gradient of a basis function, the Jacobi matrix, the inverse Jacobi matrix, the determinant of the Jacobi matrix, and a physical position inside the element, to name a few examples. 
The library also supports both Gauss-Legendre and Gauss-Lobatto quadrature rules up to 8 quadrature points in each coordinate direction. 
The examples and discussions in this paper will focus on Lagrangian solid mechanics and dynamics, but ELEMENTS can be used for many other applications.  

On top of supporting multiple types of discretization, ELEMENTS contains a rather extensive library for creating mesh structures to map between multiple different index spaces that are commonly required for solving problems in computational physics.

The example code in the examples folder does a simple average from cells to nodes and then back.  

## Building & Installation

To build ELEMENTS in place, simply enter
```
cmake ..; make
```
at the command line.

To build ELEMENTS in a separate build directory, create the build directory somewhere in your filesystem.
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

To install the ELEMENTS libraries and header files (with an in-place build or otherwise), enter
```
make install
```
at the command line.
This will copy the ELEMENTS libraries to the directory specified by the `CMAKE_INSTALL_PREFIX` variable.
In particular, the libraries will be copied to a `lib/` subdirectory and the header files will be copied to a `include/` subdirectory.
If the install prefix is not specified in a configuration file, CMake will use the default, which on Linux is `/usr/local/` and the copy will require administrative privileges.

Warning: if linking against the `libelements` library, you must define a global variable `elem` that specifies the class of element.

# Citation
```
@article{MOORE2019100257,
title = {ELEMENTS: A high-order finite element library in C++},
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
