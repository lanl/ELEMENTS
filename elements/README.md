# elements

A mesh is a collection of non-overlapping elements, and the mesh can be unstructured.  Each reference element is defined in a reference coordinate system with a spatial mapping (e.g., an interpolation polynomial) to the physical coordinates.  An element in the physical coordinates can have edges that are defined by straight lines (i.e., linear) or have edges defined by a high-order polynomial. The **elements** sub-library contains the mathematical functions for a suite of element types.    The **geometry** sub-library combines the mesh datastructures and mesh index spaces in **SWAGE** with the reference element defined in **elements**. 


## geometry definition
The position inside an element is given by a polynomial.  **elements** contains a suite of spatial polynomials and also a space-time polynomial.  The Jacobi matrix is equal to the gradient of this spatial polynomial.  The volume of an element can be calculated using the determinate of the Jacobi matrix.  The inverse of the Jacobi matrix is used to transform a nabla opperator from the physical coordinates to the reference coordinates.


## index conventions
The reference element contains nodes that are collocated with the Labatto quadrature points.  The local indexing of the quadrature points and nodes within an element follows an (i,j,k) index convenction. An index map is supplied to convert the serendipty local index convention to the local (i,j,k) index convention.
