.. _elements:

elements
========
A mesh is composed of non-overlapping elements, where the **elements** sub-library contains the mathematical functions to support a very broad range of element types including: 

* linear, quadratic, and cubic serendipity elements in 2D and 3D; 
* arbitrary-order tensor product elements in 2D and 3D;
* arbitrary-order spectral elements; and 
* a linear 4D element. 

The **elements** sub-library has functions to calculate quantities that are commonly used in finite element methods (both continuous and discontinous) such as a basis function, gradient of a basis function, the Jacobian matrix, the inverse Jacobian matrix, the determinant of the Jacobian matrix, and a physical position inside the element, to name a few examples. 
The **elements** sub-library also supports both Gauss-Legendre and Gauss-Lobatto quadrature rules up to 8 quadrature points in each coordinate direction. 

.. doxygennamespace:: elements
   :members:
