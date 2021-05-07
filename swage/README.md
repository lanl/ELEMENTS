The **sw**ift L**ag**rangian to **E**ulerian (SWAGE) library is quit general and allows users to implement a large range of methods on unstructured arbitrary-order 3D meshes.  This library supports connectivity data structures and index spaces needed for implementing low-order and high-order methods.  SWIFT is designed to work with the elements library that contains e.g., quadrature sets, basis functions, Jacobian matrices, etc.  The SWAGE library has unique index spaces to support arbitrary-order Lagrangian methods that solve the governing equations for motion on an arbitrary-order mesh that moves with the material.  The SWAGE library also allows code developers to create high-order Eulerian fluid-dynamic codes using high-order meshes that are conformal to a curved boundary (e.g., a wing). The index spaces in SWIFT are:

* elements
* vertices (kinematic degrees for freedom of an element)
* nodes (a continuous space that coincides with the Lobatto point locations)
* gauss points (a discontinuous space that are the Lobatto point locations)
* sub-cells (a decomposition of an element using the nodes)
* sub-zones (a decomposition of an element using the vertices)
* surface (the surface of the element)
* patch (a portion of the surface that is the surface of a sub-cell)
* facet (a portion of a patch)
* corner (a corner of a sub-cell)

Connectivity data structures exist to map from from an index to another index space (e.g., all nodes in an element) and to walk over neighboring mesh entities (e.g., all sub-cells around a sub-cell).
