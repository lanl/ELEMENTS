# SWAGE
The **sw**ift L**ag**rangian to **E**ulerian (SWAGE) library is quit general and allows users to implement a large range of methods on unstructured arbitrary-order 3D meshes.  This library supports connectivity data structures and index spaces needed to implement low-order or high-order numerical methods.  SWAGE is designed to work with the elements library that contains e.g., quadrature sets, basis functions, Jacobian matrices, etc.  The SWAGE library has unique index spaces to support arbitrary-order Lagrangian material dynamic methods that solve the governing equations for motion on an arbitrary-order mesh that moves with the material.  The SWAGE library also allows code developers to create high-order Eulerian fluid-dynamic codes using high-order meshes that are conformal to a curved boundary (e.g., a wing). 

## Descriptions
The index spaces in SWAGE are:

* elements (the computational mesh is decomposed into non-overlapping elements)
* vertices (kinematic degrees of freedom for an element)
* nodes (a continuous space that coincides with the Lobatto point locations)
* gauss points (a discontinuous space that are the Lobatto point locations)
* sub-cells (a decomposition of an element using the nodes), abbreviated as cells
* sub-zones (a decomposition of an element using the vertices), abbreviated as zones
* surface (a surface of the element)
* patch (a portion of the element surface that coincides with the surface of a sub-cell)
* facet (a portion of a patch)
* corner (a corner of a sub-cell)

Connectivity data structures exist to map from from an index to another index space (e.g., all nodes in an element) and to walk over neighboring mesh entities (e.g., all sub-cells around a sub-cell).  The SWAGE library is stitched together with the elements library in the geometry library.  

### Index naming conventions
The index spaces for the mesh are denoted with a _gid_ and the index spaces for local mesh entities are denoted with a _lid_.  For example, walking over all the elements in the mesh and then over all the sub-cells in the element would be, 
```
for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
   for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){ 
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);  // access the cell mesh index 
      // ...
   }
}
```

Then, to walk over all the corners in a sub-cell or all the nodes in a sub-cell would be, 

```
for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
   for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
      
      for(int node_lid = 0; node_lid < 8; node_lid++){
         int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);  // access the node mesh index
         int corner_gid = mesh.corners_in_cell(cell_gid, node_lid);  // access the corner mesh index
         // ...
      }
      
   }   
}
```

Then, to walk over all the surface facets in a corner would be,

```
for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
   for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
      
      for(int node_lid = 0; node_lid < 8; node_lid++){
         int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);  // access the node mesh index
         int corner_gid = mesh.corners_in_cell(cell_gid, node_lid);  // access the corner mesh index
         
         for(int facet_lid = 0; facet_lid < 3; facet_lid++){
            // ...
         }
         
      }
      
   }   
}
```

There are many ways to walk over or access an index space. To just walk over all the nodes in the mesh, the coding is,

```
for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
   // ...
} // end for loop over nodes
```

Then, to walk over all the corners in a node would be,

```
for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
   for(int corn_lid = 0; corn_lid < mesh.num_corners_in_node(node_gid); corn_lid++){
      int corner_gid = mesh.corners_in_node(node_gid, corn_lid);
      // ...
   }   
} // end for loop over nodes
```

The index spaces in a reference element, which comes from the elements library and are not in SWAGE, is denoted with a _rid_ and the indices local that index are denoted with a _rlid_. 


