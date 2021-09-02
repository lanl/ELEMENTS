#include "swage.h"

namespace swage {
  void evaluate_jacobian_determinants(mesh_t &mesh) {
    // Mesh parameters
    int num_dim      = 3;
    int num_elems    = mesh.num_elems();
    int elem_order   = mesh.elem_order();
    int num_verts_1d = elem_order + 1;
    int num_verts    = pow(num_verts_1d, num_dim);
    int num_nodes_1d = 2*elem_order + 1;
    int num_nodes    = pow(num_nodes_1d, num_dim);

    // Create vertex to node map
    int vert_id = 0;
    const SizeType num_rad = 3;
    SizeType radices[num_rad] = {SizeType(num_nodes_1d), 
        SizeType(num_nodes_1d), SizeType(num_nodes_1d)};
    CArray<int> vert_node_map(num_verts);
    for (int k = 0; k < num_nodes_1d; k += 2) {
      for (int j = 0; j < num_nodes_1d; j += 2) {
        for (int i = 0; i < num_nodes_1d; i += 2) {
          SizeType IJK[num_rad] = {SizeType(i), SizeType(j), SizeType(k)};
          vert_node_map(vert_id) = int(common::mixed_radix_to_base_10(num_rad, 
                radices, IJK));
          vert_id++;
        }   
      }
    }

    // Assume Lobatto nodes for quadrature
    CArray<NumType> lobatto_nodes(num_nodes_1d);
    lobatto_nodes_1D_tmp(lobatto_nodes, num_nodes_1d);

    // Extract vertex coordinates from quadrature points and evaluate
    // barycentric weights of vertices
    vert_id = 0;
    CArray<NumType> vert_coords_1d(num_verts_1d);
    for (int node_id = 0; node_id < num_nodes_1d; node_id += 2) {
      vert_coords_1d(vert_id) = lobatto_nodes(node_id);
      vert_id++;
    }

    CArray<NumType> vert_weights(num_nodes_1d);
    lagrange::compute_barycentric_weights(num_verts_1d, 
        vert_coords_1d.get_pointer(), vert_weights.get_pointer());

    // Loop over elements and...
    CArray<NumType> vert_x_coords(num_verts);
    CArray<NumType> vert_y_coords(num_verts);
    CArray<NumType> vert_z_coords(num_verts);
    CArray<NumType> workspace(3*num_verts_1d+1);
    for (int elem_id = 0; elem_id < num_elems; elem_id++) {
      // ...get spatial coordinates of vertices from mesh
      for (int vert_id = 0; vert_id < num_verts; vert_id++) {
        int node_id = vert_node_map(vert_id);
        vert_x_coords(vert_id) = mesh.node_coords(mesh.nodes_in_elem(
            elem_id, node_id), 0);
        vert_y_coords(vert_id) = mesh.node_coords(mesh.nodes_in_elem(
              elem_id, node_id), 1);
        vert_z_coords(vert_id) = mesh.node_coords(mesh.nodes_in_elem(
              elem_id, node_id), 2);
      }

      // Loop over quadrature points...
      for (int node_lid = 0; node_lid < num_nodes; node_lid++) {
        // ...and compute Jacobian determinant at each quadrature point

        // Get IJK coordinates of quadrature point
        const SizeType num_rad = 3;
        SizeType IJK[num_rad];
        common::base_10_to_mixed_radix(num_rad, radices, node_lid, IJK);

        // Get reference coordinates of quadrature point
        NumType node_coords[3];
        node_coords[0] = lobatto_nodes(IJK[0]);
        node_coords[1] = lobatto_nodes(IJK[1]);
        node_coords[2] = lobatto_nodes(IJK[2]);

        int node_gid = mesh.nodes_in_elem(elem_id, node_lid);

        SizeType num_verts_1d_unsigned = SizeType(num_verts_1d);
        lagrange::compute_jacobian_determinant(num_verts_1d_unsigned, 
            vert_coords_1d.get_pointer(), vert_weights.get_pointer(), 
            vert_x_coords.get_pointer(), vert_y_coords.get_pointer(), 
            vert_z_coords.get_pointer(), workspace.get_pointer(), node_coords, 
            mesh.gauss_pt_det_j(node_gid));
      }
    }
  }
}
