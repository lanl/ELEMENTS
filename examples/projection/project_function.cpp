#include "legendre_polynomials.h"
#include "elements.h"
#include "matar_blas_lapack_interface.h"
#include "vtk_io.h"
#include "error.h"
#include "common.h"

#include <functional>

void project_function_legendre(
    std::function<NumType(RealNumber, RealNumber, RealNumber)> func, 
    SwageMesh &mesh, CArray<NumType> &projection) {
  // Mesh parameters
  const int num_dim      = 3;
  int num_elems    = mesh.num_elems();
  int elem_order   = mesh.elem_order();
  int num_verts_1d = elem_order + 1;

  // Basis function parameters
  int num_basis = pow(num_verts_1d, num_dim);
  const SizeType num_rad = SizeType(num_dim);
  SizeType vert_radices[num_rad] = {SizeType(num_verts_1d), 
    SizeType(num_verts_1d), SizeType(num_verts_1d)};

  // Quadrature points and weights, assuming Lobatto quadrature and assuming
  // nodes are quadrature points
  int num_nodes_1d = 2*elem_order + 1;
  CArray<RealNumber> lobatto_points(num_nodes_1d);
  CArray<RealNumber> lobatto_weights(num_nodes_1d);

  elements::lobatto_nodes_1D(lobatto_points, num_nodes_1d);
  elements::lobatto_weights_1D(lobatto_weights, num_nodes_1d);

  SizeType node_radices[num_rad] = {SizeType(num_nodes_1d), 
    SizeType(num_nodes_1d), SizeType(num_nodes_1d)};

  // Allocate memory for mass matrix, right-hand side, and solution
  CArray<RealNumber> m(num_basis, num_basis);
  CArray<RealNumber> minv(num_basis, num_basis);
  CArray<RealNumber> f(num_basis);
  CArray<RealNumber> u(num_basis);

  // Loop over elements
  for (int elem_id = 0; elem_id < num_elems; elem_id++) {
    for (SizeType j = 0; j < num_basis; j++) {
      // Decompose index of 3D tensor product basis function into indices of 1D
      // basis functions
      SizeType rst[num_rad];
      common::base_10_to_mixed_radix(num_rad, vert_radices, j, rst);

      // Build mass matrix
      for (SizeType i = 0; i < num_basis; i++) {
        // Decompose index of 3D tensor product basis function into indices of 1D
        // basis functions
        SizeType uvw[num_rad];
        common::base_10_to_mixed_radix(num_rad, vert_radices, i, uvw);

        m(i,j) = 0.0;

        int num_nodes_in_elem = mesh.num_nodes_in_elem();
        for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
          // Extract reference coordinates of quadrature point
          SizeType lmn[num_rad];
          common::base_10_to_mixed_radix(num_rad, node_radices, 
              SizeType(node_lid), lmn);
          NumType X = lobatto_points(lmn[0]);
          NumType Y = lobatto_points(lmn[1]);
          NumType Z = lobatto_points(lmn[2]);

          // Evaluate 1D basis functions
          NumType psi_r = legendre::evaluate(rst[0], X);
          NumType psi_s = legendre::evaluate(rst[1], Y);
          NumType psi_t = legendre::evaluate(rst[2], Z);

          NumType psi_u = legendre::evaluate(uvw[0], X);
          NumType psi_v = legendre::evaluate(uvw[1], Y);
          NumType psi_w = legendre::evaluate(uvw[2], Z);

          // Evaluate 3D basis functions
          NumType phi_i = psi_u*psi_v*psi_w;
          NumType phi_j = psi_r*psi_s*psi_t;

          //std::cout << rst[0] << " " << rst[1] << " " << rst[2] << std::endl;
          //std::cout << uvw[0] << " " << uvw[1] << " " << uvw[2] << std::endl;
          //std::cout << lmn[0] << " " << lmn[1] << " " << lmn[2] << std::endl;
          //std::cout << X << " " << Y << " " << Z << std::endl;
          //std::cout << phi_i << " " << phi_j << std::endl;
          //std::cout << std::endl;

          // Evaluate quadrature weight times Jacobian determinant
          NumType wl = lobatto_weights(lmn[0]);
          NumType wm = lobatto_weights(lmn[1]);
          NumType wn = lobatto_weights(lmn[2]);
          int node_gid = mesh.nodes_in_elem(elem_id, node_lid);
          NumType wxj = wl*wm*wn*mesh.gauss_pt_det_j(node_gid); 

          // Accumulate contribution of quadrature point in entry of mass
          // matrix
          m(i,j) += phi_i*phi_j*wxj;
        }
      }

      // Build right-hand side
      f(j) = 0.0;
      int num_nodes_in_elem = mesh.num_nodes_in_elem();
      for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
        // Extract reference coordinates of quadrature point
        SizeType lmn[num_rad];
        common::base_10_to_mixed_radix(num_rad, node_radices, 
            SizeType(node_lid), lmn);
        NumType X = lobatto_points(lmn[0]);
        NumType Y = lobatto_points(lmn[1]);
        NumType Z = lobatto_points(lmn[2]);

        // Extract spatial coordinates of quadrature point
        int node_gid = mesh.nodes_in_elem(elem_id, node_lid);
        NumType x = mesh.node_coords(node_gid, 0);
        NumType y = mesh.node_coords(node_gid, 1);
        NumType z = mesh.node_coords(node_gid, 2);

        // Evaluate 1D basis functions
        NumType psi_r = legendre::evaluate(rst[0], X);
        NumType psi_s = legendre::evaluate(rst[1], Y);
        NumType psi_t = legendre::evaluate(rst[2], Z);

        // Evaluate 3D basis functions
        NumType phi_j = psi_r*psi_s*psi_t;

        // Evaluate quadrature weight times Jacobian determinant
        NumType wl = lobatto_weights(lmn[0]);
        NumType wm = lobatto_weights(lmn[1]);
        NumType wn = lobatto_weights(lmn[2]);
        NumType wxj = wl*wm*wn*mesh.gauss_pt_det_j(node_gid); 

        // Accumulate contribution of quadrature point in component of
        // right-hand side
        f(j) += func(x, y, z)*phi_j*wxj;
      }
    }

    // Solve for moments over element
    matar2lapack::invert(m, minv);
    matar2blas::matvec(minv, f, u);

    //std::cout << "minv = " << std::endl;
    //for (SizeType j = 0; j < num_basis; j++) {
    //  for (SizeType i = 0; i < num_basis; i++) {
    //    std::cout << minv(i,j) << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;

    // Copy local element moments into global array
    for (int n = 0; n < num_basis; n++) projection(elem_id, n) = u(n);
  }
}
