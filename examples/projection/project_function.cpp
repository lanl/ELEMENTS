#include "element_types/lagrange_element.h"
#include "element_types/legendre_element.h"
#include "element_types/legendre_polynomials.h"
#include "element_types/elements.h"
#include "io/vtk_io.h"
#include "common/matar_blas_lapack_interface.h"
#include "common/error.h"
#include "common/common.h"

#include <functional>

void project_function_lagrange(
    std::function<Real(Real, Real, Real)> func, 
    SwageMesh &mesh, CArray<Real> &projection) {
  // Mesh parameters
  const int num_dim      = 3;
  int num_elems    = mesh.num_elems();
  int elem_order   = mesh.elem_order();
  int num_verts_1d = elem_order + 1;

  // Create reference element
  CArray<Real> vert_coords_1d(num_verts_1d);
  elements::lobatto_nodes_1D(vert_coords_1d, num_verts_1d);

  LagrangeElement<Real> elem(elem_order, vert_coords_1d.pointer());

  // Quadrature points and weights, assuming Lobatto quadrature and assuming
  // nodes are quadrature points
  int num_nodes_1d = 2*elem_order + 1;
  CArray<Real> lobatto_points(num_nodes_1d);
  CArray<Real> lobatto_weights(num_nodes_1d);

  elements::lobatto_nodes_1D(lobatto_points, num_nodes_1d);
  elements::lobatto_weights_1D(lobatto_weights, num_nodes_1d);

  SizeType node_radices[num_dim] = {SizeType(num_nodes_1d), 
    SizeType(num_nodes_1d), SizeType(num_nodes_1d)};

  // Allocate memory for mass matrix, right-hand side, and solution
  CArray<Real> m(elem.Ne, elem.Ne);
  CArray<Real> minv(elem.Ne, elem.Ne);
  CArray<Real> f(elem.Ne);
  CArray<Real> u(elem.Ne);

  // Loop over elements
  for (int elem_id = 0; elem_id < num_elems; elem_id++) {
    int num_nodes_in_elem = mesh.num_nodes_in_elem();

    for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
      // Get quadrature point coordinates and weights
      SizeType lmn[num_dim];
      common::base_10_to_mixed_radix(num_dim, node_radices, 
          SizeType(node_lid), lmn);
      Real X[3] = {lobatto_points(lmn[0]), lobatto_points(lmn[1]), 
          lobatto_points(lmn[2])};

      int node_gid = mesh.nodes_in_elem(elem_id, node_lid);
      Real x = mesh.node_coords(node_gid, 0);
      Real y = mesh.node_coords(node_gid, 1);
      Real z = mesh.node_coords(node_gid, 2);

      Real wl = lobatto_weights(lmn[0]);
      Real wm = lobatto_weights(lmn[1]);
      Real wn = lobatto_weights(lmn[2]);

      // Evaluate quadrature weight times Jacobian determinant
      Real wxj = wl*wm*wn*mesh.gauss_pt_det_j(node_gid); 

      for (SizeType j = 0; j < elem.Ne; j++) {
        Real phi_j = elem.eval_basis(j, X);

        // Build mass matrix
        for (SizeType i = 0; i < elem.Ne; i++) {
          Real phi_i = elem.eval_basis(i, X);
          m(i,j) += phi_i*phi_j*wxj;
        }

        // Build right-hand side
        std::cout << func(x, y, z) << std::endl;
        f(j) += func(x, y, z)*phi_j*wxj;
      }
    }

    // Solve for moments over element
    matar2lapack::invert(m, minv);
    matar2blas::matvec(minv, f, u);

    //std::cout << "minv = " << std::endl;
    //for (SizeType j = 0; j < elem.Ne; j++) {
    //  for (SizeType i = 0; i < elem.Ne; i++) {
    //    std::cout << minv(i,j) << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;

    // Copy local element moments into global array
    for (int n = 0; n < elem.Ne; n++) projection(elem_id, n) = u(n);
  }
}

void project_function_legendre(
    std::function<Real(Real, Real, Real)> func, 
    SwageMesh &mesh, CArray<Real> &projection) {
  // Mesh parameters
  const int num_dim      = 3;
  int num_elems    = mesh.num_elems();
  int elem_order   = mesh.elem_order();
  int num_verts_1d = elem_order + 1;

  // Create reference element
  LegendreElement<Real> elem(elem_order);

  // Basis function parameters
  int num_basis = pow(num_verts_1d, num_dim);
  const SizeType num_rad = SizeType(num_dim);
  SizeType vert_radices[num_rad] = {SizeType(num_verts_1d), 
    SizeType(num_verts_1d), SizeType(num_verts_1d)};

  // Quadrature points and weights, assuming Lobatto quadrature and assuming
  // nodes are quadrature points
  int num_nodes_1d = 2*elem_order + 1;
  CArray<Real> lobatto_points(num_nodes_1d);
  CArray<Real> lobatto_weights(num_nodes_1d);

  elements::lobatto_nodes_1D(lobatto_points, num_nodes_1d);
  elements::lobatto_weights_1D(lobatto_weights, num_nodes_1d);

  SizeType node_radices[num_rad] = {SizeType(num_nodes_1d), 
    SizeType(num_nodes_1d), SizeType(num_nodes_1d)};

  // Allocate memory for mass matrix, right-hand side, and solution
  CArray<Real> m(num_basis, num_basis);
  CArray<Real> minv(num_basis, num_basis);
  CArray<Real> f(num_basis);
  CArray<Real> u(num_basis);

  // Loop over elements
  for (int elem_id = 0; elem_id < num_elems; elem_id++) {
    int num_nodes_in_elem = mesh.num_nodes_in_elem();

    for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
      // Get quadrature point coordinates and weights
      SizeType lmn[num_dim];
      common::base_10_to_mixed_radix(num_dim, node_radices, 
          SizeType(node_lid), lmn);
      Real X[3] = {lobatto_points(lmn[0]), lobatto_points(lmn[1]), 
          lobatto_points(lmn[2])};

      int node_gid = mesh.nodes_in_elem(elem_id, node_lid);
      Real x = mesh.node_coords(node_gid, 0);
      Real y = mesh.node_coords(node_gid, 1);
      Real z = mesh.node_coords(node_gid, 2);

      Real wl = lobatto_weights(lmn[0]);
      Real wm = lobatto_weights(lmn[1]);
      Real wn = lobatto_weights(lmn[2]);

      // Evaluate quadrature weight times Jacobian determinant
      Real wxj = wl*wm*wn*mesh.gauss_pt_det_j(node_gid); 

      for (SizeType j = 0; j < elem.Ne; j++) {
        Real phi_j = elem.eval_basis(j, X);

        // Build mass matrix
        for (SizeType i = 0; i < elem.Ne; i++) {
          Real phi_i = elem.eval_basis(i, X);
          m(i,j) += phi_i*phi_j*wxj;
        }

        // Build right-hand side
        std::cout << func(x, y, z) << std::endl;
        f(j) += func(x, y, z)*phi_j*wxj;
      }
    }

    // Solve for moments over element
    matar2lapack::invert(m, minv);
    matar2blas::matvec(minv, f, u);

    //std::cout << "minv = " << std::endl;
    //for (SizeType j = 0; j < elem.Ne; j++) {
    //  for (SizeType i = 0; i < elem.Ne; i++) {
    //    std::cout << minv(i,j) << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << std::endl;

    // Copy local element moments into global array
    for (int n = 0; n < elem.Ne; n++) projection(elem_id, n) = u(n);
  }
}
