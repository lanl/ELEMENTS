#include "lagrange_polynomials.h"
#include "vtk_io.h"
#include "common.h"

void create_annular_mesh(int elem_order, SwageMesh &mesh) {
  int Nd  = 3;
  int Ne  = 2;
  int Np  = elem_order;
  int Nvd = Np + 1;
  int Nve = std::pow(Nvd, Nd);
  int Nv  = Nve*Ne - Nvd*Nvd;

  NumType ri = 0.5; // 0.5*C/M_PI - 0.5*dr;
  NumType ro = 2.5; // 0.5*C/M_PI + 0.5*dr;

  mesh.init_nodes(Nv);
  mesh.init_element(Np, Nd, Ne);

  // Add vertices and connectivity for first element
  int vert_cnt = 0;
  NumType theta_lb = 0.0;
  NumType theta_ub = 0.5*M_PI;
  for (int k = 0; k < Nvd; k++) {
    for (int j = 0; j < Nvd; j++) {
      for (int i = 0; i < Nvd; i++) {
        NumType theta = theta_lb + i*(theta_ub - theta_lb)/(Nvd - 1);
        NumType r = ri + j*(ro - ri)/(Nvd - 1);

        NumType x = -r*std::cos(theta);
        NumType y = r*std::sin(theta);
        NumType z = -1.0 + 2.0*k/(Nvd - 1);

        mesh.node_coords(vert_cnt, 0) = x;
        mesh.node_coords(vert_cnt, 1) = y;
        mesh.node_coords(vert_cnt, 2) = z;

        int ii = i + j*Nvd + k*Nvd*Nvd;
        mesh.nodes_in_elem(0, ii) = vert_cnt;

        if (i == Nvd - 1) {
          mesh.nodes_in_elem(1, ii - i) = vert_cnt;
        }

        vert_cnt++;
      }
    }
  }

  // Add vertices and connectivity for second element
  theta_lb = 0.5*M_PI;
  theta_ub = M_PI;
  for (int k = 0; k < Nvd; k++) {
    for (int j = 0; j < Nvd; j++) {
      for (int i = 1; i < Nvd; i++) {
        NumType theta = theta_lb + i*(theta_ub - theta_lb)/(Nvd - 1);
        NumType r = ri + j*(ro - ri)/(Nvd - 1);

        NumType x = -r*std::cos(theta);
        NumType y = r*std::sin(theta);
        NumType z = -1.0 + 2.0*k/(Nvd - 1);

        mesh.node_coords(vert_cnt, 0) = x;
        mesh.node_coords(vert_cnt, 1) = y;
        mesh.node_coords(vert_cnt, 2) = z;

        int ii = i + j*Nvd + k*Nvd*Nvd;
        mesh.nodes_in_elem(1, ii) = vert_cnt;

        vert_cnt++;
      }
    }
  }
}

CArray<NumType> compute_jacobian_determinants(SwageMesh &mesh) {
  int Nd  = 3;                  // number of dimensions
  int Ne  = mesh.num_elems();   // number of elements
  int Np  = mesh.elem_order();  // polynomial order
  int Nvd = Np + 1;             // number of vertices in 1D
  int Nve = std::pow(Nvd, Nd);  // number of vertices per element

  CArray<NumType> jacobian_determinants(mesh.num_elems(), Nve);

  NumType ze[Nve]; 
  NumType we[Nve]; 
  NumType cx[Nve]; 
  NumType cy[Nve]; 
  NumType cz[Nve]; 
  NumType g[3*Nvd+1];
  NumType zl = -1.0, zr = 1.0;
  lagrange::equispaced_points(Nvd, zl, zr, ze);
  lagrange::compute_barycentric_weights(Nvd, ze, we);

  for (int elem_id = 0; elem_id < Ne; elem_id++) {
    // Extract the vertex coordinates into the cx, cy, cz arrays
    for (int vert_lid = 0; vert_lid < Nve; vert_lid++) {
      int vert_gid = mesh.nodes_in_elem(elem_id, vert_lid);
      cx[vert_lid] = mesh.node_coords(vert_gid, 0);
      cy[vert_lid] = mesh.node_coords(vert_gid, 1);
      cz[vert_lid] = mesh.node_coords(vert_gid, 2);
    }

    // Compute Jacobian at each node
    for (int k = 0; k < Nvd; k++) {
      for (int j = 0; j < Nvd; j++) {
        for (int i = 0; i < Nvd; i++) {
          int ii = i + j*Nvd + k*Nvd*Nvd;
          NumType X[3] = {ze[i], ze[j], ze[k]};
          NumType J;

          SizeType N = Nvd;
          lagrange::compute_jacobian_determinant(
              N, ze, we, cx, cy, cz, g, X, J);
          jacobian_determinants(elem_id, ii) = J;
        }
      }
    }
  }

  return jacobian_determinants;
}

int main() {
  // Create a 2-element half-annulus mesh 
  int elem_order = 2;
  SwageMesh mesh;
  create_annular_mesh(elem_order, mesh);

  // Compute the Jacobians on the vertices of the input mesh
  CArray<NumType> jac_dets = compute_jacobian_determinants(mesh);

  // Initialize a VTK grid file from the  the mesh and output
  std::string sol_name("j");
  VtkGrid grid = swage2vtk::init_vtk_grid(mesh, sol_name, jac_dets);

  std::string output_name("output.vtk");
  swage2vtk::write_grid(grid, output_name);

  return 0;
}
