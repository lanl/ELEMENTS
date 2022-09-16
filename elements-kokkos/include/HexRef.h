#pragma once


// MATAR includes
#include "MatarTypeDefs.h"


enum BasisType {
  LAGRANGE,
  LEGENDRE,
  BERNSTEIN
};


// Forward declarations
template <typename NumType>
struct TensorProductBasis;


struct HexRef {
  HexRef(int order_basis, BasisType type_basis);
  ~HexRef();


  // Cells (subelements used in volume integration)
  int num_ref_cells_1D_;

  int num_ref_cells_in_elem_;
  int num_ref_cells_in_elem() const;

  int cell_rid(int i, int j, int k) const;

  MatarRealCArray ref_cell_positions_;
  Real ref_cell_positions(int cell_rid, int dim) const;

  MatarRealCArray ref_cell_g_weights_;
  Real ref_cell_g_weights(int cell_rid) const;


  // Patches (subelement faces used in surface integration)
  int num_patches_in_elem_;
  int num_sides_in_elem_;

  const int patch_rlid_cell_neighbor_[6] = {
    1, // side in neighbor cell mated with -xi side in this cell
    0, // side in neighbor cell mated with +xi side in this cell
    3, // side in neighbor cell mated with -eta side in this cell
    2, // side in neighbor cell mated with +eta side in this cell
    5, // side in neighbor cell mated with -zeta side in this cell
    4  // side in neighbor cell mated with +zeta side in this cell
  };
  int patch_rlid_in_cell_neighbor(int patch_rlid) const;

  MatarIntCArray ref_patches_in_cell_list_;
  int ref_patches_in_cell(int cell_lid, int patch_rlid) const;

  MatarRealCArray ref_patch_positions_;
  Real ref_patch_positions(int patch_rid, int dim) const;

  MatarRealCArray ref_patch_g_weights_;
  Real ref_patch_g_weights(int patch_rid) const;

  const Real cell_side_unit_normals_[18] = {
    -1,  0,  0, // -xi 
     1,  0,  0, // +xi 
     0, -1,  0, // -eta 
     0,  1,  0, // +eta 
     0,  0, -1, // -zeta 
     0,  0,  1  // +zeta 
  };
  Real cell_side_unit_normals(int side_rlid, int dim) const;


  // Nodes (quadrature points) and vertices (degrees of freedom)
  int num_ref_nodes_1D_;
  int num_ref_nodes_in_elem_;
  int num_ref_nodes() const;

  int num_ref_verts_1D_;
  int num_ref_verts_in_elem_;

  Real *ref_verts_1D_;


  // Evaluations of basis functions and their gradients at quadrature points
  // in cells and on patches
  static const int num_dim_ = 3;
  int num_dim() const;

  int num_basis_;
  int num_basis() const;
  
  MatarRealCArray ref_cell_basis_;
  Real ref_cell_basis(int cell_rid, int basis_id) const;

  MatarRealCArray ref_cell_gradient_;
  Real ref_cell_gradient(int cell_rid, int basis_id, int dim) const;

  MatarRealCArray ref_patch_basis_;
  Real ref_patch_basis(int patch_rid, int basis_id) const;

  MatarRealCArray ref_patch_gradient_;
  Real ref_patch_gradient(int patch_rid, int basis_id, int dim) const;

  TensorProductBasis<Real> *basis;
};
