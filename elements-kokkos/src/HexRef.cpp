#include "HexRef.h"
#include "LagrangeProductBasis.h"
#include "LegendreProductBasis.h"
#include "BernsteinProductBasis.h"
#include "points_and_weights.h"
#include "error.h"


HexRef::HexRef(int order_basis, BasisType type_basis=LAGRANGE) {
  // Determine sizes of data structures based on the order of the basis
  assert(order_basis >= 0 && "Error: order of basis must be positive");
  if (order_basis == 0){       
    num_ref_nodes_1D_ = 2;
    num_ref_verts_1D_ = 2;
  } else if (order_basis > 0) {
    num_ref_nodes_1D_ = 2*order_basis + 1;
    num_ref_verts_1D_ = order_basis + 1;
  }

  num_ref_cells_1D_ = num_ref_nodes_1D_ - 1;
  num_ref_cells_in_elem_ = (num_ref_nodes_1D_ - 1)*(num_ref_nodes_1D_ - 1)
      *(num_ref_nodes_1D_ - 1);

  num_patches_in_elem_ = (num_ref_cells_1D_ + 1)
      *(num_ref_cells_1D_*num_ref_cells_1D_)*3;
  num_sides_in_elem_ = num_ref_cells_in_elem_*6;

  num_ref_verts_in_elem_ = num_ref_verts_1D_*num_ref_verts_1D_
      *num_ref_verts_1D_;
  num_basis_ = num_ref_verts_in_elem_;


  // Allocate memory for data structures
  ref_cell_positions_ = MatarRealCArray(num_ref_cells_in_elem_, 3);
  ref_cell_g_weights_ = MatarRealCArray(num_ref_cells_in_elem_);

  ref_patches_in_cell_list_ = MatarIntCArray(num_sides_in_elem_);
  ref_patch_positions_ = MatarRealCArray(num_patches_in_elem_, 3);
  ref_patch_g_weights_ = MatarRealCArray(num_patches_in_elem_);

  ref_cell_basis_ = MatarRealCArray(num_ref_cells_in_elem_, num_basis_);
  ref_cell_gradient_ = MatarRealCArray(num_ref_cells_in_elem_, num_basis_, num_dim_);

  ref_patch_basis_ =  MatarRealCArray(num_patches_in_elem_, num_basis_);
  ref_patch_gradient_ = MatarRealCArray(num_patches_in_elem_, num_basis_, num_dim_);

  ref_verts_1D_ = new Real[num_ref_verts_1D_];


  // Retrieve Gauss-Legendre and Gauss-Lobatto quadrature points and weights
  auto leg_points_1D = MatarRealCArray(num_ref_cells_1D_);
  legendre_points(leg_points_1D, num_ref_cells_1D_);

  auto leg_weights_1D = MatarRealCArray(num_ref_cells_1D_);
  legendre_weights(leg_weights_1D, num_ref_cells_1D_);

  auto lob_points_1D = MatarRealCArray(num_ref_nodes_1D_);
  lobatto_points(lob_points_1D, num_ref_nodes_1D_);
  
  auto lob_weights_1D = MatarRealCArray(num_ref_nodes_1D_);
  lobatto_weights(lob_weights_1D, num_ref_nodes_1D_);


  // Store every other node (quadrature point) coordinate as a vertex coordinate
  int vert_id = 0, node_id = 0;
  while (node_id < num_ref_nodes_1D_) {
    ref_verts_1D_[vert_id] = lob_points_1D(node_id);
    vert_id += 1;
    node_id += 2;
  }


  // Fill in cell data structures
  for (int k = 0; k < num_ref_cells_1D_; k++) {
    for (int j = 0; j < num_ref_cells_1D_; j++) {
      for (int i = 0; i < num_ref_cells_1D_; i++) {
        int this_cell_rid = cell_rid(i, j, k);
        
        ref_cell_positions_(this_cell_rid, 0) = leg_points_1D(i);
        ref_cell_positions_(this_cell_rid, 1) = leg_points_1D(j);
        ref_cell_positions_(this_cell_rid, 2) = leg_points_1D(k);
        
        ref_cell_g_weights_(this_cell_rid) = leg_weights_1D(i)
                *leg_weights_1D(j)*leg_weights_1D(k);
      }
    }
  }


  // Fill in patch data structures
  for (int side_rid = 0; side_rid < num_sides_in_elem_; side_rid++) {
    ref_patches_in_cell_list_(side_rid) = -1;
  }

  int patch_rid = 0; 
  for (int k = 0; k < num_ref_cells_1D_; k++) {
    for (int j = 0; j < num_ref_cells_1D_; j++) {
      for (int i = 0; i < num_ref_cells_1D_; i++) {
        int this_cell_rid = cell_rid(i, j, k);

        for (int side_lid = 0; side_lid < 6; side_lid++) {
          int side_rid = side_lid + this_cell_rid*6;

          // Check to see that this patch has not been seen already
          if (ref_patches_in_cell_list_(side_rid) == -1) {
            // Store the index for this patch in the cell
            ref_patches_in_cell_list_(side_rid) = patch_rid;
            
            // Store the quadrature point area values on the patch
            switch (side_lid) {
              case 0:  // patch with -xi normal
                ref_patch_positions_(patch_rid,0) = lob_points_1D(i); 
                ref_patch_positions_(patch_rid,1) = leg_points_1D(j);
                ref_patch_positions_(patch_rid,2) = leg_points_1D(k);
                
                ref_patch_g_weights_(patch_rid) = leg_weights_1D(j)*leg_weights_1D(k);
                break;
              case 1:  // patch with +xi normal
                ref_patch_positions_(patch_rid,0) = lob_points_1D(i+1);
                ref_patch_positions_(patch_rid,1) = leg_points_1D(j);
                ref_patch_positions_(patch_rid,2) = leg_points_1D(k);
                
                ref_patch_g_weights_(patch_rid) = leg_weights_1D(j)*leg_weights_1D(k);
                break;
              case 2:  // patch with -eta normal
                ref_patch_positions_(patch_rid,0) = leg_points_1D(i);
                ref_patch_positions_(patch_rid,1) = lob_points_1D(j); 
                ref_patch_positions_(patch_rid,2) = leg_points_1D(k);
                
                ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(k);
                break;
              case 3:  // patch with +eta normal
                ref_patch_positions_(patch_rid,0) = leg_points_1D(i); 
                ref_patch_positions_(patch_rid,1) = lob_points_1D(j+1);
                ref_patch_positions_(patch_rid,2) = leg_points_1D(k); 
                
                ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(k);
                break;
              case 4:  // patch with -zeta normal
                ref_patch_positions_(patch_rid,0) = leg_points_1D(i);
                ref_patch_positions_(patch_rid,1) = leg_points_1D(j);
                ref_patch_positions_(patch_rid,2) = lob_points_1D(k);
                
                ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(j);
                break;
              case 5:  // patch with +zeta normal
                ref_patch_positions_(patch_rid,0) = leg_points_1D(i);
                ref_patch_positions_(patch_rid,1) = leg_points_1D(j);
                ref_patch_positions_(patch_rid,2) = lob_points_1D(k+1);
                
                ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(j);
                break;
            }
                     
            // Get the rID for the neighboring cell in the element
            // that communicates with this cell across this patch
            int neighbor_cell_rid;
            if (side_lid == 0 && 0 < i) {
              neighbor_cell_rid = cell_rid(i-1, j, k);
            } else if (side_lid == 1 && i < num_ref_cells_1D_ - 1) {
              neighbor_cell_rid = cell_rid(i+1, j, k);
            } else if (side_lid == 2 && 0 < j) {
              neighbor_cell_rid = cell_rid(i, j-1, k);
            } else if (side_lid == 3 && j < num_ref_cells_1D_ - 1) {
              neighbor_cell_rid = cell_rid(i, j+1, k);
            } else if (side_lid == 4 && 0 < k) {  
              neighbor_cell_rid = cell_rid(i, j, k-1);
            } else if (side_lid == 5 && k < num_ref_cells_1D_ - 1) {
              neighbor_cell_rid = cell_rid(i, j, k+1);
            } else { 
              neighbor_cell_rid = -1;  // no neighboring cell in element
            } 
            
            // If there is a neighboring cell, set this patch rID in its list
            if (0 <= neighbor_cell_rid && neighbor_cell_rid < num_ref_cells_in_elem_) {
              // Get the side lid for this patch in the neighboring cell
              int neighbor_side_lid = patch_rlid_in_cell_neighbor(side_lid);
              
              // Get side rID from that side lid
              int neighbor_side_rid = neighbor_side_lid + (neighbor_cell_rid)*6;
              
              // Store the patch rID for this patch in the neighboring cell
              ref_patches_in_cell_list_(neighbor_side_rid) = patch_rid;
            }
            
            // Increment the patch index
            patch_rid ++;
          }
        }
      }
    }
  }


  // Initialize basis object based on type specified
  switch (type_basis) {
    case BasisType::LAGRANGE:
      basis = new LagrangeProductBasis<Real>(order_basis, ref_verts_1D_);
      break;
    case BasisType::LEGENDRE:
      basis = new LegendreProductBasis<Real>(order_basis);
      break;
    case BasisType::BERNSTEIN:
      basis = new BernsteinProductBasis<Real>(order_basis);
      break;
    default:
      basis = nullptr;
      std::string err_msg("Error: basis class not implemented for this type.");
      throw NotImplementedError(err_msg);
  }


  // Fill in evaluations of basis functions and their gradients in cells 
  for (int cell_rlid = 0; cell_rlid < num_ref_cells_in_elem_; cell_rlid++) {
    // Get the cell coordinates
    Real point[num_dim_];
    for (int dim = 0; dim < num_dim_; dim++)
      point[dim] = ref_cell_positions(cell_rlid, dim);
    
    // Compute and store all the basis values at the cell
    for (int basis_id = 0; basis_id < num_basis_; basis_id++)
      ref_cell_basis_(cell_rlid, basis_id) = basis->eval_basis(basis_id, point);
    
    // Compute and store all the basis gradient values at the cell
    for (int basis_id = 0; basis_id < num_basis_; basis_id++) {
      Real grad_basis[num_dim_];
      basis->eval_grad_basis(basis_id, point, grad_basis);
      ref_cell_gradient_(cell_rlid, basis_id, 0) = grad_basis[0];
      ref_cell_gradient_(cell_rlid, basis_id, 1) = grad_basis[1];
      ref_cell_gradient_(cell_rlid, basis_id, 2) = grad_basis[2];
    } 
  } 


  // Fill in evaluations of basis functions and their gradients on patches
  for (int patch_rlid = 0; patch_rlid < num_patches_in_elem_; patch_rlid++) {
    // Get the cell coordinates
    // Get the patch coordinates
    Real point[num_dim_];
    for (int dim = 0; dim < num_dim_; dim++)
      point[dim] = ref_patch_positions(patch_rlid, dim);
    
    // Compute and store all the basis values at the patch
    for (int basis_id = 0; basis_id < num_basis_; basis_id++)
      ref_patch_basis_(patch_rlid, basis_id) = basis->eval_basis(basis_id, 
          point);
    
    // loop over the basis polynomials where there is one from each vertex
    for (int basis_id = 0; basis_id < num_basis_; basis_id++) {
      Real grad_basis[num_dim_];
      basis->eval_grad_basis(basis_id, point, grad_basis);
      ref_patch_gradient_(patch_rlid, basis_id, 0) = grad_basis[0];
      ref_patch_gradient_(patch_rlid, basis_id, 1) = grad_basis[1];
      ref_patch_gradient_(patch_rlid, basis_id, 2) = grad_basis[2];
    } 
  } 
}


HexRef::~HexRef() {
  delete[] ref_verts_1D_;
  delete basis;
}


int HexRef::num_ref_cells_in_elem() const {
    return num_ref_cells_in_elem_;
};


int HexRef::cell_rid(int i, int j, int k) const {
  return i + j*num_ref_cells_1D_ + k*num_ref_cells_1D_*num_ref_cells_1D_;
};


Real HexRef::ref_cell_positions(int cell_rid, int dim) const {
  return ref_cell_positions_(cell_rid, dim);
}


Real HexRef::ref_cell_g_weights(int cell_rid) const {
  return ref_cell_g_weights_(cell_rid);
}


int HexRef::patch_rlid_in_cell_neighbor(int patch_rlid) const {
  return patch_rlid_cell_neighbor_[patch_rlid];
}


int HexRef::ref_patches_in_cell(int cell_rid, int patch_rlid) const {
  int index = patch_rlid + cell_rid*6;
  return ref_patches_in_cell_list_(index);
}


Real HexRef::ref_patch_positions(int patch_rid, int dim) const {
  return ref_patch_positions_(patch_rid, dim);
}


Real HexRef::ref_patch_g_weights(int patch_rid) const {
  return ref_patch_g_weights_(patch_rid);
}


Real HexRef::cell_side_unit_normals(int side_rlid, int dim) const {
  int index = side_rlid*3 + dim;
  return cell_side_unit_normals_[index];
}


int HexRef::num_ref_nodes() const {
  return num_ref_nodes_in_elem_;
}


int HexRef::num_dim() const {
  return num_dim_;
}


int HexRef::num_basis() const {
  return num_basis_;
}


Real HexRef::ref_cell_basis(int cell_rid, int basis_id) const {
    return ref_cell_basis_(cell_rid, basis_id);
};


Real HexRef::ref_cell_gradient(int cell_rid, int basis_id, 
    int dim) const {
    return ref_cell_gradient_(cell_rid, basis_id, dim);
};


Real HexRef::ref_patch_basis(int patch_rid, int basis_id) const {
    return ref_patch_basis_(patch_rid, basis_id);
};


Real HexRef::ref_patch_gradient(int patch_rid, int basis_id, 
    int dim) const {
    return ref_patch_gradient_(patch_rid, basis_id, dim);
};
