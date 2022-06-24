#pragma once

struct HexRef {

  /*
   * BEGIN FORMER REF_ELEMENT MEMBERS
   */

  // this assumes i,j,k structured indexing of nodes
  const int node_rlid_in_patch_in_cell_[24] = {
      // patches in xi-dir
      0,2,6,4,
      1,3,7,5,
      // patches in eta-dir
      0,4,5,1,
      2,6,7,3,
      // patches in mu-dir
      0,1,3,2,
      4,5,7,6
  };
  
  // return the local ref id for the same patch in the cell-patch neighbor
  const int patch_rlid_cell_neighbor_[6] = {
      1, // xi left side of my cell
      0, // xi right side of my cell
      3, // eta front side of my cell
      2, // eta bact side of my cell
      5, // mu bottom side of my cell
      4  // mu top side of my cell
  };
  
  // return the reference coordinate unit normal of the sides in a cell
  const real_t cell_side_unit_normals_[18] = {
      -1, 0, 0, // xi left side of my cell
      1, 0, 0, // xi right side of my cell
      0,-1, 0, // eta minus side of my cell
      0, 1, 0, // eta plus side of my cell
      0, 0,-1, // mu minus side of my cell
      0, 0, 1  // mu plus side of my cell
  };
  
  int num_ref_nodes_1D_;
  int num_ref_cells_1D_;
  int num_ref_corners_1D_;
  
  int num_ref_surface_nodes_in_elem_;
  int num_ref_inside_nodes_in_elem_;
  
  // Zones
  int num_zones_1d_;
  int num_zones_in_elem_;
  
  // cells
  int num_ref_cells_in_elem_;
  
  CArray <int>cells_in_zone_list_;
  
  // nodes
  int num_ref_nodes_in_elem_;
  int num_ref_nodes_in_cell_;
  
  CArray <int> ref_nodes_in_cell_;
  CArray <int> ref_surface_nodes_in_elem_;
  CArray <int> ref_inside_nodes_in_elem_;
  
  CArray <real_t> ref_node_positions_;
  CArray <real_t> ref_node_g_weights_;
  
  CArray <int> cell_nodes_in_elem_list_;
  
  // Vertices
  int num_ref_verts_1d_;
  int num_ref_verts_in_elem_;
  
  // corners
  int num_ref_corners_in_cell_;
  int num_ref_corners_in_elem_;
  
  CArray <int> ref_corners_in_cell_;
  
  CArray <real_t> ref_corner_surf_normals_;
  CArray <real_t> ref_corner_g_weights_;
  CArray <real_t> ref_corner_surf_g_weights_;
  
  // patches and sides
  int num_patches_in_elem_;
  int num_sides_in_elem_;
  
  CArray <real_t> ref_cell_positions_;
  CArray <real_t> ref_cell_g_weights_;
  CArray <real_t> ref_patch_positions_;
  CArray <real_t> ref_patch_g_weights_;
  
  CArray <int> ref_patches_in_cell_list_;
  
  CArray <real_t> ref_patch_basis_;
  CArray <real_t> ref_patch_gradient_;  // grad basis at patch
  
  CArray <real_t> ref_cell_basis_;
  CArray <real_t> ref_cell_gradient_;   // grad basis at cell
  
  // Num basis functions
  int num_basis_;
  
  // Basis evaluation at nodes
  CArray <real_t> ref_nodal_basis_;
  
  // Pointer to HexN object
  HexN *elem_ptr;

  // Gradient of basis
  CArray <real_t> ref_nodal_gradient_;
  //real_t * ref_nodal_gradient_;
  
  // Function Declarations
  
  // Default constructor
  ref_element() : elem_ptr(nullptr) {};
  
  // Initialize reference element information
  void init(int elem_order, HexN &elem);
  
  int num_dim() const;
  
  int num_basis() const;
  
  int num_ref_nodes() const;
  
  int num_ref_cells_in_elem() const;
  int num_ref_corners_in_cell() const;
  
  // local i,j,k indexing
  int node_rid(int i, int j, int k) const;
  int cell_rid(int i, int j, int k) const;
  int corner_rid(int i, int j, int k) const;
  
  int ref_corners_in_cell(int cell_rid, int corner_rlid) const;
  int ref_nodes_in_cell(int cell_rid, int node_rlid) const;
  
  int ref_surface_nodes_in_elem(int node_rlid) const;
  int ref_inside_nodes_in_elem(int node_rlid) const;
  int num_ref_surface_nodes_in_elem() const;
  int num_ref_inside_nodes_in_elem() const;
  
  real_t ref_node_positions(int node_rid, int dim) const;
  
  real_t ref_corner_surface_normals(int corner_rid, int surf_rlid, int dim) const;
  
  real_t ref_corner_g_surface_weights(int corner_rid, int surf_rlid) const;
  
  real_t ref_node_g_weights(int node_rid) const;
  
  real_t ref_corner_g_weights(int corner_rid) const;
  
  real_t &ref_nodal_gradient(int node_rid, int basis_id, int dim) const;
  
  real_t &ref_nodal_basis(int node_rid, int basis_id) const;
  
  int& cell_lid_in_zone(int zone_lid, int cell_lid) const;
  
  real_t ref_cell_positions(int cell_rid, int dim) const;
  real_t ref_cell_g_weights(int cell_rid) const;
  real_t ref_cell_basis(int cell_rid, int basis_id) const;
  real_t ref_cell_gradient(int cell_rid, int basis_id, int dim) const;
  
  real_t cell_side_unit_normals(int side_rlid, int dim) const;
  
  int& cell_nodes_in_elem(int cell_lid, int node_lid) const;
  
  int node_in_patch_in_cell(int patch_lid, int node_lid) const;
  
  int patch_rlid_in_cell_neighbor(int patch_rlid) const;
  
  int ref_patches_in_cell(int cell_lid, int patch_rlid) const;
  
  real_t ref_patch_positions(int patch_rid, int dim) const;
  real_t ref_patch_g_weights(int patch_rid) const;
  real_t ref_patch_basis(int patch_rid, int basis_id) const;
  real_t ref_patch_gradient(int patch_rid, int basis_id, int dim) const;
  
  int vert_node_map(int vert_lid);

  /*
   * END FORMER REF_ELEMENT MEMBERS
   */


  /*
   * BEGIN FORMER HEXN MEMBERS
   */
            
  const static int num_dim_ = 3;

  // Nodes
  int num_nodes_1d_;
  int num_nodes_;

  MatarRealCArray HexN_Nodes_1d_;
  MatarRealCArray HexN_Nodes_;

  // Vertices
  int num_verts_1d_;
  int num_verts_;
  int num_basis_;

  MatarRealCArray HexN_Verts_1d_;
  MatarRealCArray HexN_Verts_;

  MatarUIntCArray Vert_Node_map_;
  
  int order_;

  void setup_HexN(int elem_order);

  int num_verts();
  int num_nodes();
  int num_basis();
  int node_rid(int i, int j, int k) const;
  int vert_rid(int i, int j, int k) const;

  // Return the noda coordinates in reference space
  real_t &node_coords(int node_rlid, int dim);


  int vert_node_map(int vert_rid) const;
  
  // Evaluate the basis at a given point
  void evaluate_basis(MatarRealCArray &basis, MatarRealCArray &point);

  void build_nodal_gradient(
      MatarRealCArray &gradient);

  // calculate the partial of the basis w.r.t xi at a given point
  void evaluate_derivative_basis_xi(
      MatarRealCArray &partial_xi, 
      MatarRealCArray &point);

  // calculate the partial of the basis w.r.t eta at a given point
  void evaluate_derivative_basis_eta(
      MatarRealCArray &partial_eta, 
      MatarRealCArray &point);

  // calculate the partial of the basis w.r.t mu at a given point
  void evaluate_derivative_basis_zeta(
      MatarRealCArray &partial_mu, 
      MatarRealCArray &point);

  void lagrange_basis_1D(
      MatarRealCArray &interp,    // interpolant
      const real_t &x_point);     // point of interest in element
  
  void lagrange_derivative_1D(
      MatarRealCArray &partials,  //derivative
      const real_t &x_point);     // point of interest in element

  void create_lobatto_nodes(int element_order);

  /*
   * END FORMER HEXN MEMBERS
   */
  

  HexRef(int elem_order);
  ~HexRef();
  

  
};
