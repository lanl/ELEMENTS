#include <HexRef.h>

HexRef::HexRef(int orderBasis, BasisType typeBasis=LAGRANGE, 
    QuadratureType typeQuad=GAUSS_LOBATTO, bool verbose=false) {

  /* Determine sizes of data structures based on the order of the basis */
  int numGaussPts1d = 2;
  if(orderBasis == 0){       
      numGaussPts1d = 2;

      num_ref_nodes_1D_ = numGaussPts1d;
      num_ref_verts_1d_ = 2;
      num_zones_1d_ = 1;
  
      cells_in_zone_list_ = MatarIntCArray(num_zones_in_elem_, 1);
  } else if (orderBasis > 0) {
      numGaussPts1d = 2*orderBasis + 1;

      num_ref_nodes_1D_ = numGaussPts1d;
      num_ref_verts_1d_ = orderBasis+1;
      num_zones_1d_ = (num_ref_nodes_1D_ - 1) / 2;

      cells_in_zone_list_ = MatarIntCArray(num_zones_in_elem_, 8);
  } else {
    std::cout << "Error: element order must be positive" << std::endl;
  }

  num_ref_cells_1D_ = num_ref_nodes_1D_ - 1;
  num_ref_corners_1D_ = 2*(num_ref_nodes_1D_ - 2) + 2;

  num_zones_in_elem_ = num_zones_1d_*num_zones_1d_*num_zones_1d_;
  num_ref_nodes_in_elem_ = num_ref_nodes_1D_*num_ref_nodes_1D_*num_ref_nodes_1D_;
  num_ref_verts_in_elem_ = num_ref_verts_1d_*num_ref_verts_1d_*num_ref_verts_1d_;
  num_ref_cells_in_elem_ = num_ref_cells_1D_*num_ref_cells_*num_ref_cells_1D;

  if (verbose) {
    std::cout << "Num cells in reference  element: "
              << num_ref_cells_in_elem_
              << std::endl;
  }
  
  cell_nodes_in_elem_list_ = MatarIntCArray (num_ref_cells_in_elem_, 8);
  num_ref_corners_in_cell_ = 8;
  num_ref_corners_in_elem_ = num_ref_corners_1D_*num_ref_corners_1D_*num_ref_corners_1D_;
  num_basis_ = num_ref_verts_in_elem_;

  num_ref_inside_nodes_in_elem_  = (num_ref_nodes_1D_ - 2)*(num_ref_nodes_1D_ - 2)*(num_ref_nodes_1D_ - 2);
  num_ref_surface_nodes_in_elem_ = num_ref_nodes_1D_*num_ref_nodes_1D_*num_ref_nodes_1D_  - num_ref_inside_nodes_in_elem_;
  

  /* Allocate arrays */

  ref_node_positions_ = CArray <real_t> (num_ref_nodes_in_elem_, num_dim_);
  ref_node_g_weights_ = CArray <real_t> (num_ref_nodes_in_elem_);

  ref_corners_in_cell_ = MatarIntCArray (num_ref_corners_in_elem_);
  ref_corner_g_weights_ = CArray <real_t> (num_ref_corners_in_elem_);
  ref_corner_surf_g_weights_ = CArray <real_t> (num_ref_corners_in_elem_, num_dim_);
  ref_corner_surf_normals_ = CArray <real_t> (num_ref_corners_in_elem_, num_dim_, num_dim_);

  ref_nodes_in_cell_ = MatarIntCArray(num_ref_corners_in_elem_);

  ref_nodal_basis_ = CArray <real_t> (num_ref_nodes_in_elem_, num_basis_);
  ref_nodal_gradient_ = CArray <real_t> (num_ref_nodes_in_elem_, num_basis_, num_dim_);

  ref_surface_nodes_in_elem_ = MatarIntCArray(num_ref_surface_nodes_in_elem_);
  ref_inside_nodes_in_elem_  = MatarIntCArray(num_ref_inside_nodes_in_elem_);


  /* Populate arrays containing quadrature points and weights */
  auto lob_nodes_1D = MatarRealCArray(num_ref_nodes_1D_);
  lobatto_nodes_1D(lob_nodes_1D, num_ref_nodes_1D_);

  auto lob_weights_1D = MatarRealCArray(num_ref_nodes_1D_);
  lobatto_weights_1D(lob_weights_1D, num_ref_nodes_1D_);

  for (int k = 0; k < num_ref_nodes_1D_; k++){
      for (int j = 0; j < num_ref_nodes_1D_; j++){
          for (int i = 0; i < num_ref_nodes_1D_; i++){
              
              int n_rid = node_rid(i,j,k);
              
              ref_node_positions_(n_rid,0) = lob_nodes_1D(i);
              ref_node_positions_(n_rid,1) = lob_nodes_1D(j);
              ref_node_positions_(n_rid,2) = lob_nodes_1D(k);
              
              ref_node_g_weights_(n_rid) = lob_weights_1D(i)*lob_weights_1D(j)*lob_weights_1D(k);
          }
      }
  }

  /* 
   * BEGIN ref_element::init source up to call to HexN::setup_HexN
   */

  // must partition the nodal gauss weights to the corners
  auto corner_lob_weights_1D = CArray <real_t>(num_ref_corners_1D_);
  auto r_corner_part_g_weights = CArray <real_t>(num_ref_corners_in_elem_, num_dim_);

  // loop over interior corners in 1D
  corner_lob_weights_1D(0) = lob_weights_1D(0);
  
  for(int i = 1; i < num_ref_nodes_1D_ - 1; i++){
      
      // get the corner_rid index in 1D for the left and right corners
      int corner_left = (2*i) - 1;
      int corner_right = 2*i;
      
      // WARNING WARNING WARNING: partitioning using an average
      corner_lob_weights_1D(corner_left)  = 0.5*lob_weights_1D(i);
      corner_lob_weights_1D(corner_right) = 0.5*lob_weights_1D(i);
  }
  
  corner_lob_weights_1D(num_ref_corners_1D_ - 1) = lob_weights_1D(num_ref_nodes_1D_ - 1);

  if(orderBasis == 0){

      cell_lid_in_zone(0, 0) = 0;
  }

  else{
      // Save the Cell_lid to the zones
      int zone_rid = 0;
      for(int k = 0; k < num_zones_1d_; k++){
          for(int j = 0; j < num_zones_1d_; j++){
              for(int i = 0; i < num_zones_1d_; i++){


                  cell_lid_in_zone(zone_rid, 0) = cell_rid(2*i  , 2*j  , 2*k  );
                  cell_lid_in_zone(zone_rid, 1) = cell_rid(2*i+1, 2*j  , 2*k  );
                  cell_lid_in_zone(zone_rid, 2) = cell_rid(2*i  , 2*j+1, 2*k  );
                  cell_lid_in_zone(zone_rid, 3) = cell_rid(2*i+1, 2*j+1, 2*k  );
                  cell_lid_in_zone(zone_rid, 4) = cell_rid(2*i  , 2*j  , 2*k+1);
                  cell_lid_in_zone(zone_rid, 5) = cell_rid(2*i+1, 2*j  , 2*k+1);
                  cell_lid_in_zone(zone_rid, 6) = cell_rid(2*i  , 2*j+1, 2*k+1);
                  cell_lid_in_zone(zone_rid, 7) = cell_rid(2*i+1, 2*j+1, 2*k+1);

                  zone_rid++;
              }
          }   
      }

  }
  
  // Save the node in cell lid to node rid in elem map
  int cell_lid = 0;
  int num_ref_cells_1d = (num_ref_nodes_1D_-1);
  for(int k = 0; k < num_ref_cells_1d; k++){
      for(int j = 0; j < num_ref_cells_1d; j++){
          for(int i = 0; i < num_ref_cells_1d; i++){
              
              cell_nodes_in_elem(cell_lid, 0) = node_rid(i  , j,   k);
              cell_nodes_in_elem(cell_lid, 1) = node_rid(i+1, j,   k);
              cell_nodes_in_elem(cell_lid, 2) = node_rid(i  , j+1, k);
              cell_nodes_in_elem(cell_lid, 3) = node_rid(i+1, j+1, k);
              cell_nodes_in_elem(cell_lid, 4) = node_rid(i  , j,   k+1);
              cell_nodes_in_elem(cell_lid, 5) = node_rid(i+1, j,   k+1);
              cell_nodes_in_elem(cell_lid, 6) = node_rid(i  , j+1, k+1);
              cell_nodes_in_elem(cell_lid, 7) = node_rid(i+1, j+1, k+1);
              
              cell_lid++;
          }
      }
  }
  
  // i,j,k indicies are for corners
  for(int k = 0; k < num_ref_corners_1D_; k++){
      for(int j = 0; j < num_ref_corners_1D_; j++){
          for(int i = 0; i < num_ref_corners_1D_; i++){
              
              int crn_rid = corner_rid(i, j, k);  // the ref space id
              
              r_corner_part_g_weights(crn_rid, 0) = corner_lob_weights_1D(i);
              r_corner_part_g_weights(crn_rid, 1) = corner_lob_weights_1D(j);
              r_corner_part_g_weights(crn_rid, 2) = corner_lob_weights_1D(k);
              
              ref_corner_g_weights_(crn_rid) =
                  corner_lob_weights_1D(i)*corner_lob_weights_1D(j)*corner_lob_weights_1D(k);
              
          }
      }
  } // end loop over i,j,k for corners in element



  // --- build corners ---
  auto unit_normals = CArray <real_t> (num_ref_corners_in_cell_, num_dim_);
  
  set_unit_normals(unit_normals);
  
  
  // loop over the cells in the elem (i,j,k are ref cell indices)
  int index = 0; // index = num_cells*num_corners_in_cell + corner_lid
  for(int k = 0; k < num_ref_cells_1D_; k++){
      for(int j = 0; j < num_ref_cells_1D_; j++){
          for(int i = 0; i < num_ref_cells_1D_; i++){
              
              int corner_rid_in_cell_0 = 2*i + 2*(num_ref_corners_1D_)*j
                                       + 2*(num_ref_corners_1D_)*(num_ref_corners_1D_)*k;
              
              // loop over the corners in the cell
              int corner_rlid = 0;
              for(int k_rlid = 0; k_rlid < 2; k_rlid++){
                  for(int j_rlid = 0; j_rlid < 2; j_rlid++){
                      for(int i_rlid = 0; i_rlid < 2; i_rlid++){
                          
                          // calculate the reference corner index from the rlid's
                          int crn_rid = corner_rid_in_cell_0
                              + i_rlid + (num_ref_corners_1D_)*j_rlid
                              + (num_ref_corners_1D_)*(num_ref_corners_1D_)*k_rlid;
                          
                          ref_corners_in_cell_(index) = crn_rid; // save the rid
                          
                          // node ref index
                          ref_nodes_in_cell_(index) = node_rid(i + i_rlid, j + j_rlid, k + k_rlid);
                          int n_rid = node_rid(i + i_rlid, j + j_rlid, k + k_rlid);
                          
                          
                          // 3 surfaces vectors in each ref corner that have 3 components
                          real_t surf_vec0[3] = {unit_normals(corner_rlid,0), 0.0, 0.0};  // s_0
                          real_t surf_vec1[3] = {0.0, unit_normals(corner_rlid,1), 0.0};  // s_1
                          real_t surf_vec2[3] = {0.0, 0.0, unit_normals(corner_rlid,2)};  // s_2
                          
                          
                          // surface unit normal 0 and surface quadrature 0
                          int surf_rlid = 0;
                          for(int dim = 0; dim < num_dim_; dim++){
                              ref_corner_surf_normals_(crn_rid, surf_rlid, dim) = surf_vec0[dim];
                          }
                          
                          // power coef is =1 for the quadrature values for the surf normal else =0
                          ref_corner_surf_g_weights_(crn_rid, surf_rlid) =
                                pow( r_corner_part_g_weights(crn_rid, 0), (1.0 - fabs(surf_vec0[0])) )
                              * pow( r_corner_part_g_weights(crn_rid, 1), (1.0 - fabs(surf_vec0[1])) )
                              * pow( r_corner_part_g_weights(crn_rid, 2), (1.0 - fabs(surf_vec0[2])) );
                          
                          
                          surf_rlid = 1;
                          
                          for(int dim = 0; dim < num_dim_; dim++){
                              ref_corner_surf_normals_(crn_rid, surf_rlid, dim) = surf_vec1[dim];
                          }

                          // power coef is =1 for the quadrature values for the surf normal else =0
                          ref_corner_surf_g_weights_(crn_rid, surf_rlid) =
                                pow( r_corner_part_g_weights(crn_rid, 0), (1.0 - fabs(surf_vec1[0])) )
                              * pow( r_corner_part_g_weights(crn_rid, 1), (1.0 - fabs(surf_vec1[1])) )
                              * pow( r_corner_part_g_weights(crn_rid, 2), (1.0 - fabs(surf_vec1[2])) );
                          
                          surf_rlid = 2;
                          
                          for (int dim = 0; dim < num_dim_; dim++){
                              ref_corner_surf_normals_(crn_rid, surf_rlid, dim) = surf_vec2[dim];
                          }
                          
                          // power coef is =1 for the quadrature values for the surf normal else =0
                          ref_corner_surf_g_weights_(crn_rid, surf_rlid) =
                                pow( r_corner_part_g_weights(crn_rid, 0), (1.0 - fabs(surf_vec2[0])) )
                              * pow( r_corner_part_g_weights(crn_rid, 1), (1.0 - fabs(surf_vec2[1])) )
                              * pow( r_corner_part_g_weights(crn_rid, 2), (1.0 - fabs(surf_vec2[2])) );
                          

                          index++;       // increment the index in the ref element
                          corner_rlid++; // increment the corner index in the ref cell
                          
                      }
                  }
              } // end of loop over corners in ref cell 
          }
      }
  } // end of loop over the ref cells in element


  /* 
   * END ref_element::init source up to call to HexN::setup_HexN
   */


  /*
   * BEGIN HexN::setup_HexN source
   */

  if (orderBasis == 0) {

      num_nodes_1d_ = 2;
      num_nodes_ = pow(num_nodes_1d_, 3);

      HexN_Nodes_1d_ = MatarRealCArray (num_nodes_);
      HexN_Nodes_ = MatarRealCArray (num_nodes_, 3);


      // Vertices
      num_verts_1d_ = 2;
      num_verts_ = pow(num_verts_1d_, 3);
      num_basis_ = pow(num_verts_1d_, 3);
      
      HexN_Verts_1d_ = MatarRealCArray (num_verts_);
      HexN_Verts_ = MatarRealCArray (num_verts_, 3);


      Vert_Node_map_ = MatarUIntCArray (num_verts_);

      order_ = orderBasis+1;


  }


  else{
      
      // Nodes
      num_nodes_1d_ = 2 * orderBasis + 1;
      num_nodes_ = pow(num_nodes_1d_, 3);

      HexN_Nodes_1d_ = MatarRealCArray (num_nodes_);
      HexN_Nodes_ = MatarRealCArray (num_nodes_, 3);


      // Vertices
      num_verts_1d_ = orderBasis + 1;
      num_verts_ = pow(num_verts_1d_, 3);
      num_basis_ = pow(num_verts_1d_, 3);
      
      HexN_Verts_1d_ = MatarRealCArray (num_verts_);
      HexN_Verts_ = MatarRealCArray (num_verts_, 3);


      Vert_Node_map_ = MatarUIntCArray (num_verts_);

      order_ = orderBasis;

  }
  
  create_lobatto_nodes(orderBasis);


  // Set the vertex to node map (every other node)
  if (orderBasis == 0) {

      int vert_rid = 0;
      for(int k = 0; k < num_nodes_1d_; k++){
          for(int j = 0; j < num_nodes_1d_; j++){
              for(int i = 0; i < num_nodes_1d_; i++){

                  int node_id = node_rid(i, j, k);
                  
                  Vert_Node_map_(vert_rid) = node_id;

                  vert_rid++;                        
              }   
          }
      }


  }


  if (orderBasis >= 1){

      int vert_rid = 0;
      for(int k = 0; k < num_nodes_1d_; k=k+2){
          for(int j = 0; j < num_nodes_1d_; j=j+2){
              for(int i = 0; i < num_nodes_1d_; i=i+2){

                  int node_id = node_rid(i, j, k);
                  Vert_Node_map_(vert_rid) = node_id;
                  vert_rid++;
              }   
          }
      }
  }


  /*
   * END HexN::setup_HexN source
   */


  /*
   * BEGIN ref_element::init source after call to HexN::setup_HexN
   */


  // --- evaluate the basis at the nodal positions
  for(int node_rlid = 0; node_rlid < num_ref_nodes_in_elem_; node_rlid++){

      auto point = MatarRealCArray(3);

      // Get the nodal coordinates
      for(int dim = 0; dim < 3; dim++){
          point(dim) = ref_node_positions_(node_rlid, dim);
      }

      auto node_basis = MatarRealCArray(num_ref_verts_in_elem_);

      evaluate_basis(node_basis, point);

      for(int vert_rlid = 0; vert_rlid < num_ref_verts_in_elem_; vert_rlid++){

          ref_nodal_basis_(node_rlid, vert_rlid) = node_basis(vert_rlid);
      }

  }


  // --- evaluate grad_basis functions at the ref nodes ---

  auto partial_xi  = MatarRealCArray(num_ref_nodes_in_elem_);
  auto partial_eta = MatarRealCArray(num_ref_nodes_in_elem_);
  auto partial_zeta  = MatarRealCArray(num_ref_nodes_in_elem_);

  

  for(int node_lid = 0; node_lid < num_ref_nodes_in_elem_; node_lid++){

      auto point = MatarRealCArray(3);
      
      // Get the nodal coordinates
      for(int dim = 0; dim < 3; dim++){
          point(dim) = ref_node_positions_(node_lid, dim);
      }

      evaluate_derivative_basis_xi(partial_xi, point);
      evaluate_derivative_basis_eta(partial_eta, point);
      evaluate_derivative_basis_zeta(partial_zeta, point);

      for(int basis_id = 0; basis_id < num_ref_verts_in_elem_; basis_id++){


          ref_nodal_gradient_(node_lid, basis_id, 0) = partial_xi(basis_id);
          ref_nodal_gradient_(node_lid, basis_id, 1) = partial_eta(basis_id);
          ref_nodal_gradient_(node_lid, basis_id, 2) = partial_zeta(basis_id);

          partial_xi(basis_id)  = 0.0;
          partial_eta(basis_id) = 0.0;
          partial_zeta(basis_id)  = 0.0;
      }
  }
  
  // build a local rlid for the interior and surface nodes of the ref element
  int count_surf = 0;
  int count_inside = 0;
  for (int k=0; k<num_ref_nodes_1D_; k++){
      for (int j=0; j<num_ref_nodes_1D_; j++){
          for (int i=0; i<num_ref_nodes_1D_; i++){
              
              int rid = node_rid(i,j,k); // follows i,j,k convention
              
              if (i>0 && i<num_ref_nodes_1D_-1 &&
                  j>0 && j<num_ref_nodes_1D_-1 &&
                  k>0 && k<num_ref_nodes_1D_-1)
              {
                  ref_inside_nodes_in_elem_(count_inside) = rid;
                  count_inside ++;
              }
              else
              {
                  ref_surface_nodes_in_elem_(count_surf) = rid;
                  count_surf ++;
              }
              
          } // end for
      } // end for
  } // end for
  // checking the data structure
  if (num_ref_inside_nodes_in_elem_ < count_inside) {
      std::cout << "*** Error in inside nodes in ref elem ***" << std::endl;
  }
  if (num_ref_surface_nodes_in_elem_ < count_surf) {
      std::cout << "*** Error in surface nodes in ref elem ***" << std::endl;
  }
  
  
  // ---- build reference patches and reference sides of a cell and assign guass values ----
  
  num_patches_in_elem_ = (num_ref_cells_1D_+1) * (num_ref_cells_1D_*num_ref_cells_1D_)*3;
  num_sides_in_elem_ = num_ref_cells_in_elem_*6;
  ref_patches_in_cell_list_ = MatarIntCArray (num_sides_in_elem_);
  
  
  // allocate memory for the 1D Gauss Legendre points and weights
  auto leg_points_1D = CArray <real_t> (num_ref_cells_1D_);
  legendre_nodes_1D(leg_points_1D, num_ref_cells_1D_);
  
  auto leg_weights_1D = CArray <real_t> (num_ref_cells_1D_);
  legendre_weights_1D(leg_weights_1D, num_ref_cells_1D_);
  
  
  // allocate memory for global variables related to Gauss Legendre points and weights
  ref_cell_positions_ = CArray <real_t> (num_ref_cells_in_elem_, 3);
  ref_cell_g_weights_ = CArray <real_t> (num_ref_cells_in_elem_); // this one
  ref_patch_positions_ = CArray <real_t> (num_patches_in_elem_, 3);
  ref_patch_g_weights_ = CArray <real_t> (num_patches_in_elem_);

  
  
  for (int side_rid = 0; side_rid < num_sides_in_elem_; side_rid++){
      ref_patches_in_cell_list_(side_rid) = -1;
  }
  
  int patch_rid = 0;  // the patch reference index
  // loop over the cells
  for (int k=0; k<num_ref_cells_1D_; k++){
      for (int j=0; j<num_ref_cells_1D_; j++){
          for (int i=0; i<num_ref_cells_1D_; i++){
              
              // get the cell_rid
              int this_cell_rid = cell_rid(i, j, k);
              
              // save the gauss point volume values at cell
              // i,j,k tensor product for gauss Legendre
              ref_cell_positions_(this_cell_rid,0) = leg_points_1D(i);
              ref_cell_positions_(this_cell_rid,1) = leg_points_1D(j);
              ref_cell_positions_(this_cell_rid,2) = leg_points_1D(k);
              
              ref_cell_g_weights_(this_cell_rid) = leg_weights_1D(i)*leg_weights_1D(j)*leg_weights_1D(k);
              
              
              // loop over the sides in this cell
              for (int side_lid = 0; side_lid < 6; side_lid++){
                  
                  // get the local side index in the cell
                  int side_rid = side_lid + this_cell_rid*6; // side index
                  
                  
                  // check to see if a patch index was already saved
                  if (ref_patches_in_cell_list_(side_rid) == -1){
                      
                      // save the patch_rid for this patch in the cell
                      ref_patches_in_cell_list_(side_rid) = patch_rid;
                      
                      // save the gauss point area vaules on the patch
                      if (side_lid == 0){
                          
                          ref_patch_positions_(patch_rid,0) = lob_nodes_1D(i);  // xi-coord is the labatto point value at i
                          ref_patch_positions_(patch_rid,1) = leg_points_1D(j); // legendre points
                          ref_patch_positions_(patch_rid,2) = leg_points_1D(k); // legendre points
                          
                          // j,k tensor product for this surface
                          ref_patch_g_weights_(patch_rid) = leg_weights_1D(j)*leg_weights_1D(k);
                      }
                      else if (side_lid == 1){
                          
                          ref_patch_positions_(patch_rid,0) = lob_nodes_1D(i+1); // xi-coord is the labatto point value at i+1
                          ref_patch_positions_(patch_rid,1) = leg_points_1D(j);  // legendre points
                          ref_patch_positions_(patch_rid,2) = leg_points_1D(k);  // legendre points
                          
                          // j,k tensor product for this surface
                          ref_patch_g_weights_(patch_rid) = leg_weights_1D(j)*leg_weights_1D(k);
                      }
                      // end xi direction patches
                      //
                      else if (side_lid == 2){
                          
                          ref_patch_positions_(patch_rid,0) = leg_points_1D(i); // legendre points
                          ref_patch_positions_(patch_rid,1) = lob_nodes_1D(j);  // eta-coord is the labatto point value at j
                          ref_patch_positions_(patch_rid,2) = leg_points_1D(k); // legendre points
                          
                          // i,k tensor product for this surface
                          ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(k);
                          
                      }
                      else if(side_lid == 3){
                          
                          ref_patch_positions_(patch_rid,0) = leg_points_1D(i);  // legendre points
                          ref_patch_positions_(patch_rid,1) = lob_nodes_1D(j+1); // eta-coord is the labatto point value at j+1
                          ref_patch_positions_(patch_rid,2) = leg_points_1D(k);  // legendre points
                          
                          // i,k tensor product for this surface
                          ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(k);
                      }
                      // end eta direction patches
                      //
                      else if (side_lid == 4){
                          
                          ref_patch_positions_(patch_rid,0) = leg_points_1D(i); // legendre points
                          ref_patch_positions_(patch_rid,1) = leg_points_1D(j); // legendre points
                          ref_patch_positions_(patch_rid,2) = lob_nodes_1D(k);  // zeta-coord is the labatto point value at k
                          
                          // i,j tensor product for this surface
                          ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(j);
                      }
                      else if (side_lid == 5){
                          
                          ref_patch_positions_(patch_rid,0) = leg_points_1D(i);  // legendre points
                          ref_patch_positions_(patch_rid,1) = leg_points_1D(j);  // legendre points
                          ref_patch_positions_(patch_rid,2) = lob_nodes_1D(k+1); // zeta-coord is the labatto point value at k
                          
                          // i,j tensor product for this surface
                          ref_patch_g_weights_(patch_rid) = leg_weights_1D(i)*leg_weights_1D(j);
                      }
                      // end zeta direction patches
                           
                               
                      // also save the patch_rid in the neighboring cell for this patch
                      // select neighbors for +/- xi, eta, and zeta directions
                      int neighbor_cell_rid;
                      if (side_lid == 0 && i>0){
                          neighbor_cell_rid = cell_rid(i-1, j, k);
                      } // end xi minus direction
                      else if (side_lid == 1 && i<num_ref_cells_1D_-1){
                          neighbor_cell_rid = cell_rid(i+1, j, k);
                      } // end xi plus direction
                      else if (side_lid == 2 && j>0){
                          neighbor_cell_rid = cell_rid(i, j-1, k);
                      } // end eta minus direction
                      else if (side_lid == 3 && j<num_ref_cells_1D_-1){
                          neighbor_cell_rid = cell_rid(i, j+1, k);
                      } // end eta plus direction
                      else if (side_lid == 4 && k>0){
                          neighbor_cell_rid = cell_rid(i, j, k-1);
                      } // end zeta minus direction
                      else if (side_lid == 5 && k<num_ref_cells_1D_-1){
                          neighbor_cell_rid = cell_rid(i, j, k+1);
                      } // end zeta plus direction
                      else{
                          neighbor_cell_rid = -1;
                      } // end no neighboring cell
                      
                      // if the neighbor is valid, then save the patch_rid
                      if (neighbor_cell_rid >= 0 && neighbor_cell_rid < num_ref_cells_in_elem_){
                          
                          // get the patch_lid in the neighboring cell for this patch_lid
                          int neighbor_side_lid = patch_rlid_in_cell_neighbor(side_lid);
                          
                          // get the side index of cell for the same patch
                          int neighbor_side_rid = neighbor_side_lid + (neighbor_cell_rid)*6; // neighboring side index
                          
                          // save the patch_rid for this patch in the neighboring cell
                          ref_patches_in_cell_list_(neighbor_side_rid) = patch_rid;
                          
                      } // end if neighboring_cell_rid is valid
                      
                      
                      // increment the patch_rid index
                      patch_rid ++;
                      
                  }  // end if the patch was saved
                  
              } // end for sides in a cell
              
          } // end for i-dir cells
      } // end for j-dir cells
  } // end for k-dir cells
  
  std::cout << "  patch_rid = " << patch_rid << "  num_patches = " << num_patches_in_elem_ << "\n";
  std::cout << "  num cells * 6 = " << num_ref_cells_in_elem_*6 << " and num_cells_1D = "  << num_ref_cells_1D_ << std::endl;
  
  
  // --- evaluate the Lagrange basis at the patch positions
  ref_patch_basis_ =  CArray <real_t> (num_patches_in_elem_, num_ref_verts_in_elem_);
  
  for (int patch_rlid = 0; patch_rlid < num_patches_in_elem_; patch_rlid++){
      
      auto point = MatarRealCArray(3);
      
      // Get the patch coordinates
      for(int dim = 0; dim < 3; dim++){
          point(dim) = ref_patch_positions(patch_rlid, dim);
      }
      
      // the basis function values at the patch for each vertex
      auto patch_basis = MatarRealCArray(num_ref_verts_in_elem_);
      
      // calculate the patch basis function values at the point, for each vertex
      evaluate_basis(patch_basis, point);
      
      // save the basis values at the patch for each vertex
      for(int vert_rlid = 0; vert_rlid < num_ref_verts_in_elem_; vert_rlid++){
          ref_patch_basis_(patch_rlid, vert_rlid) = patch_basis(vert_rlid);
      }
      
  } // end for patch_rlid
  
  
  
  // --- evaluate grad_basis functions at the ref patches ---
  ref_patch_gradient_ = CArray <real_t> (num_patches_in_elem_, num_basis_, num_dim_);
  
  // loop over the patches
  for (int patch_rlid = 0; patch_rlid < num_patches_in_elem_; patch_rlid++){
      
      auto point = MatarRealCArray(3);
      
      // Get the patch coordinates
      for(int dim = 0; dim < 3; dim++){
          point(dim) = ref_patch_positions(patch_rlid, dim);
      }
      
      // calculate the partials at the patch location for each vertex
      evaluate_derivative_basis_xi(partial_xi, point);
      evaluate_derivative_basis_eta(partial_eta, point);
      evaluate_derivative_basis_zeta(partial_zeta, point);
      
      // loop over the basis polynomials where there is one from each vertex
      for(int basis_id = 0; basis_id < num_ref_verts_in_elem_; basis_id++){
          
          ref_patch_gradient_(patch_rlid, basis_id, 0) = partial_xi(basis_id);
          ref_patch_gradient_(patch_rlid, basis_id, 1) = partial_eta(basis_id);
          ref_patch_gradient_(patch_rlid, basis_id, 2) = partial_zeta(basis_id);
          
          partial_xi(basis_id)  = 0.0;
          partial_eta(basis_id) = 0.0;
          partial_zeta(basis_id)  = 0.0;
      } // end for basis_id
      
  } // end for patch_rlid
  
  
  
  // --- evaluate the Lagrange basis at the cell positions
  ref_cell_basis_ = CArray <real_t> (num_ref_cells_in_elem_, num_basis_, num_dim_);
  
  for (int cell_rlid = 0; cell_rlid < num_ref_cells_in_elem_; cell_rlid++){
      
      auto point = MatarRealCArray(3);
      
      // Get the patch coordinates
      for(int dim = 0; dim < 3; dim++){
          point(dim) = ref_cell_positions(cell_rlid, dim);
      }
      
      // the basis function values at the patch for each vertex
      auto cell_basis = MatarRealCArray(num_ref_verts_in_elem_);
      
      // calculate the cell basis function values at the point, for each vertex
      evaluate_basis(cell_basis, point);
      
      // save the basis values at the patch for each vertex
      for(int vert_rlid = 0; vert_rlid < num_ref_verts_in_elem_; vert_rlid++){
          ref_cell_basis_(cell_rlid, vert_rlid) = cell_basis(vert_rlid);
      }
      
  } // end for cell_rlid
  
  
  // --- evaluate grad_basis functions at the reference cell ---
  ref_cell_gradient_ = CArray <real_t> (num_ref_cells_in_elem_, num_basis_, num_dim_);
  
  // loop over the patches
  for (int cell_rlid = 0; cell_rlid < num_ref_cells_in_elem_; cell_rlid++){
      
      auto point = MatarRealCArray(3);
      
      // Get the cell coordinates
      for(int dim = 0; dim < 3; dim++){
          point(dim) = ref_cell_positions_(cell_rlid, dim);
      }
      
      // calculate the partials at the cell location for each vertex
      evaluate_derivative_basis_xi(partial_xi, point);
      evaluate_derivative_basis_eta(partial_eta, point);
      evaluate_derivative_basis_zeta(partial_zeta, point);
      
      // loop over the basis polynomials where there is one from each vertex
      for(int basis_id = 0; basis_id < num_ref_verts_in_elem_; basis_id++){
          
          ref_cell_gradient_(cell_rlid, basis_id, 0) = partial_xi(basis_id);
          ref_cell_gradient_(cell_rlid, basis_id, 1) = partial_eta(basis_id);
          ref_cell_gradient_(cell_rlid, basis_id, 2) = partial_zeta(basis_id);
          
          partial_xi(basis_id)  = 0.0;
          partial_eta(basis_id) = 0.0;
          partial_zeta(basis_id)  = 0.0;
      } // end for basis_id
      
  } // end for patch_rlid


  /*
   * END ref_element::init source after call to HexN::setup_HexN
   */
        
}


HexRef::evaluate_basis(MatarRealCArray &basis, MatarRealCArray &point) {
  auto val_1d = MatarRealCArray (num_verts_1d_);
  auto val_3d = MatarRealCArray (num_verts_1d_, 3);

  // Calculate 1D basis for the X coordinate of the point
  lagrange_basis_1D(val_1d, point(0));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      val_3d(i, 0) = val_1d(i);
      val_1d(i) = 0.0;
  }

  // Calculate 1D basis for the Y coordinate of the point
  lagrange_basis_1D(val_1d, point(1));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      val_3d(i, 1) = val_1d(i);
      val_1d(i) = 0.0;
  }

  // Calculate 1D basis for the Z coordinate of the point
  lagrange_basis_1D(val_1d, point(2));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      val_3d(i, 2) = val_1d(i);
      val_1d(i) = 0.0;
  }
  
  // Multiply the i, j, k components of the basis from each node
  // to get the tensor product basis for the node
  for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
          for(int i = 0; i < num_verts_1d_; i++){

              int vert_rlid = vert_rid(i,j,k);
              basis(vert_rlid) = val_3d(i, 0)*val_3d(j, 1)*val_3d(k, 2);
          }
      }
  }
}


HexRef::evaluate_derivative_basis_xi(MatarRealCArray &partial_xi, 
    MatarRealCArray &point) {
  auto val_1d = MatarRealCArray (num_verts_1d_);
  auto val_3d = MatarRealCArray (num_verts_1d_, 3);

  auto Dval_1d = MatarRealCArray (num_verts_1d_);
  auto Dval_3d = MatarRealCArray (num_verts_1d_, 3);

  // Calculate 1D partial w.r.t. xi for the X coordinate of the point
  lagrange_derivative_1D(Dval_1d, point(0));


  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      Dval_3d(i,0) = Dval_1d(i);
      Dval_1d(i) = 0.0;
  }


  // Calculate 1D basis for the Y coordinate of the point
  lagrange_basis_1D(val_1d, point(1));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      val_3d(i,1) = val_1d(i);
      val_1d(i) = 0.0;
  }


  // Calculate 1D basis for the Z coordinate of the point
  lagrange_basis_1D(val_1d, point(2));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      val_3d(i,2) = val_1d(i);
      val_1d(i) = 0.0;
  }

  // Multiply the i, j, k components of the basis and partial_xi from each node
  // to get the tensor product partial derivatives of the basis at each node
  for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
          for(int i = 0; i < num_verts_1d_; i++){
              
              int vert_rlid = vert_rid(i,j,k);

              // Partial w.r.t xi
              partial_xi(vert_rlid) = Dval_3d(i, 0)*val_3d(j, 1)*val_3d(k, 2);

          }
      }
  }
}


HexRef::evaluate_derivative_basis_eta(MatarRealCArray &partial_xi, 
    MatarRealCArray &point) {   
  auto val_1d = MatarRealCArray (num_verts_1d_);
  auto val_3d = MatarRealCArray (num_verts_1d_, 3);

  auto Dval_1d = MatarRealCArray (num_verts_1d_);
  auto Dval_3d = MatarRealCArray (num_verts_1d_, 3);

  // Calculate 1D basis for the Y coordinate of the point
  lagrange_basis_1D(val_1d, point(0));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      val_3d(i,0) = val_1d(i);
      val_1d(i) = 0.0;
  }


  // Calculate 1D partial w.r.t. eta for the Y coordinate of the point
  lagrange_derivative_1D(Dval_1d, point(1));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      Dval_3d(i,1) = Dval_1d(i);
      Dval_1d(i) = 0.0;
  }


  // Calculate 1D basis for the Z coordinate of the point
  lagrange_basis_1D(val_1d, point(2));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      val_3d(i,2) = val_1d(i);
      val_1d(i) = 0.0;
  }

  // Multiply the i, j, k components of the basis and partial_eta from each node
  // to get the tensor product partial derivatives of the basis at each node
  for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
          for(int i = 0; i < num_verts_1d_; i++){
              
              int vert_rlid = vert_rid(i,j,k);

              // Partial w.r.t xi
              partial_eta(vert_rlid) = val_3d(i, 0)*Dval_3d(j, 1)*val_3d(k, 2);

          }
      }
  }
}


HexRef::evaluate_derivative_basis_zeta(MatarRealCArray &partial_zeta, 
    MatarRealCArray &point) {
  auto val_1d = MatarRealCArray (num_verts_1d_);
  auto val_3d = MatarRealCArray (num_verts_1d_, 3);

  auto Dval_1d = MatarRealCArray (num_verts_1d_);
  auto Dval_3d = MatarRealCArray (num_verts_1d_, 3);

  // Calculate 1D basis for the X coordinate of the point
  lagrange_basis_1D(val_1d, point(0));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      val_3d(i,0) = val_1d(i);
      val_1d(i) = 0.0;
  }


  // Calculate 1D basis for the Y coordinate of the point
  lagrange_basis_1D(val_1d, point(1));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      val_3d(i,1) = val_1d(i);
      val_1d(i) = 0.0;
  }


  // Calculate 1D partial w.r.t. zeta for the Z coordinate of the point
  lagrange_derivative_1D(Dval_1d, point(2));
  
  // Save the basis value at the point to a temp array and zero out the temp array
  for(int i = 0; i < num_verts_1d_; i++){
      
      Dval_3d(i,2) = Dval_1d(i);
      val_1d(i) = 0.0;
  }

  // Multiply the i, j, k components of the basis and partial_xi from each node
  // to get the tensor product partial derivatives of the basis at each node
  for(int k = 0; k < num_verts_1d_; k++){
      for(int j = 0; j < num_verts_1d_; j++){
          for(int i = 0; i < num_verts_1d_; i++){
              
              int vert_rlid = vert_rid(i,j,k);

              // Partial w.r.t zeta
              partial_zeta(vert_rlid) = val_3d(i, 0)*val_3d(j, 1)*Dval_3d(k, 2);

          }
      }
  }
}
