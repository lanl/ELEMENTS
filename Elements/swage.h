#ifndef SWAGE_H
#define SWAGE_H 


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "utilities.h"
#include "matar.h"
#include "Simulation_Parameters.h"


using namespace utils;

// REMOVE RK STUFF


namespace swage{

/*
==============================
          Mesh Class
==============================

==========================
Representative Local Cell 
==========================

              K
              ^         J
              |        /
              |       /
              |      /
      6------------------7
     /|                 /|
    / |                / |
   /  |               /  |
  /   |              /   | 
 /    |             /    |
4------------------5     |
|     |            |     | ----> I
|     |            |     |  
|     |            |     |
|     |            |     |
|     2------------|-----3
|    /             |    /
|   /              |   /
|  /               |  /
| /                | /         
|/                 |/
0------------------1


patch 0: [0,1,3,2]
patch 1: [4,5,7,6]
patch 2: [0,1,5,4]
patch 3: [2,3,7,6]
patch 4: [0,2,6,4]
patch 6: [1,3,7,5]

*/


class mesh_t {
public:
    // constructor
    mesh_t (Simulation_Parameters *simparam);
    // deconstructor
    ~mesh_t ( );
    
    //flags
    int local_index_set_flag; //determines whether to store local indices along with node to cell connectivity
    
    CArray <int> nodes_in_cell_list_;
    CArray <int> num_cells_in_node_;
    CArray <int> cells_in_node_list_start_;
    CArray <int> cells_in_node_list_;
    CArray <int> local_index_list_;         //stores the element-local indices of all nodes (for each element the node connects to) for ease of access later.

private:
    
    // ---- GENERAL INFORMATION ---- //
    int num_dim_;
    const int num_nodes_hex_ = 8;
    const int num_patches_hex_ = 6;
    
    const int num_nodes_patch_ = 4;
    const int this_node_in_patch_in_cell_[24] = // this assumes i,j,k structures nodes
        {0,1,3,2,
         4,5,7,6,
         0,1,5,4,
         2,3,7,6,
         0,2,6,4,
         1,3,7,5};

    int indx_; //useful for returning from internal function


// ---- ELEMENT POLYNOMIAL ORDER ---- //
    int   elem_order_;


// ---- INDEX SPACES AND MAPS ---- //

    // ---- ELEMENT ---- //
    int   num_elem_;
    int   num_g_pts_in_elem_;
    int   num_cells_in_elem_;
    int   num_nodes_in_elem_;
    int   num_mat_pts_in_elem_;

    CArray <int> cells_in_elem_;

    CArray <int> num_elems_in_elem_;
    CArray <int> elems_in_elem_list_start_;
    RaggedRightArray <int> elems_in_elem_list_;
    CArray <int> nodes_in_elem_list_;


    // ---- CELLS ---- //
    int   num_cells_;
    
    CArray <int> num_cells_in_cell_;
    CArray <int> cells_in_cell_list_start_;
    CArray <int> cells_in_cell_list_;
    CArray <int> elems_in_cell_list_;


    // ---- VERTICES ---- //



    // ---- NODES ---- //
    int   num_nodes_;

    CArray <int> num_elems_in_node_;
    CArray <int> elems_in_node_list_start_;
    CArray <int> elems_in_node_list_;

    // ---- GAUSS POINTS ---- //
    int   num_g_pts_;

    CArray <int> node_in_gauss_list_;


    // ---- CORNERS ---- //
    int   num_corners_;

    CArray <int> num_corners_in_node_;
    CArray <int> corners_in_cell_list_;
    CArray <int> corners_in_node_list_start_;
    CArray <int> corners_in_node_list_;

    // int * corner_bdy_count_ = NULL;


    // ---- PATCHES ---- //
    int   num_patches_;

    CArray <int> patch_nodes_list_;
    CArray <int> patch_local_node_list;
    CArray <int> cells_in_patch_list_;
    CArray <int> cells_in_patch_local_id_list_;  //stores the local patch id for the corresponding cell


    // ---- BOUNDARY ---- //
    int num_bdy_patches_;
    int num_bdy_sets_;
    
    CArray <int> bdy_patches_;
    CArray <int> bdy_set_list_;
    CArray <int> start_index_bdy_set_;
    CArray <int> num_bdy_patches_set_;


// ---- MESH GEOMETRIC STATE ---- //

    // ---- ELEMENT ---- //
    CArray <real_t> elem_vol_;


    // ---- CELLS ---- //
    CArray <real_t> cell_vol_;
    CArray <real_t> cell_coords_;


    // ---- NODES ---- //
    CArray <real_t> node_coords_;
    CArray <real_t> node_jacobian_;
    CArray <real_t> node_jacobian_inverse_;
    CArray <real_t> node_det_j_; 

    // ---- GHOST NODE INFO ---- //
    CArray <real_t> ghost_node_coords_;

    // ---- QUADRATURE POINTS ---- //
    CArray <real_t> jacobians_;
    CArray <real_t> jacobian_determinant_;


public:

    //**********************************//
    // Mesh class function definitions  //
    //**********************************//

    void init_element (int e_order, int dim, int num_elem);
    void init_cells (int ncells);
    void init_nodes (int num_nodes);
    void init_ghost_nodes (int num_ghost_nodes);
    void init_gauss_pts ();
    void init_bdy_sets (int num_sets);

    // ==== MESH CONSTANTS ==== // 

    // returns the number of dimensions in the mesh
    int num_dim () const;

    // returns the polynomial order of the element
    int elem_order () const;


    // ==== INDEX SPACE ACCESSORS ==== //

    // ---- ELEMENT ---- //

    // returns the number of elements
    int num_elems () const;

    // returns the number of elements
    int num_elems_in_elem (int elem_gid) const;

    // returns the number of elements (WARNING: currently assumes constant size)
    int num_cells_in_elem () const;

    int num_nodes_in_elem () const;
    // returns the nodes in an element
    int& nodes_in_elem (int elem_gid, int node_lid);

    // return array of elements connected to element (corners+patches)
    int& elems_in_elem (int elem_gid, int elem_lid);

    // return the the global cell id from local element cell id
    int& cells_in_elem (int elem_gid, int cell_lid);

    // return number of gauss points in an element (currently assumes Gauss-Lobatto)
    int& num_gauss_in_elem ();

    // return number of material points in an element
    int& num_mat_pt_in_elem ();




    // ---- CELLS ---- //

    // returns the number of cells
    int num_cells () const;

    // return the node ids local to the cell
    int num_nodes_in_cell () const;

    // return the node ids local to the cell
    int& nodes_in_cell (int cell_gid, int node_lid) const;

    // return the number of cells around the cell
    int& num_cells_in_cell (int cell_gid) const;

    // return the the cells around a cell
    int& cells_in_cell (int cell_gid, int cell_lid) const;

    // return corners connected to a cell
    int& corners_in_cell (int cell_gid, int corner_lid) const;

    // return the element this cell belongs to
    int& elems_in_cell (int cell_gid) const;


    // ---- VERTICES ---- //


    // ---- NODES ---- //

    // returns the number of nodes
    int num_nodes ();

    // returns number of cells around a node
    int& num_cells_in_node (int node_gid) const;

    // returns number of elements around a node
    int& num_elems_in_node (int node_gid) const;

    // return the cells around a node
    int& cells_in_node (int node_gid, int cell_lid) const;

    // return the elements around a node
    int& elems_in_node (int node_gid, int elem_lid) const;

    // return the Jacobian at a node
    real_t & node_jacobian (int node_gid, int dim_i, int dim_j) const;

    // return the Jacobian inverse at a node
    real_t & node_jacobian_inv (int node_gid, int dim_i, int dim_j) const;


    // return the determinant of the Jacobian at a node
    real_t & node_det_j (int node_gid) const;

    // ---- GAUSS POINTS ---- //

    // return number of gauss points in mesh
    int num_gauss () const;

    // return gauss to node map
    int& node_in_gauss (int gauss_gid) const;

    // return gauss in element map (internal structured grid)
    int& gauss_in_elem (int elem_gid, int gauss_lid); 


    // ---- CORNERS ---- //
        
    // returns the number of corners
    int num_corners () const;

    // return number of corners connected to a node
    int num_corners_in_node (int node_gid) const;

    // return corner to node map
    int corners_in_node (int node_gid, int corner_lid) const;


    // ---- PATCHES ---- //

    // returns the number of elements
    int num_patches () const;

    // returns the global node id given a cell_id, local_patch_indx(0:5), local_patchnode_indx(0:3)
    int node_in_patch_in_cell(int cell_id, int this_patch, int patchnode_lid) const;

    // returns the global id for a cell that is connected to the patch
    int cells_in_patch(int patch_gid, int this_cell) const;

    // returns the local patch id for the cell and global patch id pair
    int cells_in_patch_local_id(int patch_gid, int this_cell) const;
          
    // returns the nodes in the patch
    int node_in_patch(int patch_gid, int patchnode_lid) const;


    // ---- Boundary ---- //

    int num_bdy_sets() const;

    int num_bdy_patches() const;

    int bdy_patches(int this_bdy_patch) const;

    // returns the number per bdy-patches in a particular set
    int num_bdy_patches_in_set (int bdy_set);

    // returns a subset of the boundary patches
    int bdy_patches_in_set (int bdy_set, int this_patch);



    // ==== MESH STATE FUNCTIONS ==== // 


    // ---- ELEMENTS ---- //
    real_t& elem_vol(int elem_gid) const;



    // ---- CELLS ---- //

    // return the cell volume
    real_t& cell_vol(int cell_gid) const;

    // return the cell coordinate position
    real_t& cell_coords(int cell_gid, int this_dim);



    // ---- VERTICES ---- //



    // ---- NODES ---- //
    // return the node coordinates
    real_t& node_coords(int node_gid, int this_dim) const;

    // return the node coordinates
    real_t& ghost_node_coords(int node_gid, int this_dim) const;


    // ---- QUADRATURE POINTS ---- //

    // return jacobian at quadrature point
    real_t& jacobian(int gauss_gid, int i, int j) const;

    // return determinant of jacobian at quadrature point
    real_t& det_j(int gauss_gid) const;




    // ---- CORNERS ---- //


    // ---- PATCHES ---- //
    // geometric average patch coordinate
    real_t patch_coords(int patch_gid, int this_dim) const;


    // ==== MESH CONNECTIVITY FUNCTIONS ==== // 
        
    // initialize array for mesh connectivity: all cells around a node
    void build_connectivity();

    void build_node_cell_connectivity();

    void build_corner_connectivity();

    void build_cell_cell_connectivity();

    void build_patch_connectivity();

    void build_element_connectivity();

    // identify the boundary patches
    void build_bdy_patches ();


    // ---- bdy sets ----

    // returns a subset of the boundary patches
    int set_bdy_patches (int bdy_set, int patch_lid);

    // set planes for tagging sub sets of boundary patches
    // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    // val = plane value, radius, radius
    void tag_bdys(int this_bc_tag, real_t val, int bdy_set);

    // compress the bdy_set_list to reduce the memory
    void compress_bdy_set();

    // routine for checking to see if a vertix is on a boundary
    // bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
    // val = plane value, radius, radius
    int check_bdy(int patch_gid, int this_bc_tag, real_t val);


}; // end of mesh_t declaration


void refine_mesh(
    mesh_t& init_mesh, 
    mesh_t& mesh, 
    const int p_order,
    const int dim);

} // end namespace swage


#endif // SWAGE_H 