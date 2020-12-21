/*****************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and 
to permit others to do so.


This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
    
    1.  Redistributions of source code must retain the above copyright notice, this list of 
        conditions and the following disclaimer.
 
    2.  Redistributions in binary form must reproduce the above copyright notice, this list of 
        conditions and the following disclaimer in the documentation and/or other materials 
        provided with the distribution.
 
    3.  Neither the name of the copyright holder nor the names of its contributors may be used 
        to endorse or promote products derived from this software without specific prior 
        written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************************/


// NOTE: Rename patch to patch
// NOTE: Fence Cells and Elements

#include <iostream>  // std::cout etc.
#include <cmath>

#include "swage.h"

namespace swage{
//******************************//
// Mesh_t function definitions  //
//******************************//


// ==== MESH CONSTANTS ==== // 
 

// returns the number of dimensions in the mesh
int mesh_t::num_dim () const
{
    return num_dim_;
}

// returns the polynomial order of the element
int mesh_t::elem_order () const
{
    return elem_order_;
}

// returns the polynomial order of the CCH reconstruction
// int& mesh_t::recon_order ()
// {
//     return recon_order_;
// }

// ==== INDEX SPACE INITIALIZATIONS ==== //

// ---- ELEMENT ---- //
void mesh_t::init_element (int e_order, int dim, int num_elem){
    
    elem_order_ = e_order;

    int num_g_pts_1d;
    int num_g_pts;
    int num_subcells_per_elem;
    
    if(e_order == 0){
        
        num_g_pts_1d = 2;
        num_g_pts  = pow(num_g_pts_1d, dim);
        num_subcells_per_elem = pow((num_g_pts_1d - 1), dim);
    }

    else{

        num_g_pts_1d = 2 * e_order + 1;
        num_g_pts    = pow(num_g_pts_1d, dim);
        num_subcells_per_elem = pow((2 * e_order), dim);
    }


    num_elem_ = num_elem;

    num_g_pts_in_elem_ = num_g_pts;

    num_mat_pts_in_elem_ = 1;


    num_cells_in_elem_ = num_subcells_per_elem;

    num_cells_ = num_elem * num_subcells_per_elem;

    //DANcells_in_elem_ = new int[num_elem * num_subcells_per_elem]();
    cells_in_elem_ = CArray <int> (num_elem, num_subcells_per_elem);

    //DANelem_vol_ = new real_t[num_elem]();
    elem_vol_ = CArray <real_t> (num_elem);

    // WARNING: FOLLOWING CODE ASSUMES LOBATTO 
    num_nodes_in_elem_ = num_g_pts;

    //DANnodes_in_elem_list_ = new int[num_elem_ * num_g_pts_in_elem_]();
    nodes_in_elem_list_ = CArray <int> (num_elem_, num_g_pts_in_elem_);

    //DANnum_elems_in_elem_ = new int[num_elem_ ]();
    num_elems_in_elem_ = CArray <int> (num_elem_);

    //DANelems_in_elem_list_start_ = new int[num_elem_ + 1]();
    elems_in_elem_list_start_ = CArray <int> (num_elem_ + 1);
}

// ---- CELLS ---- //
void mesh_t::init_cells (int ncells){

    
    num_cells_ = ncells;

    //DANcell_vol_  = new real_t[num_cells_]();
    cell_vol_ = CArray <real_t> (num_cells_);
    //DANcell_coords_ = new real_t[num_cells_*num_dim_]();
    cell_coords_ = CArray <real_t> (num_cells_, num_dim_);


    //DANnodes_in_cell_list_   = new int[num_cells_*num_nodes_hex_]();
    nodes_in_cell_list_ = CArray <int> (num_cells_, num_nodes_hex_);
    //DANcorners_in_cell_list_ = new int[num_cells_*num_nodes_hex_]();
    corners_in_cell_list_ = CArray <int> (num_cells_, num_nodes_hex_);
    
    //DANnum_cells_in_cell_       = new int[num_cells_](); 
    num_cells_in_cell_ = CArray <int> (num_cells_);
    //DANcells_in_cell_list_start_ = new int[num_cells_+1]();
    cells_in_cell_list_start_ = CArray <int> (num_cells_ + 1);

    //DANelems_in_cell_list_ = new int[num_cells_]();
    elems_in_cell_list_ = CArray <int> (num_cells_);
}


// ---- VERTICES ---- //

// ---- NODES ---- //
void mesh_t::init_nodes (int num_nodes) {
  
    num_nodes_ = num_nodes;
    
    //DANnode_coords_   = new real_t[num_nodes_*num_dim_]();
    node_coords_   = CArray <real_t> (num_nodes_, num_dim_);

    //DANnum_cells_in_node_        = new int[num_nodes_]();
    num_cells_in_node_ = CArray <int> (num_nodes_);
    //DANcells_in_node_list_start_ = new int[num_nodes_+1]();
    cells_in_node_list_start_ = CArray <int> (num_nodes_ + 1);

    //DANnum_corners_in_node_      = new int[num_nodes_]();
    num_corners_in_node_ = CArray <int> (num_nodes_);

    //DANnum_elems_in_node_        = new int[num_nodes_]();
    num_elems_in_node_   = CArray <int> (num_nodes_);
    //DANelems_in_node_list_start_ = new int[num_nodes_+1]();
    elems_in_node_list_start_ = CArray <int> (num_nodes_ + 1);

    //DANnode_jacobian_  = new real_t[num_nodes_ * num_dim_ * num_dim_]();
    node_jacobian_ = CArray <real_t> (num_nodes_, num_dim_, num_dim_);
    
    node_jacobian_inverse_ = CArray <real_t> (num_nodes_, num_dim_, num_dim_);

    //DANnode_det_j_  = new real_t[num_nodes_]();
    node_det_j_ = CArray <real_t> (num_nodes_);
}


// ---- GAUSS LOBATTO POINTS ---- //
void mesh_t::init_gauss_pts (){

    // Index maps
    num_g_pts_ = num_elem_ * num_g_pts_in_elem_;
    //DANnode_in_gauss_list_ = new int[num_g_pts_];
    node_in_gauss_list_ = CArray <int> (num_g_pts_);

    // geometric state
    //DANjacobians_ = new real_t[num_g_pts_*num_dim_*num_dim_];
    jacobians_ = CArray <real_t> (num_g_pts_, num_dim_, num_dim_);
    //DANjacobian_determinant_ = new real_t[num_g_pts_];
    jacobian_determinant_ = CArray <real_t> (num_g_pts_);
}


// ---- CORNERS ---- //


// ---- PATCHES ---- //


// ---- BOUNDARY ---- //

// initializes the number of bdy sets
void mesh_t::init_bdy_sets (int num_sets){
    
    // A check
    if(num_sets == 0){
        std::cout << " ERROR: number of boundary sets = 0, setting it = 1";
        num_sets = 1;
    }
    num_bdy_sets_ = num_sets;
    //DANnum_bdy_patches_set_ = new int [num_sets];
    num_bdy_patches_set_ = CArray <int> (num_sets);
    //DANstart_index_bdy_set_ = new int [num_sets+1];
    start_index_bdy_set_ = CArray <int> (num_sets + 1);
    //DANbdy_set_list_   = new int [num_sets*num_bdy_patches_]; // largest size possible
    bdy_set_list_ = CArray <int> (num_sets, num_bdy_patches_);
}
    

// ==== INDEX SPACE ACCESSORS ==== //

// ---- ELEMENT ---- //

// returns the number of elements
int mesh_t::num_elems () const
{
    return num_elem_;
}

// returns the number of elements
int mesh_t::num_elems_in_elem (int elem_gid) const
{
    return num_elems_in_elem_(elem_gid);
}

// returns the number of nodes in elements
int mesh_t::num_nodes_in_elem () const
{
    return num_nodes_in_elem_;
}

// returns the number of elements (WARNING: currently assumes constant size)
int mesh_t::num_cells_in_elem () const
{
    return num_cells_in_elem_;
}

// returns the nodes in an element
int& mesh_t::nodes_in_elem (int elem_gid, int node_lid)
{
    //DANreturn nodes_in_elem_list_[elem_gid * num_g_pts_in_elem_ + node_lid];
    return nodes_in_elem_list_(elem_gid, node_lid);
} 

// return array of elements connected to element (corners+patches)
int& mesh_t::elems_in_elem (int elem_gid, int elem_lid) 
{
    // shift index by 1 so that it is consistent with matrix syntax
    //DANint start_indx = elems_in_elem_list_start_(elem_gid);
    
    // get the index in the global 1D array
    //DANint index = start_indx + elem_lid;
    
    //DANreturn elems_in_elem_list_(index);
    return elems_in_elem_list_(elem_gid, elem_lid);
}

// return the the global cell id from local element cell id
int& mesh_t::cells_in_elem (int elem_gid, int cell_lid) 
{
    //DANreturn cells_in_elem_[cell_lid + num_cells_in_elem_*(elem_gid)];
    return cells_in_elem_(elem_gid, cell_lid);
}

// return number of gauss points in an element (currently assumes Gauss-Lobatto)
int& mesh_t::num_gauss_in_elem () 
{
    return num_g_pts_in_elem_;
}

// return number of material points in an element
int& mesh_t::num_mat_pt_in_elem () 
{
    return num_mat_pts_in_elem_;
}

// ---- CELLS ---- //

// returns the number of cells
int mesh_t::num_cells () const
{
    return num_cells_;
}

// return the node ids local to the cell
int mesh_t::num_nodes_in_cell () const
{
    return num_nodes_hex_;
}

// return the node ids local to the cell
int& mesh_t::nodes_in_cell (int cell_gid, int node_lid) const
{
    //DANreturn nodes_in_cell_list_[node_lid + cell_gid*num_nodes_hex_];
    return nodes_in_cell_list_(cell_gid, node_lid);
}


// return the number of cells around the cell
int& mesh_t::num_cells_in_cell (int cell_gid) const
{
    //DANreturn num_cells_in_cell_[cell_gid];
    return num_cells_in_cell_(cell_gid);
}

// return the the cells around a cell
int& mesh_t::cells_in_cell (int cell_gid, int cell_lid) const
{
    // shift index by 1 so that it is consistent with matrix syntax
    int start_indx = cells_in_cell_list_start_(cell_gid);
    
    // get the index in the global 1D array
    int index = start_indx + cell_lid;
    
    return cells_in_cell_list_(index);
}

// return corners connected to a cell
int& mesh_t::corners_in_cell (int cell_gid, int corner_lid) const
{
    //DANreturn corners_in_cell_list_[corner_lid + cell_gid*num_nodes_hex_];
    return corners_in_cell_list_(cell_gid, corner_lid);
}


int& mesh_t::elems_in_cell (int cell_gid) const
{
    //DANreturn elems_in_cell_list_[cell_gid];
    return elems_in_cell_list_(cell_gid);
}
    
    

// ---- VERTICES ---- //


// ---- NODES ---- //

// returns the number of nodes
int mesh_t::num_nodes ()
{
    return num_nodes_;
}

// returns number of cells around a node
int& mesh_t::num_cells_in_node (int node_gid) const
{
    //DANreturn num_cells_in_node_[node_gid];
    return num_cells_in_node_(node_gid);
}

// returns number of elements around a node
int& mesh_t::num_elems_in_node (int node_gid) const
{
    //DANreturn num_elems_in_node_[node_gid];
    return num_elems_in_node_(node_gid);
}


// return the cells around a node
int& mesh_t::cells_in_node (int node_gid, int cell_lid) const
{
    // shift index by 1 so that it is consistent with matrix syntax
    int start_indx = cells_in_node_list_start_(node_gid);
    
    // get the index in the global 1D array
    int index = start_indx + cell_lid;
    
    return cells_in_node_list_(index);
}

// return the elements around a node
int& mesh_t::elems_in_node (int node_gid, int elem_lid) const
{
    // shift index by 1 so that it is consistent with matrix syntax
    int start_indx = elems_in_node_list_start_(node_gid);
    
    // get the index in the global 1D array
    int index = start_indx + elem_lid;
    
    return elems_in_node_list_(index);
}
    
// return the Jacobian at a node
real_t & mesh_t::node_jacobian (int node_gid, int dim_i, int dim_j) const
{

    //DANint index = node_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j;

    //DANreturn node_jacobian_[index];
    return node_jacobian_(node_gid, dim_i, dim_j);

}

real_t & mesh_t::node_jacobian_inv (int node_gid, int dim_i, int dim_j) const
{

    return node_jacobian_inverse_(node_gid, dim_i, dim_j);
}



// return the determinant of the Jacobian at a node
real_t & mesh_t::node_det_j (int node_gid) const
{

    return node_det_j_(node_gid);
}



// ---- GAUSS POINTS ---- //


// return number of gauss points in mesh
int mesh_t::num_gauss () const
{
    return num_g_pts_;
}

// return gauss to node map
int& mesh_t::node_in_gauss (int gauss_gid) const
{            
    //DANreturn node_in_gauss_list_[gauss_gid];
    return node_in_gauss_list_(gauss_gid);
}

// return gauss in element map (internal structured grid)
int& mesh_t::gauss_in_elem (int elem_gid, int gauss_lid) 
{   
    indx_ = elem_gid*num_g_pts_in_elem_ + gauss_lid;
    
    return indx_;
}


// ---- CORNERS ---- //

// returns the number of corners
int mesh_t::num_corners () const
{
    return num_corners_;
}


// return number of corners connected to a node
int mesh_t::num_corners_in_node (int node_gid) const
{
    //DANreturn num_corners_in_node_[node_gid];
    return num_corners_in_node_(node_gid);
}

// return corner to node map
int mesh_t::corners_in_node (int node_gid, int corner_lid) const
{
    // Note: cell_in_node_list_start_ is the exact same as 
    //       corner_in_node_list_start_
    int start_indx = cells_in_node_list_start_(node_gid);

    // get the index in the global 1D array
    int index = start_indx + corner_lid;

    return corners_in_node_list_(index);
}

// ---- patches ---- //

// returns the number of elements
int mesh_t::num_patches () const
{
    return num_patches_;
}

// returns the global node id given a cell_id, local_patch_indx(0:5), local_patchnode_indx(0:3)
int mesh_t::node_in_patch_in_cell(int cell_id, int this_patch, int patchnode_lid) const
{
    // get the local node index in the cell
    int this_node = this_node_in_patch_in_cell_[patchnode_lid + this_patch*num_nodes_patch_];
    
    // return the global id for the local node index
    return nodes_in_cell(cell_id, this_node);
} // end of method


// returns the global id for a cell that is connected to the patch
// DANIELLOOK
int mesh_t::cells_in_patch(int patch_gid, int this_cell) const
{
    // get the 1D index
    int this_index = patch_gid*2 + this_cell;  // this_cell = 0 or 1

    // return the global id for the cell
    return cells_in_patch_list_(this_index);
}
      
// returns the nodes in the patch
// DANIELLOOK
int mesh_t::node_in_patch(int patch_gid, int patchnode_lid) const
{
    // get the 1D index
    int this_index = patch_gid*4;
    
    // patchnode_lid is in the range of 0:3
    return patch_nodes_list_(this_index + patchnode_lid);
}


// ---- Boundary ---- //

int mesh_t::num_bdy_sets() const
{
    return num_bdy_sets_;
}

int mesh_t::num_bdy_patches() const
{
    return num_bdy_patches_;
}

int mesh_t::bdy_patches(int this_bdy_patch) const
{
    return bdy_patches_(this_bdy_patch);
}

// returns the number per bdy-patches in a particular set
int mesh_t::num_bdy_patches_in_set (int bdy_set){
    //DANreturn num_bdy_patches_set_[bdy_set];
    return num_bdy_patches_set_(bdy_set);
}


// returns a subset of the boundary patches
int mesh_t::bdy_patches_in_set (int bdy_set, int this_patch){
    
    int start = start_index_bdy_set_(bdy_set);
    
    //DANreturn bdy_set_list_[start+this_patch];
    return bdy_set_list_(start+this_patch);
}


// ==== MESH STATE FUNCTIONS ==== // 

// ---- ELEMENTS ---- //
real_t& mesh_t::elem_vol(int elem_gid) const
{
    //DANreturn elem_vol_[elem_gid];
    return elem_vol_(elem_gid);
}


// ---- CELLS ---- //

// return the cell volume
real_t& mesh_t::cell_vol(int cell_gid) const
{
    return cell_vol_(cell_gid);
}

// return the cell coordinate position
real_t& mesh_t::cell_coords(int cell_gid, int this_dim)
{
    
    //DANcell_coords_[this_dim + cell_gid*num_dim_] = 0;
    cell_coords_(cell_gid, this_dim) = 0;
    
    // #pragma omp simd
    for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
        
        int node_gid = nodes_in_cell(cell_gid, node_lid); // get the global_id
        
        //DANcell_coords_[this_dim + cell_gid*num_dim_] += node_coords(node_gid, this_dim);
        cell_coords_(cell_gid, this_dim) += node_coords(node_gid, this_dim);
        
    } // end for loop over vertices in the cell
    
    // divide by the number of vertices
    //DANcell_coords_[this_dim + cell_gid*num_dim_] /= ( (real_t)num_nodes_hex_ );
    cell_coords_(cell_gid, this_dim) /= ( (real_t)num_nodes_hex_ );
    
    //DANreturn cell_coords_[cell_gid*num_dim_ + this_dim];
    return cell_coords_(cell_gid, this_dim);
}


// ---- VERTICES ---- //



// ---- NODES ---- //
// return the node coordinates
real_t& mesh_t::node_coords(int node_gid, int this_dim) const
{
    //DANreturn node_coords_[node_gid*num_dim_ + this_dim];
    return node_coords_(node_gid, this_dim);
}




// ---- QUADRATURE POINTS ---- //

// return jacobian at quadrature point
real_t& mesh_t::jacobian(int gauss_gid, int i, int j) const
{
    //int index = elem_gid*num_g_pts_in_elem_*num_dim_*num_dim_
    //            + gauss_lid*num_dim_*num_dim_
    //            + i*num_dim_
    //            + j;
    
    //DANreturn jacobians_[index];
    return jacobians_(gauss_gid, i, j);
}


// return determinant of jacobian at quadrature point
real_t& mesh_t::det_j(int gauss_gid) const
{
    //int index = elem_gid*num_g_pts_in_elem_
    //            + gauss_lid;
    
    //DANreturn jacobian_determinant_[index];
    return jacobian_determinant_(gauss_gid);
}


// ---- CORNERS ---- //


// ---- patches ---- //
// geometric average patch coordinate
real_t mesh_t::patch_coords(int patch_id, int this_dim) const
{
    
    real_t this_patch_coord = 0.0;

    // loop over all the vertices on the this patch
    for (int this_patchvert = 0; this_patchvert < num_nodes_patch_; this_patchvert++){
        
        // get the global vertex id
        int vert_gid = node_in_patch(patch_id, this_patchvert);
        
        // calc the coord
        this_patch_coord += node_coords(vert_gid, this_dim)/((real_t)num_nodes_patch_);

    } // end for this_patchvert
    
    return this_patch_coord;
    
} // end of patch_coords



// ---- BOUNDARY ---- //


// ==== MESH CONNECTIVITY FUNCTIONS ==== // 
    
// initialize array for mesh connectivity: all cells around a node
void mesh_t::build_connectivity(){
    
    // -- NODE TO CELL CONNECTIVITY -- //
    build_node_cell_connectivity(); 

    // -- CORNER CONNECTIVITY -- //
    build_corner_connectivity(); 

    // -- CELL TO CELL CONNECTIVITY -- //
    build_cell_cell_connectivity(); 

    // -- patches -- //
    build_patch_connectivity(); 

    // -- ELEMENTS -- //
    build_element_connectivity();

} // end of build connectivity


void mesh_t::build_node_cell_connectivity(){

    // initializing the number of corners (node-cell pair) to be zero
    int num_corners_saved[num_nodes_]; // local variable
    for (int node_gid = 0; node_gid < num_nodes_; node_gid++){
        
        num_cells_in_node_(node_gid) = 0;
        num_corners_saved[node_gid] = 0;
    }
    
    
    // count the number of corners (cell-node pair) in the mesh and in each node
    int num_corners = 0;
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            // each point-cell pair makes a corner
            num_corners ++;  // total number of corners in the entire mesh
            
            // increment the number of corners attached to this point
            int node_gid = nodes_in_cell(cell_gid, node_lid); // get the global_id
            
            num_cells_in_node_(node_gid) ++;
            
        }  // end for this_point
    } // end for cell_gid
    
    // Save number of corners in mesh
    num_corners_ = num_corners;

    // create memory for a list for all cell-node pairs
    //DANcells_in_node_list_ = new int [num_corners];
    cells_in_node_list_ = CArray <int> (num_corners);
    

    // Loop over nodes to set the start point of the ragged right array indices
    cells_in_node_list_start_(0) = 0;
    for (int node_gid = 0; node_gid < num_nodes_; node_gid++){

        // This is the start of the indices for the corners connected to this node
        cells_in_node_list_start_(node_gid+1) = cells_in_node_list_start_(node_gid)
                                                + num_cells_in_node_(node_gid);
    }
    
    
    // populate the cells connected to a node list
    int corner_gid = 0;
    
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            // increment the number of corners attached to this point
            int node_gid = nodes_in_cell(cell_gid, node_lid); // get the global_id of the node

            // Assign global index values to nodes_in_cell_list_
            //DANnodes_in_cell_list_[node_lid + cell_gid*num_nodes_hex_] = node_gid;
            // nodes_in_cell_list_(cell_gid, node_lid) = node_gid;
            
            // the global index in the cells_in_node_list_
            int index = cells_in_node_list_start_(node_gid) + num_corners_saved[node_gid];
            
            // save the global cell index to the list
            cells_in_node_list_(index) = cell_gid;
            
            // each point-cell pair makes a corner
            num_corners_saved[node_gid] ++;  //number of corners saved to this node
            
        }  // end for this_point
    } // end for cell_gid
} // end of build_node_cell_connectivity


void mesh_t::build_corner_connectivity(){

    // initializing the number of corners (node-cell pair) to be zero
    int num_corners_saved[num_nodes_]; // local variable
    for (int node_gid = 0; node_gid < num_nodes_; node_gid++){
        num_corners_saved[node_gid] = 0;
    }

    // count the number of corners (cell-node pair) in the mesh and in each node
    int num_corners = 0;
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            // each point-cell pair makes a corner
            num_corners ++;  // total number of corners in the entire mesh

        }  // end for this_point
    } // end for cell_gid
    
    // Save number of corners in mesh
    num_corners_ = num_corners;

    // create memory for a list of all corners in node
    //DANcorners_in_node_list_ = new int [num_corners];
    corners_in_node_list_ = CArray <int> (num_corners);

    // populate the cells connected to a node list and corners in a node
    int corner_gid = 0;
    
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            // increment the number of corners attached to this point
            int node_gid = nodes_in_cell(cell_gid, node_lid); // get the global_id of the node
            
            // the global index in the cells_in_node_list_
            int index = cells_in_node_list_start_(node_gid) + num_corners_saved[node_gid];

            // each point-cell pair makes a corner
            num_corners_saved[node_gid] ++;  //number of corners saved to this node

            // Save index for corner to global cell index
            corners_in_node_list_(index) = corner_gid;

            int corner_lid = node_lid;
            corners_in_cell(cell_gid, corner_lid) = corner_gid;

            corner_gid ++;
            
        }  // end for this_point
    } // end for cell_gid

    for (int node_gid = 0; node_gid < num_nodes_; node_gid++ ){

        num_corners_in_node_(node_gid) = num_corners_saved[node_gid]; 
    }
} // end of build_corner_connectivity


void mesh_t::build_cell_cell_connectivity(){
    
    // initializing the number of cell-cell pairs to be zero
    
    int num_cell_cell_saved[num_cells_]; // local variable
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        num_cells_in_cell_(cell_gid) = 0;
        num_cell_cell_saved[cell_gid] = 0;
    }
    
    int num_c_c_pairs = 0;
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            // get the global node id
            int node_gid = nodes_in_cell(cell_gid, node_lid);
            
            // loop over all cells connected to node_gid
            for (int cell_lid = 0; cell_lid < num_cells_in_node_(node_gid); cell_lid++){
                
                int neighbor_cell_id = cells_in_node(node_gid, cell_lid);
                
                // a true neighbor_cell_id is not equal to cell_gid
                if (neighbor_cell_id != cell_gid){
                    
                    // increment the number of cell-cell pairs in the mesh
                    num_c_c_pairs ++;
                    
                    // increment the number of cell-cell pairs for this cell
                    num_cells_in_cell_(cell_gid) ++;
                    
                } // end if neighbor_cell_id
            } // end for cell_lid
            
        }  // end for node_lid
    } // end for cell_gid
    
    
    // create memory for the list of cells around a cell (num_c_c_pairs is ~2x larger than needed)
    //int * temp_cell_in_cell_list = new int [num_c_c_pairs];
    auto temp_cell_in_cell_list = CArray <int> (num_c_c_pairs);

    cells_in_cell_list_start_(0) = 0;
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){

        // This is the start of the indices for the corners connected to this node
        cells_in_cell_list_start_(cell_gid+1) = cells_in_cell_list_start_(cell_gid)
                                               + num_cells_in_cell_(cell_gid);
    }


    int actual_size = 0; // the actual size of array of cell-neighbors
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            // get the global_id node id
            int node_id = nodes_in_cell(cell_gid, node_lid);
            
            // loop over all cells connected to node_id
            for (int cell_lid = 0; cell_lid < num_cells_in_node_(node_id); cell_lid++){
                
                // get the global id for the neighboring cell
                int neighbor_cell_id = cells_in_node(node_id, cell_lid);
                
                // the global index in the cells_in_cell_list_
                int index = cells_in_cell_list_start_(cell_gid) + num_cell_cell_saved[cell_gid];
                
                int save = 1; // a flag to save (=1) or not (=0)
                
                // a true neighbor_cell_id is not equal to cell_gid
                if (neighbor_cell_id == cell_gid ){
                    save = 0;  // don't save
                } // end if neighbor_cell_id
                
                // check to see if neighbor_id has been saved already
                for (int i=cells_in_cell_list_start_(cell_gid); i<index; i++){
                    
                    if (neighbor_cell_id == temp_cell_in_cell_list(i)){
                        save=0;   // don't save, it has been saved already
                    } // end if
                    
                } // end for i
                
                
                if (save==1){
                    // save the neighboring cell_id
                    temp_cell_in_cell_list(index) = neighbor_cell_id;
                    
                    // increment the number of neighboring cells saved
                    num_cell_cell_saved[cell_gid]++;
                    
                    // increment the actual number of neighboring cells saved
                    actual_size++;
                } // end if save
            } // end for cell_lid

        }  // end for node_lid
    } // end for cell_gid
    
    
    // create memory for the list of cells around a cell
    //DANcells_in_cell_list_ = new int [actual_size];
    cells_in_cell_list_ = CArray <int> (actual_size);
    
    // update the number of cells around a cell (because the estimate had duplicates)
    int index = 0;
    cells_in_cell_list_(0) = 0;
    
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        
        num_cells_in_cell_(cell_gid) = num_cell_cell_saved[cell_gid];
        
        // note: there is a buffer in the temp_cell_in_cell_list that is not
        // in num_cells_in_cell_ (it is smaller) so I will copying the
        // values from the temp to the actual array and resetting the start
        // array to use num_cell_cell_saved[cell_gid]
        
        // the global index range in the temp_cell_in_cell_list
        int start_cell = cells_in_cell_list_start_(cell_gid);
        int stop_cell  = start_cell + num_cell_cell_saved[cell_gid];
        
        cells_in_cell_list_start_(cell_gid) = index; // update the index
        
        // save the neighbors to the list
        for (int i = start_cell; i < stop_cell; i++){
            
            // save neighboring cell_gid to the final list
            cells_in_cell_list_(index) = temp_cell_in_cell_list(i);
            
            // increment the global index
            index++;
            
        } // end for i
    }// end for cell_gid
    
    // delete the temporary list of cells around a cell
    //DANdelete[] temp_cell_in_cell_list;

} // end of build_cell_cell_connectivity


void mesh_t::build_patch_connectivity(){
    
    int patch_gid = 0;
    real_t node_hash_delta = 1.0e64;
    
    real_t coord_min[num_dim_];
    real_t coord_max[num_dim_];

    // initialize to large values
    for (int dim = 0; dim < num_dim_; dim++){
        
        coord_min[dim] = 1.0e64;
        coord_max[dim] = -1.0e64;

    } // end for dim


    // Get min and max points in the mesh
    real_t coords[num_dim_];
    for(int node_gid = 0; node_gid < num_nodes_; node_gid++){
        for (int dim = 0; dim < num_dim_; dim++){

            coords[dim] = node_coords(node_gid, dim);

            coord_max[dim] = fmax(coord_max[dim], coords[dim]);
            coord_min[dim] = fmin(coord_min[dim], coords[dim]);

        }
    }

    // get minimum distance between any two points (WARNING: ONLY WORKS IN 3D)
    real_t dist_min;
    real_t dist_max;
    real_t cell_nodes[24];
    
    auto vert1 = ViewCArray <real_t> (cell_nodes, num_nodes_hex_, 3);
    
    real_t distance[28]; 
    auto dist = ViewCArray <real_t> (distance, 28);

    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        
        // Getting the coordinates of the element
        for(int node = 0; node < num_nodes_hex_; node++){
            for (int dim = 0; dim < 3; dim++)
                vert1(node, dim) = node_coords(nodes_in_cell(cell_gid, node), dim);
        }

        // loop conditions needed for distance calculation
        int countA = 0;
        int countB = 1;
        int a;
        int b;
        int loop = 0;
        
        
        // Solving for the magnitude of distance between each node
        for (int i = 0; i < 28; i++){
            
            a = countA;
            b = countB;
            
            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt(( pow((vert1(b,0) - vert1(a,0)), 2.0)
                                + pow((vert1(b,1) - vert1(a,1)), 2.0)
                                + pow((vert1(b,2) - vert1(a,2)), 2.0))));

            countB++;
            countA++;
            
            //tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        }

        dist_min = dist(0);
        dist_max = 0.0;
        
        for(int i = 0; i < 28; ++i){
            dist_min = fmin(dist(i),dist_min);
            dist_max = fmax(dist(i),dist_max);
        }
    }

    node_hash_delta = fmin(node_hash_delta, dist_min);
    node_hash_delta = node_hash_delta/2.0;

    // calculate the 1d array length of the hash table for nodes
    real_t num_bins[num_dim_];
    
    for (int dim = 0; dim < num_dim_; dim++){
        num_bins[dim] = (coord_max[dim] - coord_min[dim]  + 1.e-12)/node_hash_delta;
    }

    // Set size of hash key array
    int hash_keys[num_cells_*num_patches_hex_];

    real_t patch_hash_idx_real[num_dim_];

    // Calculate the hash keys and find max key
    int hash_count = 0;
    int max_key = 0;
    
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int patch_lid = 0; patch_lid < num_patches_hex_; patch_lid++){
            
            // the patch coordinates
            real_t patch_coords[num_dim_];
            
            // initialize to zero
            for (int dim = 0; dim < num_dim_; dim++){
                patch_coords[dim] = 0.0;
            } // end for dim
            
            
            // loop over all the nodes on the this patch
            for (int patchnode_lid = 0; patchnode_lid < num_nodes_patch_; patchnode_lid++){
                
                // get the global node id
                int node_gid = node_in_patch_in_cell(cell_gid, patch_lid, patchnode_lid);
                
                for (int dim = 0; dim < num_dim_; dim++){
                    patch_coords[dim] += node_coords(node_gid, dim)/num_nodes_patch_;
                } // end for dim
                    
            } // end for patchnode_lid
            
            // calculate the patch hash index for these patch_coordinates
            for (int dim = 0; dim < num_dim_; dim++){
                patch_hash_idx_real[dim] = fmax(1e-16, (patch_coords[dim]-coord_min[dim] + 1e-12)/node_hash_delta);
            } // end for dim

            // the 1D index
            real_t patch_hash_idx_real_1d;
            
            // i + j*num_x + k*num_x*num_y
            if (num_dim_ == 2){
                patch_hash_idx_real_1d = patch_hash_idx_real[0] + patch_hash_idx_real[1]*num_bins[0];
            }

            else{
                patch_hash_idx_real_1d =
                      patch_hash_idx_real[0]
                    + patch_hash_idx_real[1]*num_bins[0]
                    + patch_hash_idx_real[2]*num_bins[0]*num_bins[1];
            }
            
            hash_keys[hash_count] = (int)patch_hash_idx_real_1d;
            max_key = std::max(max_key, hash_keys[hash_count]);
            hash_count++;

        } // end for patch_lid
    } // end for cell_gid


    // Allocate hash table and initialize key locations to -1
    //DANint * hash_table = new int [max_key + 10];
    auto hash_table = CArray <int> (max_key + 10);
    
    for (int key_idx = 0; key_idx < hash_count; key_idx++){
        hash_table(hash_keys[key_idx]) = -1;
    } // end for cell_gid


    // count the number of cells around each patch
    patch_gid = 0;
    hash_count = 0;
    
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int patch_lid = 0; patch_lid < num_patches_hex_; patch_lid++){
            
            // count the number of patches
            if (hash_table(hash_keys[hash_count]) == -1){
                // save the patch_id to the hash table
                hash_table(hash_keys[hash_count]) = patch_gid;
                patch_gid++;
            }

            hash_count++;
        } // end for patch_lid
    } // end for cell_gid

    // set the number of patches in the mesh
    num_patches_ = patch_gid;

    // create memory for the patch structures
    //DANcells_in_patch_list_ = new int [num_patches_*2];   // two cells per patch
    cells_in_patch_list_ = CArray <int> (num_patches_ * 2);
    //DANpatch_nodes_list_    = new int [num_patches_*4];   // four nodes per patch in hex
    patch_nodes_list_    = CArray <int> (num_patches_ * 4);

    // initialize the cells_in_patch_list to -1
    for (int patch_gid = 0; patch_gid < 2*num_patches_; patch_gid++){
        cells_in_patch_list_(patch_gid) = -1;
    }   

    for (int patch_gid = 0; patch_gid < 4*num_patches_; patch_gid++){
        patch_nodes_list_(patch_gid) = -1;
    }


    hash_count = 0;
    
    // save the cell_gid to the cells_in_patch_list
    for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
        for (int patch_lid = 0; patch_lid < num_patches_hex_; patch_lid++){
            
            // get the patch_id
            patch_gid = hash_table(hash_keys[hash_count]);

            // a temp variable for storing the node global ids on the patch
            int these_nodes[num_nodes_patch_];

            // loop over all the vertices on the this patch
            for (int patchnode_lid = 0; patchnode_lid < num_nodes_patch_; patchnode_lid++){
                
                // get the global node id
                int node_gid = node_in_patch_in_cell(cell_gid, patch_lid, patchnode_lid);
                
                // save the global node id
                these_nodes[patchnode_lid] = node_gid;
                
            } // end for patchnode_lid

            // save the cell_gid to the cells_in_patch_list
            if (cells_in_patch_list_(patch_gid*2) == -1){
                
                // save the cell index for this patch
                int this_index = patch_gid*2;  // index in the list
                
                cells_in_patch_list_(this_index) = cell_gid;
                
                // save the nodes for this patch
                this_index = patch_gid*4;
                
                for (int patchnode_lid = 0; patchnode_lid < num_nodes_patch_; patchnode_lid++){
                    patch_nodes_list_(this_index + patchnode_lid) = these_nodes[patchnode_lid];
                }
            }

            else{
                // it is the other cell connected to this patch

                int this_index = patch_gid*2 + 1; // + num_patches_; // index in the list

                cells_in_patch_list_(this_index) = cell_gid;
            }

            hash_count++;

        } // end for patch_lid
    } // end for cell_gid

    // Delete memeory for the hash table
    //DANdelete[] hash_table;

} // end of build patches


void mesh_t::build_element_connectivity(){

    // Initialize list to count number of times a node has been 
    // touched through the node_in_gauss list

    int times_hit[num_nodes_];
    
    for(int node_gid = 0; node_gid < num_nodes_; node_gid++){
        num_elems_in_node_(node_gid) = 0;
        times_hit[node_gid] = 0;
    }

    for(int elem_gid = 0; elem_gid < num_elem_; elem_gid++){
        for(int gauss_lid = 0; gauss_lid < num_g_pts_in_elem_; gauss_lid++){

            // get gauss global id and use to get node global id
            int gauss_gid = gauss_in_elem(elem_gid, gauss_lid);
            int node_gid  = node_in_gauss(gauss_gid);

            // every time the node gid is hit add one to the hit count
            num_elems_in_node_(node_gid) += 1;

            // add nodes in element here
            nodes_in_elem(elem_gid, gauss_lid) = node_gid;
        }
    }


    // Walk over each node, if the node was touched more than once throught the node_in_gauss list add 
    // the number of times to the elem_in_node array size

    // base size is the number of nodes, one is added for each double counted one
    int elem_in_node_size = 0;
    for(int node_gid = 0; node_gid < num_nodes_; node_gid++){
        elem_in_node_size += num_elems_in_node_(node_gid);
    }

    // get access pattern and total size of ragged right for elems_in_node
    
    elems_in_node_list_start_(0) = 0;
    // starts at one because zero index has been saved
    for(int node_gid = 0; node_gid < num_nodes_; node_gid++){
        elems_in_node_list_start_(node_gid+1) = elems_in_node_list_start_(node_gid) + num_elems_in_node_(node_gid);
    }

    // std::cout<<"Before getting size of elems in node list"<<std::endl;

    // create memory for elems_in_node_list_
    //DANelems_in_node_list_ = new int [elem_in_node_size];
    elems_in_node_list_ = CArray <int> (elem_in_node_size);

    for(int elem_gid = 0; elem_gid < num_elem_; elem_gid++){
        for(int node_lid = 0; node_lid < num_g_pts_in_elem_; node_lid++){

            // get gauss global id and use to get node global id
            int gauss_gid = gauss_in_elem(elem_gid, node_lid);
            int node_gid  = node_in_gauss(gauss_gid);

            int indx = elems_in_node_list_start_(node_gid) + times_hit[node_gid];

            elems_in_node_list_(indx) = elem_gid;

            times_hit[node_gid]++;
        }
    }

    // verify that things makes sense
    // for(int node_gid = 0; node_gid < num_nodes_; node_gid++){

    //     int test = num_elems_in_node_[node_gid] - times_hit[node_gid];

    //     if (test != 0){
    //         std::cout<<"ERROR IN ELEMENTS IN NODE"<<std::endl;
    //     }
    // }

    // Find all elements connected to an element (same coding as cells_in_cell)

    // initializing the number of element-element pairs to be zero
    int num_elem_elem_saved[num_cells_]; // local variable
    for (int elem_gid = 0; elem_gid < num_elem_; elem_gid++){
        
        num_elems_in_elem_(elem_gid) = 0;
        num_elem_elem_saved[elem_gid] = 0;
    }
    
    int num_e_e_pairs = 0;
    
    for (int elem_gid = 0; elem_gid < num_elem_; elem_gid++){
        for (int node_lid = 0; node_lid < num_nodes_in_elem_; node_lid++){
            
            // get the global node id
            int node_gid = nodes_in_elem(elem_gid, node_lid);
            
            // loop over all cells connected to node_gid
            for (int elem_lid = 0; elem_lid < num_elems_in_node_(node_gid); elem_lid++){
                
                int neighbor_elem_gid = elems_in_node(node_gid, elem_lid);

                // a true neighbor_elem_gid is not equal to elem_gid
                if (neighbor_elem_gid != elem_gid){
                    
                    // increment the number of elem_elem pairs in the mesh
                    num_e_e_pairs ++;
                    
                    // increment the number of elem_elem pairs for this element
                    num_elems_in_elem_(elem_gid) ++;
                    
                } // end if neighbor_cell_id
            } // end for cell_lid
            
        }  // end for node_lid
    } // end for cell_gid
    
    
    // create memory for the list of cells around a cell (num_e_e_pairs is 2x larger than needed)
    //DANint * temp_elems_in_elem_list = new int [num_e_e_pairs];
    auto temp_elems_in_elem_list = CArray <int> (num_e_e_pairs);
    
    elems_in_elem_list_start_(0) = 0;
    for (int elem_gid = 0; elem_gid < num_elem_; elem_gid++){

        // This is the start of the indices for the corners connected to this node
        elems_in_elem_list_start_(elem_gid+1) = elems_in_elem_list_start_(elem_gid)
                                               + num_elems_in_elem_(elem_gid);
    }
    

    int actual_size = 0; // the actual size of array of cell-neighbors
    for (int elem_gid = 0; elem_gid < num_elem_; elem_gid++){
        for (int node_lid = 0; node_lid < num_nodes_in_elem_; node_lid++){
            
            // get the global node id
            int node_gid = nodes_in_elem(elem_gid, node_lid);
            
            // loop over all cells connected to node_gid
            for (int elem_lid = 0; elem_lid < num_elems_in_node_(node_gid); elem_lid++){
                
                int neighbor_elem_gid = elems_in_node(node_gid, elem_lid);
                
                // the global index in the cells_in_cell_list_
                int index = elems_in_elem_list_start_(elem_gid) + num_elem_elem_saved[elem_gid];
                
                int save = 1; // a flag to save (=1) or not (=0)
                
                // a true neighbor_elem_gid is not equal to elem_gid
                if (neighbor_elem_gid == elem_gid ){
                    save = 0;  // don't save
                } // end if neighbor_elem_gid
                
                // check to see if neighbor_id has been saved already
                for (int i = elems_in_elem_list_start_(elem_gid); i < index; i++){
                    
                    if (neighbor_elem_gid == temp_elems_in_elem_list(i)){
                        save=0;   // don't save, it has been saved already
                    } // end if
                    
                } // end for i
                
                
                if (save==1){
                    // save the neighboring cell_id
                    temp_elems_in_elem_list(index) = neighbor_elem_gid;
                    
                    // increment the number of neighboring cells saved
                    num_elem_elem_saved[elem_gid]++;
                    
                    // increment the actual number of neighboring cells saved
                    actual_size++;
                } // end if save
                
            } // end for elem_lid

        }  // end for node_lid
    } // end for elem_gid
    
    
    // create memory for the list of cells around a cell
    //DANelems_in_elem_list_ = new int [actual_size];
    //elems_in_elem_list_ = CArray <int> (actual_size);
    auto eiels = CArray <size_t> (num_elem_+1);
    for (int tempi = 0; tempi < num_elem_+1; tempi++) {
        eiels(tempi) = elems_in_elem_list_start_(tempi);
    }

    //printf("My start %d, also %d\n", num_elem_, actual_size);
    for (int shift = 0; shift < num_elem_; shift++) {
        // change start list to stride
        //printf("%d) %d,%d -> ", shift, eiels(shift), eiels(shift+1));
        eiels(shift) = eiels(shift+1) - eiels(shift);
        //printf("%d\n", eiels(shift));
        //elems_in_elem_list_start_(shift) = elems_in_elem_list_start_(shift+1) - elems_in_elem_list_start(shift);
    }
    elems_in_elem_list_ = RaggedRightArray <int> (eiels);
    for (int check = 0; check < num_elem_; check++) {
       // printf("stride %d\n", elems_in_elem_list_.stride(check));
    }

    // update the number of cells around a cell (because the estimate had duplicates)
    int index = 0;
    //DANelems_in_elem_list_(0) = 0;
    
    for (int elem_gid = 0; elem_gid < num_elem_; elem_gid++){
        
        num_elems_in_elem_(elem_gid) = num_elem_elem_saved[elem_gid];
        
        // note: there is a buffer in the temp_cell_in_cell_list that is not
        // in num_elems_in_elem_ (it is smaller) so I will copying the
        // values from the temp to the actual array and resetting the start
        // array to use num_elem_elem_saved[elem_gid]
        
        // the global index range in the temp_cell_in_cell_list
        int list_start = elems_in_elem_list_start_(elem_gid);
        int list_stop  = list_start + num_elem_elem_saved[elem_gid];
        
        //DANelems_in_elem_list_start_(elem_gid) = index; // update the index
        
        // save the neighbors to the list
        //DANfor (int i = list_start; i < list_stop; i++){
        //DANfor (int i = 0; i < elems_in_elem_list_.stride(elem_gid); i++){
        for (int i = 0, j = list_start; i < num_elem_elem_saved[elem_gid]; i++, j++) {
            
            // save neighboring elem_gid to the final list
            //DANelems_in_elem_list_(index) = temp_elems_in_elem_list(i);
            elems_in_elem_list_(elem_gid, i) = temp_elems_in_elem_list(j);
            
            // increment the global index
            index++;
            
        } // end for i
    }// end for elem_gid


    printf("ragged check\n");
    for (int elem_gid = 0; elem_gid < num_elem_; elem_gid++) {
        for (int i = 0; i < num_elem_elem_saved[elem_gid]; i++) {
            printf("%d ", elems_in_elem_list_(elem_gid, i));
        }
    }
    printf("\n");
    
    
    // delete the temporary list of cells around a cell
    //DANdelete[] temp_elems_in_elem_list;

} // end build element connectivity


// identify the boundary patches
void mesh_t::build_bdy_patches (){
    
    
    int bdy_patch_gid = 0;
    
    
    // loop over the patches in the mesh
    for(int patch_gid = 0; patch_gid < num_patches_; patch_gid++){
        
        // loop over the two cells on this patch
        for (int cell_lid = 0; cell_lid < 2; cell_lid++){
            
            // check to see if a cell has index of -1
            if (cells_in_patch(patch_gid, cell_lid) == -1){
               bdy_patch_gid ++;
            } // end if
                
        } // end for cell_lid
    } // end for patch_gid
    
    
    // save the number of boundary patches in the mesh
    num_bdy_patches_ = bdy_patch_gid;
    
    // allocate the memory for the boundary patches array
    //DANbdy_patches_ = new int[num_bdy_patches_];
    bdy_patches_ = CArray <int> (num_bdy_patches_);
    
    
    // save the global indices for the boundary patches
    
    // loop over the patches in the mesh
    bdy_patch_gid = 0;

    for(int patch_gid = 0; patch_gid < num_patches_; patch_gid++){
        
        // loop over the two cells on this patch
        for (int cell_lid = 0; cell_lid < 2; cell_lid++){
            
            // check to see if a cell has index of -1
            if (cells_in_patch(patch_gid, cell_lid) == -1){
                
                // save the patch index
                bdy_patches_(bdy_patch_gid) = patch_gid;
                
                // increment the counter
                bdy_patch_gid++;
                
            } // end if  
        } // end for cell_lid
    } // end for patch_gid
} // end of function



// ---- bdy sets ----

// returns a subset of the boundary patches
int mesh_t::set_bdy_patches (int bdy_set, int patch_lid){
    
    int start = start_index_bdy_set_(bdy_set);
    
    return bdy_set_list_(start+patch_lid);
}



// set planes for tagging sub sets of boundary patches
// bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
// val = plane value, radius, radius
void mesh_t::tag_bdys(int this_bc_tag, real_t val, int bdy_set){
    
    if (bdy_set == num_bdy_sets_){
        std::cout << " ERROR: number of boundary sets must be increased by "
                  << bdy_set-num_bdy_sets_+1 << std::endl;
        exit(0);
    }
    
    // the start index for the first list is zero
    start_index_bdy_set_(0) = 0;
    
    
    // save the boundary vertices to this set that are on the plane
    int counter = 0;
    int start = start_index_bdy_set_(bdy_set);
    
    for (int this_bdy_patch = 0; this_bdy_patch < num_bdy_patches_; this_bdy_patch++) {
        
        // save the patch index
        int bdy_patch_gid = bdy_patches(this_bdy_patch);
        
        // check to see if this patch is on the specified plane
        int is_on_bdy = check_bdy(bdy_patch_gid, this_bc_tag, val); // no=0, yes=1
        
        if (is_on_bdy == 1){
            bdy_set_list_(start+counter) = bdy_patch_gid;
            counter ++;
        }
    } // end for bdy_patch
    
    // save the number of bdy patches in the set
    num_bdy_patches_set_(bdy_set) = counter;
    
    // save the starting index for the next bdy_set
    start_index_bdy_set_(bdy_set+1) = start_index_bdy_set_(bdy_set) + counter;
    
    
    // compress the list to reduce the memory if it is the last set
    //DANif (bdy_set == num_bdy_sets_-1){
        //compress_bdy_set();
    //}
    
    std::cout << " tagged boundary patches " << std::endl;
    
} // end of method


// compress the bdy_set_list to reduce the memory
/*//DAN
void mesh_t::compress_bdy_set(){
    
    // the actual size of the bdy set list
    int length = start_index_bdy_set_(num_bdy_sets_);
    
    // create a temp array of correct size
    int * temp_bdy_list = new int [length];
    
    // save the values to the temp array
    for (int i = 0; i < length; i++){
        temp_bdy_list[i] = bdy_set_list_(i);
    }
    
    // delete original array and make a new one of correct size
    DANdelete[] bdy_set_list_;

    bdy_set_list_ = new int [length];
    
    // save the values to the bdy_set_list
    for (int i = 0; i < length; i++){
        bdy_set_list_[i] = temp_bdy_list[i];
    }
    
    // delete the temp array
    delete[] temp_bdy_list;
    
} // end of compress_bdy_set
*/


// routine for checking to see if a vertix is on a boundary
// bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
// val = plane value, radius, radius
int mesh_t::check_bdy(int patch_gid, int this_bc_tag, real_t val){
    
    // default bool is not on the boundary
    int is_on_bdy = 0;
    
    // the patch coordinates
    real_t these_patch_coords[num_dim_];
    
    for (int dim = 0; dim < num_dim_; dim++){
        these_patch_coords[dim] = patch_coords(patch_gid, dim);
    } // end for dim
    
    
    // a x-plane
    if (this_bc_tag == 0){
        
        if ( fabs(these_patch_coords[0] - val) <= 1.0e-8 ) is_on_bdy = 1;
        
    }// end if on type
    
    // a y-plane
    else if (this_bc_tag == 1){
        
        if ( fabs(these_patch_coords[1] - val) <= 1.0e-8 ) is_on_bdy = 1;
        
    }// end if on type
    
    // a z-plane
    else if (this_bc_tag == 2){
        
        if ( fabs(these_patch_coords[2] - val) <= 1.0e-8 ) is_on_bdy = 1;
        
    }// end if on type
    
    
    // cylinderical shell where radius = sqrt(x^2 + y^2)
    else if (this_bc_tag == 3){
        
        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1]);
        
        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;

        
    }// end if on type
    
    // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
    else if (this_bc_tag == 4){
        
        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1] +
                        these_patch_coords[2]*these_patch_coords[2]);
        
        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;
        
    } // end if on type
    
    return is_on_bdy;
    
} // end method to check bdy


// deconstructor
mesh_t::~mesh_t () {
    
    // ---- ELEMENTS ---- //
    //DANdelete[] cells_in_elem_;
    //DANdelete[] num_elems_in_elem_;
    //DANdelete[] elems_in_elem_list_start_;
    //DANdelete[] elems_in_elem_list_;
    //DANdelete[] nodes_in_elem_list_;


    // ---- CELLS ---- //
    //DANdelete[] nodes_in_cell_list_;
    //DANdelete[] num_cells_in_cell_;
    //DANdelete[] cells_in_cell_list_start_;
    //DANdelete[] cells_in_cell_list_;
    //DANdelete[] elems_in_cell_list_;


    // ---- VERTICES ---- //

    
    // ---- NODES ---- //
    //DANdelete[] num_cells_in_node_;
    //DANdelete[] cells_in_node_list_start_;
    //DANdelete[] cells_in_node_list_;

    //DANdelete[] num_elems_in_node_;
    //DANdelete[] elems_in_node_list_start_;
    //DANdelete[] elems_in_node_list_;

    // ---- GAUSS POINTS ---- //
    //DANdelete[] node_in_gauss_list_;


    // ---- CORNERS ---- //
    //DANdelete[] num_corners_in_node_;
    //DANdelete[] corners_in_cell_list_;
    //DANdelete[] corners_in_node_list_start_;
    //DANdelete[] corners_in_node_list_;


    // ---- PATCHES ---- //
    //DANdelete[] patch_nodes_list_;       
    //DANdelete[] cells_in_patch_list_;  

    // ---- BOUNDARY ---- //
    //DANdelete[] bdy_patches_;
    //DANdelete[] bdy_set_list_;
    //DANdelete[] start_index_bdy_set_;
    //DANdelete[] num_bdy_patches_set_; 

    // ---- MESH STATE ---- //
    // ---- ELEMENT ---- //
    //DANdelete[] elem_vol_; 
    

    // ---- CELLS ---- //
    //DANdelete[] cell_vol_;
    //DANdelete[] cell_coords_;


    // ---- NODES ---- //
    //DANdelete[] node_coords_;


    // ---- QUADRATURE POINTS ---- //
    //DANdelete[] jacobians_;            // size of num_g_pts_*num_dim_*num_dim_
    //DANdelete[] jacobian_determinant_; // size of num_g_pts_

} // end of mesh deconstructor


void refine_mesh(
    mesh_t& init_mesh, 
    mesh_t& mesh, 
    const int p_order, 
    const int dim){

    // High order mesh parameters
    int num_sub_1d;         // num subcells in 1d
    int num_g_pts_1d;       // num gauss points in 1d

    int num_g_pts;            // number of gauss points
    int num_subcells_per_elem; // number subcells in an element

    int num_elem = init_mesh.num_cells();

    if(p_order == 0){
                
        num_sub_1d   = 1;
        num_g_pts_1d = 2;
        
        num_g_pts  = pow(num_g_pts_1d, dim);
        num_subcells_per_elem = pow((num_g_pts_1d-1), dim);
    }

    else{

        num_sub_1d = p_order*2;         // num subcells in 1d
        num_g_pts_1d = 2 * p_order + 1; // num gauss points in 1d

        num_g_pts  = pow(num_g_pts_1d, dim);            // number of gauss points
        num_subcells_per_elem = pow((num_sub_1d), dim); // number subcells in an element
    }

    //  ---------------------------------------------------------------------------
    //  Initailize Element and cell information in on high order mesh
    //  ---------------------------------------------------------------------------

    // PLACEHOLDER CCH RECONSTRUCTION ORDER
    int r_order = 0; 

    mesh.init_element(p_order, dim, num_elem);
    mesh.init_cells(num_elem*num_subcells_per_elem);
    mesh.init_gauss_pts();

    //  ---------------------------------------------------------------------------
    //  Generate point positiont in reference space to map onto initial mesh
    //  ---------------------------------------------------------------------------

    auto temp_pts = CArray<real_t> (num_g_pts_1d, num_g_pts_1d, num_g_pts_1d, 3);

    double dx = 2.0/((double)num_g_pts_1d - 1.0);  // len/(num_nodes-1)
    double dy = 2.0/((double)num_g_pts_1d - 1.0);  // len/(num_nodes-1)
    double dz = 2.0/((double)num_g_pts_1d - 1.0);  // len/(num_nodes-1)

    for(int k = 0; k < num_g_pts_1d; k++){
        for(int j = 0; j < num_g_pts_1d; j++){
            for(int i = 0; i < num_g_pts_1d; i++){
                temp_pts(i, j, k, 0) = -1.0 + (double)i*dx;
                temp_pts(i, j, k, 1) = -1.0 + (double)j*dy;
                temp_pts(i, j, k, 2) = -1.0 + (double)k*dz;
            }
        }
    }


    //  ---------------------------------------------------------------------------
    //  Map new points to real space using basis functions
    //  ---------------------------------------------------------------------------

    // temp array to hold positions
    CArray <real_t> temp_gauss_point_coords;
    //DANtemp_gauss_point_coords = CArray <real_t> (num_elem*num_g_pts, dim);
    auto g_points_in_mesh = CArray <real_t> (num_elem*num_g_pts, dim);
    //DANauto g_points_in_mesh = ViewCArray <real_t> (temp_gauss_point_coords, num_elem*num_g_pts, dim);


    // Reference node positions for element (currently p1, replace with element library)
    real_t ref_vert[8][3] = // listed as {Xi, Eta, Mu}
        {
        // Bottom Nodes
        {-1.0, -1.0, -1.0},// 0
        {+1.0, -1.0, -1.0},// 1
        {-1.0, +1.0, -1.0},// 2
        {+1.0, +1.0, -1.0},// 3
        // Top Nodes
        {-1.0, -1.0, +1.0},// 4
        {+1.0, -1.0, +1.0},// 5
        {-1.0, +1.0, +1.0},// 6
        {+1.0, +1.0, +1.0},// 7
        };


    real_t basis[8];    // basis function evaluation value


    // Inital mesh node positions
    real_t x_init_a[8];
    real_t y_init_a[8];
    real_t z_init_a[8];

    auto x_init = ViewCArray <real_t> (x_init_a, 8);
    auto y_init = ViewCArray <real_t> (y_init_a, 8);
    auto z_init = ViewCArray <real_t> (z_init_a, 8);

    int g_point_count = 0; 

    // mapping points using the basis functions
    for(int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        // Assign initial positions from initial mesh
        for(int node_lid = 0; node_lid < 8; node_lid++){
            x_init(node_lid) = init_mesh.node_coords(init_mesh.nodes_in_cell(elem_gid, node_lid), 0);
            y_init(node_lid) = init_mesh.node_coords(init_mesh.nodes_in_cell(elem_gid, node_lid), 1);
            z_init(node_lid) = init_mesh.node_coords(init_mesh.nodes_in_cell(elem_gid, node_lid), 2);
        }

        // Walk over gauss points as i,j,k mesh to calculate basis 
        for(int k = 0; k < num_g_pts_1d; k++){
            for(int j = 0; j < num_g_pts_1d; j++){
                for(int i = 0; i < num_g_pts_1d; i++){
                    
                    for (int vert_lid = 0; vert_lid < 8; vert_lid++ ){
                        basis[vert_lid] = 1.0/8.0
                                        * (1.0 + temp_pts(i, j, k, 0)*ref_vert[vert_lid][0])
                                        * (1.0 + temp_pts(i, j, k, 1)*ref_vert[vert_lid][1])
                                        * (1.0 + temp_pts(i, j, k, 2)*ref_vert[vert_lid][2]);
                    }

                    g_points_in_mesh(g_point_count, 0) = 0.0;
                    g_points_in_mesh(g_point_count, 1) = 0.0;
                    g_points_in_mesh(g_point_count, 2) = 0.0;

                    // Walk over vertices to map new points onto mesh
                    for (int vert_lid = 0; vert_lid < 8; vert_lid++ ){
                        g_points_in_mesh(g_point_count, 0) += basis[vert_lid]*x_init(vert_lid);
                        g_points_in_mesh(g_point_count, 1) += basis[vert_lid]*y_init(vert_lid);
                        g_points_in_mesh(g_point_count, 2) += basis[vert_lid]*z_init(vert_lid);
                    } // end for vert_lid

                    g_point_count++;  
                }
            }
        }
    }

    //  ---------------------------------------------------------------------------
    //  Hash x, y, and x coordinates to eliminate double counted points for 
    //  node index. 
    //  ---------------------------------------------------------------------------


    real_t pos_max[dim];
    real_t pos_min[dim];

    for(int i = 0; i < dim; i++){
        pos_max[i] = -1.0E16;
        pos_min[i] =  1.0E16;
    }

    real_t position[3]; 

    // Get min and max points in the mesh
    for(int point = 0; point< init_mesh.num_nodes(); point++){
        for(int i = 0; i < dim; i++){
            
            position[i] = init_mesh.node_coords(point, i);
            pos_max[i] = fmax(pos_max[i], position[i]);
            pos_min[i] = fmin(pos_min[i], position[i]);
        }
    }

    // get minimum distance between any two points (WARNING: ONLY WORKS IN 3D)
    real_t dist_min;
    real_t dist_max;
    real_t cell_nodes[24];
    
    auto vert1 = ViewCArray <real_t> (cell_nodes, 8, 3);
    
    real_t distance[28]; 
    auto dist = ViewCArray <real_t> (distance, 28);

    for (int cell_gid = 0; cell_gid < init_mesh.num_cells(); cell_gid++){
        
        // Getting the coordinates of the element
        for(int node = 0; node < 8; node++){
            for (int dim = 0; dim < 3; dim++)
                vert1(node, dim) = init_mesh.node_coords(init_mesh.nodes_in_cell(cell_gid, node), dim);
        }

        // loop conditions needed for distance calculation
        int countA = 0;
        int countB = 1;
        int a;
        int b;
        int loop = 0;
        
        
        // Solving for the magnitude of distance between each node
        for (int i = 0; i < 28; i++){
            
            a = countA;
            b = countB;
            
            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt(( pow((vert1(b,0) - vert1(a,0)), 2.0)
                                + pow((vert1(b,1) - vert1(a,1)), 2.0)
                                + pow((vert1(b,2) - vert1(a,2)), 2.0))));

            countB++;
            countA++;
            
            //tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        }

        dist_min = 1e64;
        dist_max = 0.0;
        
        for(int i = 0; i < 28; ++i){
            dist_min = fmin(dist(i),dist_min);
            dist_max = fmax(dist(i),dist_max);
        }
    }

    std::cout<<"Min_dist = "<<dist_min<<std::endl;

    // Number of subcells per dimension to be created
    real_t sub;
    
    if (p_order == 0) sub = 1.0;
    else sub = p_order;

    real_t h = dist_min/(7.0*(sub)); // number of subdivisions between any two points

    std::cout<<"Hash dist = "<<h<<std::endl;

    // Define number of bins in each direction
    real_t float_bins[3];
    for(int i = 0; i < dim; i++){

        float_bins[i] = fmax(1e-16, (pos_max[i] - pos_min[i] + 1.0 + 1e-14)/h); //1e-14
    }

    // Convert # of bins to ints
    int int_bins[3];
    for(int i = 0; i < dim; i++){

        int_bins[i] = (int)float_bins[i];
    }

    real_t float_idx[3];   // float values for index
    int int_idx[3];        // int values for index
    
    int key;            // hash key 
    int max_key = 0;    // larges hash key value


    // Getting hash keys from x,y,z positions and largest key value
    int h_keys[num_g_pts*num_elem];
    for(int g_pt = 0; g_pt < num_g_pts*num_elem; g_pt++){
        
        real_t coords[3];
        for(int i = 0; i < dim; i++){

            coords[i]    = g_points_in_mesh(g_pt, i);
            float_idx[i] = fmax(1e-16, (coords[i] - pos_min[i] + 1e-14)/(h));
            int_idx[i] = (int)float_idx[i];
        }

        
        // i + j*num_x + k*num_x*num_y
        if (dim == 2){
            key = int_idx[0] + int_idx[1]*int_bins[0];
        }

        else{
            key = int_idx[0] 
                + int_idx[1]*int_bins[0] 
                + int_idx[2]*int_bins[0]*int_bins[1];
        }

        h_keys[g_pt] = key;
        max_key = std::max(max_key, key);
    }

    // Allocating array for hash table
    //DANint * hash; 
    CArray <int> hash;
    hash = CArray <int> (max_key+10); 

    // Initializing values at key positions to zero
    for(int g_pt = 0; g_pt < num_g_pts*num_elem; g_pt++){
        hash(h_keys[g_pt]) = 0;
    }

    // Temporary array for gauss->node map and node->gauss map

    //DANint * node_to_gauss_map;
    CArray <int> node_to_gauss_map;

    node_to_gauss_map = CArray <int> (num_g_pts*num_elem);

    // counters
    int num_nodes = 0;
    int node_gid = 0;

    // walk over all gauss points 
    for(int g_pt = 0; g_pt < num_g_pts*num_elem; g_pt++){
        
        // Subtract 1 every time the index is touched
        if(hash(h_keys[g_pt]) <= 0){
            hash(h_keys[g_pt]) += -1;
        }
        
        // If this is the first time the index is touched add to 
        // node_to_gauss_map (WARNING: ONLY THE FIRST TIME IS COUNTED)
        // and index the number of nodes 
        if(hash(h_keys[g_pt]) == -1){

            node_to_gauss_map(num_nodes) = g_pt;
            num_nodes++;
        }

        // If this index has been touched before, replace hash value with
        // node id
        if(hash(h_keys[g_pt]) <= -1){
            
            hash(h_keys[g_pt]) = node_gid;
            
            // gauss_node_map[g_pt] = node_gid;
            mesh.node_in_gauss(g_pt) = node_gid;

            node_gid++;
        }

        // If hash value is positive, then the value is the index
        // for the single node assiciated with this g_point
        else{
            // gauss_node_map[g_pt] = hash[h_keys[g_pt]];
            mesh.node_in_gauss(g_pt) = hash(h_keys[g_pt]);
        }
    }

    // remove hash table
    //DANdelete[] hash;

    // Initialize nodes on sub_mesh
    mesh.init_nodes(num_nodes);

    //  ---------------------------------------------------------------------------
    //  Write position to nodes 
    //  ---------------------------------------------------------------------------
    
    for(int node_gid=0; node_gid<num_nodes; node_gid++){
        for(int i = 0; i < dim; i++){

            mesh.node_coords(node_gid, i) = g_points_in_mesh(node_to_gauss_map(node_gid), i);
        }
    }

    // delete unneeded arrays
    //DANdelete[] temp_gauss_point_coords;
    //DANdelete[] node_to_gauss_map;



    //  ---------------------------------------------------------------------------
    //  Get gauss points and nodes associated with each cell, 
    //  as well as the cells associated with each element
    //  ---------------------------------------------------------------------------

    // auto gauss_id_in_cell = CArray<int> (sub_mesh.num_cells(), num_sub_1d*num_sub_1d*num_sub_1d, 8);
    int sub_in_elem = num_sub_1d*num_sub_1d*num_sub_1d;
    //int * gauss_id_in_cell;

    //DANgauss_id_in_cell = new int[mesh.num_cells()*sub_in_elem*8];
    auto gauss_in_cell = CArray<int> (mesh.num_cells(), sub_in_elem, 8);

    int p0, p1, p2, p3, p4, p5, p6, p7;
    p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0;
    
    int num_1d = num_g_pts_1d;
    int cell_index = 0;
    int cell_mesh_index = 0;

    for(int elem_gid = 0; elem_gid < num_elem; elem_gid++){
        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    // The p# point to a global gauss point index before double counting 
                    p0 = (i)     + (j)*num_1d   + (k)*num_1d*num_1d;
                    p1 = (i+1)   + (j)*num_1d   + (k)*num_1d*num_1d;
                    p2 = (i)     + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p3 = (i+1)   + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p4 = (i)     + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p5 = (i+1)   + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p6 = (i)     + (j+1)*num_1d + (k+1)*num_1d*num_1d;
                    p7 = (i+1)   + (j+1)*num_1d + (k+1)*num_1d*num_1d;

                    p0 += num_1d*num_1d*num_1d*(elem_gid); 
                    p1 += num_1d*num_1d*num_1d*(elem_gid); 
                    p2 += num_1d*num_1d*num_1d*(elem_gid); 
                    p3 += num_1d*num_1d*num_1d*(elem_gid); 
                    p4 += num_1d*num_1d*num_1d*(elem_gid); 
                    p5 += num_1d*num_1d*num_1d*(elem_gid); 
                    p6 += num_1d*num_1d*num_1d*(elem_gid); 
                    p7 += num_1d*num_1d*num_1d*(elem_gid); 

                    cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;

                    cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);

                    // if(cell_mesh_index != elem_gid) std::cout<<"ERROR IN REFINE MESH"<<std::endl;

                    gauss_in_cell(elem_gid, cell_index, 0) = p0;
                    gauss_in_cell(elem_gid, cell_index, 1) = p1;
                    gauss_in_cell(elem_gid, cell_index, 2) = p2;
                    gauss_in_cell(elem_gid, cell_index, 3) = p3;
                    gauss_in_cell(elem_gid, cell_index, 4) = p4;
                    gauss_in_cell(elem_gid, cell_index, 5) = p5;
                    gauss_in_cell(elem_gid, cell_index, 6) = p6;
                    gauss_in_cell(elem_gid, cell_index, 7) = p7;

                    mesh.cells_in_elem(elem_gid, cell_index) = cell_mesh_index;

                    mesh.elems_in_cell(cell_mesh_index) = elem_gid;

                }
            }
        }
    }



    int cell_gid = 0;
    int p[8];
    // for each cell read the list of associated nodes
    for(int elem = 0; elem < num_elem; elem++){
        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){


                    // The p# point to a global gauss point index before double counting 
                    p[0] = (i)     + (j)*num_1d   + (k)*num_1d*num_1d;
                    p[1] = (i+1)   + (j)*num_1d   + (k)*num_1d*num_1d;
                    p[2] = (i)     + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p[3] = (i+1)   + (j+1)*num_1d + (k)*num_1d*num_1d;
                    p[4] = (i)     + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p[5] = (i+1)   + (j)*num_1d   + (k+1)*num_1d*num_1d;
                    p[6] = (i)     + (j+1)*num_1d + (k+1)*num_1d*num_1d;
                    p[7] = (i+1)   + (j+1)*num_1d + (k+1)*num_1d*num_1d;


                    for (int idx = 0; idx < 8; idx++){
                        p[idx] += num_1d*num_1d*num_1d*(elem); 
                    }

                    for (int node_lid = 0; node_lid < 8; node_lid++){
                        mesh.nodes_in_cell(cell_gid, node_lid) = mesh.node_in_gauss(p[node_lid]);
                    }
                    // incriment global index for cell
                    cell_gid++;
                }
            }
        }
    }

    //DANdelete [] gauss_id_in_cell;

    mesh.build_connectivity();
}



} // end namespace swage


swage::mesh_t init_mesh;
swage::mesh_t mesh;
