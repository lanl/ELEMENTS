/*********************************************************************************************
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

/*********************************************************************************************



**********************************************************************************************/


#ifndef ELEMENTS_H
#define ELEMENTS_H 

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utilities.h"
#include "matar.h"

using namespace utils;



namespace elements{

//==============================================================================
//   Function Declaration
//==============================================================================

    // Used by Lobatto 1D/2D to set Lobatto quadrature points
    void labatto_nodes_1D(
        c_array_t <real_t> &lab_nodes_1D,
        const int &num);

    void labatto_weights_1D(
        c_array_t <real_t> &lab_weights_1D,  // Labbatto weights
        const int &num);

    void length_weights(
        c_array_t <real_t> &len_weights_1D,  // Labbatto weights
        c_array_t <real_t> &lab_weights_1D,  // Labbatto weights
        c_array_t <real_t> &lab_nodes_1D,
        const int &order);

    void sub_weights(
        c_array_t <real_t> &sub_weights_1D,  // Labbatto weights
        c_array_t <real_t> &lab_weights_1D,  // Labbatto weights
        c_array_t <real_t> &lab_nodes_1D,
        const int &order);

    void mat_inverse(
        c_array_t <real_t> &mat_inv,
        c_array_t <real_t> &matrix);

    void mat_mult(
        c_array_t <real_t> &result,
        c_array_t <real_t> &A,
        c_array_t <real_t> &B);

    void mat_trans(
        c_array_t <real_t> &trans,
        c_array_t <real_t> &mat);

    void set_nodes_wgts(
        c_array_t <real_t> &lab_nodes_1D,
        c_array_t <real_t> &lab_weights_1D,
        c_array_t <real_t> &len_weights_1D,
        c_array_t <real_t> &sub_weights_1D, 
        int p_order);

    void sub_cells(
        c_array_t <real_t> &lab_nodes_1D,
        int &p_order, 
        int &dim);





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
                       /
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


   face 0: [0,1,3,2]
   face 1: [4,5,7,6]
   face 2: [0,1,5,4]
   face 3: [2,3,7,6]
   face 4: [0,2,6,4]
   face 6; [1,3,7,5]


*/


class mesh_t {
    
private:
    
    // ---- GENERAL INFORMATION ---- //
    const int num_dim_ = 3;
    const int num_nodes_hex_ = 8;
    const int num_faces_hex_ = 6;
    
    const int num_nodes_face_ = 4;
    const int this_node_in_face_in_cell_[24] = // this assumes i,j,k structures nodes
        {0,1,3,2,
         4,5,7,6,
         0,1,5,4,
         2,3,7,6,
         0,2,6,4,
         1,3,7,5};

    int rk_storage_;


// ---- INDEX SPACES AND MAPS ---- //
    
    // ---- ELEMENT ---- //
    int   poly_order_;
    int   num_elem_;
    int   num_g_pts_in_elem_;
    int   num_cells_in_elem_;
    int   elems_in_elem_list_size_;

    int * cells_in_elem_ = NULL;
    int * elems_in_elem_ = NULL;


    // add these_elems_in_elem


    // ---- CELLS ---- //
    int   num_cells_;
    
    int * cell_nodes_list_ = NULL;         // size of num_cells*8
    int * num_cells_in_cell_ = NULL;       // size of num_cells
    int * cell_in_cell_list_start_ = NULL; // size of num_cells+1
    int * cell_in_cell_list_ = NULL;       // size depends on mesh connectivity
    int * elems_in_cell_list_ = NULL;       // size depends on mesh connectivity


    // ---- VERTICES ---- //
    int * num_verts_;


    // ---- NODES ---- //
    int   num_nodes_;

    int * num_cells_in_node_ = NULL;        // size of num_nodes
    int * cells_in_node_list_start_ = NULL; // size of num_nodes+1
    int * cells_in_node_list_ = NULL;       // size depends on mesh connectivity
    

    // ---- GAUSS POINTS ---- //
    int   num_g_pts_;

    int * node_in_gauss_list_ = NULL;


    // ---- CORNERS ---- //
    int   num_corners_;

    int * num_corners_in_node_ = NULL;

    // int * corners_in_cell_list_start;
    int * corners_in_cell_list_ = NULL;

    int * corners_in_node_list_start_ = NULL;
    int * corners_in_node_list_ = NULL;

    // int * corner_bdy_count_ = NULL;


    // ---- FACES ---- //
    int   num_faces_;

    int * face_nodes_list_ = NULL;       // size of num_faces*4
    int * cells_in_face_list_ = NULL;    // size of num_faces*2


    // ---- BOUNDARY ---- //
    int num_bdy_faces_;
    int num_bdy_sets_;
    
    int * bdy_faces_;          // size depends on mesh
    int * bdy_set_list_;
    int * start_index_bdy_set_;
    int * num_bdy_faces_set_;


// ---- MESH GEOMETRIC STATE ---- //

    // ---- ELEMENT ---- //
    real_t * elem_vol_ = NULL;  // size of num_elem


    // ---- CELLS ---- //
    real_t * cell_vol_ = NULL;      // size of num_cells
    real_t * cell_coords_ = NULL;  // size of num_cells


    // ---- NODES ---- //
    real_t * node_coords_ = NULL;      // size of rk_storage_*num_nodes*num_dim_

    

    // ---- QUADRATURE POINTS ---- //
    real_t * jacobians_ = NULL;            // size of num_g_pts_*num_dim_*num_dim_
    real_t * jacobian_determinant_ = NULL; // size of num_g_pts_


public:


// ==== INDEX SPACE FUNCTIONS ==== // 
    
    // returns the number of rk_stages
    inline int num_rk () const
    {
        return rk_storage_;
    }

    // returns the number of dimensions
    inline int num_dim () const
    {
        return num_dim_;
    }


    // ---- ELEMENT ---- //
    void init_element (int p_order, int dim, int num_elem, int num_rk){
        
        poly_order_ = p_order;

        rk_storage_ = num_rk;
        
        int num_g_pts_1d;
        int num_g_pts;
        int num_subcells_per_elem;
        
        if(p_order == 0){
            
            num_g_pts_1d = 2;
            num_g_pts  = pow(num_g_pts_1d, dim);
            num_subcells_per_elem = pow((num_g_pts_1d-1), dim);
        }

        else{

            num_g_pts_1d = 2 * p_order + 1;
            num_g_pts    = pow(num_g_pts_1d, dim);
            num_subcells_per_elem = pow((2 * p_order), dim);
        }


        num_elem_ = num_elem;

        num_g_pts_in_elem_ = num_g_pts;
        num_cells_in_elem_ = num_subcells_per_elem;

        num_cells_ = num_elem * num_subcells_per_elem;

        cells_in_elem_ = new int[num_elem * num_subcells_per_elem];

        elem_vol_ = new real_t[num_rk*num_elem];

    }

    // returns the number of elements
    inline int num_elements () const
    {
        return num_elem_;
    }

    // returns the number of elements
    inline int num_cells_in_element () const
    {
        return num_cells_in_elem_;
    }


    // returns the polynomial order of the element
    inline int poly_order () const
    {
        return poly_order_;
    }

    // return element geometric center coords


    // return number of elements connected to element

    // return array of element id's connected to element

    // return the the global cell id from local element cell id
    inline int& cells_in_element(int elem_gid, int cell_lid) 
    {
        return cells_in_elem_[cell_lid + num_cells_in_elem_*(elem_gid)];
    }

    inline int& num_gauss_in_element() 
    {
        return num_g_pts_in_elem_;
    }

    inline int& num_mat_pt_in_element() 
    {
        return num_g_pts_in_elem_;
    }

    // ---- CELLS ---- //

    // initialize basic cell data structures
    void init_cells (int ncells, int num_rk)
    {

        rk_storage_ = num_rk;
        
        num_cells_ = ncells;

        cell_vol_  = new real_t[num_cells_];
        cell_coords_ = new real_t[num_cells_*num_dim_];


        cell_nodes_list_ = new int[num_cells_*num_nodes_hex_];

        cell_nodes_list_      = new int[num_cells_*num_nodes_hex_];
        corners_in_cell_list_ = new int[num_cells_*num_nodes_hex_];
        
        num_cells_in_cell_       = new int[num_cells_]; 
        cell_in_cell_list_start_ = new int[num_cells_+1];

        elems_in_cell_list_ = new int[num_cells_];

    }
    
    // returns the number of cells
    inline int num_cells () const
    {
        return num_cells_;
    }


    // return the node ids local to the cell
    inline int& cell_nodes_id(int cell_gid, int node_lid) const
    {
        return cell_nodes_list_[node_lid + cell_gid*num_nodes_hex_];
    }


    inline auto cell_nodes_ids(int cell_gid)
    {
        auto give = view_c_array <int> (&cell_nodes_list_[cell_gid*num_nodes_hex_], num_nodes_hex_);
        return give;
    }
    


    // return the number of cells around the cell
    inline int& num_cells_in_cell(int cell_gid) const
    {
        return num_cells_in_cell_[cell_gid];
    }

    // return the the cells around a cell
    inline int& cell_in_cell(int cell_gid, int cell_lid) const
    {
        // shift index by 1 so that it is consistent with matrix syntax
        int start_indx = cell_in_cell_list_start_[cell_gid];
        
        // get the index in the global 1D array
        int index = start_indx + cell_lid;
        
        return cell_in_cell_list_[index];
    }

    inline int& cell_corners(int cell_gid, int corner_lid) const
    {
        return corners_in_cell_list_[corner_lid + cell_gid*num_nodes_hex_];
    }

    inline int& elems_in_cell(int cell_gid) const
    {
        return elems_in_cell_list_[cell_gid];
    }
        
        

    // ---- VERTICES ---- //
    // returns the number of vertices in mesh
    // inline int& num_verts () const
    // {
    //     return num_verts_;
    // }

    // ---- NODES ---- //

    // initialize the nodes
    void init_nodes (int num_nodes, int num_rk) {

        rk_storage_ = num_rk;
        
        num_nodes_ = num_nodes;
        node_coords_   = new real_t[num_rk*num_nodes_*num_dim_];
        // node_coords__rk = new real_t[num_rk*num_nodes_*num_dim_];
        
        num_cells_in_node_        = new int[num_nodes_];
        cells_in_node_list_start_ = new int[num_nodes_+1];

        num_corners_in_node_      = new int[num_nodes_];

    }

    // returns the number of vertices
    inline int num_nodes ()
    {
        return num_nodes_;
    }

    // returns number of cells around a node
    inline int& num_cells_in_node(int node_gid) const
    {
        return num_cells_in_node_[node_gid];
    }
    
    // return the cells around a node
    inline int& cells_in_node (int node_gid, int cell_lid) const
    {
        // shift index by 1 so that it is consistent with matrix syntax
        int start_indx = cells_in_node_list_start_[node_gid];
        
        // get the index in the global 1D array
        int index = start_indx + cell_lid;
        
        return cells_in_node_list_[index];
    }
        

    // ---- GAUSS POINTS ---- //

    // Initialize gauss point data
    void init_gauss_pts (){

        // Index maps
        num_g_pts_ = num_elem_ * num_g_pts_in_elem_;
        node_in_gauss_list_ = new int[num_g_pts_];

        // geometric state
        jacobians_ = new real_t[rk_storage_*num_g_pts_*num_dim_*num_dim_];
        jacobian_determinant_ = new real_t[rk_storage_*num_g_pts_];

    }

    // return number of gauss points in mesh
    inline int num_gauss () const
    {
        return num_g_pts_;
    }

    // return gauss to node map
    inline int& node_in_gauss (int gauss_gid) const
    {            
        return node_in_gauss_list_[gauss_gid];
    }


    // ---- CORNERS ---- //

    // initialize corner data

    
    // returns the number of ccorners
    inline int num_corners () const
    {
        return num_corners_;
    }

    // return corner to cell map

    inline int corners_in_cell (int cell_gid, int corner_lid) const
    {  
        return corners_in_cell_list_[cell_gid*num_nodes_hex_ + corner_lid];
    }

    // return corner to node map
    inline int corners_in_node (int node_gid, int corner_lid) const
    {
        // Note: cell_in_node_list_start_ is the exact same as 
        //       corner_in_node_list_start_
        int start_indx = cells_in_node_list_start_[node_gid];

        // get the index in the global 1D array
        int index = start_indx + corner_lid;

        return corners_in_node_list_[index];
    }

        // return corner to node map
    inline int num_corners_in_node (int node_gid) const
    {

        return num_corners_in_node_[node_gid];
    }


    // inline int corner_bdy_count(int corner_gid) const
    // {
    //     return corner_bdy_count_[node_gid];
    // }


    // ---- FACES ---- //


    // returns the number of elements
    inline int num_faces () const
    {
        return num_faces_;
    }

    // returns the global node id given a cell_id, local_face_indx(0:5), local_facenode_indx(0:3)
    inline int node_in_face_in_cell(int cell_id, int this_face, int facenode_lid) const
    {
        // get the local node index in the cell
        int this_node = this_node_in_face_in_cell_[facenode_lid + this_face*num_nodes_face_];
        
        // return the global id for the local node index
        return cell_nodes_id(cell_id, this_node);
    } // end of method
    
    
    // returns the global id for a cell that is connected to the face
    inline int cells_in_face(int face_gid, int this_cell) const
    {
        // get the 1D index
        int this_index = face_gid*2 + this_cell;  // this_cell = 0 or 1

        // std::cout <<"cells in face index = " << this_index << std::endl;
        
        // return the global id for the cell
        return cells_in_face_list_[this_index];

    }
          
    // returns the nodes in the face
    inline int node_in_face(int face_gid, int facenode_lid) const
    {
        
        // get the 1D index
        int this_index = face_gid*4;
        
        // facenode_lid is in the range of 0:3
        return face_nodes_list_[this_index + facenode_lid];
    }




    // ---- Boundary ---- //
    // initializes the number of bdy sets
    void init_bdy_sets (int num_sets){
        
        // A check
        if(num_sets == 0){
            std::cout << " ERROR: number of boundary sets = 0, setting it = 1";
            num_sets = 1;
        }
        num_bdy_sets_ = num_sets;
        num_bdy_faces_set_ = new int [num_sets];
        start_index_bdy_set_ = new int [num_sets+1];
        bdy_set_list_   = new int [num_sets*num_bdy_faces_]; // largest size possible
    }
    
    inline int num_bdy_sets() const
    {
        return num_bdy_sets_;
    }

    inline int bdy_faces(int this_bdy_face) const
    {
        return bdy_faces_[this_bdy_face];
    }
    
    
    inline int num_bdy_faces() const
    {
        return num_bdy_faces_;
    }


    // returns the number per bdy-faces in a particular set
    int num_bdy_faces_in_set (int bdy_set){
        return num_bdy_faces_set_[bdy_set];
    }
    
    
    // returns a subset of the boundary faces
    int bdy_faces_in_set (int bdy_set, int this_face){
        
        int start = start_index_bdy_set_[bdy_set];
        
        return bdy_set_list_[start+this_face];
    }



// ==== MESH STATE FUNCTIONS ==== // 

    // ---- ELEMENTS ---- //
    inline real_t& elem_vol(int rk_stage, int elem_gid) const
    {
        return elem_vol_[rk_stage*num_elem_ + elem_gid];
    }
    

    // ---- CELLS ---- //
    
    // return the cell volume
    inline real_t& cell_vol(int cell_gid) const
    {
        return cell_vol_[cell_gid];
    }
    
    // return the cell coordinate position
    inline real_t& cell_coords(int rk_stage, int cell_gid, int this_dim)
    {
        
        cell_coords_[this_dim + cell_gid*num_dim_] = 0;
        
        // #pragma omp simd
        for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
            
            int node_gid = cell_nodes_id(cell_gid, node_lid); // get the global_id
            
            cell_coords_[this_dim + cell_gid*num_dim_] += node_coords(rk_stage, node_gid, this_dim);
            
        } // end for loop over vertices in the cell
        
        // divide by the number of vertices
        cell_coords_[this_dim + cell_gid*num_dim_] /= ( (real_t)num_nodes_hex_ );
        
        return cell_coords_[cell_gid*num_dim_ + this_dim];
    }



    // ---- NODES ---- //
    // return the node coordinates
    inline real_t& node_coords(int rk_stage, int node_gid, int this_dim) const
    {
        return node_coords_[rk_stage*num_nodes_*num_dim_ + node_gid*num_dim_ + this_dim];
    }


    // // return the node coordinates at t=n
    // inline real_t& node_coords_rk(int rk_step, int node_gid, int this_dim) const
    // {
    //     return node_coords__rk[rk_step*num_nodes_*num_dim_ + node_gid*num_dim_ + this_dim];
    // }


    // ---- FACES ---- //
    // geometric average face coordinate
    inline real_t face_coords(int rk_stage, int face_id, int this_dim) const
    {
        
        real_t this_face_coord = 0.0;

        // loop over all the vertices on the this face
        for (int this_facevert = 0; this_facevert < num_nodes_face_; this_facevert++){
            
            // get the global vertex id
            int vert_gid = node_in_face(face_id, this_facevert);
            
            // calc the coord
            this_face_coord += node_coords(rk_stage, vert_gid, this_dim)/((real_t)num_nodes_face_);

        } // end for this_facevert
        
        return this_face_coord;
        
    } // end of face_coords

    // ---- QUADRATURE POINTS ---- //

    // return jacobian at quadrature point
    inline real_t& jacobian(int rk_stage, int elem_gid, int gauss_lid, int i, int j) const
    {
        int index = rk_stage*num_elem_*num_g_pts_in_elem_*num_dim_*num_dim_
                    + elem_gid*num_g_pts_in_elem_*num_dim_*num_dim_
                    + gauss_lid*num_dim_*num_dim_
                    + i*num_dim_
                    + j;
        
        return jacobians_[index];
    }


    // return determinant of jacobian at quadrature point
    inline real_t& det_j(int rk_stage, int elem_gid, int gauss_lid) const
    {
        int index = rk_stage*num_elem_*num_g_pts_in_elem_
                    + elem_gid*num_g_pts_in_elem_
                    + gauss_lid;
        
        return jacobian_determinant_[index];
    }





// ==== MESH CONNECTIVITY FUNCTIONS ==== // 
    
    // initialize array for mesh connectivity: all cells around a node
    void build_connectity(){
        
        
        // -- NODE TO CELL CONNECTIVITY -- //
        build_node_cell_connectivity(); 

        // -- CORNER CONNECTIVITY -- //
        build_corner_connectivity(); 

        // -- CELL TO CELL CONNECTIVITY -- //
        build_cell_cell_connectivity(); 

        // -- FACES -- //
        build_face_connectivity(); 
    
    } // end of build connectivity


    void build_node_cell_connectivity(){

        // initializing the number of corners (node-cell pair) to be zero
        int num_corners_saved[num_nodes_]; // local variable
        for (int node_gid = 0; node_gid < num_nodes_; node_gid++){
            
            num_cells_in_node_[node_gid] = 0;
            num_corners_saved[node_gid] = 0;
        }
        
        
        // count the number of corners (cell-node pair) in the mesh and in each node
        int num_corners = 0;
        for (int cell_id = 0; cell_id < num_cells_; cell_id++){
            for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
                
                // each point-cell pair makes a corner
                num_corners ++;  // total number of corners in the entire mesh
                
                // increment the number of corners attached to this point
                int node_gid = cell_nodes_id(cell_id, node_lid); // get the global_id
                
                num_cells_in_node_[node_gid] ++;
                
            }  // end for this_point
        } // end for cell_id
        
        // Save number of corners in mesh
        num_corners_ = num_corners;

        // create memory for a list for all cell-node pairs
        cells_in_node_list_ = new int [num_corners];
        
        cells_in_node_list_start_[0] = 0;
        for (int node_gid = 1; node_gid < num_nodes_; node_gid++){
            // note, for loops starts at 1 becuase we just saved the point 0
            
            // This is the start of the indices for the corners connected to this node
            cells_in_node_list_start_[node_gid] = cells_in_node_list_start_[node_gid-1]
                                                    + num_cells_in_node_[node_gid-1];
        }
        
        
        // populate the cells connected to a node list
        int corner_gid = 0;
        
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
                
                // increment the number of corners attached to this point
                int node_gid = cell_nodes_id(cell_gid, node_lid); // get the global_id of the node

                // Assign global index values to cell_nodes_list_
                cell_nodes_list_[node_lid + cell_gid*num_nodes_hex_] = node_gid;
                
                // the global index in the cells_in_node_list_
                int index = cells_in_node_list_start_[node_gid] + num_corners_saved[node_gid];
                
                // save the global cell index to the list
                cells_in_node_list_[index] = cell_gid;
                
                // each point-cell pair makes a corner
                num_corners_saved[node_gid] ++;  //number of corners saved to this node
                
            }  // end for this_point
        } // end for cell_gid
    } // end of build_node_cell_connectivity


    void build_corner_connectivity(){

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
        corners_in_node_list_ = new int [num_corners];

        // populate the cells connected to a node list and corners in a node
        int corner_gid = 0;
        
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
                
                // increment the number of corners attached to this point
                int node_gid = cell_nodes_id(cell_gid, node_lid); // get the global_id of the node
                
                // the global index in the cells_in_node_list_
                int index = cells_in_node_list_start_[node_gid] + num_corners_saved[node_gid];

                
                // each point-cell pair makes a corner
                num_corners_saved[node_gid] ++;  //number of corners saved to this node

                // Save index for corner to global cell index
                corners_in_node_list_[index] = corner_gid;

                int corner_lid = node_lid;
                cell_corners(cell_gid, corner_lid) = corner_gid;

                corner_gid ++;
                
            }  // end for this_point
        } // end for cell_gid

        for (int node_gid = 0; node_gid < num_nodes_; node_gid++ ){

            num_corners_in_node_[node_gid] = num_corners_saved[node_gid];
        
        }

    } // end of build_corner_connectivity


    void build_face_connectivity(){
        
        int face_gid = 0;
        real_t node_hash_delta = 1.0e64;
        
        real_t coord_min[num_dim_];
        real_t coord_max[num_dim_];

        // initialize to large values
        for (int this_dim = 0; this_dim < num_dim_; this_dim++){
            
            coord_min[this_dim] = 1.0e64;
            coord_max[this_dim] = -1.0e64;

        } // end for dim


        // Get min and max points in the mesh
        real_t coords[num_dim_];
        for(int point = 0; point < num_nodes_; point++){
            for (int this_dim = 0; this_dim < num_dim_; this_dim++){

                coords[this_dim] = node_coords(0, point, this_dim);

                coord_max[this_dim] = fmax(coord_max[this_dim], coords[this_dim]);
                coord_min[this_dim] = fmin(coord_min[this_dim], coords[this_dim]);

            }
        }

        // get minimum distance between any two points (WARNING: ONLY WORKS IN 3D)
        real_t dist_min;
        real_t dist_max;
        real_t cell_nodes[24];
        
        auto vert1 = view_c_array <real_t> (cell_nodes, num_nodes_hex_, 3);
        
        real_t distance[28]; 
        auto dist = view_c_array <real_t> (distance, 28);

        for (int cell_id=0; cell_id<num_cells_; cell_id++){
            
            // Getting the coordinates of the element
            for(int node = 0; node < num_nodes_hex_; node++){
                for (int dim = 0; dim < 3; dim++)
                    vert1(node, dim) = node_coords(0, cell_nodes_id(cell_id, node), dim);
            }

            // loop conditions needed for distance calculation
            int countA = 0;
            int countB = 1;
            int a;
            int b;
            int loop = 0;
            
            
            // Solving for the magnitude of distance between each node
            for (int i=0; i<28; i++){
                
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

        node_hash_delta = fmin(node_hash_delta, dist_min); // 0.9 used to offset face centers from nodes
        node_hash_delta = node_hash_delta/2.0;

        // calculate the 1d array length of the hash table for nodes
        real_t num_bins[num_dim_];
        
        for (int this_dim = 0; this_dim < num_dim_; this_dim++){
            num_bins[this_dim] = (coord_max[this_dim] - coord_min[this_dim]  + 1.e-12)/node_hash_delta;
        }


        int hash_keys[num_cells_*num_faces_hex_];

        real_t face_hash_idx_real[num_dim_];

        
        // Calculate the hash keys and find max key
        int hash_count = 0;
        int max_key = 0;
        
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int face_lid = 0; face_lid < num_faces_hex_; face_lid++){
                
                // the face coordinates
                real_t face_coords[num_dim_];
                
                for (int this_dim=0; this_dim<num_dim_; this_dim++){
                    face_coords[this_dim] = 0.0;
                } // end for dim
                
                
                // loop over all the nodes on the this face
                for (int facenode_lid = 0; facenode_lid < num_nodes_face_; facenode_lid++){
                    
                    // get the global node id
                    int node_gid = node_in_face_in_cell(cell_gid, face_lid, facenode_lid);
                    
                    for (int this_dim = 0; this_dim < num_dim_; this_dim++){
                        face_coords[this_dim] += node_coords(0, node_gid, this_dim)/num_nodes_face_;
                    } // end for dim
                        
                } // end for facenode_lid
                
                // calculate the face hash index for these face_coordinates
                for (int this_dim = 0; this_dim < num_dim_; this_dim++){
                    face_hash_idx_real[this_dim] = fmax(1e-16, (face_coords[this_dim]-coord_min[this_dim] + 1e-12)/node_hash_delta);
                } // end for dim

                // the 1D index
                real_t face_hash_idx_real_1d;
                
                // i + j*num_x + k*num_x*num_y
                if (num_dim_ == 2){
                    face_hash_idx_real_1d = face_hash_idx_real[0] + face_hash_idx_real[1]*num_bins[0];
                }

                else{
                    face_hash_idx_real_1d =
                          face_hash_idx_real[0]
                        + face_hash_idx_real[1]*num_bins[0]
                        + face_hash_idx_real[2]*num_bins[0]*num_bins[1];
                }
                
                hash_keys[hash_count] = (int)face_hash_idx_real_1d;
                max_key = std::max(max_key, hash_keys[hash_count]);
                hash_count++;

            } // end for face_lid
        } // end for cell_gid


        // Allocate hash table and initialize key locations to -1
        int * hash_table = new int [max_key + 10];
        
        for (int key_idx = 0; key_idx < hash_count; key_idx++){
            hash_table[hash_keys[key_idx]] = -1;
        } // end for cell_gid


        // count the number of cells around each face
        face_gid = 0;
        hash_count = 0;
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int face_lid = 0; face_lid < num_faces_hex_; face_lid++){
                
                // count the number of faces
                if (hash_table[hash_keys[hash_count]] == -1){
                    // save the face_id to the hash table
                    hash_table[hash_keys[hash_count]] = face_gid;
                    face_gid++;
                }

                hash_count++;
            } // end for face_lid
        } // end for cell_gid

        // set the number of faces in the mesh
        num_faces_ = face_gid;

        // create memory for the face structures
        cells_in_face_list_ = new int [num_faces_*2];   // two cells per face
        face_nodes_list_    = new int [num_faces_*4];   // four nodes per face in hex

        // initialize the cells_in_face_list to -1
        for (int face_gid = 0; face_gid < 2*num_faces_; face_gid++){
            cells_in_face_list_[face_gid] = -1;
        }   

        for (int face_gid = 0; face_gid < 4*num_faces_; face_gid++){
            face_nodes_list_[face_gid] = -1;
        }


        hash_count = 0;
        // save the cell_gid to the cells_in_face_list
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int face_lid = 0; face_lid < num_faces_hex_; face_lid++){
                
                // get the face_id
                face_gid = hash_table[hash_keys[hash_count]];

                // a temp variable for storing the node global ids on the face
                int these_nodes[num_nodes_face_];

                // loop over all the vertices on the this face
                for (int facenode_lid = 0; facenode_lid < num_nodes_face_; facenode_lid++){
                    
                    // get the global node id
                    int node_gid = node_in_face_in_cell(cell_gid, face_lid, facenode_lid);
                    
                    // save the global node id
                    these_nodes[facenode_lid] = node_gid;
                    
                } // end for facenode_lid

                // save the cell_gid to the cells_in_face_list
                if (cells_in_face_list_[face_gid*2] == -1){
                    
                    // save the cell index for this face
                    int this_index = face_gid*2;  // index in the list
                    
                    cells_in_face_list_[this_index] = cell_gid;
                    
                    // save the nodes for this face
                    this_index = face_gid*4;
                    
                    for (int facenode_lid = 0; facenode_lid < num_nodes_face_; facenode_lid++){
                        face_nodes_list_[this_index + facenode_lid] = these_nodes[facenode_lid];
                    }
                }

                else{
                    // it is the other cell connected to this face

                    int this_index = face_gid*2 + 1; // + num_faces_; // index in the list

                    cells_in_face_list_[this_index] = cell_gid;
                }

                hash_count++;

            } // end for face_lid
        } // end for cell_gid


        delete[] hash_table;
    } // end of build faces


    void build_cell_cell_connectivity(){
        
        // initializing the number of cell-cell pairs to be zero
        int num_cell_cell_saved[num_cells_]; // local variable
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            
            num_cells_in_cell_[cell_gid] = 0;
            num_cell_cell_saved[cell_gid] = 0;
        }
        
        int num_c_c_pairs = 0;
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int node_lid = 0; node_lid < num_nodes_hex_; node_lid++){
                
                // get the global node id
                int node_gid = cell_nodes_id(cell_gid, node_lid);
                
                // loop over all cells connected to node_gid
                for (int cell_lid = 0; cell_lid < num_cells_in_node_[node_gid]; cell_lid++){
                    
                    int neighbor_cell_id = cells_in_node(node_gid, cell_lid);
                    
                    // a true neighbor_cell_id is not equal to cell_gid
                    if (neighbor_cell_id != cell_gid){
                        
                        // increment the number of cell-cell pairs in the mesh
                        num_c_c_pairs ++;
                        
                        // increment the number of cell-cell pairs for this cell
                        num_cells_in_cell_[cell_gid] ++;
                        
                    } // end if neighbor_cell_id
                } // end for cell_lid
                
                
            }  // end for node_lid
        } // end for cell_gid
        
        
        // create memory for the list of cells around a cell (num_c_c_pairs is 2x larger than needed)
        int * temp_cell_in_cell_list = new int [num_c_c_pairs];

        
        cell_in_cell_list_start_[0] = 0;
        for (int cell_gid = 1; cell_gid < num_cells_; cell_gid++){
            // note, for loops starts at 1 becuase we just saved the cell 0
            
            // This is the start of the indices for the corners connected to this node
            cell_in_cell_list_start_[cell_gid] = cell_in_cell_list_start_[cell_gid-1]
                                                   + num_cells_in_cell_[cell_gid-1];
        }
        
        
        
        
        int actual_size = 0; // the actual size of array of cell-neighbors
        for (int cell_gid = 0; cell_gid < num_cells_; cell_gid++){
            for (int this_node = 0; this_node < num_nodes_hex_; this_node++){
                
                // get the global_id node id
                int node_id = cell_nodes_id(cell_gid, this_node);
                
                // loop over all cells connected to node_id
                for (int this_cell = 0; this_cell < num_cells_in_node_[node_id]; this_cell++){
                    
                    // get the global id for the neighboring cell
                    int neighbor_cell_id = cells_in_node(node_id, this_cell);
                    
                    // the global index in the cell_in_cell_list_
                    int index = cell_in_cell_list_start_[cell_gid] + num_cell_cell_saved[cell_gid];
                    
                    int save = 1; // a flag to save (=1) or not (=0)
                    
                    // a true neighbor_cell_id is not equal to cell_gid
                    if (neighbor_cell_id == cell_gid ){
                        save = 0;  // don't save
                    } // end if neighbor_cell_id
                    
                    // check to see if neighbor_id has been saved already
                    for (int i=cell_in_cell_list_start_[cell_gid]; i<index; i++){
                        
                        if (neighbor_cell_id == temp_cell_in_cell_list[i]){
                            save=0;   // don't save, it has been saved already
                        } // end if
                        
                    } // end for i
                    
                    
                    if (save==1){
                        // save the neighboring cell_id
                        temp_cell_in_cell_list[index] = neighbor_cell_id;
                        
                        // increment the number of neighboring cells saved
                        num_cell_cell_saved[cell_gid]++;
                        
                        // increment the actual number of neighboring cells saved
                        actual_size++;
                    } // end if save
                    
                } // end for this_cell
                
                
            }  // end for this_node
        } // end for cell_gid
        
        
        // create memory for the list of cells around a cell
        cell_in_cell_list_ = new int [actual_size];

        elems_in_elem_list_size_ = actual_size;
        
        // update the number of cells around a cell (because the estimate had duplicates)
        int index=0;
        cell_in_cell_list_[0] = 0;
        
        for (int cell_id=0; cell_id<num_cells_; cell_id++){
            
            num_cells_in_cell_[cell_id] = num_cell_cell_saved[cell_id];
            
            // note: there is a buffer in the temp_cell_in_cell_list that is not
            // in num_cells_in_cell_ (it is smaller) so I will copying the
            // values from the temp to the actual array and resetting the start
            // array to use num_cell_cell_saved[cell_id]
            
            // the global index range in the temp_cell_in_cell_list
            int list_start = cell_in_cell_list_start_[cell_id];
            int list_stop  = list_start + num_cell_cell_saved[cell_id];
            
            cell_in_cell_list_start_[cell_id] = index; // update the index
            
            // save the neighbors to the list
            for (int i = list_start; i < list_stop; i++){
                
                // save neighboring cell_id to the final list
                cell_in_cell_list_[index] = temp_cell_in_cell_list[i];
                
                // increment the global index
                index++;
                
            } // end for i
            
            
        }// end for cell_id
        
        
        // delete the temporary list of cells around a cell
        delete[] temp_cell_in_cell_list;
    } // end of build_cell_cell_connectivity
    

    // identify the boundary faces
    void build_bdy_faces (){
        
        
        int bdy_face_gid = 0;
        
        
        // loop over the faces in the mesh
        for(int face_gid = 0; face_gid < num_faces_; face_gid++){
            
            // loop over the two cells on this face
            for (int this_cell = 0; this_cell < 2; this_cell++){
                
                // check to see if a cell has index of -1
                if (cells_in_face(face_gid, this_cell) == -1){
                   bdy_face_gid ++;
                } // end if
                    
            } // end for this_cell
        } // end for face_gid
        
        
        // save the number of boundary faces in the mesh
        num_bdy_faces_ = bdy_face_gid;
        
        // allocate the memory for the boundary faces array
        bdy_faces_ = new int[num_bdy_faces_];
        
        
        // save the global indices for the boundary faces
        
        // loop over the faces in the mesh
        bdy_face_gid = 0;

        for(int face_gid = 0; face_gid < num_faces_; face_gid++){
            
            // loop over the two cells on this face
            for (int this_cell = 0; this_cell < 2; this_cell++){
                
                // check to see if a cell has index of -1
                if (cells_in_face(face_gid, this_cell) == -1){
                    
                    // save the face index
                    bdy_faces_[bdy_face_gid] = face_gid;
                    
                    // increment the counter
                    bdy_face_gid++;
                    
                } // end if  
            } // end for this_cell
        } // end for face_gid
    } // end of function













 // ---- bdy sets ----
    
    // returns a subset of the boundary faces
    int set_bdy_faces (int bdy_set, int this_face){
        
        int start = start_index_bdy_set_[bdy_set];
        
        return bdy_set_list_[start+this_face];
    }
    
    

    // set planes for tagging sub sets of boundary faces
    // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    // val = plane value, radius, radius
    void tag_bdys(int this_bc_tag, real_t val, int bdy_set){
        
        if (bdy_set == num_bdy_sets_){
            std::cout << " ERROR: number of boundary sets must be increased by "
                      << bdy_set-num_bdy_sets_+1 << std::endl;
            exit(0);
        }
        
        // the start index for the first list is zero
        start_index_bdy_set_[0] = 0;
        
        
        // save the boundary vertices to this set that are on the plane
        int counter = 0;
        int start = start_index_bdy_set_[bdy_set];
        for (int this_bdy_face = 0; this_bdy_face < num_bdy_faces_; this_bdy_face++) {
            
            // save the face index
            int bdy_face_gid = bdy_faces(this_bdy_face);
            
            // check to see if this face is on the specified plane
            int is_on_bdy = check_bdy(bdy_face_gid, this_bc_tag, val); // no=0, yes=1
            
            if (is_on_bdy == 1){
                bdy_set_list_[start+counter] = bdy_face_gid;
                counter ++;
            }
            
        } // end for bdy_face
        
        
        // save the number of bdy faces in the set
        num_bdy_faces_set_[bdy_set] = counter;
        
        // save the starting index for the next bdy_set
        start_index_bdy_set_[bdy_set+1] = start_index_bdy_set_[bdy_set] + counter;
        
        
        // compress the list to reduce the memory if it is the last set
        if (bdy_set == num_bdy_sets_-1){
            compress_bdy_set();
        }
        
        std::cout << " tag boundary faces " << std::endl;
        
    } // end of method
    
    
    
    
    // compress the bdy_set_list to reduce the memory
    void compress_bdy_set(){
        
        // the actual size of the bdy set list
        int length = start_index_bdy_set_[num_bdy_sets_];
        
        // create a temp array of correct size
        int * temp_bdy_list = new int [length];
        
        // save the values to the temp array
        for (int i=0; i<length; i++){
            temp_bdy_list[i] = bdy_set_list_[i];
        }
        
        // delete original array and make a new one of correct size
        delete[] bdy_set_list_;
        bdy_set_list_ = new int [length];
        
        // save the values to the bdy_set_list
        for (int i=0; i<length; i++){
            bdy_set_list_[i] = temp_bdy_list[i];
        }
        
        // delete the temp array
        delete[] temp_bdy_list;
        
    } // end of method
    

    // routine for checking to see if a vertix is on a boundary
    // bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
    // val = plane value, radius, radius
    int check_bdy(int face_gid, int this_bc_tag, real_t val){
        
        // default bool is not on the boundary
        int is_on_bdy = 0;
        
        
        // the face coordinates
        real_t these_face_coords[num_dim_];
        
        for (int this_dim = 0; this_dim < num_dim_; this_dim++){
            these_face_coords[this_dim] = face_coords(0, face_gid, this_dim);
        } // end for dim
        
        
        // a x-plane
        if (this_bc_tag == 0){
            
            if ( fabs(these_face_coords[0] - val) <= 1.0e-8 ) is_on_bdy = 1;
            
        }// end if on type
        
        // a y-plane
        else if (this_bc_tag == 1){
            
            if ( fabs(these_face_coords[1] - val) <= 1.0e-8 ) is_on_bdy = 1;
            
        }// end if on type
        
        // a z-plane
        else if (this_bc_tag == 2){
            
            if ( fabs(these_face_coords[2] - val) <= 1.0e-8 ) is_on_bdy = 1;
            
        }// end if on type
        
        
        // cylinderical shell where radius = sqrt(x^2 + y^2)
        else if (this_bc_tag == 3){
            
            real_t R = sqrt(these_face_coords[0]*these_face_coords[0] +
                            these_face_coords[1]*these_face_coords[1]);
            
            if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;

            
        }// end if on type
        
        // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
        else if (this_bc_tag == 4){
            
            real_t R = sqrt(these_face_coords[0]*these_face_coords[0] +
                            these_face_coords[1]*these_face_coords[1] +
                            these_face_coords[2]*these_face_coords[2]);
            
            if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;
            
        } // end if on type
        
        return is_on_bdy;
        
    } // end method to check bdy

















    // deconstructor
    ~mesh_t ( ) {
        
        // ---- ELEMENTS ---- //
        delete[] cells_in_elem_;


        // ---- CELLS ---- //

        delete[] cell_nodes_list_;
        delete[] num_cells_in_cell_;
        delete[] cell_in_cell_list_start_;
        delete[] cell_in_cell_list_;


        // ---- VERTICES ---- //

        
        // ---- NODES ---- //
        delete[] num_cells_in_node_;
        delete[] cells_in_node_list_start_;
        delete[] cells_in_node_list_;


        // ---- GAUSS POINTS ---- //
        delete[] node_in_gauss_list_;


        // ---- CORNERS ---- //
        delete[] corners_in_cell_list_;
        delete[] corners_in_node_list_start_;
        delete[] corners_in_node_list_;
        // delete[] corner_bdy_count_;


        // ---- FACES ---- //
        delete[] face_nodes_list_;       
        delete[] cells_in_face_list_;  

        // ---- BOUNDARY ---- //
        delete[] bdy_faces_;
        delete[] bdy_set_list_;
        delete[] start_index_bdy_set_;
        delete[] num_bdy_faces_set_; 

    // ---- MESH STATE ---- //
        // ---- ELEMENT ---- //
        delete[] elem_vol_; 


        // ---- CELLS ---- //
        delete[] cell_vol_;
        delete[] cell_coords_;


        // ---- NODES ---- //
        delete[] node_coords_;


        // ---- QUADRATURE POINTS ---- //
        delete[] jacobians_;            // size of rk_storage_*num_g_pts_*num_dim_*num_dim_
        delete[] jacobian_determinant_; // size of rk_storage_*num_g_pts_

    } // end of mesh deconstructor

}; // end of mesh_t



void refine_mesh(
    mesh_t& init_mesh, 
    mesh_t& mesh, 
    const int p_order, 
    const int rk_num_stages,
    const int dim);



namespace swage{

    // Used by Gauss2/3D to set quadrature points
    void line_gauss_info(
        real_t &x, real_t &w, 
        int &m,  int &p);

    // Used by Lovatto 1D/2D to set Lobatto quadrature points
    void line_lobatto_info(
        real_t &x, 
        real_t &w, 
        int &m, 
        int &p);

    // setting gauss quadrature points for 2D elements
    void gauss_2d(
        view_c_array <real_t> &these_g_pts,     // gauss points
        view_c_array <real_t> &these_weights,   // gauss weights
        view_c_array <real_t> &tot_g_weight,    // 2D product of gauss weights
        int &quad_order);                       // quadrature order (n)

    // setting gauss quadrature points for 2D elements
    void gauss_3d(
        view_c_array <real_t> &these_g_pts,   // gauss points
        view_c_array <real_t> &these_weights, // gauss weights
        view_c_array <real_t> &tot_g_weight,            // 3D product of gauss weights
        int &quad_order);                        // quadrature order (n)

    // setting gauss quadrature points for 4D elements
    void gauss_4d(
        view_c_array <real_t> &these_g_pts, // gauss points
        view_c_array <real_t> &these_weights, // gauss weights
        int &quad_order, // quadrature order (n)
        const int &dim);

    // setting Gauss-Lobatto quadrature points for 2D elements
    void lobatto_2d(
        view_c_array <real_t> &these_L_pts, // gauss points
        view_c_array <real_t> &these_weights, // gauss weights
        int &quad_order); // quadrature order (n)

    // setting Gauss-Lobatto quadrature points for 3D elements
    void lobatto_3d(
        view_c_array <real_t> &these_L_pts, // gauss points
        view_c_array <real_t> &these_weights, // gauss weights
        int &quad_order); 

    // setting gauss quadrature points for 4D elements
    void lobatto_4d(
        view_c_array <real_t> &these_L_pts, // gauss points
        view_c_array <real_t> &these_weights, // gauss weights
        int &quad_order, // quadrature order (n)
        const int &dim);

    //defining the jacobian for 2d elements
    void jacobian_2d(
        view_c_array <real_t> &J_matrix, 
        real_t &det_J,
        const view_c_array <real_t> &vertices, 
        const view_c_array <real_t> &this_partial,
        const int &num_nodes);


    void jacobian_3d(
        view_c_array <real_t> &J_matrix, 
        real_t &det_J,
        const view_c_array <real_t> &vertices, 
        const view_c_array <real_t> &this_partial,
        const int &num_nodes);

    //defining the jacobian for 4D elements
    void jacobian_4d(
        view_c_array <real_t> &J_matrix, 
        real_t &det_J,
        const view_c_array <real_t> &vertices, 
        const view_c_array <real_t> &this_partial,
        const int &num_nodes,
        const int &dim);

    //defining the inverse jacobian for 2D element    
    void jacobian_inverse_2d(
        view_c_array <real_t> &J_inverse, 
        const view_c_array <real_t> &jacobian);

    //defining the inverse jacobian for 2D element    
    void jacobian_inverse_3d(
        view_c_array <real_t> &J_inverse_matrix,
        const view_c_array <real_t> &jacobian);

    // defining  the inverse of the Jacobian for 4D elements
    void jacobian_inverse_4d(
        view_c_array <real_t> &J_inverse_matrix,
        const view_c_array <real_t> &jacobian,
        const real_t &det_J);

    // creates nodal positions with Chebyshev spacing
    void chebyshev_nodes_1D(
        view_c_array <real_t> &cheb_nodes_1D,   // Chebyshev nodes
        const int &order);                      // Interpolation order


    class Element2D {
        protected:
            const static int num_dim = 2;
      
        public:

        // calculate a physical position in an element for a given xi,eta
        virtual void physical_position(
            view_c_array <real_t>  &x_point,
            const view_c_array <real_t>  &xi_point,
            const view_c_array <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta
        virtual void basis(
            view_c_array <real_t>  &basis,
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi
        virtual void  partial_xi_basis(
            view_c_array <real_t>  &partial_xi, 
            const view_c_array <real_t>  &xi_point) = 0;


        // Partial derivative of shape functions with respect to Xi
        virtual void  partial_eta_basis(
            view_c_array <real_t> &partial_eta, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Map from vertex to node
        virtual inline int vert_node_map( const int vert_lid) = 0;

    }; // end of 2D element class

    class Element3D {
        protected:
            const static int num_dim = 3;

        public:

        // calculate a physical position in an element for a given xi,eta,mu
        virtual void physical_position(
            view_c_array <real_t>  &x_point,
            const view_c_array <real_t>  &xi_point,
            const view_c_array <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta, mu
        virtual void basis(
            view_c_array <real_t>  &basis,
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi at Xi_point
        virtual void partial_xi_basis(
            view_c_array <real_t>  &partial_xi, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Eta
        virtual void partial_eta_basis(
            view_c_array <real_t> &partial_eta, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Mu
        virtual void partial_mu_basis(
            view_c_array <real_t> &partial_mu, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Map from vertex to node
        virtual inline int vert_node_map( const int vert_lid) = 0;

    }; // end of 3D parent class

    class Element4D {

        public:

        // calculate a physical position in an element for a given xi,eta,mu,tau
        virtual void physical_position(
            view_c_array <real_t>  &x_point,
            const view_c_array <real_t>  &xi_point,
            const view_c_array <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta,mu,tau
        virtual void basis(
            view_c_array <real_t>  &basis,
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi at Xi_point
        virtual void partial_xi_basis(
            view_c_array <real_t>  &partial_xi, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Eta
        virtual void partial_eta_basis(
            view_c_array <real_t> &partial_eta, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Mu
        virtual void partial_mu_basis(
            view_c_array <real_t> &partial_mu, 
            const view_c_array <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Tau
        virtual void partial_tau_basis(
            view_c_array <real_t> &partial_tau, 
            const view_c_array <real_t>  &xi_point) = 0;


    }; // end of 3D parent class


    /*
     .-------------------------------. 
    | .----------------------------. |
    | |    _____       ________    | |
    | |   / ___ `.    |_   ___ `.  | |
    | |  |_/___) |      | |   `. \ | |
    | |   .'____.'      | |    | | | |
    | |  / /____       _| |___.' / | |
    | |  |_______|    |________.'  | |
    | |                            | |
    | '----------------------------' |
     '-------------------------------' 
    */
    /*
    ===========================
    2D Quad 4 Elements
    ===========================


     The finite element local point numbering for a 4 node Hexahedral is
     as follows

            Eta
             ^
             |
      3------+-----2
      |      |     |
      |      |     |
      |      |     |
      |      ------+------> Xi
      |            |
      |            |
      0------------1

    */
    class Quad4: public Element2D {
        public:
            const static int num_verts = 4;

        protected:
            static real_t ref_vert[num_verts*num_dim]; // listed as {Xi, Eta}
            static const int vert_to_node[num_verts];
        
        public:
            
            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                view_c_array <real_t>  &x_point, 
                const view_c_array <real_t>  &xi_point, 
                const view_c_array <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point); 

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid) = 0;

    }; // end of quad_4_2D class


    /*
    ===========================
    2D Quad 8 Elements
    ===========================


     The finite element local point numbering for a 8 node Hexahedral is
     as follows

             Eta
              ^
              |
      3-------6------2
      |       |      |
      |       |      |
      |       |      |
      |       |      |
      7       +------5-----> Xi   
      |              |
      |              |
      |              |
      0------4-------1

    */

    class Quad8: public Element2D {
        public:
            const static int num_verts = 8;

        protected:
            static real_t ref_vert[num_verts*num_dim]; // listed as {Xi, Eta}
            static const int vert_to_node[num_verts];

        public:

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                view_c_array <real_t>  &x_point, 
                const view_c_array <real_t>  &xi_point, 
                const view_c_array <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid) = 0;
    }; // end of quad_8_2D class

    /*
    ===========================
    2D Quad 12 Elements
    ===========================


    The finite element local point numbering for a 8 node Hexahedral is
    as follows (NEED TO DEFINE)

             Eta
              ^
              |
      3---7------6---2
      |       |      |
      |       |      |
     11       |      10
      |       |      |
      |       +------|-----> Xi   
      |              |
      8              9
      |              |
      0----4-----5---1

    */

    class Quad12: public Element2D {
        public:
            const static int num_verts = 12;

        protected:
            static real_t ref_vert[num_verts*num_dim]; // listed as {Xi, Eta} 
            static const int vert_to_node[num_verts];

        public:

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                view_c_array <real_t>  &x_point, 
                const view_c_array <real_t>  &xi_point, 
                const view_c_array <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid) = 0;

    }; // end of quad_8_2D class

    /*
    ==========================
     Arbitrary Order Elements
    ==========================

       __                   _ _   _
     / __ \                | | \ | |
    | |  | |_   _  __ _  __| |  \| |
    | |  | | | | |/ _` |/ _` | . ` |
    | |__| | |_| | (_| | (_| | |\  |
     \___\_\\__,_|\__,_|\__,_|_| \_| 

    representative linear element for visualization
     
           Eta (j)
            ^
            |
    3--------------2
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       +------|-----> Xi (i) 
    |              |
    |              |
    |              |
    0--------------1
    */


    class QuadN{
        public:

            const static int num_dim = 2;

            int num_verts;

            // calculates the basis values and derivatives in 1D
            // used in the basis_partials functiosn to build the 3D element
            void lagrange_1D(
                view_c_array <real_t> &interp,          // interpolant
                view_c_array <real_t> &Dinterp,         // derivative of function
                const real_t &x_point,                  // point of interest in element
                const view_c_array <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
                const int &orderN);                     // order of element

            void corners (
                view_c_array <real_t> &lag_nodes,   // Nodes of Lagrange elements 
                view_c_array <real_t> &lag_corner,  // corner nodes of HexN element
                const int &orderN);                 // Element order

            void physical_position (
                view_c_array <real_t> &x_point,             // location in real space
                const view_c_array <real_t> &lag_nodes,     // Nodes of Lagrange elements 
                const view_c_array <real_t> &lag_basis_2d,  // 2D basis values 
                const int &orderN);                         // order of the element

            void basis_partials (
                view_c_array <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
                view_c_array <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
                view_c_array <real_t> &val_1d,          // Interpolant Value in 1D
                view_c_array <real_t> &DVal_1d,         // Derivateive of basis in 1D
                view_c_array <real_t> &val_2d,          // for holding the interpolant in each direction
                view_c_array <real_t> &DVal_2d,         // for holding the derivatives in each direction
                view_c_array <real_t> &lag_basis_2d,    // 2D basis values 
                view_c_array <real_t> &lag_partial,     // Partial of basis 
                const view_c_array <real_t> &xi_point,  // point of interest
                const int &orderN);                     // Element order
    };



    /* 
     .-------------------------------. 
    | .----------------------------. |
    | |    ______      ________    | |
    | |   / ____ `.   |_   ___ `.  | |
    | |   `'  __) |     | |   `. \ | |
    | |   _  |__ '.     | |    | | | |
    | |  | \____) |    _| |___.' / | |
    | |   \______.'   |________.'  | |
    | |                            | |
    | '----------------------------' |
     '------------------------------' 
    */

    /*
    ==========================
      Hex 8
    ==========================

     The finite element local point numbering for a 8 node Hexahedral is
     as follows

             Mu (k)
             |     Eta (j)    
             |    /
             |   /
         7---+----6
        /|   |   /|
       / |   |  / |
      4--------5  |
      |  |    -|--+---> Xi (i)
      |  |     |  |
      |  3-----|--2
      | /      | /       
      |/       |/
      0--------1
     
    */

    class Hex8: public Element3D {
        
        public:
            
            const static int num_verts = 8;

        protected:
            static real_t ref_vert[num_verts*num_dim];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts];

        public:

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                view_c_array <real_t>  &x_point, 
                const view_c_array <real_t>  &xi_point, 
                const view_c_array <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                view_c_array <real_t> &partial_mu, 
                const view_c_array <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid);

    }; // end of hex 8 class


    /*
    ==========================
      Hex 20
    ==========================

    The finite element local point numbering for a 20 node Hexahedral is 
    as follows

               Mu (k)
               |     Eta (j)
               |    /
               |   /

          7----14----6
         /|         /|
       15 |       13 |
       / 19       /  18
      4----12----5   |
      |   |      |   |  --> Xi (i)
      |   |      |   |
      |   3---10-|---2
     16  /      17  /
      | 11       | 9         
      |/         |/
      0-----8----1

    */

    class Hex20: public Element3D {
        public:
            const static int num_verts = 20;

        protected:
            static real_t ref_vert[num_verts*num_dim];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts];

        public:

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                view_c_array <real_t>  &x_point, 
                const view_c_array <real_t>  &xi_point, 
                const view_c_array <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                view_c_array <real_t> &partial_mu, 
                const view_c_array <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid);

    }; // end of hex20 class

    /* 
    ==========================
      Hex 32
    ==========================

    The finite element local point numbering for a 32 node Hexahedral is 
    shown below


                   Mu (k)
                    ^         Eta (j)
                    |        /
                    |       /
                           /
            7----23------22----6
           /|                 /|
         15 |               14 |
         /  |               /  |
       12  31             13   30 
       /    |             /    |
      4-----20-----21----5     |
      |     |            |     |   ----> Xi (i)
      |    27            |     26  
      |     |            |     |
     28     |           29     |
      |     3----19------|18---2
      |    /             |    /
      |  11              |   10
     24  /              25  /
      | 8                | 9         
      |/                 |/
      0----16------17----1
    */

    class Hex32: public Element3D {
        public:
            const static int num_verts = 32;

        protected:
            static real_t ref_vert[num_verts*num_dim];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts];
            
        public:

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                view_c_array <real_t>  &x_point, 
                const view_c_array <real_t>  &xi_point, 
                const view_c_array <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                view_c_array <real_t> &partial_mu, 
                const view_c_array <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid);

    }; // end of hex32 class

    /*
     ==========================
      Arbitrary Order Elements
     ==========================
      _   _           _   _ 
     | | | | _____  _| \ | |
     | |_| |/ _ \ \/ /  \| |
     |  _  |  __/>  <| |\  |
     |_| |_|\___/_/\_\_| \_|
                            
    representative linear element for visualization
       
                k
                |     j    
                |    /
                |   /
            6---+----7
           /|   |   /|
          / |   |  / |
         2--------3  |
         |  |    -|--+---> i
         |  |     |  |
         |  4-----|--5
         | /      | /       
         |/       |/
         0--------1
        
       Note: left hand coordinate coordinates
    */


    class HexN{
        public:

            const static int num_dim = 3;


            int num_verts;

            // creates nodal positions with Chebyshev spacing
            void chebyshev_nodes_1D(
                view_c_array <real_t> &cheb_nodes_1D,   // Chebyshev nodes
                const int &order);                      // Interpolation order

            // calculates the basis values and derivatives in 1D
            // used in the basis_partials functiosn to build the 3D element
            void lagrange_1D(
                view_c_array <real_t> &interp,          // interpolant
                view_c_array <real_t> &Dinterp,         // derivative of function
                const real_t &x_point,                  // point of interest in element
                const view_c_array <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
                const int &orderN);                     // order of element

            void corners (
                view_c_array <real_t> &lag_nodes,   // Nodes of Lagrange elements 
                view_c_array <real_t> &lag_corner,  // corner nodes of HexN element
                const int &orderN);                 // Element order

            void physical_position (
                view_c_array <real_t> &x_point,             // location in real space
                const view_c_array <real_t> &lag_nodes,     // Nodes of Lagrange elements 
                const view_c_array <real_t> &lag_basis_3d,  // 3D basis values 
                const int &orderN);                         // order of the element

            void basis_partials (
                view_c_array <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
                view_c_array <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
                view_c_array <real_t> &val_1d,          // Interpolant Value in 1D
                view_c_array <real_t> &DVal_1d,         // Derivateive of basis in 1D
                view_c_array <real_t> &val_3d,          // for holding the interpolant in each direction
                view_c_array <real_t> &DVal_3d,         // for holding the derivatives in each direction
                view_c_array <real_t> &lag_basis_3d,    // 3D basis values 
                view_c_array <real_t> &lag_partial,     // Partial of basis 
                const view_c_array <real_t> &xi_point,  // point of interest
                const int &orderN);                     // Element order
    };



    /*
    ==========================
     4D Tesseract element
    ==========================
     
    The finite element local point numbering for a 16 node Tesseract is
    based on the 3D Hex8 Ensight element
     

                     _.15--------------------------------------14
                _.+<    |\                              . >-"/ |
          _ .+>         | \                         .>"" ./    |
      .>""              |  \                     .""    /      |
    12----------------------+------------------13    ./        |
    | )<=               |    \               / | _/""          |
    |     )\+           |     \            / /"|               |
    |         (\=       |   _. 7---------+--6  |               |
    |             \>   .|+<    |       / . "|  |               |
    |               '4--+------+------5'    |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |   _ .3------+-----2_ |               |
    |                |  | "   /       |   /'  '| \= _          |
    |                0--+---+---------1*"      |     '\        |
    |             ./'   |  "           \       |       ""\     |
    |            /      |/              \      |           '". |
    |         /        11----------------+-----+---------------10
    |      ./    .+<""                    )    |           .</
    |    /(   /(                            \  |     _.+</
    | ./  /"                                 \ |  >(
    8------------------------------------------9'

                      j
                      ^        k
                      |      /
                      |    / 
                      |  /
                      |/
                      +---------->i

       i = Xi
       j = Eta
       k = Mu
       t = Tau 


    Note: left hand coordinate coordinates
    */

    class Tess16: public Element4D {
        public:
            const static int num_verts = 16;
            const static int num_dim = 4;

        protected:
            static real_t ref_vert[num_verts*num_dim];  // listed as {Xi, Eta, Mu, Tau}


        public:

            // calculate a physical position in an element for a given xi,eta,mu
            void physical_position(
                view_c_array <real_t> &x_point,
                const view_c_array <real_t> &xi_point,
                const view_c_array <real_t> &vertices);
            
            // calculate the value for the basis at each node for a given xi,eta,mu,tau
            void basis(
                view_c_array <real_t>  &basis,
                const view_c_array <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi at Xi_point
            void partial_xi_basis(
                view_c_array <real_t> &partial_xi, 
                const view_c_array <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Eta
            void partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Mu
            void partial_mu_basis(
                view_c_array <real_t> &partial_mu, 
                const view_c_array <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Tau
            void partial_tau_basis(
                view_c_array <real_t> &partial_tau, 
                const view_c_array <real_t> &xi_point);                                          
    }; // End of Tess16 Element Class


} // end namespace elements::swage

} // end namespace elements


// Element choice
// extern elements::swage::Hex8       elem;
// extern elements::swage::Element3D  *element;

// Mesh definitions
// extern elements::mesh_t init_mesh;
// extern elements::mesh_t mesh;


#endif //ELEMENTS_H


