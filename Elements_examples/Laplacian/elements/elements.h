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

    void set_unit_normals(view_c_array <real_t> &unit_normals);





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

    int indx_; //useful for returning from internal function


// ---- ELEMENT POLYNOMIAL ORDER ---- //
    int   elem_order_;


// ---- CCH RECONSTRUCTION POLYNOMIAL ORDER ---- //
    int   recon_order_;

// ---- INDEX SPACES AND MAPS ---- //

    // ---- ELEMENT ---- //
    int   num_elem_;
    int   num_g_pts_in_elem_;        //number of quadrature points in an element
    int   num_cells_in_elem_;        //number of finite volume cells defined in an element
    int   num_nodes_in_elem_;        //nodal degrees of freedom per element
    int   num_mat_pts_in_elem_;

    int * cells_in_elem_ = NULL;

    int * num_elems_in_elem_ = NULL;
    int * elems_in_elem_list_start_ = NULL;
    int * elems_in_elem_list_ = NULL;
    int * nodes_in_elem_list_ = NULL;


    // ---- CELLS ---- //
    int   num_cells_;
    
    int * nodes_in_cell_list_ = NULL;      // size of num_cells*8
    int * num_cells_in_cell_ = NULL;       // size of num_cells; stores the number of adjacent elements for each element
    int * cells_in_cell_list_start_ = NULL; // size of num_nodes+1; used for ragged right 1D array index storage
    int * cells_in_cell_list_ = NULL;       // size depends on mesh connectivity; stores adjacent element indices for each element
    int * elems_in_cell_list_ = NULL;       // size depends on mesh connectivity


    // ---- VERTICES ---- //



    // ---- NODES ---- //
    int   num_nodes_;

    int * num_cells_in_node_ = NULL;        // size of num_nodes
    int * cells_in_node_list_start_ = NULL; // size of num_nodes+1; used for ragged right 1D array index storage
    int * cells_in_node_list_ = NULL;       // size depends on mesh connectivity
    
    int * num_elems_in_node_ = NULL;        // number of elements a node belongs to in the mesh
    int * elems_in_node_list_start_ = NULL; // used for ragged right 1D array storage
    int * elems_in_node_list_ = NULL;       // list of elements a node is shared with

    // ---- GAUSS POINTS ---- //
    int   num_g_pts_;

    int * node_in_gauss_list_ = NULL;


    // ---- CORNERS ---- //
    int   num_corners_;                          //the number of countable nodes if the elements in the mesh did not share nodes.

    int * num_corners_in_node_ = NULL;           //number of element/cell geometry corners a node is representing

    int * corners_in_cell_list_ = NULL;          //stores the global corner indices for each cell/element

    int * corners_in_node_list_start_ = NULL;    //used for ragged right 1D array storage of the corner indices each node represents
    int * corners_in_node_list_ = NULL;          //stores the global corner indices for each corner a node is representing

    // int * corner_bdy_count_ = NULL;

    /*not sure if we should call them faces once we move into the realm of non-linear element surfaces.
      The formal geometry definition of face is restricted to planar.*/
    // ---- FACES ---- //
    int   num_faces_;                    //total number of smooth element/cell surface segments

    int * face_nodes_list_ = NULL;       // size of num_faces*4; stores nodes that belong to the surface segment
    int * cells_in_face_list_ = NULL;    // size of num_faces*2; stores 1 or 2 cell ids that the surface belongs to


    // ---- BOUNDARY ---- //
    int num_bdy_faces_;          //total number of element/cell surface segments that coincide with the model boundary
    int num_bdy_sets_;           //number of independent boundary conditions applied to the model
    
    int * bdy_faces_;            /* size depends on mesh; stores the global element/cell surface indices for all
                                    such surfaces that coincide with the model boundary */
    int * bdy_set_list_;         /* stores the set of global element/cell surface indices that pertain to each
                                    boundary constraint. */
    int * start_index_bdy_set_;  /* stores the 1D array indices at which each boundary
                                    constraint's surface index storage commences */
    int * num_bdy_faces_set_;    //stores the number of surface segments that pertain to a boundary constraint


// ---- MESH GEOMETRIC STATE ---- //

    // ---- ELEMENT ---- //
    real_t * elem_vol_ = NULL;  // size of num_elem


    // ---- CELLS ---- //
    real_t * cell_vol_ = NULL;      // size of num_cells
    real_t * cell_coords_ = NULL;  // size of num_cells


    // ---- NODES ---- //
    real_t * node_coords_ = NULL;      // size of rk_storage_*num_nodes*num_dim_
    real_t * node_jacobian_ = NULL;    // Jacobian matrix at each node
    

    // ---- QUADRATURE POINTS ---- //
    real_t * jacobians_ = NULL;            // size of num_g_pts_*num_dim_*num_dim_
    real_t * jacobian_determinant_ = NULL; // size of num_g_pts_


public:

    //**********************************//
    // Mesh class function definitions  //
    //**********************************//
    
    /*The variable dim refers to the problem dimension. The variable num_rk refers to
      the number of vectors (of size dim) stored by each node. For example, position and velocity
      implies num_rk = 2. e_order and r_order are the element polynomial order and CCH reconstruction
      polynomial orders respectively.  */
      
    void init_element (int e_order, int r_order, int dim, int num_elem, int num_rk);
    void init_cells (int ncells, int num_rk);
    void init_nodes (int num_nodes, int num_rk);
    void init_gauss_pts ();
    void init_bdy_sets (int num_sets);

    // ==== MESH CONSTANTS ==== // 

    // returns the number of rk_storage bins
    int num_rk () const;

    // returns the number of dimensions in the mesh
    int num_dim () const;

    // returns the polynomial order of the element
    int elem_order () const;

    // returns the polynomial order of the CCH reconstruction
    int& recon_order ();



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

    // return array of elements connected to element (corners+faces)
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


    // ---- FACES ---- //

    // returns the number of elements
    int num_faces () const;

    // returns the global node id given a cell_id, local_face_indx(0:5), local_facenode_indx(0:3)
    int node_in_face_in_cell(int cell_id, int this_face, int facenode_lid) const;

    // returns the global id for a cell that is connected to the face
    int cells_in_face(int face_gid, int this_cell) const;
          
    // returns the nodes in the face
    int node_in_face(int face_gid, int facenode_lid) const;


    // ---- Boundary ---- //
    
    // returns the number of independent boundary constraints applied to the system
    int num_bdy_sets() const;
    
    //returns the total number of element/cell surfaces that coincide with the boundary of the model
    int num_bdy_faces() const;

    //returns the global id of a cell surface that coincides with the boundary of the model
    int bdy_faces(int this_bdy_face) const;

    // returns the number per bdy-faces in a particular set
    int num_bdy_faces_in_set (int bdy_set);

    // returns a subset of the boundary faces
    int bdy_faces_in_set (int bdy_set, int this_face);



    // ==== MESH STATE FUNCTIONS ==== // 


    // ---- ELEMENTS ---- //
    real_t& elem_vol(int rk_bin, int elem_gid) const;



    // ---- CELLS ---- //

    // return the cell volume
    real_t& cell_vol(int cell_gid) const;

    // return the cell coordinate position
    real_t& cell_coords(int rk_bin, int cell_gid, int this_dim);



    // ---- VERTICES ---- //



    // ---- NODES ---- //
    // return the node coordinates
    real_t& node_coords(int rk_bin, int node_gid, int this_dim) const;



    // ---- QUADRATURE POINTS ---- //

    // return jacobian at quadrature point
    real_t& jacobian(int elem_gid, int gauss_lid, int i, int j) const;

    // return determinant of jacobian at quadrature point
    real_t& det_j(int elem_gid, int gauss_lid) const;




    // ---- CORNERS ---- //


    // ---- FACES ---- //
    // geometric average face coordinate
    real_t face_coords(int rk_bin, int face_id, int this_dim) const;


    // ==== MESH CONNECTIVITY FUNCTIONS ==== // 
        
    // initialize array for mesh connectivity: all cells around a node
    void build_connectivity();            //invoke all the following connectivity routines

    void build_node_cell_connectivity();  //build connectivity of nodes to elements/cells

    void build_corner_connectivity();     /*build connectivity of globally indexed geometric corners
                                            (belonging to elements/cells) to elements/cells*/
    void build_cell_cell_connectivity();  //build connectivity of elements/cells to adjacent elements/cells

    void build_face_connectivity();       /*build connectivity of globally indexed element/cell surfaces to
                                            the global element/cell indices they belong to*/
    void build_element_connectivity();    //build connectivity of elements to nodes

    // identify the boundary faces
    void build_bdy_faces ();


    // ---- bdy sets ----

    // returns a subset of the boundary faces
    int set_bdy_faces (int bdy_set, int face_lid);

    // set planes for tagging sub sets of boundary faces
    // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    // val = plane value, radius, radius
    void tag_bdys(int this_bc_tag, real_t val, int bdy_set);

    // compress the bdy_set_list to reduce the memory
    void compress_bdy_set();

    // routine for checking to see if a vertix is on a boundary
    // bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
    // val = plane value, radius, radius
    int check_bdy(int face_gid, int this_bc_tag, real_t val);


    // deconstructor
    ~mesh_t ( );

}; // end of mesh_t declaration



void refine_mesh(
    mesh_t& init_mesh, 
    mesh_t& mesh, 
    const int p_order, 
    const int rk_num_bins,
    const int dim);


//what does swage stand for or mean?
namespace swage{

    // Used by Gauss2/3D to set quadrature points
    void line_gauss_info(
        real_t &x, 
        real_t &w, 
        int &m,  
        int &p);

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
        view_c_array <real_t> &tot_g_weight,  // 3D product of gauss weights
        int &quad_order);                     // quadrature order (n)

    // setting gauss quadrature points for 4D elements
    void gauss_4d(
        view_c_array <real_t> &these_g_pts,     // gauss points
        view_c_array <real_t> &these_weights,   // gauss weights
        int &quad_order,                        // quadrature order (n)
        const int &dim);

    // setting Gauss-Lobatto quadrature points for 2D elements
    void lobatto_2d(
        view_c_array <real_t> &these_L_pts,     // gauss points
        view_c_array <real_t> &these_weights,   // gauss weights
        int &quad_order);                       // quadrature order (n)

    // setting Gauss-Lobatto quadrature points for 3D elements
    void lobatto_3d(
        view_c_array <real_t> &these_L_pts,     // gauss points
        view_c_array <real_t> &these_weights,   // gauss weights
        int &quad_order); 

    // setting gauss quadrature points for 4D elements
    void lobatto_4d(
        view_c_array <real_t> &these_L_pts,     // gauss points
        view_c_array <real_t> &these_weights,   // gauss weights
        int &quad_order,                        // quadrature order (n)
        const int &dim);

    //defining the jacobian for 2D elements
    void jacobian_2d(
        view_c_array <real_t> &J_matrix, 
        real_t &det_J,
        const view_c_array <real_t> &vertices, 
        const view_c_array <real_t> &this_partial,
        const int &num_nodes);

    //defining the jacobian for 3D elements
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

// Reference Element Informations


class ref_element{
    private:
        
        int num_dim_;
        
        int num_ref_nodes_1D_;
        int num_ref_cells_1D_;
        int num_ref_corners_1D_;
        
        // cells
        int num_ref_cells_in_elem_;
        
        // nodes
        int num_ref_nodes_in_elem_;
        int num_ref_nodes_in_cell_;
        int *ref_nodes_in_cell_ = NULL;
        
        real_t *ref_node_positions_ = NULL;
        real_t *ref_node_g_weights_ = NULL;
        
        // corners
        int num_ref_corners_in_cell_;
        int num_ref_corners_in_elem_;
        int* ref_corners_in_cell_ = NULL;
        
        real_t *ref_corner_surf_normals_ = NULL;
        real_t *ref_corner_g_weights_ = NULL;
        real_t *ref_corner_surf_g_weights_ = NULL;
    
    
    public:
        // Function Declarations

        // Initialize reference element information
        void init(int poly_order, int num_dim);

        int num_ref_cells_in_elem() const;
        int num_ref_corners_in_cell() const;
        
        int node_rid(int i, int j, int k) const;
        int cell_rid(int i, int j, int k) const;
        int corner_rid(int i, int j, int k) const;
        
        int ref_corners_in_cell(int cell_rid, int corner_rlid) const;
        int ref_nodes_in_cell(int cell_rid, int node_rlid) const;

        real_t ref_node_positions(int node_rid, int dim) const;

        real_t ref_corner_surface_normals(int corner_rid, int surf_rlid, int dim) const;
        
        real_t ref_corner_g_surface_weights(int corner_rid, int surf_rlid) const;
        
        real_t ref_node_g_weights(int node_rid) const;

        real_t ref_corner_g_weights(int corner_rid) const;

        real_t ref_nodal_gradient(int node_rid, int basis_id, int dim) const;


        // Nodal jacobian and determinant should be in mesh class, I think....
        //real_t ref_nodal_jacobian(int node_rid, int basis_id, int dim) const;
        


        // Deconstructor
        ~ref_element();

    };


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
            // const static int num_verts;

            // const static int ref_verts[num_verts*num_dim];

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

        // Reference vertices location
        virtual real_t& ref_locs(const int vert_lid, const int dim) = 0;

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


    }; // end of 4D parent class


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


    The finite element local point numbering for a 4 node quadralateral is
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

    Representative linear element for visualization
     
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
            void partial_xi_basis(
                view_c_array <real_t>  &partial_xi, 
                const view_c_array <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void partial_eta_basis(
                view_c_array <real_t> &partial_eta, 
                const view_c_array <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                view_c_array <real_t> &partial_mu, 
                const view_c_array <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid);

            real_t& ref_locs(const int vert_lid, const int dim);

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
    16 /      17  /
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
    28    |           29     |
    |     3----19------|18---2
    |    /             |    /
    |  11              |   10
    24 /              25  /
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
                            
    Representative linear element for visualization
       
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