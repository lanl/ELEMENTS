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
#ifndef ELEMENTS_2D_H
#define ELEMENTS_2D_H

namespace elements {

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
    class Element2D {

        protected:
            const static int num_dim_ = 2;

        public:
        
        virtual int num_verts() = 0;
        virtual int num_nodes() = 0;
        virtual int num_basis() = 0;
        
        //list of local ids to basis functions needed to interpolate throughout a given element surface
        RaggedRightArray<int> surface_to_dof_lid;
        int nsurfaces;

        // calculate a physical position in an element for a given xi,eta
        virtual void physical_position(
            ViewCArray <real_t>  &x_point,
            const ViewCArray <real_t>  &xi_point,
            const ViewCArray <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta
        virtual void basis(
            ViewCArray <real_t>  &basis,
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi
        virtual void  partial_xi_basis(
            ViewCArray <real_t>  &partial_xi, 
            const ViewCArray <real_t>  &xi_point) = 0;


        // Partial derivative of shape functions with respect to Xi
        virtual void  partial_eta_basis(
            ViewCArray <real_t> &partial_eta, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Map from vertex to node
        virtual int vert_node_map( const int vert_lid) = 0;

    }; 


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
        
        protected:

            static const int num_verts_ = 4;
            static const int num_nodes_ = 4;
            static const int num_basis_ = 4;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:

            Quad4();
            ~Quad4();

            int num_verts();
            int num_nodes();
            int num_basis();
            
            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point); 

            // Map from vertex to node
            int vert_node_map( const int vert_lid);

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
        
        protected:

            static const int num_verts_ = 8;
            static const int num_nodes_ = 8;
            static const int num_basis_ = 8;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:
        
            Quad8();
            ~Quad8();

            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            int vert_node_map( const int vert_lid);
    
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
        
        protected:

            static const int num_verts_ = 12;
            static const int num_nodes_ = 12;
            static const int num_basis_ = 12;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:
        
            Quad12();
            ~Quad12();

            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            int vert_node_map( const int vert_lid);

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
    class QuadN {
        public:

            const static int num_dim = 2;

            int num_basis;
            int num_verts;

            // calculates the basis values and derivatives in 1D
            // used in the basis_partials functiosn to build the 3D element
            void lagrange_1D(
                ViewCArray <real_t> &interp,          // interpolant
                ViewCArray <real_t> &Dinterp,         // derivative of function
                const real_t &x_point,                  // point of interest in element
                const ViewCArray <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
                const int &orderN);                     // order of element

            void corners (
                ViewCArray <real_t> &lag_nodes,   // Nodes of Lagrange elements 
                ViewCArray <real_t> &lag_corner,  // corner nodes of QuadN element
                const int &orderN);                 // Element order

            void physical_position (
                ViewCArray <real_t> &x_point,             // location in real space
                const ViewCArray <real_t> &lag_nodes,     // Nodes of Lagrange elements 
                const ViewCArray <real_t> &lag_basis_2d,  // 2D basis values 
                const int &orderN);                         // order of the element

            void basis_partials (
                ViewCArray <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
                ViewCArray <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
                ViewCArray <real_t> &val_1d,          // Interpolant Value in 1D
                ViewCArray <real_t> &DVal_1d,         // Derivateive of basis in 1D
                ViewCArray <real_t> &val_2d,          // for holding the interpolant in each direction
                ViewCArray <real_t> &DVal_2d,         // for holding the derivatives in each direction
                ViewCArray <real_t> &lag_basis_2d,    // 2D basis values 
                ViewCArray <real_t> &lag_partial,     // Partial of basis 
                const ViewCArray <real_t> &xi_point,  // point of interest
                const int &orderN);                     // Element order
    };

}
#endif // ELEMENTS_2D_H
