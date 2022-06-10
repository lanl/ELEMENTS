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
#ifndef ELEMENTS_3D_H
#define ELEMENTS_3D_H 

namespace elements {

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
    class Element3D {
                
        protected:
            const static int num_dim_ = 3;

        public:
        
        virtual int num_verts() = 0;
        virtual int num_nodes() = 0;
        virtual int num_basis() = 0;

        //list of local ids to basis functions needed to interpolate throughout a given element surface
        RaggedRightArray<int> surface_to_dof_lid;
        int nsurfaces;

        // calculate a physical position in an element for a given xi,eta,mu
        virtual void physical_position(
            ViewCArray <real_t>  &x_point,
            const ViewCArray <real_t>  &xi_point,
            const ViewCArray <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta, mu
        virtual void basis(
            ViewCArray <real_t>  &basis,
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi at Xi_point
        virtual void partial_xi_basis(
            ViewCArray <real_t>  &partial_xi, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Eta
        virtual void partial_eta_basis(
            ViewCArray <real_t> &partial_eta, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Mu
        virtual void partial_mu_basis(
            ViewCArray <real_t> &partial_mu, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Map from vertex to node
        virtual inline int vert_node_map(const int vert_lid) = 0;

        // Reference vertices location
        virtual real_t& ref_locs(const int vert_lid, const int dim) = 0;

    }; // end of 3D parent class


    /*
    ==========================
      Hex 8
    ==========================

    The finite element local vertex numbering for a 8 node Hexahedral is
    as follows

             Mu (k)
             |     Eta (j)    
             |    /
             |   /
         6---+----7
        /|   |   /|
       / |   |  / |
      4--------5  |
      |  |    -|--+---> Xi (i)
      |  |     |  |
      |  2-----|--3
      | /      | /       
      |/       |/
      0----*----1

    */
    class Hex8: public Element3D {
        
        protected:

            static const int num_verts_ = 8;
            static const int num_nodes_ = 8;
            static const int num_basis_ = 8;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:
        
            Hex8();
            ~Hex8();

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
            void partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map(const int vert_lid);

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
        protected:

            static const int num_verts_ = 20;
            static const int num_nodes_ = 20;
            static const int num_basis_ = 20;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:
        
            Hex20();
            ~Hex20();

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

            // with repsect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map(const int vert_lid);

            real_t& ref_locs(const int vert_lid, const int dim);

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
        
        protected:

            static const int num_verts_ = 32;
            static const int num_nodes_ = 32;
            static const int num_basis_ = 32;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:
        
            Hex32();
            ~Hex32();

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

            // with repsect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            int vert_node_map(const int vert_lid);

            real_t& ref_locs(const int vert_lid, const int dim);

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
    4--------5  |
    |  |    -|--+---> i
    |  |     |  |
    |  2-----|--3
    | /      | /       
    |/       |/
    0--------1
    */
    class HexN {
        
        protected:
            
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


        public:

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
            void basis(
                MatarRealCArray &basis,
                MatarRealCArray &point);

            void build_nodal_gradient(
                MatarRealCArray &gradient);

            // calculate the partial of the basis w.r.t xi at a given point
            void partial_xi_basis(
                MatarRealCArray &partial_xi, 
                MatarRealCArray &point);

            // calculate the partial of the basis w.r.t eta at a given point
            void partial_eta_basis(
                MatarRealCArray &partial_eta, 
                MatarRealCArray &point);

            // calculate the partial of the basis w.r.t mu at a given point
            void partial_mu_basis(
                MatarRealCArray &partial_mu, 
                MatarRealCArray &point);

            void lagrange_basis_1D(
                MatarRealCArray &interp,    // interpolant
                const real_t &x_point);     // point of interest in element
            
            void lagrange_derivative_1D(
                MatarRealCArray &partials,  //derivative
                const real_t &x_point);     // point of interest in element

            void create_lobatto_nodes(int element_order);

    };

}

#endif // ELEMENTS_3D_H 
