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
#ifndef ELEMENTS_4D_H
#define ELEMENTS_4D_H 

namespace elements {

    class Element4D {

        public:

        //list of local ids to basis functions needed to interpolate throughout a given element surface
        CArray<int> surface_to_dof_lid;

        // calculate a physical position in an element for a given xi,eta,mu,tau
        virtual void physical_position(
            ViewCArray <real_t>  &x_point,
            const ViewCArray <real_t>  &xi_point,
            const ViewCArray <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta,mu,tau
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

        // Partial derivative of shape functions with respect to Tau
        virtual void partial_tau_basis(
            ViewCArray <real_t> &partial_tau, 
            const ViewCArray <real_t>  &xi_point) = 0;


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
            const static int num_basis = 16;

        protected:
            static real_t ref_vert[num_verts*num_dim];  // listed as {Xi, Eta, Mu, Tau}


        public:

            // calculate a physical position in an element for a given xi,eta,mu
            void physical_position(
                ViewCArray <real_t> &x_point,
                const ViewCArray <real_t> &xi_point,
                const ViewCArray <real_t> &vertices);
            
            // calculate the value for the basis at each node for a given xi,eta,mu,tau
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi at Xi_point
            void partial_xi_basis(
                ViewCArray <real_t> &partial_xi, 
                const ViewCArray <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Eta
            void partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Tau
            void partial_tau_basis(
                ViewCArray <real_t> &partial_tau, 
                const ViewCArray <real_t> &xi_point);                                          
    }; // End of Tess16 Element Class

}

#endif // ELEMENTS_4D_H 
