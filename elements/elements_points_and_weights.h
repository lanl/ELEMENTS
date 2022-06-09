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
#ifndef ELEMENTS_POINTS_AND_WEIGHTS_H
#define ELEMENTS_POINTS_AND_WEIGHTS_H 

namespace elements {

    // Used by Lobatto 1D/2D to set Lobatto quadrature points
    void lobatto_nodes_1D(MatarRealCArray &lob_nodes_1D, const int &num);
    void lobatto_weights_1D(MatarRealCArray &lob_weights_1D, const int &num);

    
    void legendre_nodes_1D(CArray <real_t> &leg_nodes_1D, const int &num);
    void legendre_weights_1D(CArray <real_t> &leg_weights_1D, const int &num);


    // creates nodal positions with Chebyshev spacing
    void chebyshev_nodes_1D(
        ViewCArray <real_t> &cheb_nodes_1D,  // Chebyshev nodes
        const int &order);                   // Interpolation order


    void length_weights(
        MatarRealCArray &len_weights_1D,  // Lobatto weights
        MatarRealCArray &lab_weights_1D,  // Lobatto weights
        MatarRealCArray &lab_nodes_1D,
        const int &order);

    void sub_weights(
        MatarRealCArray &sub_weights_1D,  // Lobatto weights
        MatarRealCArray &lab_weights_1D,  // Lobatto weights
        MatarRealCArray &lab_nodes_1D,
        const int &order);


    void set_nodes_wgts(
        CArray <real_t> &lab_nodes_1D,
        CArray <real_t> &lab_weights_1D,
        CArray <real_t> &len_weights_1D,
        CArray <real_t> &sub_weights_1D, 
        int p_order);

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
        ViewCArray <real_t> &these_g_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        ViewCArray <real_t> &tot_g_weight,    // 2D product of gauss weights
        int &quad_order);                     // quadrature order (n)

    // setting gauss quadrature points for 2D elements
    void gauss_3d(
        ViewCArray <real_t> &these_g_pts,   // gauss points
        ViewCArray <real_t> &these_weights, // gauss weights
        ViewCArray <real_t> &tot_g_weight,  // 3D product of gauss weights
        int &quad_order);                     // quadrature order (n)

    // setting gauss quadrature points for 4D elements
    void gauss_4d(
        ViewCArray <real_t> &these_g_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order,                        // quadrature order (n)
        const int &dim);

    // setting Gauss-Lobatto quadrature points for 2D elements
    void lobatto_2d(
        ViewCArray <real_t> &these_L_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order);                       // quadrature order (n)

    // setting Gauss-Lobatto quadrature points for 3D elements
    void lobatto_3d(
        ViewCArray <real_t> &these_L_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order); 

    // setting gauss quadrature points for 4D elements
    void lobatto_4d(
        ViewCArray <real_t> &these_L_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order,                        // quadrature order (n)
        const int &dim);
    
}

#endif // ELEMENTS_POINTS_AND_WEIGHTS_H
