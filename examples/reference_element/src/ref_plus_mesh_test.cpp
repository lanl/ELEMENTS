/**********************************************************************************************
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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// This pulls in kokkos, matar, mesh, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"
#include "cramers_rule.hpp" // det and solvers

using namespace mtr;
using namespace swage;    // unstructured mesh and point cloud
using namespace elements; // reference element space




// Convert from ensight FE elem to IJK elem
const size_t convert_ensight_to_ijk[8] =
        {0,
         1,
         3,
         2,
         4,
         5,
         7,
         6};

///////////////////////////////////////////////////////////////////////////////
/// Test 1: 3D Affine Transformation
/// Maps reference cube [-1,1]³ to a skewed, rotated, stretched parallelepiped
///////////////////////////////////////////////////////////////////////////////
void test_affine_transformation() {
    printf("\n=== TEST: Affine Transformation ===\n");
    
    // Reference element: unit cube [-1,1]^3
    // Physical element: defined by affine map x = A*xi + b
    // where A is the transformation matrix

    Quadrature_t Quad;
    ReferenceElement_t RefElem; 

    // create quadrature, it is a single point at this time
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               1,
                               3);

    // create ref element
    RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLobatto,
                                  Quad,
                                  1);   
    
    // Transformation matrix (rotation + stretch + shear)
    real_t A[3][3] = {
        { 0.5,  0.3,  0.2},   // dx/dxi, dx/deta, dx/dzeta
        { 0.2,  0.6,  0.1},   // dy/dxi, dy/deta, dy/dzeta
        { 0.1,  0.2,  0.4}    // dz/dxi, dz/deta, dz/dzeta
    };
    
    real_t b[3] = {1.0, 2.0, 3.0};  // Translation
    
    // Generate 8 corner nodes of transformed hexahedron
    real_t node_coords[8][3];
    real_t ref_corners[8][3] = {
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},  // bottom
        {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1}   // top
    };

    DCArrayKokkos <double> node_coords_dual(8,3);
    DCArrayKokkos <double> node_in_elem_dual(8);
    
    for(size_t n = 0; n < 8; n++) {
        for(size_t i = 0; i < 3; i++) {
            node_coords[n][i] = b[i];
            for(size_t j = 0; j < 3; j++) {
                node_coords[n][i] += A[i][j] * ref_corners[n][j];
            }

            // get i,j,k node lid for this node
            size_t node_lid = convert_ensight_to_ijk[n];
            node_coords_dual.host(node_lid,i) = node_coords[n][i];
            node_in_elem_dual.host(node_lid) = node_lid;
        }
    }
    node_coords_dual.update_device();
    node_in_elem_dual.update_device();

    
    printf("Node coordinates (Ensight FE element ordering):\n");
    for(size_t n = 0; n < 8; n++) {
        printf("  Node %zu: [%8.4f, %8.4f, %8.4f]\n", 
               n, node_coords[n][0], node_coords[n][1], node_coords[n][2]);
    }
    
    // For an AFFINE transformation, Jacobian is CONSTANT everywhere
    // and equals the transformation matrix A
    printf("\nExpected Jacobian (constant, equals A):\n");
    for(int i = 0; i < 3; i++) {
        printf("  ");
        for(int j = 0; j < 3; j++) {
            printf("%10.6f ", A[i][j]);
        }
        printf("\n");
    }
    
    real_t expected_det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
                        - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
                        + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    
    printf("Expected determinant: %12.8f\n", expected_det);
    
    // Now compute Jacobian at element center (xi=0, eta=0, zeta=0)
    // and verify it matches A

    // Call jacobian function and then compare results
    DCArrayKokkos <double> jac(3, 3);
    RUN({
        // exact the basis and grad basis for this quadrature point
        //ViewCArrayKokkos a_basis(&RefElem.qpt_basis(0,0),RefElem.num_dofs_in_elem);
        ViewCArrayKokkos a_grad_basis(&RefElem.qpt_grad_basis(0,0,0),RefElem.num_dofs_in_elem,3);

        jacobian(jac, 
                 node_coords_dual, 
                 node_in_elem_dual,
                 a_grad_basis);

        printf("\n");
        printf("calculated jacobian matrix = \n");
        for(size_t i = 0; i < 3; i++){  // looping over dimension
            for(size_t j = 0; j < 3; j++){ // looping over dimension
                printf("%10.6f, ", jac(i, j));
            }
            printf("\n");
        } // end for 

        double det = det_3x3(jac);
        printf("Calculated determinant: %12.8f\n", det);
    }); // end RUN on device
    jac.update_host();

    for(size_t i = 0; i < 3; i++)  // looping over dimension
    for(size_t j = 0; j < 3; j++){ // looping over dimension
        if(fabs(jac.host(i, j)-A[i][j])>1e-12){
            Kokkos::abort("ERROR: Jacobian calculation in affine test failed\n");
        }
    } // end for 

    printf("\nTest criteria:\n");
    printf("  1. J should equal A within tolerance (~1e-12)\n");
    printf("  2. det(J) should equal %12.8f\n", expected_det);
    printf("  3. J should be same at ALL quadrature points (affine property)\n");

    printf("\nAffine test: Passes\n\n");
} // end function


///////////////////////////////////////////////////////////////////////////////
/// Test 2: Heavily Skewed Hexahedron
/// Creates a non-trivial, fully-populated Jacobian
///////////////////////////////////////////////////////////////////////////////
void test_skewed_hexahedron() {
    printf("\n=== TEST: Heavily Skewed Hexahedron ===\n");
    
    // Define a highly skewed but valid element
    real_t node_coords[8][3] = {
        // Bottom face (z=0 plane, but skewed)
        {0.0, 0.0, 0.0},     // node 0
        {1.0, 0.2, 0.1},     // node 1
        {1.1, 1.0, 0.15},    // node 2
        {0.1, 0.9, 0.05},    // node 3
        // Top face (z~1 plane, but twisted and stretched)
        {0.1, 0.1, 0.9},     // node 4
        {1.2, 0.1, 1.0},     // node 5
        {1.3, 1.1, 1.1},     // node 6
        {0.0, 1.0, 0.95}     // node 7
    };
    
    printf("Node coordinates (Ensight FE element ordering):\n");
    for(int n = 0; n < 8; n++) {
        printf("  Node %d: [%8.4f, %8.4f, %8.4f]\n", 
               n, node_coords[n][0], node_coords[n][1], node_coords[n][2]);
    }
    
    // Compute Jacobian at center
    printf("\nJacobian at element center (xi=eta=zeta=0):\n");
    printf("Expected properties:\n");
    printf("  - All 9 entries should be non-zero\n");
    printf("  - Determinant should be positive\n");
    printf("  - Off-diagonal terms significant (non-trivial mapping)\n");
    
    // Analytical approximation (for verification):
    // At center, with trilinear shape functions:
    // J ~ (1/8) * sum of differences in each direction
    
    real_t J_approx[3][3];
    for(int i = 0; i < 3; i++) {
        J_approx[i][0] = 0.125 * ((node_coords[1][i] + node_coords[2][i] + 
                                   node_coords[5][i] + node_coords[6][i]) -
                                  (node_coords[0][i] + node_coords[3][i] + 
                                   node_coords[4][i] + node_coords[7][i]));
        
        J_approx[i][1] = 0.125 * ((node_coords[2][i] + node_coords[3][i] + 
                                   node_coords[6][i] + node_coords[7][i]) -
                                  (node_coords[0][i] + node_coords[1][i] + 
                                   node_coords[4][i] + node_coords[5][i]));
        
        J_approx[i][2] = 0.125 * ((node_coords[4][i] + node_coords[5][i] + 
                                   node_coords[6][i] + node_coords[7][i]) -
                                  (node_coords[0][i] + node_coords[1][i] + 
                                   node_coords[2][i] + node_coords[3][i]));
    }
    
    printf("\nApproximate Jacobian (analytical):\n");
    for(int i = 0; i < 3; i++) {
        printf("  ");
        for(int j = 0; j < 3; j++) {
            printf("%10.6f ", J_approx[i][j]);
        }
        printf("\n");
    }
    
    real_t det_approx = J_approx[0][0]*(J_approx[1][1]*J_approx[2][2] - 
                                        J_approx[1][2]*J_approx[2][1])
                      - J_approx[0][1]*(J_approx[1][0]*J_approx[2][2] - 
                                        J_approx[1][2]*J_approx[2][0])
                      + J_approx[0][2]*(J_approx[1][0]*J_approx[2][1] - 
                                        J_approx[1][1]*J_approx[2][0]);
    
    printf("Approximate determinant: %12.8f\n", det_approx);

    // =======================
    // using ELEMENTS
    // =======================
    DCArrayKokkos <double> node_coords_dual(8,3);
    DCArrayKokkos <double> node_in_elem_dual(8);
    
    for(size_t n = 0; n < 8; n++) {
        for(size_t i = 0; i < 3; i++) {
            // get i,j,k node lid for this node
            size_t node_lid = convert_ensight_to_ijk[n];
            node_coords_dual.host(node_lid,i) = node_coords[n][i];
            node_in_elem_dual.host(node_lid) = node_lid;
        }
    }
    node_coords_dual.update_device();
    node_in_elem_dual.update_device();

    Quadrature_t Quad;
    ReferenceElement_t RefElem; 

    // create quadrature, it is a single point at this time
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               1,
                               3);

    // create ref element
    RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLobatto,
                                  Quad,
                                  1);   

    // Call jacobian function and then compare results
    DCArrayKokkos <double> jac(3, 3);
    RUN({
        // exact the basis and grad basis for this quadrature point
        //ViewCArrayKokkos a_basis(&RefElem.qpt_basis(0,0),RefElem.num_dofs_in_elem);
        ViewCArrayKokkos a_grad_basis(&RefElem.qpt_grad_basis(0,0,0),RefElem.num_dofs_in_elem,3);

        jacobian(jac, 
                 node_coords_dual, 
                 node_in_elem_dual,
                 a_grad_basis);

        printf("\n");
        printf("calculated jacobian matrix = \n");
        for(size_t i = 0; i < 3; i++){  // looping over dimension
            for(size_t j = 0; j < 3; j++){ // looping over dimension
                printf("%10.6f, ", jac(i, j));
            }
            printf("\n");
        } // end for 

        double det = det_3x3(jac);
        printf("Calculated determinant: %12.8f\n", det);
    }); // end RUN on device
    jac.update_host();

    for(size_t i = 0; i < 3; i++)  // looping over dimension
    for(size_t j = 0; j < 3; j++){ // looping over dimension
        if(fabs(jac.host(i, j)-J_approx[i][j])>1.e-10){
            Kokkos::abort("ERROR: Jacobian calculation in skewed hex test failed\n");
        }
    } // end for 

    printf("\nSkewed hex test: Passes\n\n");

} // end function



///////////////////////////////////////////////////////////////////////////////
/// Test 3: Rotated and Scaled Cube
/// Apply known rotation matrix + scaling
///////////////////////////////////////////////////////////////////////////////
void test_rotated_scaled_cube() {
    printf("\n=== TEST: Rotated and Scaled Cube ===\n");
    
    // Rotation angle
    const real_t theta = M_PI / 6.0;  // 30 degrees
    const real_t phi   = M_PI / 4.0;  // 45 degrees
    
    // Scale factors
    const real_t sx = 0.5, sy = 0.8, sz = 1.2;
    
    // Combined transformation matrix: R_z(phi) * R_y(theta) * Scale
    real_t cos_t = cos(theta), sin_t = sin(theta);
    real_t cos_p = cos(phi),   sin_p = sin(phi);
    
    real_t T[3][3] = {
        {sx * cos_t * cos_p,  sx * (-sin_p),  sx * sin_t * cos_p},
        {sy * cos_t * sin_p,  sy * cos_p,     sy * sin_t * sin_p},
        {sz * (-sin_t),       sz * 0.0,       sz * cos_t}
    };
    
    printf("Transformation matrix (Rotation + Scale):\n");
    for(int i = 0; i < 3; i++) {
        printf("  ");
        for(int j = 0; j < 3; j++) {
            printf("%10.6f ", T[i][j]);
        }
        printf("\n");
    }
    
    // Reference cube vertices
    real_t ref_verts[8][3] = {
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
        {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1}
    };
    
    // Transform vertices
    real_t node_coords[8][3];
    for(int n = 0; n < 8; n++) {
        for(int i = 0; i < 3; i++) {
            node_coords[n][i] = 0.0;
            for(int j = 0; j < 3; j++) {
                node_coords[n][i] += T[i][j] * ref_verts[n][j];
            }
        }
    }
    
    printf("\nTransformed node coordinates:\n");
    for(int n = 0; n < 8; n++) {
        printf("  Node %d: [%8.4f, %8.4f, %8.4f]\n", 
               n, node_coords[n][0], node_coords[n][1], node_coords[n][2]);
    }
    
    // Expected Jacobian = T (constant for affine transformation)
    printf("\nExpected Jacobian (equals T):\n");
    for(int i = 0; i < 3; i++) {
        printf("  ");
        for(int j = 0; j < 3; j++) {
            printf("%10.6f ", T[i][j]);
        }
        printf("\n");
    }
    
    // Expected determinant
    real_t expected_det = sx * sy * sz * 
                         (cos_t * cos_t * cos_p * cos_p + 
                          cos_t * cos_t * sin_p * sin_p + 
                          sin_t * sin_t);
    
    // Simplified: det(R) = 1, det(Scale) = sx*sy*sz
    expected_det = sx * sy * sz;
    
    printf("\nExpected determinant: %12.8f\n", expected_det);
    printf("(Should equal sx*sy*sz = %f * %f * %f = 0.48)\n", sx, sy, sz);


    // =======================
    // using ELEMENTS
    // =======================
    DCArrayKokkos <double> node_coords_dual(8,3);
    DCArrayKokkos <double> node_in_elem_dual(8);
    
    for(size_t n = 0; n < 8; n++) {
        for(size_t i = 0; i < 3; i++) {
            // get i,j,k node lid for this node
            size_t node_lid = convert_ensight_to_ijk[n];
            node_coords_dual.host(node_lid,i) = node_coords[n][i];
            node_in_elem_dual.host(node_lid) = node_lid;
        }
    }
    node_coords_dual.update_device();
    node_in_elem_dual.update_device();

    Quadrature_t Quad;
    ReferenceElement_t RefElem; 

    // create quadrature, it is a single point at this time
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               1,
                               3);

    // create ref element
    RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLobatto,
                                  Quad,
                                  1);   

    // Call jacobian function and then compare results
    DCArrayKokkos <double> jac(3, 3);
    RUN({
        // exact the basis and grad basis for this quadrature point
        //ViewCArrayKokkos a_basis(&RefElem.qpt_basis(0,0),RefElem.num_dofs_in_elem);
        ViewCArrayKokkos a_grad_basis(&RefElem.qpt_grad_basis(0,0,0),RefElem.num_dofs_in_elem,3);

        jacobian(jac, 
                 node_coords_dual, 
                 node_in_elem_dual,
                 a_grad_basis);

        printf("\n");
        printf("calculated jacobian matrix = \n");
        for(size_t i = 0; i < 3; i++){  // looping over dimension
            for(size_t j = 0; j < 3; j++){ // looping over dimension
                printf("%10.6f, ", jac(i, j));
            }
            printf("\n");
        } // end for 

        double det = det_3x3(jac);
        printf("Calculated determinant: %12.8f\n", det);
    }); // end RUN on device
    jac.update_host();

    for(size_t i = 0; i < 3; i++)  // looping over dimension
    for(size_t j = 0; j < 3; j++){ // looping over dimension
        if(fabs(jac.host(i, j)-T[i][j])>1.e-10){
            Kokkos::abort("ERROR: Jacobian calculation in rotated and scaled test failed\n");
        }
    } // end for 

    printf("\nRotated and scaled hex test: Passes\n\n");

} // end function


void test_manufactured_solution() {
    printf("\n=== TEST: Method of Manufactured Solutions ===\n");
    
    // Define physical node positions in ENSIGHT order (your input)
    real_t node_coords_ensight[8][3] = {
        {0.0,  0.0,  0.0},      // node 0: ref(-1,-1,-1)
        {1.0,  0.1,  0.05},     // node 1: ref( 1,-1,-1)
        {0.9,  1.0,  0.1},      // node 2: ref( 1, 1,-1)
        {0.05, 0.95, 0.05},     // node 3: ref(-1, 1,-1)
        {0.1,  0.05, 1.0},      // node 4: ref(-1,-1, 1)
        {1.05, 0.15, 0.95},     // node 5: ref( 1,-1, 1)
        {1.0,  1.05, 1.05},     // node 6: ref( 1, 1, 1)
        {0.1,  1.0,  0.95}      // node 7: ref(-1, 1, 1)
    };
    
    // Convert to IJK ordering for BOTH analytical and numerical
    real_t node_coords_ijk[8][3];
    for(size_t n = 0; n < 8; n++) {
        size_t ijk_lid = convert_ensight_to_ijk[n];
        for(size_t i = 0; i < 3; i++) {
            node_coords_ijk[ijk_lid][i] = node_coords_ensight[n][i];
        }
    }
    
    printf("Physical node coordinates (IJK ordering):\n");
    for(int n = 0; n < 8; n++) {
        printf("  Node %d: [%9.4f, %9.4f, %9.4f]\n", 
               n, node_coords_ijk[n][0], node_coords_ijk[n][1], node_coords_ijk[n][2]);
    }
    
    // Compute basis function gradients in IJK ordering
    auto compute_hex8_basis_gradients_ijk = [](real_t xi, real_t eta, real_t zeta, 
                                               real_t grad_basis[8][3]) {
        // IJK reference coordinates
        real_t ref_coords[8][3] = {
            {-1, -1, -1},  // 0: i=0, j=0, k=0
            { 1, -1, -1},  // 1: i=1, j=0, k=0
            {-1,  1, -1},  // 2: i=0, j=1, k=0
            { 1,  1, -1},  // 3: i=1, j=1, k=0
            {-1, -1,  1},  // 4: i=0, j=0, k=1
            { 1, -1,  1},  // 5: i=1, j=0, k=1
            {-1,  1,  1},  // 6: i=0, j=1, k=1
            { 1,  1,  1}   // 7: i=1, j=1, k=1
        };
        
        for(int i = 0; i < 8; i++) {
            real_t xi_i   = ref_coords[i][0];
            real_t eta_i  = ref_coords[i][1];
            real_t zeta_i = ref_coords[i][2];
            
            grad_basis[i][0] = 0.125 * xi_i   * (1.0 + eta_i*eta)  * (1.0 + zeta_i*zeta);
            grad_basis[i][1] = 0.125 * (1.0 + xi_i*xi) * eta_i  * (1.0 + zeta_i*zeta);
            grad_basis[i][2] = 0.125 * (1.0 + xi_i*xi) * (1.0 + eta_i*eta)  * zeta_i;
        }
    };
    
    // Compute ANALYTICAL Jacobian using IJK ordering
    auto compute_analytical_jacobian = [&](real_t xi, real_t eta, real_t zeta, 
                                           real_t J[3][3]) {
        real_t grad_basis[8][3];
        compute_hex8_basis_gradients_ijk(xi, eta, zeta, grad_basis);
        
        // J_ij = Sum_k x_k^i * (partial N_k/partial xi_j)
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                J[i][j] = 0.0;
                for(int k = 0; k < 8; k++) {
                    J[i][j] += node_coords_ijk[k][i] * grad_basis[k][j];
                }
            }
        }
    };
    
    // Setup Kokkos arrays with IJK ordering
    DCArrayKokkos<double> node_coords_dual(8,3);
    DCArrayKokkos<double> node_in_elem_dual(8);
    DCArrayKokkos<double> jac(3, 3);

    for(size_t n = 0; n < 8; n++) {
        node_in_elem_dual.host(n) = n;
        for(size_t i = 0; i < 3; i++) {
            node_coords_dual.host(n,i) = node_coords_ijk[n][i];  // Already in IJK order!
        }
    }
    node_coords_dual.update_device();
    node_in_elem_dual.update_device();

    // ---- reference volume ----
    Quadrature_t Quad;
    ReferenceElement_t RefElem; 

    Quad.initialize_quadrature(reference_space::GaussLegendre, 3, 3); //(type,num_qpts_1d,elem_dims)

    RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                reference_space::LagrangeLobatto,
                                Quad,
                                1);

    // ---- reference surface ----
    SurfaceQuadrature_t SurfQuad;
    ReferenceSurface_t RefSurf;

    SurfQuad.initialize_quadrature(reference_space::GaussLegendre, 3, 3); //(type,num_qpts_1d,elem_dims)

    RefSurf.initialize_ref_surf(SurfQuad,
                                RefElem);

    // Test points
    const double qpt = 0.774596669241483377035853079956;  // sqrt(3/5)
    real_t test_points[3][3] = {
        {-qpt, -qpt, -qpt},  // qpt_id = 0
        {0.0, 0.0, 0.0},     // qpt_id = 13 
        {qpt, qpt, qpt},     // qpt_id = 26
    };
    
    printf("\nTesting Jacobian at quadrature points:\n");
    printf("%s\n", std::string(80, '=').c_str());
    
    size_t qpt_id = 0;
    bool all_passed = true;
    
    for(int p = 0; p < 3; p++) {
        real_t xi   = test_points[p][0];
        real_t eta  = test_points[p][1];
        real_t zeta = test_points[p][2];
        
        printf("\n  Point %d: (xi=%7.4f, eta=%7.4f, zeta=%7.4f)\n", p+1, xi, eta, zeta);
        printf("  %s\n", std::string(70, '-').c_str());
        
        // Analytical
        real_t J_analytical[3][3];
        compute_analytical_jacobian(xi, eta, zeta, J_analytical);
        
        printf("  Analytical Jacobian:\n");
        for(int i = 0; i < 3; i++) {
            printf("    ");
            for(int j = 0; j < 3; j++) {
                printf("%10.6f ", J_analytical[i][j]);
            }
            printf("\n");
        }

        // Numerical
        RUN({
            ViewCArrayKokkos<double> a_grad_basis(&RefElem.qpt_grad_basis(qpt_id,0,0),
                                                   RefElem.num_dofs_in_elem, 3);

            jacobian(jac, 
                    node_coords_dual, 
                    node_in_elem_dual,
                    a_grad_basis);
        });
        Kokkos::fence();
        jac.update_host();

        printf("\n  Numerical Jacobian:\n");
        for(size_t i = 0; i < 3; i++){
            printf("    ");
            for(size_t j = 0; j < 3; j++){ 
                printf("%10.6f ", jac.host(i, j));
            }
            printf("\n");
        }
        
        // Compare
        double max_error = 0.0;
        printf("\n  Absolute Error:\n");
        for(int i = 0; i < 3; i++) {
            printf("    ");
            for(int j = 0; j < 3; j++) {
                double error = fabs(J_analytical[i][j] - jac.host(i,j));
                max_error = fmax(max_error, error);
                printf("%10.2e ", error);
            }
            printf("\n");
        }
        
        // Determinant
        auto det_analytical = J_analytical[0][0]*(J_analytical[1][1]*J_analytical[2][2] - 
                                                  J_analytical[1][2]*J_analytical[2][1])
                            - J_analytical[0][1]*(J_analytical[1][0]*J_analytical[2][2] - 
                                                  J_analytical[1][2]*J_analytical[2][0])
                            + J_analytical[0][2]*(J_analytical[1][0]*J_analytical[2][1] - 
                                                  J_analytical[1][1]*J_analytical[2][0]);
        
        double det_numerical = det_3x3(jac);
        double det_error = fabs(det_analytical - det_numerical);
        
        printf("\n  Determinant:\n");
        printf("    Analytical: %12.8f\n", det_analytical);
        printf("    Numerical:  %12.8f\n", det_numerical);
        printf("    Error:      %12.2e\n", det_error);
        
        const double tol = 1e-10;
        bool passed = (max_error < tol) && (det_error < tol);
        
        printf("\n  Max Error: %12.2e\n", max_error);
        printf("  Result: %s\n", passed ? " PASSED" : "X FAILED");

        if(passed==false)Kokkos::abort("test failed \n");
        
        if(!passed) all_passed = false;
        
        qpt_id += 13;
    } // end for qpt loop

    // ===============================
    // Now testing a surface element
    // ===============================

    real_t surf_test_points[18][3] = {
        // side 0 (xi=-1)
        {-1, -qpt, -qpt},  // qpt_id = 0
        {-1, 0.0, 0.0},    // qpt_id = 4 
        {-1, qpt, qpt},    // qpt_id = 8
        // side 1 (xi=1)
        {1, -qpt, -qpt},   // qpt_id = 0
        {1, 0.0, 0.0},     // qpt_id = 4 
        {1, qpt, qpt},     // qpt_id = 8
        // side 2 (eta=-1)
        {-qpt, -1, -qpt},  // qpt_id = 0
        {0.0, -1, 0.0},    // qpt_id = 4 
        {qpt, -1, qpt},    // qpt_id = 8
        // side 3 (eta=1)
        {-qpt, 1, -qpt},   // qpt_id = 0
        {0.0, 1, 0.0},     // qpt_id = 4 
        {qpt, 1, qpt},     // qpt_id = 8
        //side 4 (mu=-1)
        {-qpt, -qpt, -1},  // qpt_id = 0
        {0.0, 0.0, -1},    // qpt_id = 4 
        {qpt, qpt, -1},    // qpt_id = 8
        //side 5 (mu=1)
        {-qpt, -qpt, 1},   // qpt_id = 0
        {0.0, 0.0, 1},     // qpt_id = 4 
        {qpt, qpt, 1}      // qpt_id = 8
    };
    
    printf("\nTesting Jacobian at quadrature points:\n");
    printf("%s\n", std::string(80, '=').c_str());
    
    // reset helper vars
    qpt_id = 0;
    all_passed = true;
    
    for(int side = 0; side<6; side++){
        printf("\n\n  side = %d\n", side);

        for(int lid = 0; lid < 3; lid++) {
            int p = lid + side*3;
            real_t xi   = surf_test_points[p][0];
            real_t eta  = surf_test_points[p][1];
            real_t zeta = surf_test_points[p][2];
            
            printf("\n  Point %d: (xi=%7.4f, eta=%7.4f, zeta=%7.4f)\n", p+1, xi, eta, zeta);
            // CPU prints only after here for verifying tests:
            //printf("  QPt # %zu: (xi=%7.4f, eta=%7.4f, zeta=%7.4f)\n", 
            //        qpt_id, 
            //        SurfQuad.qpt_positions(side,qpt_id,0),
            //        SurfQuad.qpt_positions(side,qpt_id,1),
            //        SurfQuad.qpt_positions(side,qpt_id,2)
            //      );
            printf("  %s\n", std::string(70, '-').c_str());
            
            // Analytical
            real_t J_analytical[3][3];
            compute_analytical_jacobian(xi, eta, zeta, J_analytical);
            
            printf("  Analytical Jacobian:\n");
            for(int i = 0; i < 3; i++) {
                printf("    ");
                for(int j = 0; j < 3; j++) {
                    printf("%10.6f ", J_analytical[i][j]);
                }
                printf("\n");
            }

            // Numerical
            RUN({
                // extract the grad_basis at a single quadrature point (num_dofs,3D)
                ViewCArrayKokkos<double> a_grad_basis(&RefSurf.qpt_grad_basis(side,qpt_id,0,0),
                                                    RefElem.num_dofs_in_elem, 3);

                jacobian(jac, 
                        node_coords_dual, 
                        node_in_elem_dual,
                        a_grad_basis);
            });
            Kokkos::fence();
            jac.update_host();

            printf("\n  Numerical Jacobian:\n");
            for(size_t i = 0; i < 3; i++){
                printf("    ");
                for(size_t j = 0; j < 3; j++){ 
                    printf("%10.6f ", jac.host(i, j));
                }
                printf("\n");
            }
            
            // Compare
            double max_error = 0.0;
            printf("\n  Absolute Error:\n");
            for(int i = 0; i < 3; i++) {
                printf("    ");
                for(int j = 0; j < 3; j++) {
                    double error = fabs(J_analytical[i][j] - jac.host(i,j));
                    max_error = fmax(max_error, error);
                    printf("%10.2e ", error);
                }
                printf("\n");
            }
            
            // Determinant
            auto det_analytical = J_analytical[0][0]*(J_analytical[1][1]*J_analytical[2][2] - 
                                                      J_analytical[1][2]*J_analytical[2][1])
                                - J_analytical[0][1]*(J_analytical[1][0]*J_analytical[2][2] - 
                                                      J_analytical[1][2]*J_analytical[2][0])
                                + J_analytical[0][2]*(J_analytical[1][0]*J_analytical[2][1] - 
                                                      J_analytical[1][1]*J_analytical[2][0]);
            
            double det_numerical = det_3x3(jac);
            double det_error = fabs(det_analytical - det_numerical);
            
            printf("\n  Determinant:\n");
            printf("    Analytical: %12.8f\n", det_analytical);
            printf("    Numerical:  %12.8f\n", det_numerical);
            printf("    Error:      %12.2e\n", det_error);
            
            const double tol = 1e-10;
            bool passed = (max_error < tol) && (det_error < tol);
            
            printf("\n  Max Error: %12.2e\n", max_error);
            printf("  Result: %s\n", passed ? " PASSED" : "X FAILED");
            
            if(passed==false)Kokkos::abort("test failed \n");

            if(!passed) all_passed = false;
            
            qpt_id += 4;
            if(qpt_id>8) qpt_id=0; // restart local counter on the surface
        } // end for qpt loop
    } // end for side of ref element

} // end function



int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Reference Element plus Mesh Example!"<<std::endl;
    
    Quadrature_t Quad;
    ReferenceElement_t FERefElem; // kinematic space
    ReferenceElement_t DGRefElem; // thermal space, it is discontinous

    Mesh_t Mesh; // unstructured mesh

    const size_t elem_dims = 3;
    const size_t elem_order = 3;  // cubic element will have 4 nodes in 1D
    const size_t num_elems_1D = 2;

    // the minimum quadrature for FE hydrodynamics based on elem order
    const size_t num_DOFs_1d = elem_order + 1; // 
    const size_t num_qpts_1d = 2*elem_order;   // using Legendre, but if using Lobatto, it requires 2*Order+1

    // ============================================
    // Create quadrature and reference element

    std::cout<<"Building reference elements and quadrature \n"<<std::endl;

    // create quadrature
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               num_qpts_1d,
                               elem_dims);

    // p_order is the basis order for the Lagrange polynomial defining the element
    FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLobatto,
                                  Quad,
                                  elem_order);    


    // for thermal DOFs, we use p_order-1 for the basis order of the Lagrange polynomial, it is DG
    DGRefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLegendre,
                                  Quad,
                                  elem_order-1); 

    // ==========================================
    // Build the unstructured mesh structure

    std::cout<<"Building unstructed mesh \n"<<std::endl;

    const size_t num_elems    = num_elems_1D*num_elems_1D*num_elems_1D;
    const size_t num_nodes_1D = elem_order*num_elems_1D + 1;  // number of nodes
    const size_t num_nodes    = num_nodes_1D*num_nodes_1D*num_nodes_1D;

    Mesh.initialize_dims(elem_dims);
    Mesh.initialize_elems_Pn(num_elems, elem_order, Quad.num_qpts_1d);
    Mesh.initialize_nodes(num_nodes);
    
    DCArrayKokkos <double> node_coords(Mesh.num_nodes, Mesh.num_dims);
    const double h = 1.0/((double)num_nodes_1D);

    // create indexing for a Pn order mesh

    // Step 1: Initialize ALL node coordinates once (no race condition)
    FOR_ALL(kc, 0, num_nodes_1D,
            jc, 0, num_nodes_1D,
            ic, 0, num_nodes_1D, {
        
        size_t node_gid = ic + (jc + kc*num_nodes_1D)*num_nodes_1D;
        
        node_coords(node_gid, 0) = (double)ic * h;
        node_coords(node_gid, 1) = (double)jc * h;
        node_coords(node_gid, 2) = (double)kc * h;
    });

    // Step 2: Build element connectivity
    FOR_ALL(i, 0, num_elems_1D,
            j, 0, num_elems_1D,
            k, 0, num_elems_1D, {
        
        size_t elem_gid = i + (j + k*num_elems_1D)*num_elems_1D;
        size_t node_lid = 0;
        
        // FIX: Multiply by elem_order to space elements correctly
        for(size_t kc = k*elem_order; kc <= k*elem_order + elem_order; kc++)
        for(size_t jc = j*elem_order; jc <= j*elem_order + elem_order; jc++)
        for(size_t ic = i*elem_order; ic <= i*elem_order + elem_order; ic++){
            size_t node_gid = ic + (jc + kc*num_nodes_1D)*num_nodes_1D;
            Mesh.nodes_in_elem(elem_gid, node_lid) = node_gid;
            node_lid++;
        }
    });

    Mesh.build_corner_connectivity();
    Mesh.build_elem_elem_connectivity();
    Mesh.build_surf_connectivity();

    // check mesh index sizes
    if(Mesh.num_nodes!=num_nodes){
        printf("num nodes = %zu and mesh.num_nodes = %zu", num_nodes, Mesh.num_nodes);
        Kokkos::abort("ERROR: wrong number of mesh nodes");
    }
    if(Mesh.num_gauss_in_elem!=Quad.num_qpts_in_elem){
        Kokkos::abort("ERROR: wrong number of Gauss points in elem");
    }

    
    printf("\n=== TEST: 2x2x2 Cubic Mesh Face Pair Test ===\n");

    // Face 1 of element 0 matches face 0 of element 1
    const size_t elem_nodes_in_face_match[16] = 
        {3, 10, 17, 24, 52, 59, 66, 73, 101, 108, 115, 122, 150, 157, 164, 171};

    // Check surfaces and patches
    FOR_ALL(elem_gid, 0, Mesh.num_elems, {
        
        for(size_t face_lid = 0; face_lid < Mesh.num_surfs_in_elem; face_lid++) {
            
            // Only test the specific faces we want to validate
            if(!((elem_gid == 0 && face_lid == 1) || (elem_gid == 1 && face_lid == 0))) {
                continue;
            }
            
            size_t count = 0;  // Count unique nodes found across all patches in this face
            bool found[16] = {false};  // Track which of the 16 face nodes we've found
            
            for(size_t patch_surf_lid = 0; patch_surf_lid < Mesh.num_patches_in_surf; patch_surf_lid++) {
                
                const size_t patch_lid = face_lid * Mesh.num_patches_in_surf + patch_surf_lid;  // FIX: num_patches_in_surf
                const size_t patch_gid = Mesh.patches_in_elem(elem_gid, patch_lid);
                const size_t surf_gid  = Mesh.surf_in_patch(patch_gid);
                
                // Validate surface ID
                if(surf_gid > 5) {
                    printf("ERROR: should be <= 5, but surf_gid = %d \n", (int)surf_gid);
                    Kokkos::abort("Surface gid is out of bounds\n");
                }
                
                // Check each node in this patch
                for(size_t node_lid = 0; node_lid < 4; node_lid++) {
                    const size_t node_gid = Mesh.nodes_in_patch(patch_gid, node_lid);
                    
                    // Check if this node is in our expected face list
                    for(size_t i = 0; i < 16; i++) {
                        if(elem_nodes_in_face_match[i] == node_gid && !found[i]) {
                            found[i] = true;
                            count++;
                            break;
                        }
                    }
                }
            } // end patch loop
            
            // Now check that we found all 16 unique nodes
            if(count != 16) {
                printf("ERROR: Elem %d, Face %zu found only %zu/16 nodes\n", 
                    elem_gid, face_lid, count);
                Kokkos::abort("Face nodes don't match!\n");
            }
            
        } // end face loop
    });

    printf("Face matching test PASSED!\n\n");



    // Test bidirectional surface-element mapping consistency
    DCArrayKokkos <size_t> bdy_surf_counter(1);
    bdy_surf_counter.set_values(0);

    FOR_ALL(surf_gid, 0, Mesh.num_surfs, {
        
        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);
        if(num_elems_in_surf==1) Kokkos::atomic_add(&bdy_surf_counter(0),1);
        
        // Test should work for both boundary (1 elem) and interior (2 elems) surfaces
        for(size_t side_lid = 0; side_lid < num_elems_in_surf; side_lid++) {
            
            const int elem_gid = Mesh.elems_in_surf(surf_gid, side_lid);
            const int face_lid = Mesh.faces_in_surf(surf_gid, side_lid);
            
            // Verify elem_gid is valid
            if(elem_gid < 0 || elem_gid >= Mesh.num_elems) {
                printf("ERROR: surf %d side_lid %zu has invalid elem_gid = %d\n", 
                    surf_gid, side_lid, elem_gid);
                Kokkos::abort("Invalid element in surface\n");
            }
            
            // Verify face_lid is valid
            if(face_lid < 0 || face_lid >= Mesh.num_surfs_in_elem) {
                printf("ERROR: surf %d side_lid %zu has invalid face_lid = %d\n", 
                    surf_gid, side_lid, face_lid);
                Kokkos::abort("Invalid face in surface\n");
            }
            
            // *** CRITICAL TEST: Verify reverse mapping ***
            const size_t reverse_surf_gid = Mesh.surfs_in_elem(elem_gid, face_lid);
            if(reverse_surf_gid != surf_gid) {
                printf("ERROR: Mapping inconsistency!\n");
                printf("  surf %d -> elem %d, face %d\n", surf_gid, elem_gid, face_lid);
                printf("  BUT elem %d, face %d -> surf %zu\n", 
                    elem_gid, face_lid, reverse_surf_gid);
                Kokkos::abort("Failed reverse mapping surf->elem->surf\n");
            }
            
        } // end side_lid loop
        
    }); // end surf loop
    Kokkos::fence();
    bdy_surf_counter.update_host();
    printf("Surface-Element bidirectional mapping test PASSED!\n");

    if(Mesh.num_bdy_surfs!=bdy_surf_counter.host(0)){
        Kokkos::abort("Failed to find correct number of boundary surfaces\n");
    } 
    if(Mesh.num_bdy_patches!=bdy_surf_counter.host(0)*Mesh.num_patches_in_surf) Kokkos::abort("Failed to find correct number of boundary patches\n");
    printf("Boundary surfaces and patches test PASSED!\n");

    // for this 2x2x2 mesh, there are:
    // 12 interior surfaces + 24 boundary surfaces
    //     3 planes of 2x2 in x-dir = 12
    //     3 planes of 2x2 in y-dir = 12
    //     3 planes of 2x2 in z-dir = 12
    // 36 total surfaces
    if(Mesh.num_surfs!=36) Kokkos::abort("Wrong number of surfaces\n");
    if(Mesh.num_bdy_surfs!=24) Kokkos::abort("Wrong number of boundary surfaces\n");
    printf("\n\n");


    // ==========================================
    // Create state on the unstructured mesh structure

    DCArrayKokkos <double> node_scalar(Mesh.num_nodes);
    DCArrayKokkos <double> gauss_scalar(Mesh.num_gauss_in_elem);


    // ===========
    RUN({

        size_t elem_gid = 0; // testing first element
        size_t qpt = 0;      // testing geometry at first quadrature point

        const size_t num_nodes_in_elem = FERefElem.num_dofs_in_elem;

        // get Jacobian
        double j_1D[9];
        ViewCArrayKokkos <double> jac(&j_1D[0], 3, 3); 

        // make views on device, extacting only element and quadrature point values
        ViewCArrayKokkos <size_t> nodes_in_an_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);
        ViewCArrayKokkos <double> a_grad_basis(&FERefElem.qpt_grad_basis(qpt, 0, 0), num_nodes_in_elem, 3);

        jacobian(jac, 
                 node_coords, 
                 nodes_in_an_elem,
                 a_grad_basis);

        printf("jacobian matrix = \n");
        for(size_t i = 0; i < elem_dims; i++){  // looping over dimension
            for(size_t j = 0; j < elem_dims; j++){ // looping over dimension
               printf("%f, ", jac(i, j));
           }
           printf("\n");
        } // end for 
    });


    printf("\nReference plus mesh test finished.\n");

    // running unit tests of Jacobian
    test_affine_transformation(); 
    test_skewed_hexahedron();
    test_rotated_scaled_cube();
    test_manufactured_solution();


} // end MATAR scope
MATAR_FINALIZE();

return 0;
}
