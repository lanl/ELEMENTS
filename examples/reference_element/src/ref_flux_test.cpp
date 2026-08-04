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


// 8x8x8 linear mesh
// num_patches = 8*8*9*3 = 1728
// bdy_patches = 8*8*6 = 384

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Reference Element Flux Example!"<<std::endl;
    
    Quadrature_t Quad;
    ReferenceElement_t FERefElem; // kinematic space
    ReferenceElement_t DGRefElem; // thermal space, it is discontinous

    SurfaceQuadrature_t SurfQuad;
    ReferenceSurface_t RefSurf;

    Mesh_t Mesh; // unstructured mesh

    const size_t elem_dims = 3;
    const size_t elem_order = 1;  
    const size_t num_elems_1D = 8;

    // the minimum quadrature for FE hydrodynamics based on elem order
    const size_t num_DOFs_1d = elem_order + 1; // 
    const size_t num_qpts_1d = 2*elem_order;   // using Legendre, but if using Lobatto, it requires 2*Order+1

    // ================================================================
    // Create quadrature along with the reference element and surface

    std::cout<<"Building reference elements and quadrature \n";

    // ---- reference element ----

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

    // ---- reference surface ----
    SurfQuad.initialize_quadrature(reference_space::GaussLegendre, 
                                   num_qpts_1d, 
                                   elem_dims); 

    RefSurf.initialize_ref_surf(SurfQuad,
                                FERefElem);


    // ==========================================
    // Build the unstructured mesh structure

    std::cout<<"Building unstructured mesh \n";

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

    std::cout<<"Building corner connectivity \n";
    Mesh.build_corner_connectivity();
    std::cout<<"Building element element connectivity \n";
    Mesh.build_elem_elem_connectivity();
    std::cout<<"Building surface connectivity \n";
    Mesh.build_surf_connectivity();

    // check mesh index sizes
    if(Mesh.num_nodes!=num_nodes){
        printf("num nodes = %zu and mesh.num_nodes = %zu", num_nodes, Mesh.num_nodes);
        Kokkos::abort("ERROR: wrong number of mesh nodes");
    }
    if(Mesh.num_gauss_in_elem!=Quad.num_qpts_in_elem){
        Kokkos::abort("ERROR: wrong number of Gauss points in elem");
    }
    if(Mesh.num_patches!=1728){
        Kokkos::abort("ERROR: wrong number of patches");
    }
    if(Mesh.num_bdy_patches!=384){
        Kokkos::abort("ERROR: wrong number of boundary patches");
    }


    // State arrays for this test
    const size_t num_qpts_in_elem = Quad.num_qpts_in_elem;

    const size_t num_surfs = Mesh.num_surfs;
    const size_t num_qpts_in_surf = SurfQuad.num_qpts_in_surf;

    const size_t num_nodes_in_elem = Mesh.num_nodes_in_elem;
    const size_t num_surfs_in_elem = Mesh.num_surfs_in_elem;

    if(num_nodes_in_elem != FERefElem.num_dofs_in_elem) Kokkos::abort("ERROR: mismatch in DOFs and num nodes in elem \n");

    DCArrayKokkos<double> elem_jac(num_qpts_in_elem, elem_dims, elem_dims, "elem_jacobian");
    DCArrayKokkos<double> surf_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_jacobian");
    DCArrayKokkos<double> surf_inv_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_inv_jacobian");
    DCArrayKokkos<double> surf_flux(num_surfs, num_qpts_in_surf, elem_dims, "surf_flux");
    DCArrayKokkos<double> field(num_elems, "elem_field"); //elem_order-1 = 1 so it is a P0 element
    DCArrayKokkos<double> mesh_velocity(num_surfs, num_qpts_in_surf, elem_dims, "surf_mesh_velocity");
    mesh_velocity.set_values(0.0);
    
    printf("\n=== TEST: Surface Integration ===\n");

    

    FOR_ALL(elem_gid, 0, num_elems, {

        double elem_vol = 0.0; 

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        // calculate volume
        for(size_t qpt_lid=0; qpt_lid<num_qpts_in_elem; qpt_lid++){
            // extract the grad_basis at a single quadrature point (qpt,dof,3D)
            ViewCArrayKokkos<double> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0),
                                                  num_nodes_in_elem, 3);

            // extract the basis at a single quadrature point (qpt,dof)    
            ViewCArrayKokkos<double> a_basis(&FERefElem.qpt_basis(qpt_lid,0),
                                             num_nodes_in_elem);
            
            // jacobian matrix
            ViewCArrayKokkos<double> jac(&elem_jac(qpt_lid,0,0),3,3);

            jacobian(jac, 
                    node_coords, 
                    nodes_in_elem,
                    a_grad_basis);

            const double det_jac = det_3x3(jac);
            
            elem_vol += det_jac*Quad.qpt_weights(qpt_lid);

        } // end for

        if(fabs(elem_vol-pow(h,elem_dims))>1.0e-12){
            printf("wrong elem volume using gauss quadrature, error = %.15f \n", fabs(elem_vol-pow(h,elem_dims)) );
            Kokkos::abort("ERROR: wrong volume \n");
        }


        double tally_div=0;    
        double tally_flux[3][3];
        double tally_normal[3];

        double eye[3][3];

        for(size_t j=0; j<elem_dims; j++){ 
            for(size_t i=0; i<elem_dims; i++){   
                tally_flux[i][j] = 0.;  
                eye[i][j] = 0.;  
            } // end i
            tally_normal[j] = 0.;
        } // end j   
        for(size_t j=0; j<elem_dims; j++){
            eye[j][j] = 1.0;
        }    

        for(size_t face_lid=0; face_lid<num_surfs_in_elem; face_lid++){

            const size_t surf_gid = Mesh.surfs_in_elem(elem_gid, face_lid);

            for(size_t qpt_lid=0; qpt_lid<num_qpts_in_surf; qpt_lid++){
                
                // extract the grad_basis at a single quadrature point (surf,qpt,dof,3D)
                ViewCArrayKokkos<double> a_grad_basis(&RefSurf.qpt_grad_basis(face_lid,qpt_lid,0,0),
                                                      num_nodes_in_elem, 3);

                // extract the basis at a single quadrature point (surf,qpt,dof)    
                ViewCArrayKokkos<double> a_basis(&RefSurf.qpt_basis(face_lid,qpt_lid,0),
                                                 num_nodes_in_elem);
                
                ViewCArrayKokkos<double> jac(&surf_jac(surf_gid,qpt_lid,0,0),3,3);
                ViewCArrayKokkos<double> inv_jac(&surf_inv_jac(surf_gid,qpt_lid,0,0),3,3);

                jacobian(jac, 
                        node_coords, 
                        nodes_in_elem,
                        a_grad_basis);

                const double det_jac = det_3x3(jac);

                invert_3x3(jac, inv_jac, det_jac);

                // Nanson's formula: s*J^-1*j*f*w
                double area_normal[3];
                area_normal[0] = 0.;
                area_normal[1] = 0.;
                area_normal[2] = 0.;
                for(size_t j=0; j<elem_dims; j++){ 
                    for(size_t i=0; i<elem_dims; i++){
                        area_normal[j] += RefSurf.outward_normal(face_lid,i)*inv_jac(i,j);
                    } // end i
                    area_normal[j] *= det_jac*SurfQuad.qpt_weights(face_lid,qpt_lid);
                } // end j


                // ============================
                //  testing a area normals

                double mag_sqrd = 0.0;
                for(size_t dim=0; dim<elem_dims; dim++){
                    mag_sqrd += area_normal[dim]*area_normal[dim];    
                }

                if(mag_sqrd < 1.e-16){
                    printf("Zero area normal = %f, %f, %f \n", area_normal[0], area_normal[1], area_normal[2]);
                    Kokkos::abort("ERROR: area normal = 0 \n");
                }

                // conservation of area normals, tally must = 0
                for(size_t dim=0; dim<elem_dims; dim++){
                    tally_normal[dim] += area_normal[dim];
                }


                // ========================================================
                //  testing a linear field, answer-> identity & trace=3

                double qpt_coords[3];
                for(size_t dim=0; dim<elem_dims; dim++){
                    qpt_coords[dim] = 0.0;
                }

                for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++)
                for(size_t dim=0; dim<elem_dims; dim++){
                    const size_t node_gid = nodes_in_elem(node_lid);
                    qpt_coords[dim] += a_basis(node_lid)*node_coords(node_gid,dim);
                } // end for

                //printf("qpt coord = (%f, %f, %f) \n", qpt_coords[0], qpt_coords[1], qpt_coords[2]);

                // divergence check, a Tr(grad x) = 3
                for(size_t dim=0; dim<elem_dims; dim++){
                    tally_div += area_normal[dim]*qpt_coords[dim];
                }                
                
                for(size_t i=0; i<elem_dims; i++)
                for(size_t j=0; j<elem_dims; j++){
                    tally_flux[i][j] += area_normal[j]*qpt_coords[i];
                }

            } // end for qpt
        } // end for face_lid


        // Test: check normals sum to zero over the element
        for(size_t dim=0; dim<elem_dims; dim++){
            if(fabs(tally_normal[dim]) >= 1.e-12){
                printf("Tally_normal of area normals in elem (should=0) = (%f, %f, %f) \n", tally_normal[0], tally_normal[1], tally_normal[2]);
                Kokkos::abort("ERROR: no conservation of surface normals \n");
            } // end check on tally_normal
        } // end for dim


        // Test: check grad x = eye
        bool passed = true;
        for(size_t i=0; i<elem_dims; i++) 
        for(size_t j=0; j<elem_dims; j++) {
            if(fabs((tally_flux[i][j]/elem_vol) - eye[i][j]) >= 1.e-12) passed = false;
        }
        if(!passed){
            printf("\n");
            printf("eye = gradx = \n ");
            for(size_t i=0; i<elem_dims; i++) {
                for(size_t j=0; j<elem_dims; j++) {
                    printf("%f, ", tally_flux[i][j]/elem_vol);
                }
                printf("\n");
            }
            printf("\n");
            Kokkos::abort("ERROR: Gradx should be equal to eye \n");
        }

        // Test: divergence is = 3, for some reason, tol needs to be 1e-10
        if(fabs((tally_div/elem_vol) - 3.0) >= 1.e-10){
            printf("Tally_div = %.15f \n",tally_div/elem_vol-3.0);
            Kokkos::abort("ERROR: divergence of Gradx should be equal to 3 \n");
        }

    }); // end parallel for elements



    FOR_ALL(surf_gid, 0, Mesh.num_surfs, {
        
        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);
        for(size_t side_lid = 0; side_lid < num_elems_in_surf; side_lid++) {
            
            const int elem_gid = Mesh.elems_in_surf(surf_gid, side_lid);
            const int face_lid = Mesh.faces_in_surf(surf_gid, side_lid);
            
            
        } // end side_lid loop
        
    }); // end surf loop
    Kokkos::fence();

    


    printf("\n Surface flux test finished.\n");


} // end MATAR scope
MATAR_FINALIZE();

return 0;
}
