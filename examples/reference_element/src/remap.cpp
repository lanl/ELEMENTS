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


inline int PointIndexFromIJK(int i, int j, int k, const int* order);
void write_vtk(const char* filename, 
               Mesh_t& Mesh,
               DCArrayKokkos<double>& node_coords,
               DCArrayKokkos<double>& elem_field,
               int elem_order);

int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Reference Element Remap Example!"<<std::endl;
    
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
    const double h = 1.0/((double)num_nodes_1D-1);

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


    // State arrays for this test
    const size_t num_qpts_in_elem = Quad.num_qpts_in_elem;

    const size_t num_surfs = Mesh.num_surfs;
    const size_t num_qpts_in_surf = SurfQuad.num_qpts_in_surf;

    const size_t num_nodes_in_elem = Mesh.num_nodes_in_elem;
    const size_t num_surfs_in_elem = Mesh.num_surfs_in_elem;

    if(num_nodes_in_elem != FERefElem.num_dofs_in_elem) Kokkos::abort("ERROR: mismatch in DOFs and num nodes in elem \n");

    DCArrayKokkos<double> elem_jac(num_elems, num_qpts_in_elem, elem_dims, elem_dims, "elem_jacobian");
    DCArrayKokkos<double> elem_det_jac(num_elems, num_qpts_in_elem, "elem_det_jacobian");
    DCArrayKokkos<double> elem_vol(num_elems, "elem_vol");
    DCArrayKokkos<double> elem_vol_n(num_elems, "elem_vol_n");
    
    
    DCArrayKokkos<double> surf_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_jacobian");
    DCArrayKokkos<double> surf_inv_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_inv_jacobian");
    DCArrayKokkos<double> surf_flux(num_surfs, "surf_flux");
    
    DCArrayKokkos<double> elem_field(num_elems, "elem_field"); //elem_order-1 = 1 so it is a P0 element
    DCArrayKokkos<double> elem_field_n(num_elems, "elem_field_n"); //elem_order-1 = 1 so it is a P0 element
    DCArrayKokkos<double> node_velocity(num_nodes, elem_dims, "node_velocity");


    const double max_vel = 1.0;
    const size_t max_cycles = 10;
    double h_cfl = h;
    double dt = 0.2*h_cfl/max_vel;  // dt from CFL at start, this time is psuedo time

    FOR_ALL(elem_gid, 0, num_elems, {

        elem_field(elem_gid) = 1.0; // constant field

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        double elem_coords[3];
        for(size_t dim=0; dim<elem_dims; dim++){elem_coords[dim] = 0.0;}

        for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
            const size_t node_gid = nodes_in_elem(node_lid);
            for(size_t dim=0; dim<elem_dims; dim++){ 
                    elem_coords[dim] += node_coords(node_gid,dim)/((double)num_nodes_in_elem);
            }
        }

        // a spatially varing field
        elem_field(elem_gid) = sin(PI*elem_coords[0]);

    });
    Kokkos::fence();

    // Step 0: calculate the element volume 
    FOR_ALL(elem_gid, 0, num_elems, {

        elem_vol(elem_gid) = 0.0; 

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
            ViewCArrayKokkos<double> jac(&elem_jac(elem_gid,qpt_lid,0,0),3,3);
            jacobian(jac, 
                    node_coords, 
                    nodes_in_elem,
                    a_grad_basis);

            // calculate det_J 
            elem_det_jac(elem_gid, qpt_lid) = det_3x3(jac);
            
            elem_vol(elem_gid) += elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

        } // end for

    });
    Kokkos::fence();


    ////--------///

    for(size_t cycle=0; cycle<max_cycles; cycle++){

        // Step 1a: store level n state
        FOR_ALL(elem_gid, 0, num_elems, {
            elem_vol_n(elem_gid)   = elem_vol(elem_gid);
            elem_field_n(elem_gid) = elem_field(elem_gid);
        });


        // Step 1b: calculate the mesh velocity
        FOR_ALL(node_gid, 0, num_nodes,{
            // new velocity, it is Taylor-Green vortex
            // PI is defined in mesh class
            node_velocity(node_gid, 0) = sin(PI*node_coords(node_gid, 0))*cos(PI*node_coords(node_gid, 1));
            node_velocity(node_gid, 1) = -cos(PI*node_coords(node_gid, 0))*sin(PI*node_coords(node_gid, 1));
            node_velocity(node_gid, 2) = 0.0;
        });


        // Step 2: get CFL pseudo time step for moving mesh
        double min_h_loc;
        FOR_REDUCE_MIN(elem_gid, 0, num_elems, 
                       min_h_loc, { 
            
            for(size_t qpt_lid=0; qpt_lid<num_qpts_in_elem; qpt_lid++){
                const double vol_qpt= Quad.qpt_weights(qpt_lid)*elem_det_jac(elem_gid, qpt_lid);
                const double h_qpt = pow(vol_qpt,0.3333333);
                if(h_qpt < min_h_loc) min_h_loc = h_qpt;
            }

        }, h_cfl);
        Kokkos::fence();

        dt = 0.2*h_cfl/max_vel; // pseudo time step


        // Step 3: calculate the surface fluxes
        surf_flux.set_values(0.0);
        FOR_ALL(surf_gid, 0, Mesh.num_surfs, {
            
            const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);

            const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0); // first element in surf
            const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0); // first element in surf

            ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

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

                const double det_jac_qpt = det_3x3(jac);

                invert_3x3(jac, inv_jac, det_jac_qpt);

                // Nanson's formula: s*J^-1*j*f*w
                double area_normal[3];
                area_normal[0] = 0.;
                area_normal[1] = 0.;
                area_normal[2] = 0.;
                for(size_t j=0; j<elem_dims; j++){ 
                    for(size_t i=0; i<elem_dims; i++){
                        area_normal[j] += RefSurf.outward_normal(face_lid,i)*inv_jac(i,j);
                    } // end i
                    area_normal[j] *= det_jac_qpt*SurfQuad.qpt_weights(face_lid,qpt_lid);
                } // end j

                double qpt_vel[3];
                for(size_t dim=0; dim<elem_dims; dim++){
                    qpt_vel[dim] = 0.0;
                }

                for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++)
                for(size_t dim=0; dim<elem_dims; dim++){
                    const size_t node_gid = nodes_in_elem(node_lid);
                    qpt_vel[dim] += a_basis(node_lid)*node_velocity(node_gid,dim);
                } // end for

                //printf("qpt normal x = %.4f \n", area_normal[0]);

                double normal_dot_vel = 0.0;
                for(size_t dim=0; dim<elem_dims; dim++){
                    normal_dot_vel += area_normal[dim]*qpt_vel[dim];
                }
                //printf("normal dot vel = %.4f \n", normal_dot_vel);
                
                size_t nbr_elem_gid = elem_gid;
                size_t nbr_face_lid = face_lid;
                if(num_elems_in_surf==2){
                    nbr_elem_gid = Mesh.elems_in_surf(surf_gid, 1); // second elem
                    nbr_face_lid = Mesh.faces_in_surf(surf_gid, 1); // second elem
                }    
                    
                surf_flux(surf_gid) += 
                    -0.5*(elem_field(elem_gid)+elem_field(nbr_elem_gid))*normal_dot_vel 
                    +0.5*fabs(normal_dot_vel)*(elem_field(elem_gid)-elem_field(nbr_elem_gid));

            } // end for qpt
            
        }); // end surf loop
        Kokkos::fence();


        // Step 4: Move the mesh to the new location
        FOR_ALL(node_gid, 0, num_nodes,{
            // new position of the mesh
            node_coords(node_gid, 0) += node_velocity(node_gid, 0)*dt; 
            node_coords(node_gid, 1) += node_velocity(node_gid, 1)*dt;
            // z-coords never change
        });
        Kokkos::fence();


        // Step 4: calculate the element volume 
        FOR_ALL(elem_gid, 0, num_elems, {

            elem_vol(elem_gid) = 0.0; 

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
                ViewCArrayKokkos<double> jac(&elem_jac(elem_gid,qpt_lid,0,0),3,3);
                jacobian(jac, 
                        node_coords, 
                        nodes_in_elem,
                        a_grad_basis);

                // calculate det_J 
                elem_det_jac(elem_gid, qpt_lid) = det_3x3(jac);
                
                elem_vol(elem_gid) += elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

            } // end for

        });
        Kokkos::fence();


        // Step 5: update the elem field
        FOR_ALL(elem_gid, 0, num_elems,{

            // Start with conservative quantity from previous time
            double conservative_quantity = elem_field_n(elem_gid) * elem_vol_n(elem_gid);
            
            // Accumulate fluxes from all faces
            for(size_t face_lid=0; face_lid<num_surfs_in_elem; face_lid++){
                const size_t surf_gid = Mesh.surfs_in_elem(elem_gid, face_lid);
                
                // Determine sign: flux is positive OUT of first element
                double flux_sign = 1.0;
                if(Mesh.elems_in_surf(surf_gid, 0) != elem_gid){
                    flux_sign = -1.0;  // We're the second element, flux goes in
                }
                
                conservative_quantity -= flux_sign * surf_flux(surf_gid) * dt;
            } // end for
        
            // Update field by dividing by new volume
            elem_field(elem_gid) = conservative_quantity / elem_vol(elem_gid);
            
        });

        // After your cycle loop ends, add:
        char filename[100];
        sprintf(filename, "output_cycle_%04zu.vtk", cycle); // Inside cycle loop
        // OR
        //sprintf(filename, "output_final.vtk"); // After all cycles

        write_vtk(filename, Mesh, node_coords, elem_field, elem_order);

    } // end loop over cycle


/*
    FOR_ALL(elem_gid, 0, num_elems,{
        printf("%.4f\n", elem_field(elem_gid));
    });
    Kokkos::fence();

    FOR_ALL(node_gid, 0, num_nodes,{
        printf("%.4f, %.4f \n", node_coords(node_gid, 0), node_coords(node_gid, 1));
    });
*/

    printf("\n Remap test finished.\n");


} // end MATAR scope
MATAR_FINALIZE();

return 0;
} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn PointIndexFromIJK
///
/// \brief Given (i,j,k) coordinates within the Lagrange hex, return an 
/// offset into the local connectivity (PointIds) array. The order parameter
/// must point to an array of 3 integers specifying the order along each 
/// axis of the hexahedron.
///
/////////////////////////////////////////////////////////////////////////////
inline int PointIndexFromIJK(int i, int j, int k, const int* order)
{
    bool ibdy = (i == 0 || i == order[0]);
    bool jbdy = (j == 0 || j == order[1]);
    bool kbdy = (k == 0 || k == order[2]);
    // How many boundaries do we lie on at once?
    int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

    if (nbdy == 3) { // Vertex DOF
        // ijk is a corner node. Return the proper index (somewhere in [0,7]):
        return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    }

    int offset = 8;
    if (nbdy == 2) { // Edge DOF
        if (!ibdy) { // On i axis
            return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) + offset;
        }
        if (!jbdy) { // On j axis
            return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) + offset;
        }
        // !kbdy, On k axis
        offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
        return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
    }

    offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
    if (nbdy == 1) { // Face DOF
        if (ibdy) { // On i-normal face
            return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
        }
        offset += 2 * (order[1] - 1) * (order[2] - 1);
        if (jbdy) { // On j-normal face
            return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
        }
        offset += 2 * (order[2] - 1) * (order[0] - 1);
        // kbdy, On k-normal face
        return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
    }

    // nbdy == 0: Body DOF
    offset += 2 * ( (order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) + (order[0] - 1) * (order[1] - 1));
    return offset + (i - 1) + (order[0] - 1) * ( (j - 1) + (order[1] - 1) * ( (k - 1)));
}

void write_vtk(const char* filename, 
               Mesh_t& Mesh,
               DCArrayKokkos<double>& node_coords,
               DCArrayKokkos<double>& elem_field,
               int elem_order)
{
    // Update host side data
    node_coords.update_host();
    elem_field.update_host();
    Mesh.nodes_in_elem.update_host();
    
    
    FILE* vtk_file = fopen(filename, "w");
    if (vtk_file == NULL) {
        printf("Error: Could not open VTK file %s\n", filename);
        return;
    }
    
    // Write VTK header
    fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
    fprintf(vtk_file, "Remap Test Output\n");
    fprintf(vtk_file, "ASCII\n");
    fprintf(vtk_file, "DATASET UNSTRUCTURED_GRID\n");
    
    // Write points (nodes)
    fprintf(vtk_file, "POINTS %lu double\n", Mesh.num_nodes);
    for (size_t node_gid = 0; node_gid < Mesh.num_nodes; node_gid++) {
        fprintf(vtk_file, "%.8e %.8e %.8e\n",
                node_coords.host(node_gid, 0),
                node_coords.host(node_gid, 1),
                node_coords.host(node_gid, 2));
    }
    
    // Write cells (elements)
    size_t num_nodes_in_elem = Mesh.num_nodes_in_elem;
    
    // For VTK, we need to specify: number of cells, and total size
    // size = num_elems * (1 + num_nodes_in_elem) where 1 is for the count
    size_t total_size = Mesh.num_elems * (1 + num_nodes_in_elem);
    fprintf(vtk_file, "CELLS %lu %lu\n", Mesh.num_elems, total_size);
    
    int order[3]   = { elem_order, elem_order, elem_order };
    for (size_t elem_gid = 0; elem_gid < Mesh.num_elems; elem_gid++) {
        fprintf(vtk_file, "%lu", num_nodes_in_elem);
        //for (int k = 0; k <= elem_order; k++) {
        //    for (int j = 0; j <= elem_order; j++) {
        //        for (int i = 0; i <= elem_order; i++) {
        //            size_t node_lid = PointIndexFromIJK(i, j, k, order);
        //            fprintf(vtk_file, "%lu ", Mesh.nodes_in_elem.host(elem_gid, node_lid));
        //        }
        //    }
        //}
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            fprintf(vtk_file, " %lu", Mesh.nodes_in_elem.host(elem_gid, node_lid));
        }
        fprintf(vtk_file, "\n");
    }
    
    // Write cell types
    fprintf(vtk_file, "CELL_TYPES %lu\n", Mesh.num_elems);
    
    // Determine VTK cell type based on element order and dimension
    int vtk_cell_type;
    //if (elem_order == 1 && Mesh.num_dims == 3) {
    //    vtk_cell_type = 12; // VTK_HEXAHEDRON (linear hex)
    //} else if (elem_order == 2 && Mesh.num_dims == 3) {
    //    vtk_cell_type = 29; // VTK_TRIQUADRATIC_HEXAHEDRON
    //} else {
        // Default to higher order hex
        vtk_cell_type = 72; // VTK_LAGRANGE_HEXAHEDRON
    //}
    vtk_cell_type = 11;
    
    for (size_t elem_gid = 0; elem_gid < Mesh.num_elems; elem_gid++) {
        fprintf(vtk_file, "%d\n", vtk_cell_type);
    }
    
    // Write cell data (element field)
    fprintf(vtk_file, "CELL_DATA %lu\n", Mesh.num_elems);
    fprintf(vtk_file, "SCALARS elem_field double 1\n");
    fprintf(vtk_file, "LOOKUP_TABLE default\n");
    for (size_t elem_gid = 0; elem_gid < Mesh.num_elems; elem_gid++) {
        fprintf(vtk_file, "%.8e\n", elem_field.host(elem_gid));
    }

    /*
    fprintf(vtk_file, "POINT_DATA %lu\n", Mesh.num_nodes);
    fprintf(vtk_file, "VECTORS node_velocity double\n");

    node_velocity.update_host();
    for (size_t node_gid = 0; node_gid < Mesh.num_nodes; node_gid++) {
        fprintf(vtk_file, "%.8e %.8e %.8e\n",
                node_velocity.host(node_gid, 0),
                node_velocity.host(node_gid, 1),
                node_velocity.host(node_gid, 2));
    }
    */
    
    fclose(vtk_file);
    printf("VTK file written to: %s\n", filename);
}