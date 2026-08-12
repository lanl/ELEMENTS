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

// for VTU writing
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

// This pulls in kokkos, matar, mesh, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"
#include "cramers_rule.hpp" // det and solvers
#include "lu_solver.hpp"

using namespace mtr;
using namespace swage;    // unstructured mesh and point cloud
using namespace elements; // reference element space


void write_lagrange_hex_mesh(
    const std::string& filename,
    const DCArrayKokkos<double>& node_coords,       // All node coordinates [num_nodes][3]
    const size_t num_nodes,
    const DCArrayKokkos<size_t>& nodes_in_elem,     // Connectivity
    const size_t num_elems,
    const size_t order,
    const DCArrayKokkos<double>& node_data,         // Nodal data
    const std::string& node_data_name,
    const DCArrayKokkos<double>& elem_data,         // Element center data (NEW)
    const std::string& elem_data_name);             // Element data name (NEW)


void write_lagrange_cells(std::ofstream& file, 
                        const DCArrayKokkos<size_t>& nodes_in_elem,
                        size_t num_elems, 
                        size_t order,
                        const DCArrayKokkos<double>& elem_data,    // Element data
                        const std::string& elem_data_name);        // Element data name

void write_points(std::ofstream& file, 
                  const DCArrayKokkos<double>& coords, 
                  size_t num_nodes);


void write_point_data(std::ofstream& file, 
                      const DCArrayKokkos<double>& data, 
                      size_t num_nodes,
                      const std::string& name);

void reorder_ijk_to_vtk_lagrange(const DCArrayKokkos<size_t>& nodes_in_elem, 
                                 CArray<size_t>& vtk_nodes,
                                 const size_t elem_gid, 
                                 const size_t order);

inline int PointIndexFromIJK(int i, int j, int k, const int* order);

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
    const size_t elem_order = 3;  
    const size_t num_elems_1D = 4;

    const size_t max_cycles = 100000;
    const double max_time   = 0.5;
    const double graphics_dt = 0.1;

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
    Kokkos::fence();
    Mesh.nodes_in_elem.update_host();


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

    const size_t num_corners = Mesh.num_corners;

    if(num_nodes_in_elem != FERefElem.num_dofs_in_elem) Kokkos::abort("ERROR: mismatch in DOFs and num nodes in elem \n");

    DCArrayKokkos<double> elem_jac(num_elems, num_qpts_in_elem, elem_dims, elem_dims, "elem_jacobian");
    DCArrayKokkos<double> elem_inv_jac(num_elems, num_qpts_in_elem, elem_dims, elem_dims, "elem_jacobian");
    DCArrayKokkos<double> elem_det_jac(num_elems, num_qpts_in_elem, "elem_det_jacobian");
    DCArrayKokkos<double> elem_vol(num_elems, "elem_vol");    
    
    DCArrayKokkos<double> surf_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_jacobian");
    DCArrayKokkos<double> surf_inv_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_inv_jacobian");
    DCArrayKokkos<double> surf_flux(num_surfs, "surf_flux");
    
    DCArrayKokkos<double> elem_field(num_elems, "elem_field");
    DCArrayKokkos<double> node_field(num_nodes, "node_field");     // for displaying field results
    DCArrayKokkos<double> node_velocity(num_nodes, elem_dims, "node_velocity");
    DCArrayKokkos<double> node_velocity_n(num_nodes, elem_dims, "node_velocity_n");
    DCArrayKokkos<double> node_coords_n(num_nodes, elem_dims, "node_coords_n");

    DCArrayKokkos<double> corner_field(num_corners, "corner_field");
    DCArrayKokkos<double> corner_field_n(num_corners, "corner_field_n");
    DCArrayKokkos<double> elem_vol_matrix(num_elems, num_nodes_in_elem, num_nodes_in_elem, "elem_vol_matrix");
    DCArrayKokkos<double> elem_inv_vol_matrix(num_elems, num_nodes_in_elem, num_nodes_in_elem, "elem_inv_vol_matrix");
    DCArrayKokkos<double> elem_vol_matrix_n(num_elems, num_nodes_in_elem, num_nodes_in_elem, "elem_vol_matrix_n");

    CArrayKokkos <size_t> perm_elem(num_elems,num_nodes_in_elem, "perm");  // permutation array
    CArrayKokkos <double> vv_elem(num_elems,num_nodes_in_elem, "vv");      // temp arrary for LU solver
    
    // Calculate RHS_surf_flux
    CArrayKokkos <double> RHS_surf_flux(num_elems, num_surfs_in_elem, num_qpts_in_surf, "RHS_surf_flux"); // used to build RHS vector 
    CArrayKokkos <double> RHS_elem(num_elems, num_nodes_in_elem, "RHS_elem"); // RHS vector 


    // ================================================================
    // Build qpt to qpt connectivity on surfaces of elements
    
    CArrayKokkos<int> surf_qpt_qpt_map(num_surfs,2,num_qpts_in_surf);
    surf_qpt_qpt_map.set_values(-1);

    build_quadrature_point_connectivity(Mesh,
                                        RefSurf,
                                        surf_qpt_qpt_map,
                                        node_coords);


    // ================================================================
    // Step 1: build the volume matrix for nodal DG at t=0

    elem_vol.set_values(0.0);
    elem_vol_matrix.set_values(0.0);
    FOR_ALL(elem_gid, 0, num_elems, {

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        // Vol_matrix = \int (phi_q \phi_p j w)
        for(size_t dof_lid=0;  dof_lid<num_nodes_in_elem;  dof_lid++)
        for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++)
        for(size_t qpt_lid=0;  qpt_lid<num_qpts_in_elem;   qpt_lid++){

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
                
                // volume contribution from qpt
                const double vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

                elem_vol(elem_gid) += vol_qpt;
                
                elem_vol_matrix(elem_gid, dof_lid, node_lid) += a_basis(dof_lid)*a_basis(node_lid)*vol_qpt;
        } // end for 

    }); // end parallel for
    Kokkos::fence();


    // -----------------------------------------------------
    const double max_vel = 1.0;             // the CFL velocity used for calculating dt
    double h_cfl = h/(double)num_nodes_1D;  // the CFL length scale for calculating dt
    double dt = 0.2*h_cfl/max_vel;          // dt from CFL at start, this time is psuedo time


    // -----------------------------------------------------
    double time = 0;                    // the time 
    double time_output = graphics_dt;   // the time for graphics outputs
    size_t output_id = 0;               // the file id for the outputs
    size_t rk_num_stages = 1;           // runge kutta time integration levels


    // ================================================================
    // Step 0: Set the initial conditions

    FOR_ALL(elem_gid, 0, num_elems, {


        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        for(size_t corner_lid=0; corner_lid<num_nodes_in_elem; corner_lid++){
            // remember: corner_lid = node_lid
            const size_t node_gid = nodes_in_elem(corner_lid);
            const size_t dim = 0;  // x-coord
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid,corner_lid);

            corner_field(corner_gid) = sin(PI*node_coords(node_gid,dim));
        }

    });  // end parallel for

    FOR_ALL(node_gid, 0, num_nodes,{
        // new velocity, it is Taylor-Green vortex
        // PI is defined in mesh class
        node_velocity(node_gid, 0) =  sin(PI*node_coords(node_gid, 0))*cos(PI*node_coords(node_gid, 1));
        node_velocity(node_gid, 1) = -cos(PI*node_coords(node_gid, 0))*sin(PI*node_coords(node_gid, 1));
        node_velocity(node_gid, 2) = 0.0;
    });
    Kokkos::fence();


    // Conservation Check
    double sum_elem = 0.0;
    double domain_mass_t0 = 0.0;
    FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_elem, {

        for(size_t dof_lid=0;  dof_lid<num_nodes_in_elem;  dof_lid++)
        for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            sum_elem += elem_vol_matrix(elem_gid, dof_lid, node_lid)*corner_field(corner_gid);
        }
    }, domain_mass_t0);

    printf("Domain Mass t=0: %f \n", domain_mass_t0);


    // export results to Paraview graphics file
    {

        elem_field.set_values(0.0);
        FOR_ALL(elem_gid,0,num_elems,{
            for(size_t node_lid=0; node_lid<Mesh.num_nodes_in_elem; node_lid++){
                const size_t corner_lid = node_lid;
                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, corner_lid);
                elem_field(elem_gid) += corner_field(corner_gid);
            } 
            elem_field(elem_gid) /= (double)Mesh.num_nodes_in_elem;
        });

        // save corner field to the nodes for graphics outputs
        node_field.set_values(0.0);
        FOR_ALL(node_gid,0,num_nodes,{
            for(size_t corner_lid=0; corner_lid<Mesh.num_corners_in_node(node_gid); corner_lid++){
                const size_t corner_gid = Mesh.corners_in_node(node_gid, corner_lid);
                node_field(node_gid) += corner_field(corner_gid);
            } 
            node_field(node_gid) /= (double)Mesh.num_corners_in_node(node_gid);
        });

        // writing the initial mesh and state
        node_coords.update_host();
        node_field.update_host();
        elem_vol.update_host();

        printf(" Writing output at time = %.4f. ", time);

        char filename[100];
        snprintf(filename, sizeof(filename), "output_time_%04zu.vtu", output_id);

        // Write the mesh state
        write_lagrange_hex_mesh(
            filename,
            node_coords,           
            Mesh.num_nodes,
            Mesh.nodes_in_elem,    
            Mesh.num_elems,
            elem_order,            
            node_field,       
            "Node_Field",
            elem_field,          // element center data
            "Elem_Field"         // element data name                  
        );
        output_id += 1;
    } // end graphics dump scope


    // --------------------------------------------------
    // Time integration loop
    for(size_t cycle=0; cycle<max_cycles; cycle++){
        
        if(cycle%10 == 0) printf(" time = %.4f \n", time);


        // --------------------------------------------------
        // Step 1a: Store time level n state

        FOR_ALL(elem_gid, 0, num_elems, {
            for(size_t dof_lid=0;  dof_lid<num_nodes_in_elem;  dof_lid++)
            for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                elem_vol_matrix_n(elem_gid, dof_lid, node_lid) = elem_vol_matrix(elem_gid, dof_lid, node_lid);
            }
        });

        FOR_ALL(corner_gid, 0, num_corners, {
            corner_field_n(corner_gid) = corner_field(corner_gid);
        });

        FOR_ALL(node_gid, 0, num_nodes, {
            for(size_t dim=0; dim<elem_dims; dim++){
                node_coords_n(node_gid, dim)   = node_coords(node_gid, dim);
                node_velocity_n(node_gid, dim) = node_velocity(node_gid, dim);
            }
        });
        Kokkos::fence();


        // Runge Kutta time integration levels
        for(size_t rk_stage=0; rk_stage<rk_num_stages; rk_stage++){

            // ---- RK coefficient ----
            const double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);


            // ------------------------------------------------------
            // Step 1b: calculate the mesh velocity at time level k

            FOR_ALL(node_gid, 0, num_nodes,{
                // new velocity, it is Taylor-Green vortex
                // PI is defined in mesh class
                node_velocity(node_gid, 0) =  sin(PI*node_coords(node_gid, 0))*cos(PI*node_coords(node_gid, 1));
                node_velocity(node_gid, 1) = -cos(PI*node_coords(node_gid, 0))*sin(PI*node_coords(node_gid, 1));
                node_velocity(node_gid, 2) = 0.0;
            });


            // ------------------------------------------------------
            // Step 2: get CFL time step for moving mesh
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

            dt = 0.1*h_cfl/max_vel; // pseudo time step used for the remap


            // ----------------------------------------------------------
            // Step 3: Calculate the surface fluxes at quadrature points

            RHS_surf_flux.set_values(0.0);
            FOR_ALL(surf_gid, 0, num_surfs, {
                
                const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);

                // get the first elem id and face in this surf
                const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0);
                const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0);

                ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

                for(size_t dof_lid=0; dof_lid<num_nodes_in_elem; dof_lid++)
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


                    double normal_dot_vel = 0.0;
                    for(size_t dim=0; dim<elem_dims; dim++){
                        normal_dot_vel += area_normal[dim]*qpt_vel[dim];
                    }
                    
                    size_t nbr_elem_gid = elem_gid;
                    size_t nbr_face_lid = face_lid;
                    if(num_elems_in_surf==2){
                        nbr_elem_gid = Mesh.elems_in_surf(surf_gid, 1); // second elem
                        nbr_face_lid = Mesh.faces_in_surf(surf_gid, 1); // second elem face
                    }    

                    const size_t nbr_qpt_lid = surf_qpt_qpt_map(surf_gid,0,qpt_lid); // matching qpt

                    ViewCArrayKokkos<double> a_nbr_basis(&RefSurf.qpt_basis(nbr_face_lid,nbr_qpt_lid,0),
                                                         num_nodes_in_elem);

                    // reconstruct the fields
                    double qpt_field     = 0.0;
                    double nbr_qpt_field = 0.0;

                    for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){

                        // Note: corner_lid = node_lid inside the element

                        const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                        qpt_field += a_basis(node_lid)*corner_field(corner_gid);

                        const size_t nbr_corner_gid = Mesh.corners_in_elem(nbr_elem_gid,node_lid);
                        nbr_qpt_field += a_nbr_basis(node_lid)*corner_field(nbr_corner_gid);  
                    } // end for

                    //
                    // Rusanov flux at the quadrature point
                    //
                    
                    // if normal_dot_vel<0 advection is out of first elem in the surf
                    const double flux_val = 0.5*(qpt_field+nbr_qpt_field)*normal_dot_vel 
                                           -0.5*fabs(normal_dot_vel)*(qpt_field-nbr_qpt_field);

                    // save flux value to the quadrature points on either side of the element
                    RHS_surf_flux(elem_gid, face_lid, qpt_lid) = flux_val;
                    if(num_elems_in_surf==2) RHS_surf_flux(nbr_elem_gid, nbr_face_lid, nbr_qpt_lid) = -flux_val;

                } // end for qpt
        
            }); // end surf loop
            Kokkos::fence();

            // -------------------------------------------------
            // Step 4: Build RHS of DG equations in the element

            RHS_elem.set_values(0.0); 
            FOR_ALL(elem_gid, 0, num_elems,{

                ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);
                
                // ----------------------------------------------
                // 4a. Initialize RHS = M*u^n
                
                for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
                    
                    // Add M * u^n term
                    for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                        const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                        RHS_elem(elem_gid, dof_lid) += 
                            elem_vol_matrix_n(elem_gid, dof_lid, node_lid) * corner_field_n(corner_gid);
                    }
                }
                
                // ----------------------------------------------
                // 4b. Subtract the VOLUME integral: \int (\nabal phi_q)J^{-1}*(v*U) dV
                
                for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++)
                for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){
                    
                    ViewCArrayKokkos<double> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0),
                                                          num_nodes_in_elem, 3);
                    ViewCArrayKokkos<double> a_basis(&FERefElem.qpt_basis(qpt_lid,0),
                                                     num_nodes_in_elem);
                    
                    // jacobian matrix and inver
                    ViewCArrayKokkos<double> jac(&elem_jac(elem_gid,qpt_lid,0,0),3,3);
                    ViewCArrayKokkos<double> inv_jac(&elem_inv_jac(elem_gid,qpt_lid,0,0),3,3);

                    jacobian(jac, 
                            node_coords, 
                            nodes_in_elem,
                            a_grad_basis);

                    elem_det_jac(elem_gid,qpt_lid) = det_3x3(jac);

                    invert_3x3(jac, inv_jac, elem_det_jac(elem_gid,qpt_lid));

                    // Reconstruct field at quadrature point
                    double qpt_field = 0.0;
                    for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                        const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                        qpt_field += a_basis(node_lid) * corner_field(corner_gid);
                    }
                    
                    // Reconstruct velocity at quadrature point 
                    double qpt_vel[3];
                    qpt_vel[0] = 0.0;
                    qpt_vel[1] = 0.0; 
                    qpt_vel[2] = 0.0;
                    
                    for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                        const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
                        for(size_t dim = 0; dim < elem_dims; dim++){
                            qpt_vel[dim] += a_basis(node_lid) * node_velocity(node_gid, dim);
                        }
                    }

                    // transform the gradient to the physical space
                    double physical_grad[3]; 
                    physical_grad[0] = 0.0;
                    physical_grad[1] = 0.0; 
                    physical_grad[2] = 0.0;
                    for(size_t i = 0; i < elem_dims; i++)
                    for(size_t j = 0; j < elem_dims; j++){
                        physical_grad[i] += a_grad_basis(dof_lid, j)*inv_jac(j,i);
                    }
                    
                    // Compute (\nabal phi_q)*J^-1*(v*U)
                    // From Anderson et. al. paper
                    double grad_dot_flux = 0.0;
                    for(size_t dim = 0; dim < elem_dims; dim++){
                        grad_dot_flux += physical_grad[dim] * qpt_vel[dim] * qpt_field;
                    }
                    
                    const double vol_qpt = elem_det_jac(elem_gid,qpt_lid) * Quad.qpt_weights(qpt_lid);
                    RHS_elem(elem_gid, dof_lid) -= rk_alpha * dt * grad_dot_flux * vol_qpt;
                } // end for dof_lid and qpt
                
                // ----------------------------------------------
                // 4c. Add SURFACE flux contribution
                
                for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++)
                for(size_t face_lid = 0; face_lid < num_surfs_in_elem; face_lid++)
                for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++){
                    
                    ViewCArrayKokkos<double> a_basis(&RefSurf.qpt_basis(face_lid, qpt_lid, 0),
                                                     num_nodes_in_elem);
                    
                    // surface flux (note: RHS_surf_flux already has correct sign)
                    RHS_elem(elem_gid, dof_lid) += 
                        rk_alpha * dt * RHS_surf_flux(elem_gid, face_lid, qpt_lid) * a_basis(dof_lid);
                }

            });
            Kokkos::fence();



            // ================================================================
            // Step 5: Move the mesh to the new location
            FOR_ALL(node_gid, 0, num_nodes,{
                // new position of the mesh
                node_coords(node_gid, 0) = node_coords_n(node_gid, 0) + 0.5*(node_velocity(node_gid, 0)+node_velocity_n(node_gid, 0)) * rk_alpha * dt; 
                node_coords(node_gid, 1) = node_coords_n(node_gid, 1) + 0.5*(node_velocity(node_gid, 1)+node_velocity_n(node_gid, 1)) * rk_alpha * dt;
                // z-coords never change
            });
            Kokkos::fence();


            // ================================================================
            // Step 6: build the volume matrix for nodal DG after the mesh moved
            elem_vol.set_values(0.0);
            elem_vol_matrix.set_values(0.0);
            FOR_ALL(elem_gid, 0, num_elems, {

                ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

                // Vol_matrix = \int (phi_q \phi_p j w)
                for(size_t dof_lid=0;  dof_lid<num_nodes_in_elem;  dof_lid++)
                for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++)
                for(size_t qpt_lid=0;  qpt_lid<num_qpts_in_elem;   qpt_lid++){

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
                        
                        // volume contribution from qpt
                        const double vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

                        elem_vol(elem_gid) += vol_qpt;
                        
                        elem_vol_matrix(elem_gid, dof_lid, node_lid) += a_basis(dof_lid)*a_basis(node_lid)*vol_qpt;
                } // end for 

            }); // end parallel for
            Kokkos::fence();


            // -----------------------------------------------------
            // 7. Solve M * u^{n+1} = RHS using LU decomposition

            FOR_ALL(elem_gid, 0, num_elems,{
                
                // populate the matrix to invert
                for(size_t dof_lid=0; dof_lid<num_nodes_in_elem; dof_lid++)
                for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                    elem_inv_vol_matrix(elem_gid, dof_lid, node_lid) = elem_vol_matrix(elem_gid, dof_lid, node_lid);
                }

                int singular = 0; 
                int parity = 0;

                ViewCArrayKokkos<size_t> perm(&perm_elem(elem_gid,0), num_nodes_in_elem);
                ViewCArrayKokkos<double> vv(&vv_elem(elem_gid,0), num_nodes_in_elem);

                ViewCArrayKokkos<double> A(&elem_inv_vol_matrix(elem_gid, 0, 0), 
                                        num_nodes_in_elem, num_nodes_in_elem); 
                ViewCArrayKokkos<double> b(&RHS_elem(elem_gid,0), num_nodes_in_elem);

                singular = LU_decompose(A, perm, vv, parity);
                if(singular == 0){
                    Kokkos::abort("ERROR: matrix is singular \n");
                }

                //printf("\nA = \n");
                //for(size_t i=0; i<elem_dims; i++){
                //    for(size_t j=0; j<elem_dims; j++){   
                //        printf(" %f, ", A(i,j));
                //    }
                //    printf("\n");
                //}
                //printf("\nb = \n");
                //for(size_t i=0; i<elem_dims; i++){
                //    printf("   %f\n", b(i));
                //}

                LU_backsub(A, perm, b); // answer sent back in b

                //printf("\nx = \n");
                //for(size_t i=0; i<elem_dims; i++){
                //    printf("   %f\n", b(i));
                //}

                // -----------------------------------------------------
                // 4e. Save the new corner DOFs
                for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
                    const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
                    corner_field(corner_gid) = b(dof_lid);
                }

            }); // end parallel for elems
            Kokkos::fence();

        } // end Runge Kutta time level loop


        // ================================================================
        // Step 7: update time
        time += dt;

        // Conservation Check
        double sum_elem = 0.0;
        double domain_mass_time = 0.0;
        FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_elem, {

            for(size_t dof_lid=0;  dof_lid<num_nodes_in_elem;  dof_lid++)
            for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                sum_elem += elem_vol_matrix(elem_gid, dof_lid, node_lid)*corner_field(corner_gid);
            }
        }, domain_mass_time);

        printf("Domain mass error= %f \n", domain_mass_time-domain_mass_t0);
        if(fabs(domain_mass_time-domain_mass_t0)>1.e-12) Kokkos::abort("ERROR: Mass is not conserved");


        // ================================================================
        // Step 8: write outputs
        if( time-time_output >= -1.e-8 ){

            elem_field.set_values(0.0);
            FOR_ALL(elem_gid,0,num_elems,{
                for(size_t node_lid=0; node_lid<Mesh.num_nodes_in_elem; node_lid++){
                    const size_t corner_lid = node_lid;
                    const size_t corner_gid = Mesh.corners_in_elem(elem_gid, corner_lid);
                    elem_field(elem_gid) += corner_field(corner_gid);
                } 
                elem_field(elem_gid) /= (double)Mesh.num_nodes_in_elem;
            });

            // save corner field to the nodes for graphics outputs
            node_field.set_values(0.0);
            FOR_ALL(node_gid,0,num_nodes,{
                for(size_t corner_lid=0; corner_lid<Mesh.num_corners_in_node(node_gid); corner_lid++){
                    const size_t corner_gid = Mesh.corners_in_node(node_gid, corner_lid);
                    node_field(node_gid) += corner_field(corner_gid);
                } 
                node_field(node_gid) /= (double)Mesh.num_corners_in_node(node_gid);
            });

            // writing the initial mesh and state
            node_coords.update_host();
            node_field.update_host();
            elem_vol.update_host();

            printf(" Writing output at time = %.4f. ", time);

            char filename[100];
            snprintf(filename, sizeof(filename), "output_time_%04zu.vtu", output_id);

            // Write the mesh state
            write_lagrange_hex_mesh(
                filename,
                node_coords,           
                Mesh.num_nodes,
                Mesh.nodes_in_elem,    
                Mesh.num_elems,
                elem_order,            
                node_field,       
                "Node_Field",
                elem_field,          // element center data
                "Elem_Field"           // element data name                  
            );
            time_output += graphics_dt;
            output_id += 1;

        } // end if

        if (time >= max_time  ){
            printf("Domain mass at time=%f: %f \n", time, domain_mass_time);
            break;
        }

    } // end loop over cycle
    printf(" time = %.4f ", time);

    printf("\n Remap test finished.\n");


} // end MATAR scope
MATAR_FINALIZE();

return 0;
} // end function



// 
void write_lagrange_hex_mesh(
    const std::string& filename,
    const DCArrayKokkos<double>& node_coords,       // All node coordinates [num_nodes][3]
    const size_t num_nodes,
    const DCArrayKokkos<size_t>& nodes_in_elem,     // Connectivity
    const size_t num_elems,
    const size_t order,
    const DCArrayKokkos<double>& node_data,         // Nodal data
    const std::string& node_data_name,
    const DCArrayKokkos<double>& elem_data,         // Element center data (NEW)
    const std::string& elem_data_name)              // Element data name (NEW)
{
    std::ofstream vtu_file(filename);
    if (!vtu_file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    vtu_file << std::fixed << std::setprecision(8);

    // Header
    vtu_file << "<?xml version=\"1.0\"?>\n";
    vtu_file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    vtu_file << "  <UnstructuredGrid>\n";
    vtu_file << "    <Piece NumberOfPoints=\"" << num_nodes 
             << "\" NumberOfCells=\"" << num_elems << "\">\n";

    // Write Points
    write_points(vtu_file, node_coords, num_nodes);

    // Write Cells (connectivity, types, AND cell data)
    write_lagrange_cells(vtu_file, nodes_in_elem, num_elems, order, 
                        elem_data, elem_data_name);  // Pass element data

    // Write Point Data
    write_point_data(vtu_file, node_data, num_nodes, node_data_name);

    // Footer
    vtu_file << "    </Piece>\n";
    vtu_file << "  </UnstructuredGrid>\n";
    vtu_file << "</VTKFile>\n";

    vtu_file.close();
    std::cout << "Wrote VTU file: " << filename << std::endl;
}

void write_points(std::ofstream& file, const DCArrayKokkos<double>& coords, size_t num_nodes)
{
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    for (size_t i = 0; i < num_nodes; i++) {
        file << "          " << coords.host(i, 0) << " " 
             << coords.host(i, 1) << " " 
             << coords.host(i, 2) << "\n";
    }
    
    file << "        </DataArray>\n";
    file << "      </Points>\n";
}

void write_lagrange_cells(std::ofstream& file, 
                          const DCArrayKokkos<size_t>& nodes_in_elem,
                          size_t num_elems, 
                          size_t order,
                          const DCArrayKokkos<double>& elem_data,    // Element data
                          const std::string& elem_data_name)         // Element data name
{
    const size_t nodes_per_elem = (order + 1) * (order + 1) * (order + 1);
    const int VTK_LAGRANGE_HEXAHEDRON = 72;

    file << "      <Cells>\n";
    
    // Connectivity
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    
    CArray<size_t> vtk_nodes(nodes_per_elem);
    
    for (size_t elem = 0; elem < num_elems; elem++) {
        // Convert to VTK ordering
        reorder_ijk_to_vtk_lagrange(nodes_in_elem, vtk_nodes, elem, order);
        
        file << "          ";
        for (size_t i = 0; i < nodes_per_elem; i++) {
            file << vtk_nodes(i) << " ";
        }
        file << "\n";
    }
    
    file << "        </DataArray>\n";

    // Offsets
    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    file << "          ";
    for (size_t elem = 0; elem < num_elems; elem++) {
        file << (elem + 1) * nodes_per_elem << " ";
    }
    file << "\n        </DataArray>\n";

    // Cell types
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    file << "          ";
    for (size_t elem = 0; elem < num_elems; elem++) {
        file << VTK_LAGRANGE_HEXAHEDRON << " ";
    }
    file << "\n        </DataArray>\n";

    file << "      </Cells>\n";

    // CellData section with HigherOrderDegrees AND user data
    file << "      <CellData Scalars=\"" << elem_data_name << "\">\n";
    
    // HigherOrderDegrees (CRITICAL for Lagrange elements!)
    file << "        <DataArray type=\"Int32\" Name=\"HigherOrderDegrees\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    file << "          ";
    for (size_t elem = 0; elem < num_elems; elem++) {
        file << order << " " << order << " " << order << " ";
    }
    file << "\n        </DataArray>\n";
    
    // User-provided element center data
    file << "        <DataArray type=\"Float64\" Name=\"" << elem_data_name << "\" format=\"ascii\">\n";
    file << "          ";
    for (size_t elem = 0; elem < num_elems; elem++) {
        file << elem_data.host(elem) << " ";
    }
    file << "\n        </DataArray>\n";
    
    file << "      </CellData>\n";
}

void write_point_data(std::ofstream& file, 
                      const DCArrayKokkos<double>& data, 
                      size_t num_nodes,
                      const std::string& name)
{
    file << "      <PointData Scalars=\"" << name << "\">\n";
    file << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
    
    // writing node field data
    for (size_t i = 0; i < num_nodes; i++) {
        file << "          " << data.host(i) << "\n";
    }
    
    file << "        </DataArray>\n";
    file << "      </PointData>\n";
}

// Keep your existing helper functions unchanged
void reorder_ijk_to_vtk_lagrange(const DCArrayKokkos<size_t>& nodes_in_elem, 
                                 CArray<size_t>& vtk_nodes,
                                 const size_t elem_gid, 
                                 const size_t order)
{
    const int n = order + 1;
    int ord[3] = {(int)order, (int)order, (int)order};
    
    std::vector<std::pair<int, size_t>> vtk_to_ijk;
    
    for(int k = 0; k < n; k++){
        for(int j = 0; j < n; j++){
            for(int i = 0; i < n; i++){
                int vtk_pos = PointIndexFromIJK(i, j, k, ord);
                size_t ijk_linear = i + j*n + k*n*n;
                vtk_to_ijk.push_back({vtk_pos, ijk_linear});
            }
        }
    }
    
    std::sort(vtk_to_ijk.begin(), vtk_to_ijk.end());
    
    for(size_t v = 0; v < vtk_to_ijk.size(); v++){
        size_t ijk_linear = vtk_to_ijk[v].second;
        vtk_nodes(v) = nodes_in_elem.host(elem_gid, ijk_linear);
    }
}

inline int PointIndexFromIJK(int i, int j, int k, const int* order)
{
    bool ibdy = (i == 0 || i == order[0]);
    bool jbdy = (j == 0 || j == order[1]);
    bool kbdy = (k == 0 || k == order[2]);
    int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

    if (nbdy == 3) { // Vertex DOF
        return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    }

    int offset = 8;
    if (nbdy == 2) { // Edge DOF
        if (!ibdy) {
            return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + 
                   (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) + offset;
        }
        if (!jbdy) {
            return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + 
                   (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) + offset;
        }
        offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
        return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
    }

    offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
    if (nbdy == 1) { // Face DOF
        if (ibdy) {
            return (j - 1) + ((order[1] - 1) * (k - 1)) + 
                   (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
        }
        offset += 2 * (order[1] - 1) * (order[2] - 1);
        if (jbdy) {
            return (i - 1) + ((order[0] - 1) * (k - 1)) + 
                   (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
        }
        offset += 2 * (order[2] - 1) * (order[0] - 1);
        return (i - 1) + ((order[0] - 1) * (j - 1)) + 
               (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
    }

    // Interior DOF
    offset += 2 * ((order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) + 
                   (order[0] - 1) * (order[1] - 1));
    return offset + (i - 1) + (order[0] - 1) * ((j - 1) + (order[1] - 1) * (k - 1));
}