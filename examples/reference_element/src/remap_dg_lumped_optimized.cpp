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
//
// GPU-tuned lumped mass discontinuous Galerkin remap.
//
// This solves the same problem as remap_dg_lumped_test.cpp and reproduces its
// results; it is organized for the GPU:
//
//   - kernels are launched over flat 1D index spaces rather than a team policy
//   - the Jacobian, its determinate and its inverse are built once per volume
//     quadrature point in their own pass, instead of once per DOF
//   - the field and velocity reconstruction is hoisted out of the DOF loop, so
//     each quadrature point is reconstructed once per element
//   - the reference element tables are replicated in the layouts the kernels
//     read, so a warp walks contiguous doubles
//   - the surface Jacobian never leaves registers, and the L1 and L2 error
//     norms are gathered in a single pass over the elements
//   - consecutive kernels share the default execution space, so they are
//     already stream ordered and the host only waits where it reads a result
//
// Optional arguments:
//   num_elems_x num_elems_y [num_elems_z [max_time [graphics_dt]]]
//

#include <cstdlib>
#include <cstdint>
#include <cstring>
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

using namespace mtr;
using namespace swage;    // unstructured mesh and point cloud
using namespace elements; // reference element space


using REAL_t = double; // REAL_t precision


#define USE_NOTCHED_CIRCLE
//#define USE_SIN_FUNCTION
// #define USE_GAUSSIAN

KOKKOS_INLINE_FUNCTION
REAL_t test_function(const REAL_t x, 
                     const REAL_t y);


void write_lagrange_hex_mesh(
    const std::string& filename,
    const DCArrayKokkos<REAL_t>& node_coords,       // All node coordinates [num_nodes][3]
    const size_t num_nodes,
    const DCArrayKokkos<size_t>& nodes_in_elem,     // Connectivity
    const size_t num_elems,
    const size_t order,
    const DCArrayKokkos<REAL_t>& node_data,         // Nodal data
    const std::string& node_data_name,
    const DCArrayKokkos<REAL_t>& elem_data,         // Element center data
    const std::string& elem_data_name);             // Element data name


void write_lagrange_cells(std::ofstream& file, 
                        const DCArrayKokkos<size_t>& nodes_in_elem,
                        size_t num_elems, 
                        size_t order,
                        const DCArrayKokkos<REAL_t>& elem_data,
                        const std::string& elem_data_name);

void write_points(std::ofstream& file, 
                  const DCArrayKokkos<REAL_t>& coords, 
                  size_t num_nodes);


void write_point_data(std::ofstream& file, 
                      const DCArrayKokkos<REAL_t>& data, 
                      size_t num_nodes,
                      const std::string& name);

void reorder_ijk_to_vtk_lagrange(const DCArrayKokkos<size_t>& nodes_in_elem, 
                                 CArray<size_t>& vtk_nodes,
                                 const size_t elem_gid, 
                                 const size_t order);

inline int PointIndexFromIJK(int i, int j, int k, const int* order);


REAL_t lagrange_basis(const REAL_t xi, const size_t i, const CArrayKokkos<REAL_t>& nodes);
void interpolate_to_uniform(const DCArrayKokkos<size_t>& nodes_in_elem,
                            const CArrayKokkos<REAL_t>& lob_nodes_1D,
                            const DCArrayKokkos<REAL_t>& node_coords_lob,  // Lobatto node positions
                            DCArrayKokkos<REAL_t>& node_coords_uniform,   // Output uniform positions 
                            const size_t num_elems); 


// ============================================================================
// Reference-element tables replicated in the layouts the GPU kernels want.
// Kernels threaded over quadrature points need qpt as the fastest index;
// kernels threaded over DOFs need dof as the fastest index.
// ============================================================================
struct BasisTables_t
{
    CArrayKokkos<REAL_t> basis_dq;       // (dof, qpt)
    CArrayKokkos<REAL_t> grad_basis_njq; // (dof, dim, qpt)
    CArrayKokkos<REAL_t> grad_basis_qjd; // (qpt, dim, dof)
    CArrayKokkos<REAL_t> basis_row_sum;  // (qpt) = sum over DOFs of the basis

    CArrayKokkos<REAL_t> surf_basis_fdq; // (face, dof, surf qpt)
    CArrayKokkos<REAL_t> surf_grad_fjdq; // (face, dim, dof, surf qpt)
};

// ============================================================================
// Geometry of the moving mesh at every volume quadrature point.
//
// The Jacobian is accumulated in registers and only its determinate and its
// inverse are written out.  The inverse is stored with the quadrature point as
// the fastest index so the kernels that read it back are coalesced.
// ============================================================================
static void build_element_geometry(const Mesh_t& Mesh,
                                   const BasisTables_t& tables,
                                   const DCArrayKokkos<REAL_t>& node_coords,
                                   const CArrayKokkos<REAL_t>& elem_det_jac,
                                   const CArrayKokkos<REAL_t>& inv_jac_ijq,
                                   const size_t num_elems,
                                   const size_t num_qpts_in_elem,
                                   const size_t num_nodes_in_elem)
{
    FOR_ALL(idx, 0, num_elems*num_qpts_in_elem, {

        const size_t elem_gid = idx / num_qpts_in_elem;
        const size_t qpt_lid  = idx % num_qpts_in_elem;

        REAL_t j00 = 0.0; REAL_t j01 = 0.0; REAL_t j02 = 0.0;
        REAL_t j10 = 0.0; REAL_t j11 = 0.0; REAL_t j12 = 0.0;
        REAL_t j20 = 0.0; REAL_t j21 = 0.0; REAL_t j22 = 0.0;

        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, dof_lid);
            const REAL_t x0 = node_coords(node_gid, 0);
            const REAL_t x1 = node_coords(node_gid, 1);
            const REAL_t x2 = node_coords(node_gid, 2);
            const REAL_t g0 = tables.grad_basis_njq(dof_lid, 0, qpt_lid);
            const REAL_t g1 = tables.grad_basis_njq(dof_lid, 1, qpt_lid);
            const REAL_t g2 = tables.grad_basis_njq(dof_lid, 2, qpt_lid);
            j00 += x0*g0; j01 += x0*g1; j02 += x0*g2;
            j10 += x1*g0; j11 += x1*g1; j12 += x1*g2;
            j20 += x2*g0; j21 += x2*g1; j22 += x2*g2;
        }

        const REAL_t det = det_3x3(j00, j01, j02, j10, j11, j12, j20, j21, j22);
        elem_det_jac(elem_gid, qpt_lid) = det;

        REAL_t i00; REAL_t i01; REAL_t i02;
        REAL_t i10; REAL_t i11; REAL_t i12;
        REAL_t i20; REAL_t i21; REAL_t i22;
        invert_3x3(det, j00, j01, j02, j10, j11, j12, j20, j21, j22,
                   i00, i01, i02, i10, i11, i12, i20, i21, i22);

        inv_jac_ijq(elem_gid, 0, 0, qpt_lid) = i00;
        inv_jac_ijq(elem_gid, 0, 1, qpt_lid) = i01;
        inv_jac_ijq(elem_gid, 0, 2, qpt_lid) = i02;
        inv_jac_ijq(elem_gid, 1, 0, qpt_lid) = i10;
        inv_jac_ijq(elem_gid, 1, 1, qpt_lid) = i11;
        inv_jac_ijq(elem_gid, 1, 2, qpt_lid) = i12;
        inv_jac_ijq(elem_gid, 2, 0, qpt_lid) = i20;
        inv_jac_ijq(elem_gid, 2, 1, qpt_lid) = i21;
        inv_jac_ijq(elem_gid, 2, 2, qpt_lid) = i22;
    });
} // end build_element_geometry


// ============================================================================
// Row-lumped volume (mass) vector, elem_corner_vol(elem, node).
// ============================================================================
static void build_lumped_volume(const ReferenceElement_t& FERefElem,
                                const Quadrature_t& Quad,
                                const BasisTables_t& tables,
                                const CArrayKokkos<REAL_t>& elem_det_jac,
                                DCArrayKokkos<REAL_t>& elem_corner_vol,
                                const size_t num_elems,
                                const size_t num_qpts_in_elem,
                                const size_t num_nodes_in_elem)
{
    // The inner DOF loop of the original only contributes the row sum of the
    // basis, which is the same at every element, so it is tabulated once.
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t node_lid = idx % num_nodes_in_elem;

        REAL_t vol = 0.0;
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){
            const REAL_t vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
            vol += tables.basis_row_sum(qpt_lid)*FERefElem.qpt_basis(qpt_lid, node_lid)*vol_qpt;
        }
        elem_corner_vol(elem_gid, node_lid) = vol;
    });
} // end build_lumped_volume


// ============================================================================
// Rusanov flux at the surface quadrature points.
// ============================================================================
static void build_surface_flux(const Mesh_t& Mesh,
                               const ReferenceSurface_t& RefSurf,
                               const SurfaceQuadrature_t& SurfQuad,
                               const BasisTables_t& tables,
                               const DCArrayKokkos<REAL_t>& node_coords,
                               const DCArrayKokkos<REAL_t>& node_velocity,
                               const DCArrayKokkos<REAL_t>& corner_field,
                               const CArrayKokkos<int>& surf_qpt_qpt_map,
                               CArrayKokkos<REAL_t>& RHS_surf_flux,
                               const size_t num_surfs,
                               const size_t num_qpts_in_surf,
                               const size_t num_nodes_in_elem)
{
    RHS_surf_flux.set_values(0.0);

    // The reference tables are indexed with the surface quadrature point as the
    // fastest axis so that a warp reads contiguous doubles, the surface
    // Jacobian never leaves registers, and the geometry, the velocity and the
    // field are all reconstructed in the same sweep over the DOFs.
    FOR_ALL(idx, 0, num_surfs*num_qpts_in_surf, {

        const size_t surf_gid = idx / num_qpts_in_surf;
        const size_t qpt_lid  = idx % num_qpts_in_surf;

        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);
        const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0);
        const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0);

        REAL_t j00 = 0.0; REAL_t j01 = 0.0; REAL_t j02 = 0.0;
        REAL_t j10 = 0.0; REAL_t j11 = 0.0; REAL_t j12 = 0.0;
        REAL_t j20 = 0.0; REAL_t j21 = 0.0; REAL_t j22 = 0.0;

        REAL_t qpt_vel0 = 0.0;
        REAL_t qpt_vel1 = 0.0;
        REAL_t qpt_vel2 = 0.0;
        REAL_t qpt_field = 0.0;

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
            const REAL_t x0 = node_coords(node_gid, 0);
            const REAL_t x1 = node_coords(node_gid, 1);
            const REAL_t x2 = node_coords(node_gid, 2);
            const REAL_t g0 = tables.surf_grad_fjdq(face_lid, 0, node_lid, qpt_lid);
            const REAL_t g1 = tables.surf_grad_fjdq(face_lid, 1, node_lid, qpt_lid);
            const REAL_t g2 = tables.surf_grad_fjdq(face_lid, 2, node_lid, qpt_lid);
            j00 += x0*g0; j01 += x0*g1; j02 += x0*g2;
            j10 += x1*g0; j11 += x1*g1; j12 += x1*g2;
            j20 += x2*g0; j21 += x2*g1; j22 += x2*g2;

            const REAL_t phi = tables.surf_basis_fdq(face_lid, node_lid, qpt_lid);
            qpt_vel0 += phi*node_velocity(node_gid, 0);
            qpt_vel1 += phi*node_velocity(node_gid, 1);
            qpt_vel2 += phi*node_velocity(node_gid, 2);

            qpt_field += phi*corner_field(Mesh.corners_in_elem(elem_gid, node_lid));
        }

        const REAL_t det_jac_qpt = det_3x3(j00, j01, j02, j10, j11, j12, j20, j21, j22);

        REAL_t i00; REAL_t i01; REAL_t i02;
        REAL_t i10; REAL_t i11; REAL_t i12;
        REAL_t i20; REAL_t i21; REAL_t i22;
        invert_3x3(det_jac_qpt, j00, j01, j02, j10, j11, j12, j20, j21, j22,
                   i00, i01, i02, i10, i11, i12, i20, i21, i22);

        const REAL_t scale = det_jac_qpt*SurfQuad.qpt_weights(face_lid, qpt_lid);

        REAL_t area_normal0; REAL_t area_normal1; REAL_t area_normal2;
        nanson_area_normal(RefSurf.outward_normal(face_lid, 0),
                           RefSurf.outward_normal(face_lid, 1),
                           RefSurf.outward_normal(face_lid, 2),
                           scale,
                           i00, i01, i02, i10, i11, i12, i20, i21, i22,
                           area_normal0, area_normal1, area_normal2);

        const REAL_t normal_dot_vel = area_normal0*qpt_vel0
                                    + area_normal1*qpt_vel1
                                    + area_normal2*qpt_vel2;

        size_t nbr_elem_gid = elem_gid;
        size_t nbr_face_lid = face_lid;
        if(num_elems_in_surf == 2){
            nbr_elem_gid = Mesh.elems_in_surf(surf_gid, 1);
            nbr_face_lid = Mesh.faces_in_surf(surf_gid, 1);
        }

        const size_t nbr_qpt_lid = surf_qpt_qpt_map(surf_gid, 0, qpt_lid);

        REAL_t nbr_qpt_field = 0.0;
        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            nbr_qpt_field += tables.surf_basis_fdq(nbr_face_lid, node_lid, nbr_qpt_lid)
                            *corner_field(Mesh.corners_in_elem(nbr_elem_gid, node_lid));
        }

        const REAL_t flux_val = rusanov_flux(qpt_field, nbr_qpt_field, normal_dot_vel);

        RHS_surf_flux(elem_gid, face_lid, qpt_lid) = flux_val;
        if(num_elems_in_surf == 2) RHS_surf_flux(nbr_elem_gid, nbr_face_lid, nbr_qpt_lid) = -flux_val;
    });
} // end build_surface_flux


// ============================================================================
// RHS of the DG equations.
// ============================================================================
static void assemble_rhs(const Mesh_t& Mesh,
                         const ReferenceSurface_t& RefSurf,
                         const Quadrature_t& Quad,
                         const BasisTables_t& tables,
                         const CArrayKokkos<REAL_t>& elem_det_jac,
                         const CArrayKokkos<REAL_t>& inv_jac_ijq,
                         const DCArrayKokkos<REAL_t>& corner_field,
                         const DCArrayKokkos<REAL_t>& corner_field_n,
                         const DCArrayKokkos<REAL_t>& node_velocity,
                         const DCArrayKokkos<REAL_t>& elem_corner_vol_n,
                         const CArrayKokkos<REAL_t>& RHS_surf_flux,
                         CArrayKokkos<REAL_t>& RHS_elem,
                         const CArrayKokkos<REAL_t>& qpt_vol_flux,
                         const REAL_t rk_alpha,
                         const REAL_t dt,
                         const size_t num_elems,
                         const size_t num_qpts_in_elem,
                         const size_t num_nodes_in_elem,
                         const size_t num_surfs_in_elem,
                         const size_t num_qpts_in_surf,
                         const size_t elem_dims)
{
    // Pass 1: the reconstruction at a quadrature point is shared by all DOFs
    // of the element, so it is formed once instead of num_dofs times.
    FOR_ALL(idx, 0, num_elems*num_qpts_in_elem, {

        const size_t elem_gid = idx / num_qpts_in_elem;
        const size_t qpt_lid  = idx % num_qpts_in_elem;

        REAL_t qpt_field = 0.0;
        REAL_t qpt_vel_0 = 0.0;
        REAL_t qpt_vel_1 = 0.0;
        REAL_t qpt_vel_2 = 0.0;

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            const REAL_t basis_val = tables.basis_dq(node_lid, qpt_lid);
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            qpt_field += basis_val*corner_field(corner_gid);
        }

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            const REAL_t basis_val = tables.basis_dq(node_lid, qpt_lid);
            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
            qpt_vel_0 += basis_val*node_velocity(node_gid, 0);
            qpt_vel_1 += basis_val*node_velocity(node_gid, 1);
            qpt_vel_2 += basis_val*node_velocity(node_gid, 2);
        }

        // grad(phi).J^-1.(v U) = sum_j dphi/dxi_j * [ sum_i Jinv(j,i) v_i U ],
        // so the DOF-independent bracket is tabulated here.
        const REAL_t vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
        for(size_t j = 0; j < elem_dims; j++){
            const REAL_t flux = inv_jac_ijq(elem_gid, j, 0, qpt_lid)*qpt_vel_0
                              + inv_jac_ijq(elem_gid, j, 1, qpt_lid)*qpt_vel_1
                              + inv_jac_ijq(elem_gid, j, 2, qpt_lid)*qpt_vel_2;
            qpt_vol_flux(elem_gid, qpt_lid, j) = flux*qpt_field*vol_qpt;
        }
    });
    // Pass 2 reads what pass 1 wrote into qpt_vol_flux.
    Kokkos::fence();

    // Pass 2: one thread per (element, DOF)
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t dof_lid  = idx % num_nodes_in_elem;

        // 4a. the M*u^n term; remember node_lid = dof_lid = corner_lid
        const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
        REAL_t rhs = elem_corner_vol_n(elem_gid, dof_lid)*corner_field_n(corner_gid);

        // 4b. subtract the VOLUME integral: \int (\nabla phi_q) J^{-1} (v U) dV
        REAL_t vol_integral = 0.0;
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){
            vol_integral += tables.grad_basis_qjd(qpt_lid, 0, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 0)
                          + tables.grad_basis_qjd(qpt_lid, 1, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 1)
                          + tables.grad_basis_qjd(qpt_lid, 2, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 2);
        }
        rhs -= rk_alpha*dt*vol_integral;

        // 4c. add the SURFACE flux contribution
        REAL_t surf_integral = 0.0;
        for(size_t face_lid = 0; face_lid < num_surfs_in_elem; face_lid++)
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++){
            surf_integral += RHS_surf_flux(elem_gid, face_lid, qpt_lid)
                           * RefSurf.qpt_basis(face_lid, qpt_lid, dof_lid);
        }
        rhs += rk_alpha*dt*surf_integral;

        RHS_elem(elem_gid, dof_lid) = rhs;
    });
} // end assemble_rhs


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
    
    const REAL_t L_x = 1.;
    const REAL_t L_y = 1.;
    const REAL_t L_z = 0.0625;  // 0.5, 0.25, 0.125, 0.0625
    // The defaults reproduce the reference case; the optional arguments only
    // exist so the problem can be swept without recompiling.
    size_t num_elems_x = 10;
    size_t num_elems_y = 10;
    size_t num_elems_z = 1;
    REAL_t max_time    = 0.2;
    REAL_t graphics_dt = 0.01;

    if (argc >= 3) {
        num_elems_x = (size_t)std::strtod(argv[1], nullptr);
        num_elems_y = (size_t)std::strtod(argv[2], nullptr);
    }
    if (argc >= 4) num_elems_z = (size_t)std::strtod(argv[3], nullptr);
    if (argc >= 5) max_time    = std::strtod(argv[4], nullptr);
    if (argc >= 6) graphics_dt = std::strtod(argv[5], nullptr);

    const size_t rk_num_stages = 2;    // number of runge kutta time integration levels
    const size_t max_cycles = 10000000;

    printf("mesh %zux%zux%zu, max_time %g, graphics_dt %g\n",
           num_elems_x, num_elems_y, num_elems_z, max_time, graphics_dt);


    // ================================================================
    // Create quadrature along with the reference element and surface

    std::cout<<"Building reference elements and quadrature \n";

    // the minimum quadrature for FE hydrodynamics based on elem order
    const size_t num_DOFs_1d = elem_order + 1; // 
    const size_t num_qpts_1d = 2*elem_order;   // using Legendre, but if using Lobatto, it requires 2*Order+1


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

    const size_t num_elems   = num_elems_x*num_elems_y*num_elems_z;
    const size_t num_nodes_x = elem_order*num_elems_x + 1;  // number of nodes
    const size_t num_nodes_y = elem_order*num_elems_y + 1;  // number of nodes
    const size_t num_nodes_z = elem_order*num_elems_z + 1;  // number of nodes
    const size_t num_nodes   = num_nodes_x*num_nodes_y*num_nodes_z;


    Mesh.initialize_dims(elem_dims);
    Mesh.initialize_elems_Pn(num_elems, elem_order, Quad.num_qpts_1d);
    Mesh.initialize_nodes(num_nodes);
    
    DCArrayKokkos <REAL_t> node_coords(Mesh.num_nodes, Mesh.num_dims);

    // Physical element sizes
    const REAL_t h_x = L_x / (REAL_t)num_elems_x;
    const REAL_t h_y = L_y / (REAL_t)num_elems_y;
    const REAL_t h_z = L_z / (REAL_t)num_elems_z;



    // create indexing for a Pn order mesh

    CArrayKokkos <REAL_t> lob_nodes_1D(num_DOFs_1d);
    RUN({
        get_lobatto_nodes_1D(lob_nodes_1D,num_DOFs_1d); 
    });
    const REAL_t ref_length = 2.0;


    // Step 1: Initialize ALL node coordinates once
    FOR_ALL(kp, 0, num_nodes_z,
            jp, 0, num_nodes_y,
            ip, 0, num_nodes_x, {
        
        size_t node_gid = ip + (jp + kp*num_nodes_y)*num_nodes_x;
        
        size_t elem_x = ip / (num_DOFs_1d - 1); // integer division floors to nearest whole id
        size_t loc_x = ip % (num_DOFs_1d - 1);  // modulo gives the local node id in elem
        
        size_t elem_y = jp / (num_DOFs_1d - 1);
        size_t loc_y = jp % (num_DOFs_1d - 1);
        
        size_t elem_z = kp / (num_DOFs_1d - 1);
        size_t loc_z = kp % (num_DOFs_1d - 1);
        
        node_coords(node_gid, 0) = elem_x*h_x + (lob_nodes_1D(loc_x) - lob_nodes_1D(0))*h_x/ref_length;
        node_coords(node_gid, 1) = elem_y*h_y + (lob_nodes_1D(loc_y) - lob_nodes_1D(0))*h_y/ref_length;
        node_coords(node_gid, 2) = elem_z*h_z + (lob_nodes_1D(loc_z) - lob_nodes_1D(0))*h_z/ref_length;
    });

    // Step 2: Build element connectivity
    FOR_ALL(i, 0, num_elems_x,
            j, 0, num_elems_y,
            k, 0, num_elems_z, {
        
        size_t elem_gid = i + (j + k*num_elems_y)*num_elems_x;
        size_t node_lid = 0;
        
        // Multiply by elem_order to space elements correctly
        for(size_t kc = k*elem_order; kc <= k*elem_order + elem_order; kc++)
        for(size_t jc = j*elem_order; jc <= j*elem_order + elem_order; jc++)
        for(size_t ic = i*elem_order; ic <= i*elem_order + elem_order; ic++){
            size_t node_gid = ic + (jc + kc*num_nodes_y)*num_nodes_x;
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


    // ==========================================
    // Build state 

    const size_t num_qpts_in_elem = Quad.num_qpts_in_elem;

    const size_t num_surfs = Mesh.num_surfs;
    const size_t num_qpts_in_surf = SurfQuad.num_qpts_in_surf;

    const size_t num_nodes_in_elem = Mesh.num_nodes_in_elem;
    const size_t num_surfs_in_elem = Mesh.num_surfs_in_elem;

    const size_t num_corners = Mesh.num_corners;

    if(num_nodes_in_elem != FERefElem.num_dofs_in_elem) Kokkos::abort("ERROR: mismatch in DOFs and num nodes in elem \n");

    // Per quadrature point geometry of the moving mesh.  Both are only ever
    // read and written on the device, so neither needs a host mirror.
    CArrayKokkos<REAL_t> elem_det_jac(num_elems, num_qpts_in_elem, "elem_det_jacobian");
    CArrayKokkos<REAL_t> inv_jac_ijq(num_elems, elem_dims, elem_dims, num_qpts_in_elem, "inv_jac_ijq");

    DCArrayKokkos<REAL_t> elem_field(num_elems, "elem_field");
    DCArrayKokkos<REAL_t> node_field(num_nodes, "node_field");     // for displaying field results
    DCArrayKokkos<REAL_t> node_velocity(num_nodes, elem_dims, "node_velocity");
    DCArrayKokkos<REAL_t> node_velocity_n(num_nodes, elem_dims, "node_velocity_n");
    DCArrayKokkos<REAL_t> node_coords_n(num_nodes, elem_dims, "node_coords_n");

    DCArrayKokkos<REAL_t> corner_field(num_corners, "corner_field");
    DCArrayKokkos<REAL_t> corner_field_n(num_corners, "corner_field_n");
    DCArrayKokkos<REAL_t> elem_corner_vol(num_elems, num_nodes_in_elem, "elem_corner_vol");
    DCArrayKokkos<REAL_t> elem_corner_vol_n(num_elems, num_nodes_in_elem, "elem_corner_vol_n");

    // Calculate RHS_surf_flux
    CArrayKokkos <REAL_t> RHS_surf_flux(num_elems, num_surfs_in_elem, num_qpts_in_surf, "RHS_surf_flux"); // used to build RHS vector 
    CArrayKokkos <REAL_t> RHS_elem(num_elems, num_nodes_in_elem, "RHS_elem"); // RHS vector 


    // ---- transposed reference tables and per quadrature point scratch ----
    BasisTables_t tables;
    CArrayKokkos<REAL_t> qpt_vol_flux(num_elems, num_qpts_in_elem, elem_dims, "qpt_vol_flux");
    CArrayKokkos<REAL_t> elem_err(num_elems, 2, "elem_error_norms");

    tables.basis_row_sum = CArrayKokkos<REAL_t>(num_qpts_in_elem, "basis_row_sum");
    FOR_ALL(qpt_lid, 0, num_qpts_in_elem, {
        REAL_t sum = 0.0;
        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            sum += FERefElem.qpt_basis(qpt_lid, dof_lid);
        }
        tables.basis_row_sum(qpt_lid) = sum;
    });
    tables.basis_dq       = CArrayKokkos<REAL_t>(num_nodes_in_elem, num_qpts_in_elem, "basis_dq");
    tables.grad_basis_njq = CArrayKokkos<REAL_t>(num_nodes_in_elem, elem_dims, num_qpts_in_elem, "grad_basis_njq");
    tables.grad_basis_qjd = CArrayKokkos<REAL_t>(num_qpts_in_elem, elem_dims, num_nodes_in_elem, "grad_basis_qjd");
    FOR_ALL(qpt_lid, 0, num_qpts_in_elem, {
        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            tables.basis_dq(dof_lid, qpt_lid) = FERefElem.qpt_basis(qpt_lid, dof_lid);
            for(size_t dim = 0; dim < elem_dims; dim++){
                const REAL_t val = FERefElem.qpt_grad_basis(qpt_lid, dof_lid, dim);
                tables.grad_basis_njq(dof_lid, dim, qpt_lid) = val;
                tables.grad_basis_qjd(qpt_lid, dim, dof_lid) = val;
            }
        }
    });
    tables.surf_basis_fdq = CArrayKokkos<REAL_t>(num_surfs_in_elem, num_nodes_in_elem,
                                                 num_qpts_in_surf, "surf_basis_fdq");
    tables.surf_grad_fjdq = CArrayKokkos<REAL_t>(num_surfs_in_elem, elem_dims, num_nodes_in_elem,
                                                 num_qpts_in_surf, "surf_grad_fjdq");
    FOR_ALL(face_lid, 0, num_surfs_in_elem, {
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++){
            for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
                tables.surf_basis_fdq(face_lid, dof_lid, qpt_lid) =
                    RefSurf.qpt_basis(face_lid, qpt_lid, dof_lid);
                for(size_t dim = 0; dim < elem_dims; dim++){
                    tables.surf_grad_fjdq(face_lid, dim, dof_lid, qpt_lid) =
                        RefSurf.qpt_grad_basis(face_lid, qpt_lid, dof_lid, dim);
                }
            }
        }
    });
    Kokkos::fence();


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

    build_element_geometry(Mesh, tables, node_coords, elem_det_jac, inv_jac_ijq,
                           num_elems, num_qpts_in_elem, num_nodes_in_elem);
    build_lumped_volume(FERefElem, Quad, tables, elem_det_jac, elem_corner_vol,
                        num_elems, num_qpts_in_elem, num_nodes_in_elem);
    Kokkos::fence();


    // -----------------------------------------------------
    const REAL_t max_vel = 1.0; // the CFL velocity used for calculating dt
    REAL_t h_cfl = 1.e-6;       // the CFL length scale for calculating dt
    REAL_t dt = 1.e-6;          // dt from CFL at start, this time is psuedo time


    // -----------------------------------------------------
    REAL_t time = 0;                    // the time 
    REAL_t time_output = graphics_dt;   // the time for graphics outputs
    size_t output_id = 0;               // the file id for the outputs
    
    std::ofstream err_file("ErrorNorms.txt");
    if (!err_file.is_open()) {
        std::cerr << "Error: Cannot open file " << "ErrorNorms.txt" << std::endl;
        return 0;
    }
    err_file <<" time       L1          L2 \n";
    err_file << std::fixed << std::setprecision(8);


    // ================================================================
    // Step 0: Set the initial conditions

    FOR_ALL(elem_gid, 0, num_elems, {


        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        for(size_t corner_lid=0; corner_lid<num_nodes_in_elem; corner_lid++){
            // remember: corner_lid = node_lid
            const size_t node_gid = nodes_in_elem(corner_lid);
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid,corner_lid);

            //the function, e.g., sin(PI*node_coords(node_gid,dim));
            corner_field(corner_gid) = test_function(node_coords(node_gid,0),node_coords(node_gid,1));
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
    REAL_t sum_elem = 0.0;
    REAL_t domain_mass_t0 = 0.0;
    FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_elem, {

        for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            sum_elem += elem_corner_vol(elem_gid, node_lid)*corner_field(corner_gid);
        }

    }, domain_mass_t0);

    printf("Domain Mass t=0: %f \n", domain_mass_t0);

    DCArrayKokkos <REAL_t> output_node_coords(num_nodes,3);
    

    // export results to Paraview graphics file
    {
        elem_field.set_values(0.0);
        FOR_ALL(elem_gid,0,num_elems,{
            for(size_t node_lid=0; node_lid<Mesh.num_nodes_in_elem; node_lid++){
                const size_t corner_lid = node_lid;
                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, corner_lid);
                elem_field(elem_gid) += corner_field(corner_gid);
            } 
            elem_field(elem_gid) /= (REAL_t)Mesh.num_nodes_in_elem;
        });

        // save corner field to the nodes for graphics outputs
        node_field.set_values(0.0);
        FOR_ALL(node_gid,0,num_nodes,{
            for(size_t corner_lid=0; corner_lid<Mesh.num_corners_in_node(node_gid); corner_lid++){
                const size_t corner_gid = Mesh.corners_in_node(node_gid, corner_lid);
                node_field(node_gid) += corner_field(corner_gid);
            } 
            node_field(node_gid) /= (REAL_t)Mesh.num_corners_in_node(node_gid);
        });

        // map nodes to uniform locations
        output_node_coords.set_values(0.0);
        interpolate_to_uniform(Mesh.nodes_in_elem,
                               lob_nodes_1D,
                               node_coords,  
                               output_node_coords,   
                               num_elems); 

        // writing the initial mesh and state
        output_node_coords.update_host();
        node_field.update_host();
        elem_field.update_host();

        printf(" Writing output at time = %.4f. ", time);

        char filename[100];
        snprintf(filename, sizeof(filename), "output_time_%04zu.vtu", output_id);

        // Write the mesh state
        write_lagrange_hex_mesh(
            filename,
            output_node_coords,           
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

        FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

            const size_t elem_gid = idx / num_nodes_in_elem;
            const size_t node_lid = idx % num_nodes_in_elem;

            elem_corner_vol_n(elem_gid, node_lid) = elem_corner_vol(elem_gid, node_lid);

            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            corner_field_n(corner_gid) = corner_field(corner_gid);
        });

        FOR_ALL(node_gid, 0, num_nodes, {
            for(size_t dim=0; dim<elem_dims; dim++){
                node_coords_n(node_gid, dim)   = node_coords(node_gid, dim);
                node_velocity_n(node_gid, dim) = node_velocity(node_gid, dim);
            }
        });


        // ------------------------------------------------------
        // Step 1b: get CFL time step for moving mesh

        // x -> x^(1/3) is monotone, so the smallest length comes from the
        // smallest quadrature volume: one transcendental instead of one per qpt.
        REAL_t min_vol_loc;
        REAL_t min_vol_qpt;
        FOR_REDUCE_MIN(idx, 0, num_elems*num_qpts_in_elem,
                       min_vol_loc, {
            const size_t elem_gid = idx / num_qpts_in_elem;
            const size_t qpt_lid  = idx % num_qpts_in_elem;

            const REAL_t vol_qpt = Quad.qpt_weights(qpt_lid)*elem_det_jac(elem_gid, qpt_lid);
            if(vol_qpt < min_vol_loc) min_vol_loc = vol_qpt;
        }, min_vol_qpt);
        h_cfl = pow(min_vol_qpt, 0.3333333);

        dt = 0.1*h_cfl/max_vel; // pseudo time step used for the remap

        // A tangled element gives a non-positive quadrature volume, so the CFL
        // length becomes NaN. NaN then defeats both the max_time test and the
        // mass conservation check below, so the run has to stop here.
        if(!(dt > 0.0)){
            printf("\n STOPPING at time = %.6f, cycle %zu: CFL length is %g,"
                   " the mesh has tangled.\n", time, cycle, h_cfl);
            break;
        }


        // Runge Kutta time integration levels
        for(size_t rk_stage=0; rk_stage<rk_num_stages; rk_stage++){

            // RK coefficient
            const REAL_t rk_alpha = 1.0 / ((REAL_t)rk_num_stages - (REAL_t)rk_stage);


            // ------------------------------------------------------
            // Step 2: calculate the mesh velocity at time level k

            FOR_ALL(node_gid, 0, num_nodes,{
                // new velocity, it is Taylor-Green vortex
                // PI is defined in mesh class
                node_velocity(node_gid, 0) =  sin(PI*node_coords(node_gid, 0))*cos(PI*node_coords(node_gid, 1));
                node_velocity(node_gid, 1) = -cos(PI*node_coords(node_gid, 0))*sin(PI*node_coords(node_gid, 1));
                node_velocity(node_gid, 2) = 0.0;
            });


            // ----------------------------------------------------------
            // Step 3: Calculate the surface fluxes at quadrature points

            build_surface_flux(Mesh, RefSurf, SurfQuad, tables, node_coords, node_velocity, corner_field,
                               surf_qpt_qpt_map, RHS_surf_flux,
                               num_surfs, num_qpts_in_surf, num_nodes_in_elem);


            // -------------------------------------------------
            // Step 4: Build RHS of DG equations in the element

            assemble_rhs(Mesh, RefSurf, Quad, tables, elem_det_jac, inv_jac_ijq,
                         corner_field, corner_field_n, node_velocity, elem_corner_vol_n,
                         RHS_surf_flux, RHS_elem, qpt_vol_flux,
                         rk_alpha, dt,
                         num_elems, num_qpts_in_elem, num_nodes_in_elem,
                         num_surfs_in_elem, num_qpts_in_surf, elem_dims);


            // ================================================================
            // Step 5: Move the mesh to the new location
            FOR_ALL(node_gid, 0, num_nodes,{
                // new position of the mesh
                node_coords(node_gid, 0) = node_coords_n(node_gid, 0) + 0.5*(node_velocity(node_gid, 0)+node_velocity_n(node_gid, 0)) * rk_alpha * dt; 
                node_coords(node_gid, 1) = node_coords_n(node_gid, 1) + 0.5*(node_velocity(node_gid, 1)+node_velocity_n(node_gid, 1)) * rk_alpha * dt;
                // z-coords never change
            });


            // ================================================================
            // Step 6: build the diagonal volume matrix for nodal DG after the mesh moved
            build_element_geometry(Mesh, tables, node_coords, elem_det_jac, inv_jac_ijq,
                                   num_elems, num_qpts_in_elem, num_nodes_in_elem);
            build_lumped_volume(FERefElem, Quad, tables, elem_det_jac, elem_corner_vol,
                                num_elems, num_qpts_in_elem, num_nodes_in_elem);


            // -----------------------------------------------------
            // 7. Solve M * u^{n+1} = RHS where M is diagonal

            FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

                const size_t elem_gid = idx / num_nodes_in_elem;
                const size_t dof_lid  = idx % num_nodes_in_elem;

                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
                corner_field(corner_gid) = RHS_elem(elem_gid, dof_lid)/elem_corner_vol(elem_gid, dof_lid);
            });


            // -----------------------------------------------------
            // 8. A slope/bound limiter would be applied to corner_field here;
            //    see limit_corner_field in remap_dg_lumped_test.cpp.

        } // end Runge Kutta time level loop


        // ================================================================
        // Step 7: update time
        time += dt;

        // Conservation Check
        REAL_t sum_elem = 0.0;
        REAL_t domain_mass_time = 0.0;
        FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_elem, {

            for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                sum_elem += elem_corner_vol(elem_gid, node_lid)*corner_field(corner_gid);
            }
        }, domain_mass_time);

        printf("Domain mass error= %f \n", domain_mass_time-domain_mass_t0);
        if(fabs(domain_mass_time-domain_mass_t0)>1.e-12) Kokkos::abort("ERROR: Mass is not conserved");


        // ================================================================
        // Step 8: write outputs
        if( time-time_output >= -1.e-8 ){

            //// L1 and L2 error norms
            // The original walks the same 216x64 reconstruction twice, once per
            // norm. One pass fills both per element, then two cheap reductions
            // over elements keep the original summation tree.
            FOR_ALL(elem_gid, 0, num_elems, {

                REAL_t l1_elem = 0.0;
                REAL_t l2_elem = 0.0;

                for(size_t qpt_lid=0; qpt_lid<num_qpts_in_elem; qpt_lid++){

                    const REAL_t vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

                    REAL_t val_qpt = 0.0;
                    REAL_t x_qpt   = 0.0;
                    REAL_t y_qpt   = 0.0;

                    for(size_t corner_lid=0; corner_lid<num_nodes_in_elem; corner_lid++) {
                        const size_t node_gid   = Mesh.nodes_in_elem(elem_gid,corner_lid);
                        const size_t corner_gid = Mesh.corners_in_elem(elem_gid, corner_lid);
                        const REAL_t phi = FERefElem.qpt_basis(qpt_lid,corner_lid);

                        val_qpt += corner_field(corner_gid)*phi;
                        x_qpt   += node_coords(node_gid,0)*phi;
                        y_qpt   += node_coords(node_gid,1)*phi;
                    }

                    const REAL_t diff = val_qpt - test_function(x_qpt,y_qpt);
                    l1_elem += fabs(diff)*vol_qpt;
                    l2_elem += diff*diff*vol_qpt;
                }

                elem_err(elem_gid, 0) = l1_elem;
                elem_err(elem_gid, 1) = l2_elem;
            });

            REAL_t L1;
            REAL_t L1_lcl;
            FOR_REDUCE_SUM(elem_gid, 0, num_elems, L1_lcl, {
                L1_lcl += elem_err(elem_gid, 0);
            }, L1);

            REAL_t L2;
            REAL_t L2_lcl;
            FOR_REDUCE_SUM(elem_gid, 0, num_elems, L2_lcl, {
                L2_lcl += elem_err(elem_gid, 1);
            }, L2);
            L2 = sqrt(L2);

            printf("=====\n");
            printf("L1 error = %f, L2 error = %f \n", L1, L2);
            err_file << time << " " << L1 << " " << L2 << "\n";

            //////



            elem_field.set_values(0.0);
            FOR_ALL(elem_gid,0,num_elems,{
                for(size_t node_lid=0; node_lid<Mesh.num_nodes_in_elem; node_lid++){
                    const size_t corner_lid = node_lid;
                    const size_t corner_gid = Mesh.corners_in_elem(elem_gid, corner_lid);
                    elem_field(elem_gid) += corner_field(corner_gid);
                } 
                elem_field(elem_gid) /= (REAL_t)Mesh.num_nodes_in_elem;
            });

            // save corner field to the nodes for graphics outputs
            node_field.set_values(0.0);
            FOR_ALL(node_gid,0,num_nodes,{
                for(size_t corner_lid=0; corner_lid<Mesh.num_corners_in_node(node_gid); corner_lid++){
                    const size_t corner_gid = Mesh.corners_in_node(node_gid, corner_lid);
                    node_field(node_gid) += corner_field(corner_gid);
                    //node_field(node_gid) = fmax(node_field(node_gid), corner_field(corner_gid));
                } 
                node_field(node_gid) /= (REAL_t)Mesh.num_corners_in_node(node_gid);
            });

            // map nodes to uniform locations
            output_node_coords.set_values(0.0);
            interpolate_to_uniform(Mesh.nodes_in_elem,
                                lob_nodes_1D,
                                node_coords,  
                                output_node_coords,   
                                num_elems); 

            // writing the initial mesh and state
            output_node_coords.update_host();
            node_field.update_host();
            elem_field.update_host();


            printf(" Writing output at time = %.4f. ", time);

            char filename[100];
            snprintf(filename, sizeof(filename), "output_time_%04zu.vtu", output_id);

            // Write the mesh state
            write_lagrange_hex_mesh(
                filename,
                output_node_coords,           
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

    err_file.close();

    printf("\n Remap test finished.\n");


} // end MATAR scope
MATAR_FINALIZE();

return 0;
} // end function



// 
void write_lagrange_hex_mesh(
    const std::string& filename,
    const DCArrayKokkos<REAL_t>& node_coords,       // All node coordinates [num_nodes][3]
    const size_t num_nodes,
    const DCArrayKokkos<size_t>& nodes_in_elem,     // Connectivity
    const size_t num_elems,
    const size_t order,
    const DCArrayKokkos<REAL_t>& node_data,         // Nodal data
    const std::string& node_data_name,
    const DCArrayKokkos<REAL_t>& elem_data,         // Element center data
    const std::string& elem_data_name)              // Element data name
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

void write_points(std::ofstream& file, const DCArrayKokkos<REAL_t>& coords, size_t num_nodes)
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
                          const DCArrayKokkos<REAL_t>& elem_data,    // Element data
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
                      const DCArrayKokkos<REAL_t>& data, 
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


// Lagrange basis function
KOKKOS_INLINE_FUNCTION
REAL_t lagrange_basis(const REAL_t xi, const size_t i, const CArrayKokkos<REAL_t>& nodes) {
    REAL_t L = 1.0;
    for (size_t j = 0; j < nodes.dims(0); j++) {
        if (j != i) {
            L *= (xi - nodes(j)) / (nodes(i) - nodes(j));
        }
    }
    return L;
}


void interpolate_to_uniform(const DCArrayKokkos<size_t>& nodes_in_elem,
                            const CArrayKokkos<REAL_t>& lob_nodes_1D,
                            const DCArrayKokkos<REAL_t>& node_coords_lob,  // Lobatto node positions
                            DCArrayKokkos<REAL_t>& node_coords_uniform,   // Output uniform positions 
                            const size_t num_elems)   
{

    const size_t num_DOFs_1d = lob_nodes_1D.dims(0);

    // One launch for the whole mesh instead of one launch plus fence per element.
    const size_t num_dofs_in_elem = num_DOFs_1d*num_DOFs_1d*num_DOFs_1d;
    FOR_ALL(idx, 0, num_elems*num_dofs_in_elem, {

        const size_t elem_gid = idx / num_dofs_in_elem;
        const size_t node_lcl = idx % num_dofs_in_elem;

        const size_t i = node_lcl % num_DOFs_1d;
        const size_t j = (node_lcl / num_DOFs_1d) % num_DOFs_1d;
        const size_t k = node_lcl / (num_DOFs_1d*num_DOFs_1d);

        // Uniform parametric coordinates in [-1, 1]
        REAL_t xi   = -1.0 + 2.0 * (REAL_t)i / ((REAL_t)(num_DOFs_1d - 1));
        REAL_t eta  = -1.0 + 2.0 * (REAL_t)j / ((REAL_t)(num_DOFs_1d - 1));
        REAL_t zeta = -1.0 + 2.0 * (REAL_t)k / ((REAL_t)(num_DOFs_1d - 1));

        // Interpolate using Lagrange basis at Lobatto nodes
        REAL_t x = 0.0;
        REAL_t y = 0.0;
        REAL_t z = 0.0;
        for (size_t kk = 0; kk < num_DOFs_1d; kk++) {
            REAL_t Lk = lagrange_basis(zeta, kk, lob_nodes_1D);
            for (size_t jj = 0; jj < num_DOFs_1d; jj++) {
                REAL_t Lj = lagrange_basis(eta, jj, lob_nodes_1D);
                for (size_t ii = 0; ii < num_DOFs_1d; ii++) {
                    REAL_t Li = lagrange_basis(xi, ii, lob_nodes_1D);

                    size_t node_lid = ii + (jj + kk*num_DOFs_1d)*num_DOFs_1d;
                    size_t node = nodes_in_elem(elem_gid, node_lid);
                    REAL_t basis = Li * Lj * Lk;

                    x += basis * node_coords_lob(node, 0);
                    y += basis * node_coords_lob(node, 1);
                    z += basis * node_coords_lob(node, 2);
                }
            }
        }

        // Store interpolated position
        size_t node_gid = nodes_in_elem(elem_gid, node_lcl);
        node_coords_uniform(node_gid,0) = x;
        node_coords_uniform(node_gid,1) = y;
        node_coords_uniform(node_gid,2) = z;
    });
    Kokkos::fence();
} // end function


// ============================================================================
// Notched CIRCLE FUNCTION IMPLEMENTATION
// ============================================================================
#ifdef USE_NOTCHED_CIRCLE 

KOKKOS_INLINE_FUNCTION
REAL_t test_function(const REAL_t x, 
                     const REAL_t y){


    const REAL_t x0 = 0.5;
    const REAL_t y0 = 0.5;
    const REAL_t radius = 0.25;
    const REAL_t notch_width = 0.15;
    const REAL_t notch_depth = 0.4-radius;  // How far the notch cuts INTO the circle
    const REAL_t smoothing = 0.01;   // Smoothing width, 
    const REAL_t eps = 1e-10;


    // Check if inside circle
    const REAL_t dx = x - x0;
    const REAL_t dy = y - y0;
    const REAL_t r = sqrt(dx*dx + dy*dy);


    // Smooth circle
    REAL_t circle = 0.5 * (1.0 - tanh((r - radius) / smoothing));
    
    // Deep notch from TOP, extending PAST center
    if(fabs(x - x0) < notch_width/2.0){  // Remove y > y0 condition!
        REAL_t notch = 0.5 * (1.0 + tanh((y - (y0 - notch_depth)) / smoothing));
        circle *= (1.0 - notch);
    }
    
    return circle;  // Inside circle, outside notch
    
} // end function

#endif // USE_NOTCHED_CIRCLE



// ============================================================================
// SIN FUNCTION IMPLEMENTATION
// ============================================================================
#ifdef USE_SIN_FUNCTION

KOKKOS_INLINE_FUNCTION
REAL_t test_function(REAL_t x, REAL_t y) {
    return sin(PI*x);
}

#endif // USE_SIN_FUNCTION



// ============================================================================
// GAUSSIAN FUNCTION IMPLEMENTATION
// ============================================================================
#ifdef USE_GAUSSIAN

KOKKOS_INLINE_FUNCTION
REAL_t test_function(REAL_t x, REAL_t y) {
    const REAL_t cx = 0.5;
    const REAL_t cy = 0.5;
    const REAL_t sigma = 0.1;
    
    REAL_t dx = x - cx;
    REAL_t dy = y - cy;
    REAL_t r2 = dx*dx + dy*dy;
    
    return exp(-r2 / (2.0 * sigma * sigma));
}

#endif // USE_GAUSSIAN
