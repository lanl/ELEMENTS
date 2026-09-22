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
// GPU-tuned variant of remap_dg_lumped_test.cpp.
//
// REMAP_OPT selects how much of the tuning is active so the contribution of
// each step can be measured against the untouched baseline:
//
//   0  kernels identical to remap_dg_lumped_test.cpp (reference)
//   1  + flat 1D index spaces instead of the FOR_FIRST/FOR_SECOND team policy
//   2  + Jacobian / determinant / inverse evaluated once per quadrature point
//   3  + field and velocity reconstruction hoisted out of the DOF loop
//   4  + coalesced basis tables, factored volume flux, single-launch output map
//   5  + coalesced surface tables, register-resident surface Jacobian, fused
//        L1/L2 error pass
//   6  + drop the host fences between stream-ordered steps (only visible in a
//        build without ENABLE_PHASE_TIMERS, whose timers fence anyway)
//
// Level 1 is kept only as a measurement point: flattening the team policy
// exposes a write race on the shared Jacobian scratch that levels >= 2 remove
// by computing the geometry in its own pass, so level 1 does not reproduce the
// ground truth.
//
// Optional arguments (any order of numbers vs the field file):
//   num_elems_x num_elems_y [num_elems_z [max_time [graphics_dt]]] [field.rfld]
//
// The optional .rfld raster is produced by
//   examples/reference_element/scripts/make_field.py
// and replaces the compiled-in test_function for both the initial condition
// and the Eulerian error-norm comparison. REMAP_FIELD=path is also accepted.
//
#ifndef REMAP_OPT
#define REMAP_OPT 5
#endif

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
#include <chrono>

// This pulls in kokkos, matar, mesh, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"
#include "cramers_rule.hpp" // det and solvers
#include "lu_solver.hpp"

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


// Eulerian field sampled by bilinear interpolation from a raster written by
// make_field.py. nx<=1 means "use the compiled-in test_function".
struct EulerianField_t
{
    int nx = 0;
    int ny = 0;
    REAL_t xmin = 0.0;
    REAL_t xmax = 1.0;
    REAL_t ymin = 0.0;
    REAL_t ymax = 1.0;
    DCArrayKokkos<REAL_t> values; // (ny, nx), j=0 at ymin

    KOKKOS_INLINE_FUNCTION
    REAL_t operator()(const REAL_t x, const REAL_t y) const
    {
        if (nx <= 1 || ny <= 1) {
            return test_function(x, y);
        }

        const REAL_t fx = (x - xmin) / (xmax - xmin) * (REAL_t)(nx - 1);
        const REAL_t fy = (y - ymin) / (ymax - ymin) * (REAL_t)(ny - 1);

        REAL_t uc = fx;
        REAL_t vc = fy;
        if (uc < 0.0) uc = 0.0;
        if (vc < 0.0) vc = 0.0;
        if (uc > (REAL_t)(nx - 1)) uc = (REAL_t)(nx - 1);
        if (vc > (REAL_t)(ny - 1)) vc = (REAL_t)(ny - 1);

        const int i0 = (int)uc;
        const int j0 = (int)vc;
        const int i1 = (i0 + 1 < nx) ? i0 + 1 : i0;
        const int j1 = (j0 + 1 < ny) ? j0 + 1 : j0;
        const REAL_t tx = uc - (REAL_t)i0;
        const REAL_t ty = vc - (REAL_t)j0;

        const REAL_t v00 = values(j0, i0);
        const REAL_t v10 = values(j0, i1);
        const REAL_t v01 = values(j1, i0);
        const REAL_t v11 = values(j1, i1);
        return (1.0 - ty) * ((1.0 - tx) * v00 + tx * v10)
             +        ty  * ((1.0 - tx) * v01 + tx * v11);
    }
};


static bool is_numeric_token(const char* s)
{
    char* end = nullptr;
    std::strtod(s, &end);
    return end != s && *end == '\0';
}


static bool load_eulerian_field(const char* path, EulerianField_t& field)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        printf("ERROR: cannot open field file %s\n", path);
        return false;
    }

    char mag[4];
    in.read(mag, 4);
    if (!in || std::memcmp(mag, "RF01", 4) != 0) {
        printf("ERROR: %s is not an RF01 field raster (run make_field.py)\n", path);
        return false;
    }

    int32_t nx = 0;
    int32_t ny = 0;
    in.read(reinterpret_cast<char*>(&nx), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&ny), sizeof(int32_t));
    double box[4] = {0.0, 0.0, 0.0, 0.0};
    in.read(reinterpret_cast<char*>(box), 4 * sizeof(double));
    if (!in || nx < 2 || ny < 2) {
        printf("ERROR: %s has a truncated or invalid header\n", path);
        return false;
    }

    const size_t n = (size_t)nx * (size_t)ny;
    std::vector<double> buf(n);
    in.read(reinterpret_cast<char*>(buf.data()),
            (std::streamsize)(n * sizeof(double)));
    if (!in) {
        printf("ERROR: %s is truncated (%zu samples expected)\n", path, n);
        return false;
    }

    field.nx   = (int)nx;
    field.ny   = (int)ny;
    field.xmin = box[0];
    field.xmax = box[1];
    field.ymin = box[2];
    field.ymax = box[3];
    field.values = DCArrayKokkos<REAL_t>((size_t)ny, (size_t)nx, "eulerian_field");
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            field.values.host(j, i) = buf[(size_t)j * (size_t)nx + (size_t)i];
        }
    }
    field.values.update_device();
    printf("loaded Eulerian field %s  (%d x %d over [%g,%g] x [%g,%g])\n",
           path, nx, ny, field.xmin, field.xmax, field.ymin, field.ymax);
    return true;
}


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

// Per quadrature point geometry of the moving mesh.
struct GeometryState_t
{
    DCArrayKokkos<REAL_t> elem_jac;     // (elem, qpt, 3, 3)
    DCArrayKokkos<REAL_t> elem_det_jac; // (elem, qpt)
    DCArrayKokkos<REAL_t> elem_inv_jac; // (elem, qpt, 3, 3)
    CArrayKokkos<REAL_t>  inv_jac_ijq;  // (elem, 3, 3, qpt) - qpt fastest
};


// Consecutive kernels all run in the default execution space, so they are
// already stream ordered; the host only has to wait where it reads a result.
// The fences between steps exist to make the phase timers meaningful.
#if REMAP_OPT >= 6
#define STEP_FENCE() ((void)0)
#else
#define STEP_FENCE() Kokkos::fence()
#endif

// ============================================================================
// Coarse phase timers. stop() fences, so they are only compiled in on request.
// ============================================================================
enum phase_id {
    PH_STORE = 0, PH_CFL, PH_VELOCITY, PH_SURF_FLUX, PH_RHS,
    PH_MOVE, PH_MASS, PH_SOLVE, PH_CONSERVE, PH_OUTPUT, PH_COUNT
};

struct PhaseTimer
{
    double total[PH_COUNT] = {};
    std::chrono::high_resolution_clock::time_point t0;

    inline void start()
    {
#ifdef ENABLE_PHASE_TIMERS
        Kokkos::fence();
        t0 = std::chrono::high_resolution_clock::now();
#endif
    }

    inline void stop(phase_id id)
    {
#ifdef ENABLE_PHASE_TIMERS
        Kokkos::fence();
        total[id] += std::chrono::duration<double>(
                         std::chrono::high_resolution_clock::now() - t0).count();
#else
        (void)id;
#endif
    }

    void report() const
    {
#ifdef ENABLE_PHASE_TIMERS
        static const char* name[PH_COUNT] = {
            "store state (1a)", "cfl dt (1b)", "mesh velocity (2)", "surface flux (3)",
            "rhs assembly (4)", "move mesh (5)", "lumped mass (6)", "dof solve (7)",
            "conservation", "graphics + norms"
        };
        double sum = 0.0;
        printf("\n---- phase breakdown (s) ----\n");
        for (int i = 0; i < PH_COUNT; i++) {
            printf("  %-20s %10.4f\n", name[i], total[i]);
            sum += total[i];
        }
        printf("  %-20s %10.4f\n", "TOTAL timed", sum);
#endif
    }
};


// ============================================================================
// Geometry of the moving mesh at every volume quadrature point.
// At REMAP_OPT < 2 this is folded into build_lumped_volume, exactly as the
// original code does (the Jacobian is then rebuilt once per DOF).
// ============================================================================
static void build_element_geometry(const Mesh_t& Mesh,
                                   const ReferenceElement_t& FERefElem,
                                   const BasisTables_t& tables,
                                   const DCArrayKokkos<REAL_t>& node_coords,
                                   const GeometryState_t& geom,
                                   const size_t num_elems,
                                   const size_t num_qpts_in_elem,
                                   const size_t num_nodes_in_elem)
{
#if REMAP_OPT >= 4
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
        geom.elem_det_jac(elem_gid, qpt_lid) = det;

        // same expression as invert_3x3 in cramers_rule.hpp
        const REAL_t den = det + 1e-16;
        geom.inv_jac_ijq(elem_gid, 0, 0, qpt_lid) = +(j11*j22 - j12*j21) / den;
        geom.inv_jac_ijq(elem_gid, 0, 1, qpt_lid) = -(j01*j22 - j02*j21) / den;
        geom.inv_jac_ijq(elem_gid, 0, 2, qpt_lid) = +(j01*j12 - j02*j11) / den;
        geom.inv_jac_ijq(elem_gid, 1, 0, qpt_lid) = -(j10*j22 - j12*j20) / den;
        geom.inv_jac_ijq(elem_gid, 1, 1, qpt_lid) = +(j00*j22 - j02*j20) / den;
        geom.inv_jac_ijq(elem_gid, 1, 2, qpt_lid) = -(j00*j12 - j02*j10) / den;
        geom.inv_jac_ijq(elem_gid, 2, 0, qpt_lid) = +(j10*j21 - j11*j20) / den;
        geom.inv_jac_ijq(elem_gid, 2, 1, qpt_lid) = -(j00*j21 - j01*j20) / den;
        geom.inv_jac_ijq(elem_gid, 2, 2, qpt_lid) = +(j00*j11 - j01*j10) / den;
    });
#elif REMAP_OPT >= 2
    FOR_ALL(idx, 0, num_elems*num_qpts_in_elem, {

        const size_t elem_gid = idx / num_qpts_in_elem;
        const size_t qpt_lid  = idx % num_qpts_in_elem;

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);
        ViewCArrayKokkos<REAL_t> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0), num_nodes_in_elem, 3);
        ViewCArrayKokkos<REAL_t> jac(&geom.elem_jac(elem_gid,qpt_lid,0,0), 3, 3);
        ViewCArrayKokkos<REAL_t> inv_jac(&geom.elem_inv_jac(elem_gid,qpt_lid,0,0), 3, 3);

        jacobian(jac, node_coords, nodes_in_elem, a_grad_basis);

        geom.elem_det_jac(elem_gid, qpt_lid) = det_3x3(jac);
        invert_3x3(jac, inv_jac, geom.elem_det_jac(elem_gid, qpt_lid));
    });
#else
    (void)Mesh; (void)FERefElem; (void)tables; (void)node_coords; (void)geom;
    (void)num_elems; (void)num_qpts_in_elem; (void)num_nodes_in_elem;
#endif
} // end build_element_geometry


// ============================================================================
// Row-lumped volume (mass) vector, elem_corner_vol(elem, node).
// ============================================================================
static void build_lumped_volume(const Mesh_t& Mesh,
                                const ReferenceElement_t& FERefElem,
                                const Quadrature_t& Quad,
                                const BasisTables_t& tables,
                                const DCArrayKokkos<REAL_t>& node_coords,
                                const GeometryState_t& geom,
                                DCArrayKokkos<REAL_t>& elem_corner_vol,
                                const size_t num_elems,
                                const size_t num_qpts_in_elem,
                                const size_t num_nodes_in_elem)
{
#if REMAP_OPT >= 2
    // The inner DOF loop of the original only contributes the row sum of the
    // basis, which is the same at every element, so it is tabulated once.
    (void)node_coords;
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t node_lid = idx % num_nodes_in_elem;

        REAL_t vol = 0.0;
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){
            const REAL_t vol_qpt = geom.elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
            vol += tables.basis_row_sum(qpt_lid)*FERefElem.qpt_basis(qpt_lid, node_lid)*vol_qpt;
        }
        elem_corner_vol(elem_gid, node_lid) = vol;
    });
#elif REMAP_OPT >= 1
    elem_corner_vol.set_values(0.0);
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t node_lid = idx % num_nodes_in_elem;

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){

            ViewCArrayKokkos<REAL_t> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0), num_nodes_in_elem, 3);
            ViewCArrayKokkos<REAL_t> a_basis(&FERefElem.qpt_basis(qpt_lid,0), num_nodes_in_elem);
            ViewCArrayKokkos<REAL_t> jac(&geom.elem_jac(elem_gid,qpt_lid,0,0), 3, 3);
            ViewCArrayKokkos<REAL_t> inv_jac(&geom.elem_inv_jac(elem_gid,qpt_lid,0,0), 3, 3);

            jacobian(jac, node_coords, nodes_in_elem, a_grad_basis);

            geom.elem_det_jac(elem_gid, qpt_lid) = det_3x3(jac);
            invert_3x3(jac, inv_jac, geom.elem_det_jac(elem_gid, qpt_lid));

            const REAL_t vol_qpt = geom.elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

            for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
                elem_corner_vol(elem_gid, node_lid) += a_basis(dof_lid)*a_basis(node_lid)*vol_qpt;
            }
        }
    });
    (void)tables;
#else
    elem_corner_vol.set_values(0.0);
    FOR_FIRST(elem_gid, 0, num_elems, {

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        FOR_SECOND(node_lid, 0, num_nodes_in_elem, {

            for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){

                ViewCArrayKokkos<REAL_t> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0), num_nodes_in_elem, 3);
                ViewCArrayKokkos<REAL_t> a_basis(&FERefElem.qpt_basis(qpt_lid,0), num_nodes_in_elem);
                ViewCArrayKokkos<REAL_t> jac(&geom.elem_jac(elem_gid,qpt_lid,0,0), 3, 3);
                ViewCArrayKokkos<REAL_t> inv_jac(&geom.elem_inv_jac(elem_gid,qpt_lid,0,0), 3, 3);

                jacobian(jac, node_coords, nodes_in_elem, a_grad_basis);

                geom.elem_det_jac(elem_gid, qpt_lid) = det_3x3(jac);
                invert_3x3(jac, inv_jac, geom.elem_det_jac(elem_gid, qpt_lid));

                const REAL_t vol_qpt = geom.elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);

                for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
                    elem_corner_vol(elem_gid, node_lid) += a_basis(dof_lid)*a_basis(node_lid)*vol_qpt;
                }
            }
        });
    });
    (void)tables;
#endif
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
                               const DCArrayKokkos<REAL_t>& surf_jac,
                               const CArrayKokkos<int>& surf_qpt_qpt_map,
                               CArrayKokkos<REAL_t>& RHS_surf_flux,
                               const size_t num_surfs,
                               const size_t num_qpts_in_surf,
                               const size_t num_nodes_in_elem,
                               const size_t elem_dims)
{
    RHS_surf_flux.set_values(0.0);

#if REMAP_OPT >= 5
    // Same math as below, but the reference tables are indexed with the surface
    // quadrature point as the fastest axis so that a warp reads contiguous
    // doubles, and the surface Jacobian never leaves registers.
    (void)surf_jac;
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

        const REAL_t den = det_jac_qpt + 1e-16;
        const REAL_t i00 = +(j11*j22 - j12*j21) / den;
        const REAL_t i01 = -(j01*j22 - j02*j21) / den;
        const REAL_t i02 = +(j01*j12 - j02*j11) / den;
        const REAL_t i10 = -(j10*j22 - j12*j20) / den;
        const REAL_t i11 = +(j00*j22 - j02*j20) / den;
        const REAL_t i12 = -(j00*j12 - j02*j10) / den;
        const REAL_t i20 = +(j10*j21 - j11*j20) / den;
        const REAL_t i21 = -(j00*j21 - j01*j20) / den;
        const REAL_t i22 = +(j00*j11 - j01*j10) / den;

        const REAL_t n0 = RefSurf.outward_normal(face_lid, 0);
        const REAL_t n1 = RefSurf.outward_normal(face_lid, 1);
        const REAL_t n2 = RefSurf.outward_normal(face_lid, 2);
        const REAL_t scale = det_jac_qpt*SurfQuad.qpt_weights(face_lid, qpt_lid);

        // Nanson's formula: s*J^-1*n
        const REAL_t area_normal0 = (n0*i00 + n1*i10 + n2*i20)*scale;
        const REAL_t area_normal1 = (n0*i01 + n1*i11 + n2*i21)*scale;
        const REAL_t area_normal2 = (n0*i02 + n1*i12 + n2*i22)*scale;

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

        const REAL_t flux_val = 0.5*(qpt_field+nbr_qpt_field)*normal_dot_vel
                               -0.5*fabs(normal_dot_vel)*(qpt_field-nbr_qpt_field);

        RHS_surf_flux(elem_gid, face_lid, qpt_lid) = flux_val;
        if(num_elems_in_surf == 2) RHS_surf_flux(nbr_elem_gid, nbr_face_lid, nbr_qpt_lid) = -flux_val;
    });
#else
    (void)tables;
#if REMAP_OPT >= 1
    FOR_ALL(idx, 0, num_surfs*num_qpts_in_surf, {

        const size_t surf_gid = idx / num_qpts_in_surf;
        const size_t qpt_lid  = idx % num_qpts_in_surf;

        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);
        const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0);
        const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0);

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);
#else
    FOR_FIRST(surf_gid, 0, num_surfs, {

        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);
        const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0);
        const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0);

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);

        FOR_SECOND(qpt_lid, 0, num_qpts_in_surf, {
#endif

            // extract the grad_basis at a single quadrature point (surf,qpt,dof,3D)
            ViewCArrayKokkos<REAL_t> a_grad_basis(&RefSurf.qpt_grad_basis(face_lid,qpt_lid,0,0),
                                                  num_nodes_in_elem, 3);

            // extract the basis at a single quadrature point (surf,qpt,dof)
            ViewCArrayKokkos<REAL_t> a_basis(&RefSurf.qpt_basis(face_lid,qpt_lid,0),
                                             num_nodes_in_elem);

            ViewCArrayKokkos<REAL_t> jac(&surf_jac(surf_gid,qpt_lid,0,0), 3, 3);

            REAL_t surf_inv_jac_1D[9];
            ViewCArrayKokkos<REAL_t> inv_jac(&surf_inv_jac_1D[0], 3, 3);

            jacobian(jac, node_coords, nodes_in_elem, a_grad_basis);

            const REAL_t det_jac_qpt = det_3x3(jac);

            invert_3x3(jac, inv_jac, det_jac_qpt);

            // Nanson's formula: s*J^-1*j*f*w
            REAL_t area_normal[3];
            area_normal[0] = 0.;
            area_normal[1] = 0.;
            area_normal[2] = 0.;
            for(size_t j = 0; j < elem_dims; j++){
                for(size_t i = 0; i < elem_dims; i++){
                    area_normal[j] += RefSurf.outward_normal(face_lid,i)*inv_jac(i,j);
                }
                area_normal[j] *= det_jac_qpt*SurfQuad.qpt_weights(face_lid,qpt_lid);
            }

            REAL_t qpt_vel[3];
            for(size_t dim = 0; dim < elem_dims; dim++){
                qpt_vel[dim] = 0.0;
            }

            for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
            for(size_t dim = 0; dim < elem_dims; dim++){
                const size_t node_gid = nodes_in_elem(node_lid);
                qpt_vel[dim] += a_basis(node_lid)*node_velocity(node_gid,dim);
            }

            REAL_t normal_dot_vel = 0.0;
            for(size_t dim = 0; dim < elem_dims; dim++){
                normal_dot_vel += area_normal[dim]*qpt_vel[dim];
            }

            size_t nbr_elem_gid = elem_gid;
            size_t nbr_face_lid = face_lid;
            if(num_elems_in_surf == 2){
                nbr_elem_gid = Mesh.elems_in_surf(surf_gid, 1); // second elem
                nbr_face_lid = Mesh.faces_in_surf(surf_gid, 1); // second elem face
            }

            const size_t nbr_qpt_lid = surf_qpt_qpt_map(surf_gid,0,qpt_lid); // matching qpt

            ViewCArrayKokkos<REAL_t> a_nbr_basis(&RefSurf.qpt_basis(nbr_face_lid,nbr_qpt_lid,0),
                                                 num_nodes_in_elem);

            // reconstruct the fields
            REAL_t qpt_field     = 0.0;
            REAL_t nbr_qpt_field = 0.0;

            for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

                // Note: corner_lid = node_lid inside the element
                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                qpt_field += a_basis(node_lid)*corner_field(corner_gid);

                const size_t nbr_corner_gid = Mesh.corners_in_elem(nbr_elem_gid, node_lid);
                nbr_qpt_field += a_nbr_basis(node_lid)*corner_field(nbr_corner_gid);
            }

            // Rusanov flux at the quadrature point.
            // if normal_dot_vel<0 advection is out of first elem in the surf
            const REAL_t flux_val = 0.5*(qpt_field+nbr_qpt_field)*normal_dot_vel
                                   -0.5*fabs(normal_dot_vel)*(qpt_field-nbr_qpt_field);

            // save flux value to the quadrature points on either side of the element
            RHS_surf_flux(elem_gid, face_lid, qpt_lid) = flux_val;
            if(num_elems_in_surf == 2) RHS_surf_flux(nbr_elem_gid, nbr_face_lid, nbr_qpt_lid) = -flux_val;

#if REMAP_OPT >= 1
    }); // end surf-qpt loop
#else
        }); // end parallel for qpt
    }); // end surf loop
#endif
#endif // REMAP_OPT >= 5
} // end build_surface_flux


// ============================================================================
// RHS of the DG equations.
// ============================================================================
static void assemble_rhs(const Mesh_t& Mesh,
                         const ReferenceElement_t& FERefElem,
                         const ReferenceSurface_t& RefSurf,
                         const Quadrature_t& Quad,
                         const BasisTables_t& tables,
                         const GeometryState_t& geom,
                         const DCArrayKokkos<REAL_t>& corner_field,
                         const DCArrayKokkos<REAL_t>& corner_field_n,
                         const DCArrayKokkos<REAL_t>& node_velocity,
                         const DCArrayKokkos<REAL_t>& elem_corner_vol_n,
                         const CArrayKokkos<REAL_t>& RHS_surf_flux,
                         CArrayKokkos<REAL_t>& RHS_elem,
                         const CArrayKokkos<REAL_t>& qpt_field_elem,
                         const CArrayKokkos<REAL_t>& qpt_vel_elem,
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
#if REMAP_OPT >= 3
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
    #if REMAP_OPT >= 4
            const REAL_t basis_val = tables.basis_dq(node_lid, qpt_lid);
    #else
            const REAL_t basis_val = FERefElem.qpt_basis(qpt_lid, node_lid);
    #endif
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            qpt_field += basis_val*corner_field(corner_gid);
        }

        for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
    #if REMAP_OPT >= 4
            const REAL_t basis_val = tables.basis_dq(node_lid, qpt_lid);
    #else
            const REAL_t basis_val = FERefElem.qpt_basis(qpt_lid, node_lid);
    #endif
            const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
            qpt_vel_0 += basis_val*node_velocity(node_gid, 0);
            qpt_vel_1 += basis_val*node_velocity(node_gid, 1);
            qpt_vel_2 += basis_val*node_velocity(node_gid, 2);
        }

    #if REMAP_OPT >= 4
        // grad(phi).J^-1.(v U) = sum_j dphi/dxi_j * [ sum_i Jinv(j,i) v_i U ],
        // so the DOF-independent bracket is tabulated here.
        const REAL_t vol_qpt = geom.elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
        for(size_t j = 0; j < elem_dims; j++){
            const REAL_t flux = geom.inv_jac_ijq(elem_gid, j, 0, qpt_lid)*qpt_vel_0
                              + geom.inv_jac_ijq(elem_gid, j, 1, qpt_lid)*qpt_vel_1
                              + geom.inv_jac_ijq(elem_gid, j, 2, qpt_lid)*qpt_vel_2;
            qpt_vol_flux(elem_gid, qpt_lid, j) = flux*qpt_field*vol_qpt;
        }
    #else
        qpt_field_elem(elem_gid, qpt_lid) = qpt_field;
        qpt_vel_elem(elem_gid, qpt_lid, 0) = qpt_vel_0;
        qpt_vel_elem(elem_gid, qpt_lid, 1) = qpt_vel_1;
        qpt_vel_elem(elem_gid, qpt_lid, 2) = qpt_vel_2;
    #endif
    });
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
    #if REMAP_OPT >= 4
            vol_integral += tables.grad_basis_qjd(qpt_lid, 0, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 0)
                          + tables.grad_basis_qjd(qpt_lid, 1, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 1)
                          + tables.grad_basis_qjd(qpt_lid, 2, dof_lid)*qpt_vol_flux(elem_gid, qpt_lid, 2);
    #else
            ViewCArrayKokkos<REAL_t> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0),
                                                  num_nodes_in_elem, 3);
            ViewCArrayKokkos<REAL_t> inv_jac(&geom.elem_inv_jac(elem_gid,qpt_lid,0,0), 3, 3);

            REAL_t physical_grad[3];
            physical_grad[0] = 0.0;
            physical_grad[1] = 0.0;
            physical_grad[2] = 0.0;
            for(size_t i = 0; i < elem_dims; i++)
            for(size_t j = 0; j < elem_dims; j++){
                physical_grad[i] += a_grad_basis(dof_lid, j)*inv_jac(j,i);
            }

            REAL_t grad_dot_flux = 0.0;
            for(size_t dim = 0; dim < elem_dims; dim++){
                grad_dot_flux += physical_grad[dim]*qpt_vel_elem(elem_gid, qpt_lid, dim)
                               * qpt_field_elem(elem_gid, qpt_lid);
            }

            vol_integral += grad_dot_flux*geom.elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
    #endif
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
#else
    (void)tables; (void)qpt_field_elem; (void)qpt_vel_elem; (void)qpt_vol_flux;

    RHS_elem.set_values(0.0);
    #if REMAP_OPT >= 1
    FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

        const size_t elem_gid = idx / num_nodes_in_elem;
        const size_t dof_lid  = idx % num_nodes_in_elem;
    #else
    FOR_FIRST(elem_gid, 0, num_elems, {

        FOR_SECOND(dof_lid, 0, num_nodes_in_elem, {
    #endif

            // 4a. First add to RHS the M*u^n term
            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
            RHS_elem(elem_gid, dof_lid) +=
                elem_corner_vol_n(elem_gid, dof_lid)*corner_field_n(corner_gid);

            // 4b. Subtract the VOLUME integral
            for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++){

                ViewCArrayKokkos<REAL_t> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid,0,0),
                                                      num_nodes_in_elem, 3);
                ViewCArrayKokkos<REAL_t> a_basis(&FERefElem.qpt_basis(qpt_lid,0),
                                                 num_nodes_in_elem);
                ViewCArrayKokkos<REAL_t> inv_jac(&geom.elem_inv_jac(elem_gid,qpt_lid,0,0), 3, 3);

                REAL_t qpt_field = 0.0;
                for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                    const size_t inner_corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                    qpt_field += a_basis(node_lid)*corner_field(inner_corner_gid);
                }

                REAL_t qpt_vel[3];
                qpt_vel[0] = 0.0;
                qpt_vel[1] = 0.0;
                qpt_vel[2] = 0.0;

                for(size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                    const size_t node_gid = Mesh.nodes_in_elem(elem_gid, node_lid);
                    for(size_t dim = 0; dim < elem_dims; dim++){
                        qpt_vel[dim] += a_basis(node_lid)*node_velocity(node_gid, dim);
                    }
                }

                REAL_t physical_grad[3];
                physical_grad[0] = 0.0;
                physical_grad[1] = 0.0;
                physical_grad[2] = 0.0;
                for(size_t i = 0; i < elem_dims; i++)
                for(size_t j = 0; j < elem_dims; j++){
                    physical_grad[i] += a_grad_basis(dof_lid, j)*inv_jac(j,i);
                }

                REAL_t grad_dot_flux = 0.0;
                for(size_t dim = 0; dim < elem_dims; dim++){
                    grad_dot_flux += physical_grad[dim]*qpt_vel[dim]*qpt_field;
                }

                const REAL_t vol_qpt = geom.elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
                RHS_elem(elem_gid, dof_lid) -= rk_alpha*dt*grad_dot_flux*vol_qpt;
            }

            // 4c. Add SURFACE flux contribution
            for(size_t face_lid = 0; face_lid < num_surfs_in_elem; face_lid++)
            for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++){

                ViewCArrayKokkos<REAL_t> a_basis(&RefSurf.qpt_basis(face_lid, qpt_lid, 0),
                                                 num_nodes_in_elem);

                RHS_elem(elem_gid, dof_lid) +=
                    rk_alpha*dt*RHS_surf_flux(elem_gid, face_lid, qpt_lid)*a_basis(dof_lid);
            }

    #if REMAP_OPT >= 1
    }); // end elem-dof loop
    #else
        }); // end parallel for over dof_lid
    });
    #endif
#endif
} // end assemble_rhs


int main(int argc, char** argv) {

MATAR_INITIALIZE(argc, argv);
{ // MATAR scope
    std::cout<<"Reference Element Remap Example! (REMAP_OPT="<<REMAP_OPT<<")"<<std::endl;

    const auto wall_t0 = std::chrono::high_resolution_clock::now();
    PhaseTimer timer;
    double io_seconds = 0.0;

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
    const char* field_path = std::getenv("REMAP_FIELD");

    int nnum = 0;
    double nums[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (int a = 1; a < argc; ++a) {
        if (is_numeric_token(argv[a])) {
            if (nnum < 5) nums[nnum++] = std::strtod(argv[a], nullptr);
        } else {
            field_path = argv[a];
        }
    }
    if (nnum >= 2) {
        num_elems_x = (size_t)nums[0];
        num_elems_y = (size_t)nums[1];
    }
    if (nnum >= 3) num_elems_z  = (size_t)nums[2];
    if (nnum >= 4) max_time     = nums[3];
    if (nnum >= 5) graphics_dt  = nums[4];

    const size_t rk_num_stages = 2;    // number of runge kutta time integration levels
    const size_t max_cycles = 10000000;

    EulerianField_t field;
    if (field_path) {
        if (!load_eulerian_field(field_path, field)) {
            return 1;
        }
    }

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

    GeometryState_t geom;
    geom.elem_jac     = DCArrayKokkos<REAL_t>(num_elems, num_qpts_in_elem, elem_dims, elem_dims, "elem_jacobian");
    geom.elem_det_jac = DCArrayKokkos<REAL_t>(num_elems, num_qpts_in_elem, "elem_det_jacobian");
    geom.elem_inv_jac = DCArrayKokkos<REAL_t>(num_elems, num_qpts_in_elem, elem_dims, elem_dims, "elem_inv_jacobian");
#if REMAP_OPT >= 4
    geom.inv_jac_ijq  = CArrayKokkos<REAL_t>(num_elems, elem_dims, elem_dims, num_qpts_in_elem, "inv_jac_ijq");
#endif
    const DCArrayKokkos<REAL_t>& elem_det_jac = geom.elem_det_jac;

    DCArrayKokkos<REAL_t> surf_jac(num_surfs, num_qpts_in_surf, elem_dims, elem_dims, "surf_jacobian");
    DCArrayKokkos<REAL_t> surf_flux(num_surfs, "surf_flux");
    
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
    CArrayKokkos<REAL_t> qpt_field_elem;
    CArrayKokkos<REAL_t> qpt_vel_elem;
    CArrayKokkos<REAL_t> qpt_vol_flux;
    CArrayKokkos<REAL_t> elem_err;
#if REMAP_OPT >= 5
    elem_err = CArrayKokkos<REAL_t>(num_elems, 2, "elem_error_norms");
#endif

#if REMAP_OPT >= 2
    tables.basis_row_sum = CArrayKokkos<REAL_t>(num_qpts_in_elem, "basis_row_sum");
    FOR_ALL(qpt_lid, 0, num_qpts_in_elem, {
        REAL_t sum = 0.0;
        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            sum += FERefElem.qpt_basis(qpt_lid, dof_lid);
        }
        tables.basis_row_sum(qpt_lid) = sum;
    });
#endif
#if REMAP_OPT == 3
    qpt_field_elem = CArrayKokkos<REAL_t>(num_elems, num_qpts_in_elem, "qpt_field_elem");
    qpt_vel_elem   = CArrayKokkos<REAL_t>(num_elems, num_qpts_in_elem, elem_dims, "qpt_vel_elem");
#endif
#if REMAP_OPT >= 4
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
    qpt_vol_flux = CArrayKokkos<REAL_t>(num_elems, num_qpts_in_elem, elem_dims, "qpt_vol_flux");
#endif
#if REMAP_OPT >= 5
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
#endif
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

    build_element_geometry(Mesh, FERefElem, tables, node_coords, geom,
                           num_elems, num_qpts_in_elem, num_nodes_in_elem);
    build_lumped_volume(Mesh, FERefElem, Quad, tables, node_coords, geom, elem_corner_vol,
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
            corner_field(corner_gid) = field(node_coords(node_gid,0),node_coords(node_gid,1));
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
        const auto io_t0 = std::chrono::high_resolution_clock::now();

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

        io_seconds += std::chrono::duration<double>(
                          std::chrono::high_resolution_clock::now() - io_t0).count();
    } // end graphics dump scope


    // --------------------------------------------------
    // Time integration loop
    for(size_t cycle=0; cycle<max_cycles; cycle++){
        
        if(cycle%10 == 0) printf(" time = %.4f \n", time);


        // --------------------------------------------------
        // Step 1a: Store time level n state

        timer.start();
#if REMAP_OPT >= 1
        FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

            const size_t elem_gid = idx / num_nodes_in_elem;
            const size_t node_lid = idx % num_nodes_in_elem;

            elem_corner_vol_n(elem_gid, node_lid) = elem_corner_vol(elem_gid, node_lid);

            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
            corner_field_n(corner_gid) = corner_field(corner_gid);
        });
#else
        FOR_FIRST(elem_gid, 0, num_elems, {

            FOR_SECOND(node_lid, 0, num_nodes_in_elem,{
                elem_corner_vol_n(elem_gid, node_lid) = elem_corner_vol(elem_gid, node_lid);

                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                corner_field_n(corner_gid) = corner_field(corner_gid);
            });
            
        });
#endif

        FOR_ALL(node_gid, 0, num_nodes, {
            for(size_t dim=0; dim<elem_dims; dim++){
                node_coords_n(node_gid, dim)   = node_coords(node_gid, dim);
                node_velocity_n(node_gid, dim) = node_velocity(node_gid, dim);
            }
        });
        STEP_FENCE();
        timer.stop(PH_STORE);


        // ------------------------------------------------------
        // Step 1b: get CFL time step for moving mesh

        timer.start();
#if REMAP_OPT >= 4
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
#else
        REAL_t min_h_loc;
        FOR_REDUCE_MIN(elem_gid, 0, num_elems, 
                        min_h_loc, { 
            
            for(size_t qpt_lid=0; qpt_lid<num_qpts_in_elem; qpt_lid++){

                // jacobian matrix was already calculated in corner volume vector
                ViewCArrayKokkos<REAL_t> jac(&geom.elem_jac(elem_gid,qpt_lid,0,0),3,3);

                // calculate det_J 
                REAL_t det = det_3x3(jac);

                const REAL_t vol_qpt= Quad.qpt_weights(qpt_lid)*det;
                const REAL_t h_qpt = pow(vol_qpt,0.3333333);
                if(h_qpt < min_h_loc) min_h_loc = h_qpt;
            }

        }, h_cfl);
        Kokkos::fence();
#endif
        timer.stop(PH_CFL);

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

            timer.start();
            FOR_ALL(node_gid, 0, num_nodes,{
                // new velocity, it is Taylor-Green vortex
                // PI is defined in mesh class
                node_velocity(node_gid, 0) =  sin(PI*node_coords(node_gid, 0))*cos(PI*node_coords(node_gid, 1));
                node_velocity(node_gid, 1) = -cos(PI*node_coords(node_gid, 0))*sin(PI*node_coords(node_gid, 1));
                node_velocity(node_gid, 2) = 0.0;
            });
            timer.stop(PH_VELOCITY);


            // ----------------------------------------------------------
            // Step 3: Calculate the surface fluxes at quadrature points

            timer.start();
            build_surface_flux(Mesh, RefSurf, SurfQuad, tables, node_coords, node_velocity, corner_field,
                               surf_jac, surf_qpt_qpt_map, RHS_surf_flux,
                               num_surfs, num_qpts_in_surf, num_nodes_in_elem, elem_dims);
            STEP_FENCE();
            timer.stop(PH_SURF_FLUX);


            // -------------------------------------------------
            // Step 4: Build RHS of DG equations in the element

            timer.start();
            assemble_rhs(Mesh, FERefElem, RefSurf, Quad, tables, geom,
                         corner_field, corner_field_n, node_velocity, elem_corner_vol_n,
                         RHS_surf_flux, RHS_elem,
                         qpt_field_elem, qpt_vel_elem, qpt_vol_flux,
                         rk_alpha, dt,
                         num_elems, num_qpts_in_elem, num_nodes_in_elem,
                         num_surfs_in_elem, num_qpts_in_surf, elem_dims);
            timer.stop(PH_RHS);


            // ================================================================
            // Step 5: Move the mesh to the new location
            timer.start();
            FOR_ALL(node_gid, 0, num_nodes,{
                // new position of the mesh
                node_coords(node_gid, 0) = node_coords_n(node_gid, 0) + 0.5*(node_velocity(node_gid, 0)+node_velocity_n(node_gid, 0)) * rk_alpha * dt; 
                node_coords(node_gid, 1) = node_coords_n(node_gid, 1) + 0.5*(node_velocity(node_gid, 1)+node_velocity_n(node_gid, 1)) * rk_alpha * dt;
                // z-coords never change
            });
            STEP_FENCE();
            timer.stop(PH_MOVE);


            // ================================================================
            // Step 6: build the diagonal volume matrix for nodal DG after the mesh moved
            timer.start();
            build_element_geometry(Mesh, FERefElem, tables, node_coords, geom,
                                   num_elems, num_qpts_in_elem, num_nodes_in_elem);
            build_lumped_volume(Mesh, FERefElem, Quad, tables, node_coords, geom, elem_corner_vol,
                                num_elems, num_qpts_in_elem, num_nodes_in_elem);
            STEP_FENCE();
            timer.stop(PH_MASS);


            // -----------------------------------------------------
            // 7. Solve M * u^{n+1} = RHS where M is diagonal

            timer.start();
#if REMAP_OPT >= 1
            FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

                const size_t elem_gid = idx / num_nodes_in_elem;
                const size_t dof_lid  = idx % num_nodes_in_elem;

                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
                corner_field(corner_gid) = RHS_elem(elem_gid, dof_lid)/elem_corner_vol(elem_gid, dof_lid);
            });
#else
            FOR_FIRST(elem_gid, 0, num_elems,{
        
                // -----------------------------------------------------
                // 4e. Save the new corner DOFs

                // for_all_second here
                FOR_SECOND(dof_lid, 0, num_nodes_in_elem, {
                    const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
                    corner_field(corner_gid) = RHS_elem(elem_gid, dof_lid)/elem_corner_vol(elem_gid, dof_lid);
                });

            }); // end parallel for elems
#endif
            STEP_FENCE();
            timer.stop(PH_SOLVE);

        } // end Runge Kutta time level loop


        // ================================================================
        // Step 7: update time
        time += dt;

        // Conservation Check
        timer.start();
        REAL_t sum_elem = 0.0;
        REAL_t domain_mass_time = 0.0;
        FOR_REDUCE_SUM(elem_gid, 0, num_elems, sum_elem, {

            for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                sum_elem += elem_corner_vol(elem_gid, node_lid)*corner_field(corner_gid);
            }
        }, domain_mass_time);
        timer.stop(PH_CONSERVE);

        printf("Domain mass error= %f \n", domain_mass_time-domain_mass_t0);
        if(fabs(domain_mass_time-domain_mass_t0)>1.e-12) Kokkos::abort("ERROR: Mass is not conserved");


        // ================================================================
        // Step 8: write outputs
        if( time-time_output >= -1.e-8 ){

            timer.start();
            const auto io_t0 = std::chrono::high_resolution_clock::now();

            //// L1 and L2 error norms
#if REMAP_OPT >= 5
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

                    const REAL_t diff = val_qpt - field(x_qpt,y_qpt);
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
#else
            REAL_t L1;
            REAL_t L1_lcl;
            FOR_REDUCE_SUM(elem_gid, 0, num_elems,  L1_lcl, {

                for(size_t qpt_lid=0; qpt_lid<num_qpts_in_elem; qpt_lid++){

                        // volume contribution from qpt
                        const REAL_t vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
                        
                        REAL_t val_qpt = 0.0;
                        REAL_t x_qpt   = 0.0;
                        REAL_t y_qpt   = 0.0;

                        for(size_t corner_lid=0; corner_lid<num_nodes_in_elem; corner_lid++) {
                            
                            // remmeber node_lid = corner_lid
                            const size_t node_gid   = Mesh.nodes_in_elem(elem_gid,corner_lid);
                            const size_t corner_gid = Mesh.corners_in_elem(elem_gid, corner_lid);
                            
                            val_qpt += corner_field(corner_gid)*FERefElem.qpt_basis(qpt_lid,corner_lid); 
                            x_qpt   += node_coords(node_gid,0)*FERefElem.qpt_basis(qpt_lid,corner_lid);
                            y_qpt   += node_coords(node_gid,1)*FERefElem.qpt_basis(qpt_lid,corner_lid);
                                
                        } // loop over corners of the element

                        L1_lcl += fabs(val_qpt - field(x_qpt,y_qpt))*vol_qpt;
                        
                } // end for dof_lid and qpt_lid

            }, L1); // end parallel for

            REAL_t L2;
            REAL_t L2_lcl;
            FOR_REDUCE_SUM(elem_gid, 0, num_elems, L2_lcl, {

                for(size_t qpt_lid=0; qpt_lid<num_qpts_in_elem; qpt_lid++){

                        // volume contribution from qpt
                        const REAL_t vol_qpt = elem_det_jac(elem_gid, qpt_lid)*Quad.qpt_weights(qpt_lid);
                        
                        REAL_t val_qpt = 0.0;
                        REAL_t x_qpt   = 0.0;
                        REAL_t y_qpt   = 0.0;

                        for(size_t corner_lid=0; corner_lid<num_nodes_in_elem; corner_lid++) {
                            
                            // remmeber node_lid = corner_lid
                            const size_t node_gid   = Mesh.nodes_in_elem(elem_gid,corner_lid);
                            const size_t corner_gid = Mesh.corners_in_elem(elem_gid,corner_lid);
                            
                            val_qpt += corner_field(corner_gid)*FERefElem.qpt_basis(qpt_lid,corner_lid); 
                            x_qpt   += node_coords(node_gid,0)*FERefElem.qpt_basis(qpt_lid,corner_lid);
                            y_qpt   += node_coords(node_gid,1)*FERefElem.qpt_basis(qpt_lid,corner_lid);

                        } // loop over corners of the element

                        L2_lcl += (val_qpt - field(x_qpt,y_qpt))*(val_qpt - field(x_qpt,y_qpt))*vol_qpt;
                        
                } // end for dof_lid and qpt_lid

            }, L2); // end parallel for
            L2 = sqrt(L2);
#endif

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

            io_seconds += std::chrono::duration<double>(
                              std::chrono::high_resolution_clock::now() - io_t0).count();
            timer.stop(PH_OUTPUT);

        } // end if

        if (time >= max_time  ){
            printf("Domain mass at time=%f: %f \n", time, domain_mass_time);
            break;
        }

    } // end loop over cycle
    printf(" time = %.4f ", time);

    err_file.close();

    printf("\n Remap test finished.\n");

    Kokkos::fence();
    const double wall_seconds = std::chrono::duration<double>(
                                    std::chrono::high_resolution_clock::now() - wall_t0).count();
    timer.report();
    printf("\n==== REMAP_OPT=%d : %.3f s inside MATAR scope (%.3f s of that is graphics/norms) ====\n",
           REMAP_OPT, wall_seconds, io_seconds);


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

#if REMAP_OPT >= 4
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
#else
    for(size_t elem_gid=0; elem_gid<num_elems; elem_gid++){
    
        // loop over structured nodes in this element
        FOR_ALL(k,0,num_DOFs_1d, 
                j,0,num_DOFs_1d, 
                i,0,num_DOFs_1d, {
                    
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
            size_t node_lcl = i + (j + k*num_DOFs_1d)*num_DOFs_1d;
            size_t node_gid = nodes_in_elem(elem_gid, node_lcl);
            node_coords_uniform(node_gid,0) = x;
            node_coords_uniform(node_gid,1) = y;
            node_coords_uniform(node_gid,2) = z;
        
        }); // end parallel for
        Kokkos::fence();
    }
#endif
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
