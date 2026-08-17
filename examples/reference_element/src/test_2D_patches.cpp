
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

using namespace mtr;
using namespace swage;    // unstructured mesh and point cloud
using namespace elements; // reference element space

void verify_patch_basics(const Mesh_t& mesh) {
    
    // Check 1: All patch nodes exist and are valid
    FOR_ALL(patch_gid, 0, mesh.num_patches, {
        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++) {
            size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);
            
            // Node ID should be valid
            assert(node_gid < mesh.num_nodes);
            
            // No duplicate nodes in same patch
            for (size_t other_lid = node_lid + 1; other_lid < mesh.num_nodes_in_patch; other_lid++) {
                assert(mesh.nodes_in_patch(patch_gid, other_lid) != node_gid);
            }
        }
    });
    
    printf("[Yes] Basic patch node validation passed\n");
}


void verify_patches_in_surfaces(const Mesh_t& mesh) {
    
    FOR_ALL(patch_gid, 0, mesh.num_patches, {
        size_t surf_gid = mesh.surf_in_patch(patch_gid);
        
        // Each patch node must be in the surface
        for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; patch_node_lid++) {
            size_t patch_node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
            
            bool found = false;
            for (size_t surf_node_lid = 0; surf_node_lid < mesh.num_nodes_in_surf; surf_node_lid++) {
                if (mesh.nodes_in_surf(surf_gid, surf_node_lid) == patch_node_gid) {
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                printf("ERROR: Patch %d node %zu (gid=%zu) not found in surface %zu\n",
                       patch_gid, patch_node_lid, patch_node_gid, surf_gid);
            }
            assert(found);
        }
    });
    
    printf("[Yes] All patch nodes belong to their parent surface\n");
}


void verify_patch_orientation(const Mesh_t& mesh, 
                              const DCArrayKokkos<double>& node_coords) {
    
    if (mesh.num_dims == 3) {

        FOR_ALL(patch_gid, 0, mesh.num_patches, {
            // Get the 4 nodes of the patch (for 3D)
            size_t n0 = mesh.nodes_in_patch(patch_gid, 0);
            size_t n1 = mesh.nodes_in_patch(patch_gid, 1);
            size_t n2 = mesh.nodes_in_patch(patch_gid, 2);
            size_t n3 = mesh.nodes_in_patch(patch_gid, 3);
            
            // Calculate two edge vectors
            double v1[3];
            double v2[3];
            for (int d = 0; d < 3; d++) {
                v1[d] = node_coords(1, n1, d) - node_coords(1, n0, d);
                v2[d] = node_coords(1, n3, d) - node_coords(1, n0, d);
            }
            
            // Cross product gives normal
            double normal[3];
            normal[0] = v1[1]*v2[2] - v1[2]*v2[1];
            normal[1] = v1[2]*v2[0] - v1[0]*v2[2];
            normal[2] = v1[0]*v2[1] - v1[1]*v2[0];
            
            // Get element center to determine inward/outward
            size_t elem_gid = mesh.elems_in_patch(patch_gid, 0);
            double elem_center[3];
            elem_center[0] = 0;
            elem_center[1] = 0; 
            elem_center[2] = 0;
            
            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                for (int d = 0; d < 3; d++) {
                    elem_center[d] += node_coords(1, node_gid, d);
                }
            }
            for (int d = 0; d < 3; d++) {
                elem_center[d] /= mesh.num_nodes_in_elem;
            }
            
            // Patch center
            double patch_center[3];
            patch_center[0]= 0;
            patch_center[1]= 0; 
            patch_center[2]= 0;
            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++) {
                size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);
                for (int d = 0; d < 3; d++) {
                    patch_center[d] += node_coords(1, node_gid, d);
                }
            }
            for (int d = 0; d < 3; d++) {
                patch_center[d] /= mesh.num_nodes_in_patch;
            }
            
            // Vector from elem center to patch center
            double outward[3];
            for (int d = 0; d < 3; d++) {
                outward[d] = patch_center[d] - elem_center[d];
            }
            
            // Dot product should be positive (outward normal)
            double dot = normal[0]*outward[0] + normal[1]*outward[1] + normal[2]*outward[2];
            
            if (dot < 0) {
                printf("ERROR: Patch %d has inward-facing normal (dot=%f)\n", patch_gid, dot);
            }
            assert(dot > -1e-12); // Allow small numerical error
        });
    }
    
    printf("[Yes] Patch orientations follow right-hand outward normal convention\n");
}



void verify_internal_patch_continuity(const Mesh_t& mesh) {
    
    DCArrayKokkos <size_t> errors(1); 
    errors.set_values(0);
    
    FOR_ALL(patch_gid, 0, mesh.num_patches, {
        size_t surf_gid = mesh.surf_in_patch(patch_gid);
        
        // Only check internal surfaces (have 2 elements)
        if (mesh.num_elems_in_surf(surf_gid) == 2) {
            
            size_t elem0 = mesh.elems_in_patch(patch_gid, 0);
            int elem1 = mesh.elems_in_patch(patch_gid, 1);
            
            // elem1 should be valid (non-negative)
            if (elem1 < 0) {
                printf("ERROR: Internal patch %d has negative neighbor element\n", patch_gid);
                Kokkos::atomic_add(&errors(0),1);
            }
            
            // Both elements should contain all patch nodes
            for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; patch_node_lid++) {
                size_t patch_node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
                
                bool found_in_elem0 = false;
                bool found_in_elem1 = false;
                
                for (size_t elem_node_lid = 0; elem_node_lid < mesh.num_nodes_in_elem; elem_node_lid++) {
                    if (mesh.nodes_in_elem(elem0, elem_node_lid) == patch_node_gid) {
                        found_in_elem0 = true;
                    }
                    if (elem1 >= 0 && mesh.nodes_in_elem(elem1, elem_node_lid) == patch_node_gid) {
                        found_in_elem1 = true;
                    }
                }
                
                if (!found_in_elem0 || !found_in_elem1) {
                    printf("ERROR: Patch %d node %zu not in both elements\n", 
                           patch_gid, patch_node_gid);
                    Kokkos::atomic_add(&errors(0),1);
                }
            }
        }
    });
    Kokkos::fence();
    
    errors.update_host();
    if (errors.host(0) == 0) {
        printf("[Yes] Internal patch continuity verified\n");
    } else {
        printf("[No] Found %zu patch continuity errors\n", errors.host(0));
    }
}



void verify_boundary_patches(const Mesh_t& mesh) {
    
    FOR_ALL(bdy_patch_gid, 0, mesh.num_bdy_patches, {
        size_t patch_gid = mesh.bdy_patches(bdy_patch_gid);
        size_t surf_gid = mesh.surf_in_patch(patch_gid);
        
        // Boundary patches should have only 1 element
        assert(mesh.num_elems_in_surf(surf_gid) == 1);
        
        // Second element should be negative
        int elem1 = mesh.elems_in_patch(patch_gid, 1);
        if (elem1 >= 0) {
            printf("ERROR: Boundary patch %zu has valid second element %d\n", 
                   patch_gid, elem1);
        }
        assert(elem1 < 0);
    });
    
    printf("[Yes] Boundary patch validation passed\n");
}


void print_patch_diagnostics(const Mesh_t& mesh) {

    printf("\n=== Patch Diagnostics ===\n");
    printf("Total patches: %zu\n", mesh.num_patches);
    printf("Boundary patches: %zu\n", mesh.num_bdy_patches);
    printf("Internal patches: %zu\n", mesh.num_patches - mesh.num_bdy_patches);
    printf("Patches per surface: %zu\n", mesh.num_patches_in_surf);
    printf("Nodes per patch: %zu\n", mesh.num_nodes_in_patch);
    printf("Pn order: %zu\n", mesh.Pn);
    
    // Sample a few patches
    printf("\nSample patches:\n");
    RUN({
        for (size_t i = 0; i < mesh.num_patches; i++) {
            printf("Patch %zu: surf_in_patch=%zu, elems_in_patch=[%d,%d], nodes=[",
                i, mesh.surf_in_patch(i),
                mesh.elems_in_patch(i,0), mesh.elems_in_patch(i,1));
            for (size_t n = 0; n < mesh.num_nodes_in_patch; n++) {
                printf("%zu%s", mesh.nodes_in_patch(i,n),
                    n < mesh.num_nodes_in_patch-1 ? "," : "");
            }
            printf("]\n");
        }
    });
}


void verify_all_patch_connectivity(const Mesh_t& mesh, 
                                   const DCArrayKokkos<double>& node_coords) {
    printf("\n=== Starting Patch Connectivity Verification ===\n\n");
    
    verify_patch_basics(mesh);
    verify_patches_in_surfaces(mesh);
    verify_patch_orientation(mesh, node_coords);
    verify_internal_patch_continuity(mesh);
    verify_boundary_patches(mesh);
    print_patch_diagnostics(mesh);
    
    printf("\n=== All Patch Verification Tests Passed! ===\n\n");
}







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

        const size_t elem_dims = 2;  // MUST be 2D!!!!
        const size_t elem_order = 1; // Must be 1!!!! 
        const size_t num_elems_1D = 1;

        const size_t max_cycles = 1000;
        const double max_time   = 0.5;
        const double graphics_dt = 0.1;


        // ================================================================
        // Create quadrature along with the reference element and surface

        std::cout<<"Building classical P1 2D element \n";

        if(elem_dims!=2) Kokkos::abort("ERROR: elem_dims must be equal to 2\n");


        // ==========================================
        // Build the unstructured mesh structure

        std::cout<<"Building unstructured mesh \n";

        const size_t num_elems    = num_elems_1D*num_elems_1D;
        const size_t num_nodes_1D = elem_order*num_elems_1D + 1;  // number of nodes
        const size_t num_nodes    = num_nodes_1D*num_nodes_1D;

        Mesh.initialize_dims(elem_dims);
        Mesh.initialize_elems(num_elems);  // NOTE: no Pn, no Guass points, etc. It's classic element
        Mesh.initialize_nodes(num_nodes);
        
        DCArrayKokkos <double> node_coords(Mesh.num_nodes, Mesh.num_dims);
        const double h = 1.0/((double)num_nodes_1D-1);

        // create indexing for a Pn order mesh

        // Step 1: Initialize ALL node coordinates once (no race condition)
        FOR_ALL(jc, 0, num_nodes_1D,
                ic, 0, num_nodes_1D, {
            
            size_t node_gid = ic + jc*num_nodes_1D;
            
            node_coords(node_gid, 0) = (double)ic * h;
            node_coords(node_gid, 1) = (double)jc * h;
        });

        const size_t ensight_ordering[4] = {0, 1, 3, 2};

        // Step 2: Build element connectivity
        FOR_ALL(i, 0, num_elems_1D,
                j, 0, num_elems_1D, {
            
            size_t elem_gid = i + j*num_elems_1D;
            size_t node_lid = 0;
            
            // Create nodes in the order: (ic,jc) pairs
            // (0,0), (1,0), (0,1), (1,1) which gives global IDs 0,1,2,3
            for(size_t jc = j; jc <= j + elem_order; jc++)
            for(size_t ic = i; ic <= i + elem_order; ic++){
                size_t node_gid = ic + jc*num_nodes_1D;
                
                // Map to element node position based on convention
                // Global layout: 2-3    Element layout: 3-2
                //                0-1                    0-1
                size_t elem_node_lid;
                if (jc == j) {
                    // Bottom row: natural order
                    elem_node_lid = (ic - i);  // 0 or 1
                } else {
                    // Top row: reversed order  
                    elem_node_lid = 3 - (ic - i);  // 3 or 2
                }
                
                Mesh.nodes_in_elem(elem_gid, elem_node_lid) = node_gid;
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
        std::cout<<"build_node_node_connectivity \n";
        Mesh.build_node_node_connectivity();  // not used in this test, though


        verify_all_patch_connectivity(Mesh, node_coords);

    } // end MATAR scope

    MATAR_FINALIZE();
    
    return 1;

} // end main