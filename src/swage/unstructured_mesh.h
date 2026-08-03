/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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
#ifndef UNSTRUCTURED_MESH_H
#define UNSTRUCTURED_MESH_H

#include "matar.h"
#include "indexing_utils.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

namespace mesh_init
{
    // element mesh types
    enum ElementNameType
    {
        linearTensorElement = 1,   // single quadrature point element
        arbitraryTensorElement = 2 // fully integrated arbitrary-order element
    };

    // other enums could go here on the mesh
} // end namespace

namespace swage
{


////////////////////////////////////////////////////////////////////////////////////
///
/// \fn Mesh_t
///
/// \brief Builds and stores mesh sizes and connectivity data structures for 
///        arbitrary-order unstructured meshes in 2D/3D
///
///  Mesh entity definitions:
///      Element: Aribtrary-order hexahedral or quadralateral volume
///      Zone:    A discretization of an element by subdividing it using the nodes 
///               The zone has 8 nodes (3D) or 4 nodes (2D) for any order mesh
///      Node:    A kinematic degree of freedom
///      Corner:  A element-node pair
///      Surface: The surface of the element, it is one dimension lower than the volume
///      Patch:   A discretization of a surface by subdividing it using the nodes
///      Face:    The local surface entity of the Element, equal to 6 (3D) or 4 (2D)
///      Side:    A element-surface pair -- not in the mesh type at this time
///
///////////////////////////////////////////////////////////////////////////////////
struct Mesh_t
{

    bool verbose = false;

    // ---- Global Mesh Definitions ---- //
    mesh_init::ElementNameType elem_kind = mesh_init::linearTensorElement; ///< The type of elements used in the mesh

    size_t Pn = 1;       ///< Polynomial order of kinematic space defining element
    size_t num_dims = 0; ///< Number of spatial dimension


    // ---- Element Data Definitions ---- //
    size_t num_elems = 0;           ///< Number of elements in the mesh
    size_t num_nodes_in_elem   = 0; ///< Number of nodes in an element
    size_t num_patches_in_elem = 0; ///< Number of patches in an element
    size_t num_surfs_in_elem   = 0; ///< Number of surfaces in an element
    size_t num_zones_in_elem   = 0; ///< Number of zones in an element

    size_t num_gauss_in_elem   = 0; ///< Number of Gauss points in an element

    DCArrayKokkos<size_t> nodes_in_elem; ///< Nodes in an element
    corners_in_elem_t     corners_in_elem;

    RaggedRightArrayKokkos<size_t> elems_in_elem; ///< Elements connected to an element
    CArrayKokkos<size_t> num_elems_in_elem;       ///< Number of elements connected to an element

    CArrayKokkos<size_t> patches_in_elem; ///< Patches in an element (including internal patches)
    CArrayKokkos<size_t> surfs_in_elem;   ///< Surfaces on an element

    zones_in_elem_t   zones_in_elem;   ///< Zones in an element
    gauss_in_elem_t   gauss_in_elem;   ///< Gauss points in an element


    // ---- Node Data Definitions ---- //
    size_t num_nodes = 0; ///< Number of nodes in the mesh

    RaggedRightArrayKokkos<size_t> corners_in_node; ///< Corners connected to a node
    CArrayKokkos<size_t> num_corners_in_node;       ///< Number of corners connected to a node
    RaggedRightArrayKokkos<size_t> elems_in_node;   ///< Elements connected to a given node
    RaggedRightArrayKokkos<size_t> nodes_in_node;   ///< Nodes connected to a node along an edge
    CArrayKokkos<size_t> num_nodes_in_node;         ///< Number of nodes connected to a node along an edge


    // ---- Surface Data Definitions ---- //
    size_t num_surfs = 0;           ///< Number of surfaces in the mesh
    size_t num_nodes_in_surf = 0;   ///< Number of nodes in a surface
    size_t num_patches_in_surf = 0; ///< Number of patches in a surface

    patches_in_surf_t    patches_in_surf;   ///< Patches in a surface
    CArrayKokkos<size_t> nodes_in_surf;     ///< Nodes in a surface
    CArrayKokkos<int>    elems_in_surf;     ///< Elements connected to a surface
    CArrayKokkos<size_t> num_elems_in_surf; ///<Number of elements in this surface
    CArrayKokkos<int>    faces_in_surf;     ///< Local face index of the element 


    // ---- Patch Data Definitions ---- //
    size_t num_patches = 0;         ///< Number of patches in the mesh
    size_t num_nodes_in_patch = 0;  ///< Number of nodes in a patch

    CArrayKokkos<size_t> nodes_in_patch; ///< Nodes connected to a patch
    CArrayKokkos<int>    elems_in_patch; ///< Elements connected to a patch
    CArrayKokkos<size_t> surf_in_patch;  ///< the surface the patch belongs to


    // ---- Corner Data Definitions ---- //
    size_t num_corners = 0; ///< Number of corners (define) in the mesh


    // ---- Zone Data Definitions ---- //
    size_t num_zones = 0;           ///< Number of zones in the mesh
    size_t num_nodes_in_zone = 0;   ///< Number of nodes in a zone

    CArrayKokkos<size_t> nodes_in_zone; ///< Nodes defining a zone


    // ---- Boundary Data Definitions ---- //

    size_t num_bdy_surfs = 0;   ///< Number of boundary surfaces
    size_t num_bdy_patches = 0; ///< Number of boundary patches
    size_t num_bdy_nodes = 0;   ///< Number of boundary nodes
    size_t num_bdy_sets = 0;    ///< Number of boundary sets
    
    CArrayKokkos<size_t> bdy_surfs;   ///< Boundary patches
    CArrayKokkos<size_t> bdy_patches; ///< Boundary patches
    CArrayKokkos<size_t> bdy_nodes;   ///< Boundary nodes

    RaggedRightArrayKokkos<size_t> bdy_patches_in_set; ///< Boundary patches in a boundary set
    DCArrayKokkos<size_t> num_bdy_patches_in_set;      ///< Number of boundary nodes in a set

    RaggedRightArrayKokkos<size_t> bdy_nodes_in_set; ///< Boundary nodes in a boundary set
    DCArrayKokkos<size_t> num_bdy_nodes_in_set; ///< Number of boundary nodes in a set


    // ---- Internal Condition Data Definitions ---- //
    size_t num_internal_sets = 0; ///< Number of internal sets

    RaggedRightArrayKokkos<size_t> internal_nodes_in_set; ///< Internal nodes in an internal set
    DCArrayKokkos<size_t> num_internal_nodes_in_set;      ///< Number of internal nodes in a set


    // MPI Decomposition Data Definitions ---- //
    DCArrayKokkos<size_t> local_to_global_node_mapping; ///< Local to global node mapping
    DCArrayKokkos<size_t> local_to_global_elem_mapping; ///< Local to global element mapping

    // Element communicaiton data definitions
    size_t num_owned_elems = 0;    ///< Number of owned elements on this rank
    size_t num_boundary_elems = 0; ///< Number of boundary elements on this rank (send data to neighboring MPI ranks)
    DCArrayKokkos<size_t> boundary_elem_local_ids; ///< Local IDs of boundary elements on this rank (send data to neighboring MPI ranks)
    size_t num_ghost_elems = 0; ///< Number of ghost elements on this rank (receive data from neighboring MPI ranks)
    
    // Node communicaiton data definitions
    size_t num_owned_nodes = 0;    ///< Number of owned nodes on this rank
    size_t num_boundary_nodes = 0; ///< Number of boundary nodes on this rank (send data to neighboring MPI ranks)
    DCArrayKokkos<bool> shared_tally_owned_nodes; ///< Owned-node mask: true where this rank is the min MPI rank among ranks that own the global node (domain tally contributor); length num_owned_nodes
    // DCArrayKokkos<size_t> boundary_node_local_ids; ///< Local IDs of boundary nodes on this rank (send data to neighboring MPI ranks)
    size_t num_ghost_nodes = 0; ///< Number of ghost nodes on this rank (receive data from neighboring MPI ranks)
    

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_dims
    ///
    /// \brief Set up mesh dimensions
    ///
    /// \param elem_dims_in The number dimensions 
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_dims(const size_t num_dims_inp)
    {
        num_dims  = num_dims_inp;
    }; // end method


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_nodes
    ///
    /// \brief Set the number of nodes in the mesh
    ///
    /// \param num_nodes_inp The number dimensions 
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_nodes(const size_t num_nodes_inp)
    { 
        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at initialize_nodes().");
        }
        num_nodes = num_nodes_inp;
        return;
    }; // end method


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_elems
    ///
    /// \brief Setup the elemenets in a mesh comprised of linear edges and one
    ///        quadrature point per element
    ///
    /// \param num_elems_inp The number elems
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_elems(const size_t num_elems_inp)
    {

        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at initialize_elems().");
        }

        // --- Basic element bookkeeping ---
        num_elems = num_elems_inp;

        // initializes a linear element with a single gauss point for saving results
        Pn = 1;

        // --- Derived sizes ---
        num_nodes_in_elem = (size_t)std::pow(2, num_dims); 
        num_nodes_in_zone = (size_t)std::pow(2, num_dims); // (4, or 8, always)
        num_gauss_in_elem = 1;  // 1 Gauss point per element
        num_zones_in_elem = 1;  // 1 zone per element
        num_surfs_in_elem = num_dims == 2 ? 4 : 6; // 4 or 6 (always)
        num_zones = num_zones_in_elem * num_elems;

        num_corners = num_nodes_in_elem*num_elems;

        
        // --- Allocations ---
        nodes_in_elem    = DCArrayKokkos<size_t>(num_elems, num_nodes_in_elem, "mesh.nodes_in_elem");
        corners_in_elem  = corners_in_elem_t(num_nodes_in_elem); 
        gauss_in_elem    = gauss_in_elem_t(num_gauss_in_elem);
        zones_in_elem    = zones_in_elem_t(num_zones_in_elem);
        surfs_in_elem    = CArrayKokkos<size_t>(num_elems, num_surfs_in_elem, "mesh.surfs_in_zone");

        return;
    }; // end method

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_elems_Pn
    ///
    /// \brief Setup the elements in an arbitrary-order mesh
    ///
    /// \param num_elems_inp The number elems
    /// \param elem_Pn_order The element order, where num_nodes_1D = Pn+1
    /// \param num_gauss_1D  The number of gauss points in each direction
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_elems_Pn(
        const size_t num_elems_inp,
        const size_t elem_Pn_order,
        const size_t num_gauss_1D)
    {

        // Note: num_gauss_1D creates an index space that can be used to register state on

        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at initialize_elems_Pn().");
        }

        // --- Set element details ---
        Pn = elem_Pn_order; // Note: element Pn_order = dofs_1D-1, where dofs are the element nodes
        if (Pn == 0) {
            Kokkos::abort("Error: Pn must be greater than 0. Exiting at initialize_elems_Pn().");
        }

        num_elems = num_elems_inp;  
        if (num_elems == 0) {
            Kokkos::abort("Error: num_elems must be greater than 0. Exiting at initialize_elems_Pn().");
        }

        elem_kind = mesh_init::arbitraryTensorElement;


        // --- Derived sizes ---
        num_gauss_in_elem = (size_t)std::pow(num_gauss_1D, num_dims); // Note: 2*Pn with Legendre is needed for solids mechanics
        num_nodes_in_elem = (size_t)std::pow(Pn + 1, num_dims);
        num_nodes_in_zone = (size_t)std::pow(2, num_dims);   // (4, or 8, always)
        num_zones_in_elem = (size_t)std::pow(Pn, num_dims);  // Pn^dim
        num_surfs_in_elem = num_dims == 2 ? 4 : 6;           // 4 or 6 (always)
        num_zones = num_zones_in_elem * num_elems;        
        
        num_corners = num_nodes_in_elem*num_elems;

        // --- Allocations ---
        nodes_in_elem    = DCArrayKokkos<size_t>(num_elems, num_nodes_in_elem, "mesh.nodes_in_elem");
        corners_in_elem  = corners_in_elem_t(num_nodes_in_elem); 
        zones_in_elem    = zones_in_elem_t(num_zones_in_elem);
        surfs_in_elem    = CArrayKokkos<size_t>(num_elems, num_surfs_in_elem, "mesh.surfs_in_zone");
        nodes_in_zone    = CArrayKokkos<size_t>(num_zones, num_nodes_in_zone, "mesh.nodes_in_zone");
        gauss_in_elem    = gauss_in_elem_t(num_gauss_in_elem);

        return;
    }; // end method


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_zones
    ///
    /// \brief Using the nodes in the element, decompose the element into zones.
    ///        
    /// The zones have 8 nodes (3D) or 4 nodes (2D) for all element orders. The
    /// zones are built after the mesh dims are set, and the nodes and elements 
    /// are initialized.
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_zones(){

        const size_t num_1D = Pn + 1;

        if(num_zones==0) Kokkos::abort("ERROR: Elements must be initialized prior to creating zones.");

        if(num_dims==3){

            FOR_ALL_CLASS(elem_gid, 0, num_elems, {
                size_t node_lids[8]; // temp storage for local node ids
                for (size_t k = 0; k < num_1D-1; k++) 
                for (size_t j = 0; j < num_1D-1; j++) 
                for (size_t i = 0; i < num_1D-1; i++) {
                    node_lids[0] = i + j * (num_1D) + k * (num_1D) * (num_1D);                 // i,j,k
                    node_lids[1] = i + 1 + j * (num_1D) + k * (num_1D) * (num_1D);             // i+1, j, k
                    node_lids[2] = i + (j + 1) * (num_1D) + k * (num_1D) * (num_1D);           // i,j+1,k
                    node_lids[3] = i + 1 + (j + 1) * (num_1D) + k * (num_1D) * (num_1D);       // i+1, j+1, k
                    node_lids[4] = i + j * (num_1D) + (k + 1) * (num_1D) * (num_1D);           // i, j , k+1
                    node_lids[5] = i + 1 + j * (num_1D) + (k + 1) * (num_1D) * (num_1D);       // i + 1, j , k+1
                    node_lids[6] = i + (j + 1) * (num_1D) + (k + 1) * (num_1D) * (num_1D);     // i,j+1,k+1
                    node_lids[7] = i + 1 + (j + 1) * (num_1D) + (k + 1) * (num_1D) * (num_1D); // i+1, j+1, k+1

                    size_t zone_lid = i + j * (num_1D - 1) + k * (num_1D - 1) * (num_1D - 1);
                    size_t zone_gid = zones_in_elem(elem_gid, zone_lid);

                    for (size_t node_lid = 0; node_lid < 8; node_lid++) {
                        // get global id for the node
                        size_t node_gid = nodes_in_elem(elem_gid, node_lids[node_lid]);
                        nodes_in_zone(zone_gid, node_lid) = node_gid;
                    }
                } // end for
            }); // end FOR_ALL elem_gid

        } 
        else if (num_dims==2){

            FOR_ALL_CLASS(elem_gid, 0, num_elems, {
                size_t node_lids[4]; // temp storage for local node ids
                for (size_t j = 0; j < num_1D-1; j++) 
                for (size_t i = 0; i < num_1D-1; i++) {
                    node_lids[0] = i + j * (num_1D);            // i,   j
                    node_lids[1] = i + 1 + j * (num_1D);        // i+1, j
                    node_lids[2] = i + (j + 1) * (num_1D);      // i,   j+1
                    node_lids[3] = i + 1 + (j + 1) * (num_1D);  // i+1, j+1

                    size_t zone_lid = i + j * (num_1D - 1);
                    size_t zone_gid = zones_in_elem(elem_gid, zone_lid);

                    for (size_t node_lid = 0; node_lid < 4; node_lid++) {
                        // get global id for the node
                        size_t node_gid = nodes_in_elem(elem_gid, node_lids[node_lid]);
                        nodes_in_zone(zone_gid, node_lid) = node_gid;
                    }
                } // end for
            }); // end FOR_ALL elem_gid
            
        }
        else 
        {
            Kokkos::abort("ERROR: incorrect mesh dimensions- only 2D and 3D are supported. Wow! How did this happen?");
        } // end if

    } // end build zones


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_corner_connectivity
    ///
    /// \brief Build the corner mesh connectivity arrays
    ///        
    /// The corner is defined as an element node pair.  For high-order elements,
    /// corners exist inside the element.
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_corner_connectivity()
    {
        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at build_corner_connectivity().");
        }
        num_corners_in_node = CArrayKokkos<size_t>(num_nodes, "mesh.num_corners_in_node"); // stride sizes

        // initializing the number of corners (node-cell pair) to be zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            num_corners_in_node(node_gid) = 0;
        });

        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // increment the number of corners attached to this point
                num_corners_in_node(node_gid) = num_corners_in_node(node_gid) + 1;
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid

        // the stride sizes are the num_corners_in_node at the node
        corners_in_node = RaggedRightArrayKokkos<size_t>(num_corners_in_node, "mesh.corners_in_node");

        CArrayKokkos<size_t> count_saved_corners_in_node(num_nodes, "count_saved_corners_in_node");

        // reset num_corners to zero
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            count_saved_corners_in_node(node_gid) = 0;
        });

        // the elems_in_elem data type
        elems_in_node = RaggedRightArrayKokkos<size_t>(num_corners_in_node, "mesh.elems_in_node");

        // populate the elements connected to a node list and corners in a node
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // the column index is the num corners saved
                size_t j = count_saved_corners_in_node(node_gid);

                // Save corner index to this node_gid
                size_t corner_gid = node_lid + elem_gid * num_nodes_in_elem;  // this can be a functor
                corners_in_node(node_gid, j) = corner_gid;

                elems_in_node(node_gid, j) = elem_gid; // save the elem_gid

                // increment the number of corners saved to this node_gid
                count_saved_corners_in_node(node_gid) = count_saved_corners_in_node(node_gid) + 1;
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid

        return;
    } // end of build_corner_connectivity


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_elem_elem_connectivity
    ///
    /// \brief Build the neighboring element to element connectivity 
    ///        
    /// The elements surrounding an element accounts for all adjacent elements
    /// sharing a common node. 
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_elem_elem_connectivity()
    {
        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at build_elem_elem_connectivity().");
        }
        if (num_corners_in_node.size() == 0) {
            Kokkos::abort("Error: build_corner_connectivity must be called first. Exiting at build_elem_elem_connectivity().");
        }
   
        // find the max number of elems around a node
        size_t max_num_elems_in_node;
        size_t max_num_lcl;
        FOR_REDUCE_MAX_CLASS(node_gid, 0, num_nodes, max_num_lcl, {
            // num_corners_in_node = num_elems_in_node
            size_t max_num = num_corners_in_node(node_gid);

            if (max_num > max_num_lcl) {
                max_num_lcl = max_num;
            }
        }, max_num_elems_in_node); // end parallel reduction on max
        Kokkos::fence();

        // a temporary ragged array to save the elems around an elem
        DynamicRaggedRightArrayKokkos<size_t> temp_elems_in_elem(num_nodes, num_nodes_in_elem * max_num_elems_in_node, "temp_elems_in_elem");

        num_elems_in_elem = CArrayKokkos<size_t>(num_elems, "mesh.num_elems_in_elem");
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            num_elems_in_elem(elem_gid) = 0;
        });
        Kokkos::fence();

        // find and save neighboring elem_gids of an elem
        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                // get the gid for the node
                size_t node_id = nodes_in_elem(elem_gid, node_lid);

                // loop over all elems connected to node_gid
                for (int elem_lid = 0; elem_lid < num_corners_in_node(node_id); elem_lid++) {
                    // get the global id for the neighboring elem
                    size_t neighbor_elem_gid = elems_in_node(node_id, elem_lid);

                    // a flag to save (=1) or not (=0)
                    size_t save = 1;

                    // a true neighbor_elem_id is not equal to elem_gid
                    if (neighbor_elem_gid == elem_gid) {
                        save = 0;  // don't save
                    } // end if

                    // check to see if the neighbor_elem_gid has been saved already
                    size_t num_saved = temp_elems_in_elem.stride(elem_gid);
                    for (size_t i = 0; i < num_saved; i++) {
                        if (neighbor_elem_gid == temp_elems_in_elem(elem_gid, i)) {
                            save = 0;   // don't save, it has been saved already
                        } // end if
                    } // end for i

                    if (save == 1) {
                        // increment the number of neighboring elements saved
                        temp_elems_in_elem.stride(elem_gid)++;

                        // save the neighboring elem_gid
                        temp_elems_in_elem(elem_gid, num_saved) = neighbor_elem_gid;
                    } // end if save
                } // end for elem_lid in a node
            }  // end for node_lid in an elem

            // save the actial stride size
            num_elems_in_elem(elem_gid) = temp_elems_in_elem.stride(elem_gid);
        }); // end FOR_ALL elems
        Kokkos::fence();

        // compress out the extra space in the temp_elems_in_elem
        elems_in_elem = RaggedRightArrayKokkos<size_t>(num_elems_in_elem, "mesh.elems_in_elem");

        FOR_ALL_CLASS(elem_gid, 0, num_elems, {
            for (size_t i = 0; i < num_elems_in_elem(elem_gid); i++) {
                elems_in_elem(elem_gid, i) = temp_elems_in_elem(elem_gid, i);
            } // end for i
        });  // end FOR_ALL elems
        Kokkos::fence();

        return;
    } // end of build_elem_elem_connectivity


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_surf_connectivity
    ///
    /// \brief Build the surfaces and decompose them into patches
    ///        
    /// A surface separates two face-adjacent elements, where the shared nodes 
    /// define the surface.  On mesh boundaries, there is no face-adjacent 
    /// element.  This rountine calculates and populates boundary index spaces.  
    //  Each surface is decomposed into patches using the nodes on the surface. 
    /// Patches have 4 nodes (3D) or 2 nodes (2D) for all element orders.  
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_surf_connectivity() {

        if(num_dims==0) Kokkos::abort("Error: mesh.num_dims is not set. Exiting at build_surf_connectivity \n");
        if(num_elems_in_elem.size()==0) Kokkos::abort("Error: build_elem_elem_connectivity must be called first.  Existing build_surf_connectivity(). \n");
        if(num_corners_in_node.size() == 0) Kokkos::abort("Error: build_corner_connectivity must be called first. Exiting at build_surf_connectivity().");
        

        //const bool mk_sides = false; // not built at this time to save memory, their need is unknown

        num_surfs_in_elem  = 2*num_dims;      // 4 (2D) or 6 (3D)
        const size_t num_1D = Pn+1;
        num_nodes_in_surf = pow(num_1D, num_dims - 1);

        // -----------------------------------------------------------------------
        // 1a. Get element face nodes and sort (small...biggest index) for building connectivity
        // -----------------------------------------------------------------------

        // get the elem nodes on the surface
        CArrayKokkos<size_t> surf_node_ordering_in_elem(num_surfs_in_elem, num_nodes_in_surf,"surf_node_ordering_in_elem");
        get_surf_node_lids(surf_node_ordering_in_elem, num_1D, num_dims);
        
        // sort the nodes on each elem face from smallest to largest, these are the hash keys 
        CArrayKokkos <size_t> face_hash_keys(num_elems,num_surfs_in_elem,num_nodes_in_surf,"face_hash_keys");

        FOR_ALL_CLASS(elem_gid, 0, num_elems,
                      face_lid, 0, num_surfs_in_elem,
                      node_lid, 0, num_nodes_in_surf, {
                    
            const size_t elem_node_lid = surf_node_ordering_in_elem(face_lid,node_lid); // elem nodes on face 
            const size_t node_gid      = nodes_in_elem(elem_gid,elem_node_lid); // all nodes in element
            face_hash_keys(elem_gid,face_lid,node_lid) = node_gid;  // save the gid
        });
        Kokkos::fence();

        FOR_ALL_CLASS(elem_gid, 0, num_elems,
                      face_lid, 0, num_surfs_in_elem, {
            ViewCArrayKokkos <size_t> sorted_face_nodes(&face_hash_keys(elem_gid,face_lid,0),num_nodes_in_surf);
            bubble_sort(sorted_face_nodes); // sort nodes from smallest to largest

            // remember that sorted_face_nodes are in order of smallest to largest now
        });
        Kokkos::fence();

        DCArrayKokkos <size_t> surf_counter(1);
        surf_counter.set_values(0);
        DCArrayKokkos <size_t> bdy_surf_counter(1);
        bdy_surf_counter.set_values(0);

        CArrayKokkos<int> face_elems_in_elem(num_elems, num_surfs_in_elem,"face_elems_in_elem");
        face_elems_in_elem.set_values(-1);

        surfs_in_elem = CArrayKokkos<size_t>(num_elems, num_surfs_in_elem,"mesh.surfs_in_elem");
        //if(mk_sides) sides_in_elem = CArrayKokkos<size_t>(num_elems*num_surfs_in_elem);

        // helper variables for temporary storage, its sized larger than num_surfs
        CArrayKokkos<int>    elems_in_surf_helper(num_elems*num_surfs_in_elem,2,"elems_in_surf_helper"); 
        CArrayKokkos<size_t> faces_in_surf_helper(num_elems*num_surfs_in_elem,2,"faces_in_surf_helper"); 
        CArrayKokkos<size_t> num_elems_in_surf_helper(num_elems*num_surfs_in_elem, "num_elems_in_surf_helper"); 
        //if(mk_sides) CArrayKokkos<size_t> sides_in_surf_helper(num_elems*num_surfs_in_elem,2); 
        CArrayKokkos<size_t> bdy_surfs_helper(num_elems*num_surfs_in_elem,"bdy_surfs_helper"); 

        // -----------------------------------------------------------------------
        // 1b. Build Surfaces
        // -----------------------------------------------------------------------

        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {

            
            FOR_ALL_CLASS(face_lid, 0, num_surfs_in_elem, {


                // only search for a matching surface if it wasn't already found
                if(face_elems_in_elem(elem_gid,face_lid)<0){

                    bool found_nbr_surf = false;

                    // loop over neighboring faces to see if they have the same nodes as this face
                    for(size_t nbr_elem_lid=0; nbr_elem_lid<num_elems_in_elem(elem_gid); nbr_elem_lid++){
                        
                        const size_t nbr_elem_gid = elems_in_elem(elem_gid,nbr_elem_lid);

                        // loop faces of nbr elem
                        for(size_t nbr_face_lid=0; nbr_face_lid<num_surfs_in_elem; nbr_face_lid++){

                            // skip over the neighboring faces if they are tagged with an ID
                            if(face_elems_in_elem(nbr_elem_gid,nbr_face_lid) >= 0) continue;


                            // loop over the nodes in the hash and compare
                            size_t tally = 0; // if tally = num_nodes_in_surf it is a match
                            for(size_t node_lid=0; node_lid<num_nodes_in_surf; node_lid++){
                                const size_t nbr_node_gid = face_hash_keys(nbr_elem_gid,nbr_face_lid,node_lid);
                                const size_t node_gid = face_hash_keys(elem_gid,face_lid,node_lid);
                                if(nbr_node_gid==node_gid) tally++;
                            }

                            if(tally==num_nodes_in_surf){

                                // its a match! Yah, this face has a soul mate!
                                const size_t surf_gid = Kokkos::atomic_fetch_add(&surf_counter(0), 1);
                                surfs_in_elem(elem_gid,face_lid)         = surf_gid;
                                surfs_in_elem(nbr_elem_gid,nbr_face_lid) = surf_gid;
                                
                                // --- must compress these later to size num_surfs ---
                                elems_in_surf_helper(surf_gid,0) = elem_gid; 
                                elems_in_surf_helper(surf_gid,1) = nbr_elem_gid; 
                                faces_in_surf_helper(surf_gid,0) = face_lid; 
                                faces_in_surf_helper(surf_gid,1) = nbr_face_lid; 

                                num_elems_in_surf_helper(surf_gid) = 2;
                                
                                //if(mk_sides){
                                //    const size_t side_gid = face_lid + elem_gid*num_surfs_in_elem;
                                //    const size_t nbr_side_gid = nbr_face_lid + nbr_elem_gid*num_surfs_in_elem;
                                //    sides_in_surf_helper(surf_gid,0) = side_gid; 
                                //    sides_in_surf_helper(surf_gid,1) = nbr_side_gid; 
                                //}

                                face_elems_in_elem(elem_gid,face_lid)         = nbr_elem_gid; 
                                face_elems_in_elem(nbr_elem_gid,nbr_face_lid) = elem_gid;

                                found_nbr_surf = true;

                                break; // exit 
                            } // else not a match to this face_lid

                        } // end for nbr faces in nbr elem

                        if(found_nbr_surf==true) break;

                    } // end for loop over nbr elems
                

                    // boundary surfaces don't have a match
                    if(found_nbr_surf==false){
                        
                        // you didn't find a match, no soulmate for you, until then, you are on the bdy
                        const size_t surf_gid = Kokkos::atomic_fetch_add(&surf_counter(0), 1);
                        surfs_in_elem(elem_gid,face_lid) = surf_gid;
                                                
                        // --- must compress these later to have size num_surfs ---
                        elems_in_surf_helper(surf_gid,0) =  elem_gid; 
                        elems_in_surf_helper(surf_gid,1) = -elem_gid-1; // negative because elem does not exist
                        faces_in_surf_helper(surf_gid,0) =  face_lid; 
                        faces_in_surf_helper(surf_gid,1) = -face_lid-1; // negative because elem does not exist

                        num_elems_in_surf_helper(surf_gid) = 1; // no neighbor, it's a boundary

                        // Important: 
                        // By design, we use a negative index for surface elem and face connectivity ;
                        // on the boundary for the second accessor.  A negative index will cause code to 
                        // stop because element and face data structures are size_t.

                        //if(mk_sides){
                        //    const size_t side_gid = face_lid + elem_gid*num_surfs_in_elem;
                        //    sides_in_surf_helper(surf_gid,0) = side_gid; 
                        //}
                        
                        // --- must compress these later to have size num_bdy_surfs ---
                        const size_t bdy_surf_gid = Kokkos::atomic_fetch_add(&bdy_surf_counter(0), 1);
                        bdy_surfs_helper(bdy_surf_gid) = surf_gid;

                    } // end if bdy surface

                } // end if this surface was already saved
                // remember if I found a face neighbor, the face_elems_in_elem has index >=0

            }); // end parallel for over faces in elem
            Kokkos::fence();  // don't go to the next elem until this one is finished
            // jumping to the next element will break check on face_elems_in_elem(elem_gid,face_lid), thus a fence is needed

        } // end for elems
        surf_counter.update_host();
        bdy_surf_counter.update_host();

        // allocate memory for surface data structures
        num_surfs = surf_counter.host(0);
        num_bdy_surfs = bdy_surf_counter.host(0);
        

        // -----------------------------------------------------------------------
        // 1c. Finish populating values in surface data structures
        // -----------------------------------------------------------------------

        nodes_in_surf     = CArrayKokkos<size_t>(num_surfs,num_nodes_in_surf, "mesh.nodes_in_surf");
        elems_in_surf     = CArrayKokkos<int>(num_surfs, 2, "mesh.elems_in_surf");
        num_elems_in_surf = CArrayKokkos<size_t>(num_surfs, "mesh.num_elems_in_surf");
        faces_in_surf     = CArrayKokkos<int>(num_surfs, 2, "mesh.elem_faces_in_surf");
        //if(mk_sides)sides_in_surf = CArrayKokkos<size_t>(num_surfs, 2, "mesh.sides_in_surf");
        
        FOR_ALL_CLASS(surf_gid, 0, num_surfs, {
            elems_in_surf(surf_gid,0) = elems_in_surf_helper(surf_gid,0); // = elem_gid 
            elems_in_surf(surf_gid,1) = elems_in_surf_helper(surf_gid,1); // = nbr_elem_gid (on bdy = -elem_gid-1)     
            faces_in_surf(surf_gid,0) = faces_in_surf_helper(surf_gid,0); // = face_lid   
            faces_in_surf(surf_gid,1) = faces_in_surf_helper(surf_gid,1); // = nbr_face_lid (on bdy = -face_lid-1)
            num_elems_in_surf(surf_gid) = num_elems_in_surf_helper(surf_gid);
            
            //if(mk_sides){
            //    sides_in_surf(surf_gid,0) = sides_in_surf_helper(surf_gid,0); // = side_gid; 
            //    if(num_elems_in_surf(surf_gid)==2)
            //       sides_in_surf(surf_gid,1) = sides_in_surf_helper(surf_gid,1); // = nbr_side_gid; 
            //}

            const size_t elem_gid = elems_in_surf(surf_gid,0);
            const size_t face_lid = faces_in_surf(surf_gid,0);
            surfs_in_elem(elem_gid, face_lid) = surf_gid;
            

            // set the neighboring element, if it exists
            if(num_elems_in_surf(surf_gid)==2){
                const size_t nbr_elem_gid = elems_in_surf(surf_gid,1);
                const size_t nbr_face_lid = faces_in_surf(surf_gid,1);
                surfs_in_elem(nbr_elem_gid, nbr_face_lid) = surf_gid;
            }

            for(size_t node_lid=0; node_lid<num_nodes_in_surf; node_lid++){
                const size_t elem_node_lid = surf_node_ordering_in_elem(face_lid,node_lid);
                nodes_in_surf(surf_gid,node_lid) = nodes_in_elem(elem_gid, elem_node_lid); // elem nodes on face 
            }    
        }); // end parallel for


        // -----------------------------------------------------------------------
        // 1b. Save boundary surfaces and build boundary nodes
        // -----------------------------------------------------------------------
        bdy_surfs = CArrayKokkos<size_t> (num_bdy_surfs, "mesh.bdy_surfs"); 
        FOR_ALL_CLASS(bdy_surf_gid, 0, num_bdy_surfs,{
            bdy_surfs(bdy_surf_gid) = bdy_surfs_helper(bdy_surf_gid);
        });

        CArrayKokkos<size_t> num_bdy_surfs_in_bdy_node(num_nodes, "mesh.num_bdy_surfs_in_bdy_node");
        num_bdy_surfs_in_bdy_node.set_values(0);

        CArrayKokkos<int> bdy_surf_node_storage_bin(num_bdy_surfs, num_nodes_in_surf, "mesh.bdy_surf_node_storage_bin");
        bdy_surf_node_storage_bin.set_values(-1);

        DCArrayKokkos<size_t> bdy_node_counter(1, "mesh.bdy_node_counter");
        CArrayKokkos<size_t> bdy_node_helper(num_bdy_surfs*num_nodes_in_surf, "mesh.bdy_node_helper");

        FOR_ALL_CLASS(bdy_surf_gid, 0, num_bdy_surfs,{

            const size_t surf_gid = bdy_surfs(bdy_surf_gid);

            const size_t elem_gid = elems_in_surf(surf_gid,0);
            const size_t face_lid = faces_in_surf(surf_gid,0);

            for(size_t node_lid=0; node_lid<num_nodes_in_surf; node_lid++){
                const size_t node_gid = nodes_in_surf(surf_gid,node_lid); // surface nodes
                bdy_surf_node_storage_bin(bdy_surf_gid,node_lid) = Kokkos::atomic_fetch_add(&num_bdy_surfs_in_bdy_node(node_gid), 1);
            }
        });
        Kokkos::fence();

        FOR_ALL_CLASS(bdy_surf_gid, 0, num_bdy_surfs,{
            
            const size_t surf_gid = bdy_surfs(bdy_surf_gid);

            const size_t elem_gid = elems_in_surf(surf_gid,0);
            const size_t face_lid = faces_in_surf(surf_gid,0);

            for(size_t node_lid=0; node_lid<num_nodes_in_surf; node_lid++){

                // when the storage_bin==0, it is the first surface to have this node
                if (bdy_surf_node_storage_bin(bdy_surf_gid,node_lid)==0){
                    const size_t bdy_node_gid = Kokkos::atomic_fetch_add(&bdy_node_counter(0), 1);
                    const size_t node_gid = face_hash_keys(elem_gid,face_lid,node_lid); // sorted nodes on elem face small...big
                    bdy_node_helper(bdy_node_gid) = node_gid;
                } 

            } // end for node_lid
        }); // end parallel for
        Kokkos::fence();
        bdy_node_counter.update_host();

        num_bdy_nodes = bdy_node_counter.host(0);

        CArrayKokkos<size_t> bdy_nodes(num_bdy_nodes, "mesh.bdy_nodes");

        // compress the storage of boundary nodes
        FOR_ALL_CLASS(bdy_node_gid, 0, num_bdy_nodes,{
            bdy_nodes(bdy_node_gid) = bdy_node_helper(bdy_node_gid);
        });


        // -----------------------------------------------------------------------
        // 2a. Build Patches, a decomposition of the surface
        // -----------------------------------------------------------------------
        
        num_patches_in_surf = pow(Pn,(num_dims-1));  // = Pn_order or = Pn_order*Pn_order
        num_nodes_in_patch = 2*(num_dims-1);         // 2 (2D) or 4 (3D)
        num_patches_in_elem = num_surfs_in_elem*num_patches_in_surf;

        patches_in_surf  = patches_in_surf_t(num_patches_in_surf); // initializing

        num_patches     = num_surfs*num_patches_in_surf;
        elems_in_patch  = CArrayKokkos<int>(num_patches, 2, "mesh.elems_in_patch");
        nodes_in_patch  = CArrayKokkos<size_t>(num_patches, num_nodes_in_patch, "mesh.nodes_in_patch");
        
        surf_in_patch   = CArrayKokkos<size_t>(num_patches, "mesh.surf_in_patch");
        patches_in_elem = CArrayKokkos<size_t>(num_elems, num_patches_in_elem,"patches_in_elem");
        

        DCArrayKokkos<size_t> patch_node_ordering_in_elem (num_surfs_in_elem, num_patches_in_surf, num_nodes_in_patch);
        get_patch_node_lids(patch_node_ordering_in_elem, num_1D, num_dims); // R-hand rule node convention for patch nodes


        // now break up the surface into patches
        FOR_ALL_CLASS(surf_gid, 0, num_surfs,{

            const size_t elem_gid     = elems_in_surf(surf_gid,0);
            const size_t face_lid     = faces_in_surf(surf_gid,0);

            const int nbr_elem_gid = elems_in_surf(surf_gid,1); // num_nbrs (if negative, does not exist)
            const int nbr_face_lid = faces_in_surf(surf_gid,1); // num_nbrs (if negative, does not exist)

            // loop patches on this surface
            for(size_t patch_lid =0; patch_lid<num_patches_in_surf; patch_lid++){
                const size_t patch_gid = patches_in_surf(surf_gid,patch_lid);

                elems_in_patch(patch_gid,0) = elem_gid;
                elems_in_patch(patch_gid,1) = nbr_elem_gid; 
 
                surf_in_patch(patch_gid) = surf_gid;

                size_t elem_patch_lid = patch_lid + face_lid*num_patches_in_surf;
                patches_in_elem(elem_gid, elem_patch_lid) = patch_gid;

 
                if(num_elems_in_surf(surf_gid)==2){
                    size_t nbr_elem_patch_lid = patch_lid + nbr_face_lid*num_patches_in_surf;
                    patches_in_elem(nbr_elem_gid, nbr_elem_patch_lid) = patch_gid;
                }
 
                // save the nodes in this patch using the first element to find the surface
                for(size_t patch_node_lid=0; patch_node_lid<num_nodes_in_patch; patch_node_lid++){
                    size_t elem_node_lid=patch_node_ordering_in_elem(face_lid, patch_lid, patch_node_lid);
                     
                    nodes_in_patch(patch_gid,patch_node_lid) = nodes_in_elem(elem_gid, elem_node_lid);
                }
            } // end loop over patches

        });


        // -----------------------------------------------------------------------
        // 2b. Build boundary patches
        // -----------------------------------------------------------------------
        
        num_bdy_patches = num_bdy_surfs*num_patches_in_surf;
        bdy_patches = CArrayKokkos<size_t> (num_bdy_patches, "mesh.bdy_patches");

        FOR_ALL_CLASS(bdy_surf_gid, 0, num_bdy_surfs,{

            const size_t surf_gid = bdy_surfs(bdy_surf_gid);

            // all patches on this surface are on the boundary
            for(size_t patch_lid =0; patch_lid<num_patches_in_surf; patch_lid++){
                const size_t patch_gid = patches_in_surf(surf_gid,patch_lid);
                const size_t bdy_patch_gid = patch_lid + bdy_surf_gid*num_patches_in_surf; 
                bdy_patches(bdy_patch_gid) = patch_gid;
            }
            
        }); // end parallel for
        Kokkos::fence();
 
    } // end of build surface connectivity

        // testing:
        // 8x8x8 linear mesh
        // num_patches = 8*8*9*3 = 1728
        // bdy_patches = 8*8*6 = 384

        // see indexing_utils for ASCII art of element
        // size_t patch_node_ordering[24] = { 0, 4, 6, 2,
        //                                    1, 3, 7, 5,
        //                                    0, 1, 5, 4,
        //                                    3, 2, 6, 7,
        //                                    0, 2, 3, 1,
        //                                    4, 5, 7, 6 };

        //   J
        //   |
        // 3---2
        // |   |  -- I
        // 0---1
        //
        //size_t patch_node_ordering[8] = { 0, 3,
        //                                  1, 2,
        //                                  0, 1,
        //                                  3, 2 };


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_node_node_connectivity
    ///
    /// \brief Build the connectivity between a node and the nodes connected by edges
    ///        
    ///
    /////////////////////////////////////////////////////////////////////////////
    // build the patches
    void build_node_node_connectivity()
    {
        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at build_node_node_connectivity().");
        }
        // find the max number of elems around a node
        size_t max_num_elems_in_node;
        size_t max_num_lcl;
        FOR_REDUCE_MAX_CLASS(node_gid, 0, num_nodes, max_num_lcl, {
            // num_corners_in_node = num_elems_in_node
            size_t max_num = num_corners_in_node(node_gid);

            if (max_num > max_num_lcl) {
                max_num_lcl = max_num;
            }
        }, max_num_elems_in_node); // end parallel reduction on max
        Kokkos::fence();

        // each elem corner will contribute 3 edges to the node. Those edges will likely be the same
        // ones from an adjacent element so it is a safe estimate to multiply by 3
        DynamicRaggedRightArrayKokkos<size_t> temp_nodes_in_nodes(num_nodes, max_num_elems_in_node * 3, "temp_nodes_in_nodes");

        num_nodes_in_node = CArrayKokkos<size_t>(num_nodes, "mesh.num_nodes_in_node");

        // walk over the patches and save the node node connectivity
        RUN_CLASS({
            if (num_dims == 3) {
                for (size_t patch_gid = 0; patch_gid < num_patches; patch_gid++) {
                    for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                        // the first node on the edge
                        size_t node_gid_0 = nodes_in_patch(patch_gid, node_lid);

                        // second node on this edge
                        size_t node_gid_1;

                        if (node_lid == num_nodes_in_patch - 1) {
                            node_gid_1 = nodes_in_patch(patch_gid, 0);
                        }
                        else {
                            node_gid_1 = nodes_in_patch(patch_gid, node_lid + 1);
                        } // end if

                        size_t num_saved_0 = temp_nodes_in_nodes.stride(node_gid_0);
                        size_t num_saved_1 = temp_nodes_in_nodes.stride(node_gid_1);

                        size_t save_0 = 1;
                        size_t save_1 = 1;

                        // check to see if the node_gid_1 was already saved
                        for (size_t contents_lid = 0; contents_lid < num_saved_0; contents_lid++) {
                            if (temp_nodes_in_nodes(node_gid_0, contents_lid) == node_gid_1) {
                                save_0 = 0; // don't save, it was already saved
                            }
                        }

                        // check to see if the node_gid_0 was already saved
                        for (size_t contents_lid = 0; contents_lid < num_saved_1; contents_lid++) {
                            if (temp_nodes_in_nodes(node_gid_1, contents_lid) == node_gid_0) {
                                save_1 = 0;  // don't save, it was already saved
                            }
                        }

                        if (save_0 == 1) {
                            // increment the number of nodes in a node saved
                            temp_nodes_in_nodes.stride(node_gid_0)++;

                            // save the second node to the first node
                            temp_nodes_in_nodes(node_gid_0, num_saved_0) = node_gid_1;
                        }

                        if (save_1 == 1) {
                            // increment the number of nodes in a node saved
                            temp_nodes_in_nodes.stride(node_gid_1)++;

                            // save the first node to the second node
                            temp_nodes_in_nodes(node_gid_1, num_saved_1) = node_gid_0;
                        }

                        // save the strides
                        num_nodes_in_node(node_gid_0) = temp_nodes_in_nodes.stride(node_gid_0);
                        num_nodes_in_node(node_gid_1) = temp_nodes_in_nodes.stride(node_gid_1);
                    } // end for node in patch
                } // end for patches
            } // end if 3D
            else {
                for (size_t patch_gid = 0; patch_gid < num_patches; patch_gid++) {
                    // the first node on the edge
                    size_t node_gid_0 = nodes_in_patch(patch_gid, 0);

                    // second node on this edge
                    size_t node_gid_1 = nodes_in_patch(patch_gid, 1);

                    size_t num_saved_0 = temp_nodes_in_nodes.stride(node_gid_0);
                    size_t num_saved_1 = temp_nodes_in_nodes.stride(node_gid_1);

                    // increment the number of nodes in a node saved
                    temp_nodes_in_nodes.stride(node_gid_0)++;
                    temp_nodes_in_nodes.stride(node_gid_1)++;

                    // save the second node to the first node
                    temp_nodes_in_nodes(node_gid_0, num_saved_0) = node_gid_1;

                    // save the first node to the second node
                    temp_nodes_in_nodes(node_gid_1, num_saved_1) = node_gid_0;

                    // save the strides
                    num_nodes_in_node(node_gid_0) = temp_nodes_in_nodes.stride(node_gid_0);
                    num_nodes_in_node(node_gid_1) = temp_nodes_in_nodes.stride(node_gid_1);
                } // end for patches
            } // end if 2D
        });  // end RUN
        Kokkos::fence();

        nodes_in_node = RaggedRightArrayKokkos<size_t>(num_nodes_in_node, "mesh.nodes_in_node");

        // save the connectivity
        FOR_ALL_CLASS(node_gid, 0, num_nodes, {
            size_t num_saved = 0;
            for (size_t node_lid = 0; node_lid < num_nodes_in_node(node_gid); node_lid++) {
                nodes_in_node(node_gid, num_saved) = temp_nodes_in_nodes(node_gid, num_saved);

                // increment the number of nodes in node saved
                num_saved++;
            } // end for node_lid
        }); // end parallel for over nodes

    } // end of node node connectivity


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_connectivity
    ///
    /// \brief Calls multiple build connectivity functions
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_connectivity()
    {
        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at build_connectivity().");
        }
        build_corner_connectivity();
        if (verbose) printf("Built corner connectivity \n");

        build_elem_elem_connectivity();
        if (verbose) printf("Built element-element connectivity \n");

        build_surf_connectivity();
        if (verbose) printf("Built surface and patch connectivity \n");

        build_node_node_connectivity();
        if (verbose) printf("Built node-node connectivity \n");
    }


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize_bdy_sets
    ///
    /// \brief Initialize memory for boundary sets
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize_bdy_sets(size_t num_bcs)
    {
        if (num_dims == 0) {
            Kokkos::abort("Error: mesh.num_dims is not set. Exiting at init_bdy_sets().");
        }
        // if (num_bcs == 0) {
        //     printf("ERROR: number of boundary sets = 0, set it = 1");
        //     num_bcs = 1;
        // }
        num_bdy_sets = num_bcs;
        num_bdy_patches_in_set = DCArrayKokkos<size_t>(num_bcs, "mesh.num_bdy_patches_in_set");

        // bdy_patches_in_set is a raggedRight array, it is allocated 
        // in tag_bdys fcn after the sparsity is known, see geometry_new.cpp

        return;
    } // end of initialize_bdy_sets method

    
}; // end Mesh_t

} // end namespace swage

#endif