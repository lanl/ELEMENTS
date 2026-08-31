#include <gtest/gtest.h>
#include <matar.h>
#include "swage/unstructured_mesh.h"

// Ports the legitimate patch-connectivity checks from
// examples/reference_element/src/test_2D_patches.cpp, converting the
// example's bare assert()s (compiled out under NDEBUG) into host-readable
// error counters checked with EXPECT_*. Skips the example's
// verify_patch_orientation helper, which is dead code guarded on
// mesh.num_dims==3 but indexes node_coords(1, n1, d) -- a 3-index access
// on what is actually a 2-index DCArrayKokkos<double>(num_nodes, num_dims)
// array; it is broken and is intentionally NOT ported.

namespace test_patches_2d_detail {

using namespace mtr;

struct PatchErrorCounts {
    size_t basics_errors;
    size_t surface_membership_errors;
    size_t internal_continuity_errors;
    size_t boundary_errors;
};

// Builds a classic P1 quad mesh (num_elems_1D x num_elems_1D elements),
// matching the element-node ordering convention used in
// examples/reference_element/src/test_2D_patches.cpp.
swage::Mesh_t build_p1_quad_mesh(size_t num_elems_1D)
{
    swage::Mesh_t mesh;

    const size_t elem_dims = 2;
    const size_t elem_order = 1;
    const size_t num_elems = num_elems_1D * num_elems_1D;
    const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
    const size_t num_nodes = num_nodes_1D * num_nodes_1D;

    mesh.initialize_dims(elem_dims);
    mesh.initialize_elems(num_elems);
    mesh.initialize_nodes(num_nodes);

    FOR_ALL(i, 0, num_elems_1D,
            j, 0, num_elems_1D, {
        size_t elem_gid = i + j * num_elems_1D;

        for (size_t jc = j; jc <= j + elem_order; jc++)
        for (size_t ic = i; ic <= i + elem_order; ic++) {
            size_t node_gid = ic + jc * num_nodes_1D;

            // Global layout: 2-3    Element layout: 3-2
            //                0-1                    0-1
            size_t elem_node_lid;
            if (jc == j) {
                elem_node_lid = (ic - i);       // 0 or 1
            } else {
                elem_node_lid = 3 - (ic - i);   // 3 or 2
            }

            mesh.nodes_in_elem(elem_gid, elem_node_lid) = node_gid;
        }
    });
    Kokkos::fence();
    mesh.nodes_in_elem.update_host();

    mesh.build_corner_connectivity();
    mesh.build_elem_elem_connectivity();
    mesh.build_surf_connectivity();
    mesh.build_node_node_connectivity();

    return mesh;
}

size_t count_patch_basics_errors(swage::Mesh_t& mesh)
{
    DCArrayKokkos<size_t> errors(1);
    errors.set_values(0);

    FOR_ALL(patch_gid, 0, mesh.num_patches, {
        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++) {
            size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);

            if (node_gid >= mesh.num_nodes) {
                Kokkos::atomic_add(&errors(0), 1);
            }

            for (size_t other_lid = node_lid + 1; other_lid < mesh.num_nodes_in_patch; other_lid++) {
                if (mesh.nodes_in_patch(patch_gid, other_lid) == node_gid) {
                    Kokkos::atomic_add(&errors(0), 1);
                }
            }
        }
    });
    Kokkos::fence();
    errors.update_host();

    return errors.host(0);
}

size_t count_patches_in_surfaces_errors(swage::Mesh_t& mesh)
{
    DCArrayKokkos<size_t> errors(1);
    errors.set_values(0);

    FOR_ALL(patch_gid, 0, mesh.num_patches, {
        size_t surf_gid = mesh.surf_in_patch(patch_gid);

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
                Kokkos::atomic_add(&errors(0), 1);
            }
        }
    });
    Kokkos::fence();
    errors.update_host();

    return errors.host(0);
}

size_t count_internal_patch_continuity_errors(swage::Mesh_t& mesh)
{
    DCArrayKokkos<size_t> errors(1);
    errors.set_values(0);

    FOR_ALL(patch_gid, 0, mesh.num_patches, {
        size_t surf_gid = mesh.surf_in_patch(patch_gid);

        if (mesh.num_elems_in_surf(surf_gid) == 2) {
            size_t elem0 = mesh.elems_in_patch(patch_gid, 0);
            int elem1 = mesh.elems_in_patch(patch_gid, 1);

            if (elem1 < 0) {
                Kokkos::atomic_add(&errors(0), 1);
            }

            for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; patch_node_lid++) {
                size_t patch_node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);

                bool found_in_elem0 = false;
                bool found_in_elem1 = false;

                for (size_t elem_node_lid = 0; elem_node_lid < mesh.num_nodes_in_elem; elem_node_lid++) {
                    if (mesh.nodes_in_elem(elem0, elem_node_lid) == patch_node_gid) {
                        found_in_elem0 = true;
                    }
                    if (elem1 >= 0 && mesh.nodes_in_elem((size_t)elem1, elem_node_lid) == patch_node_gid) {
                        found_in_elem1 = true;
                    }
                }

                if (!found_in_elem0 || !found_in_elem1) {
                    Kokkos::atomic_add(&errors(0), 1);
                }
            }
        }
    });
    Kokkos::fence();
    errors.update_host();

    return errors.host(0);
}

size_t count_boundary_patch_errors(swage::Mesh_t& mesh)
{
    DCArrayKokkos<size_t> errors(1);
    errors.set_values(0);

    FOR_ALL(bdy_patch_gid, 0, mesh.num_bdy_patches, {
        size_t patch_gid = mesh.bdy_patches(bdy_patch_gid);
        size_t surf_gid = mesh.surf_in_patch(patch_gid);

        if (mesh.num_elems_in_surf(surf_gid) != 1) {
            Kokkos::atomic_add(&errors(0), 1);
        }

        int elem1 = mesh.elems_in_patch(patch_gid, 1);
        if (elem1 >= 0) {
            Kokkos::atomic_add(&errors(0), 1);
        }
    });
    Kokkos::fence();
    errors.update_host();

    return errors.host(0);
}

PatchErrorCounts verify_all_patch_connectivity(swage::Mesh_t& mesh)
{
    PatchErrorCounts counts;
    counts.basics_errors = count_patch_basics_errors(mesh);
    counts.surface_membership_errors = count_patches_in_surfaces_errors(mesh);
    counts.internal_continuity_errors = count_internal_patch_continuity_errors(mesh);
    counts.boundary_errors = count_boundary_patch_errors(mesh);
    return counts;
}

} // namespace test_patches_2d_detail

TEST(Patches2D, SingleElementMeshPatchConnectivity)
{
    // Matches the example's single P1 quad element (elem_dims=2, elem_order=1,
    // num_elems_1D=1): every surface is on the boundary.
    swage::Mesh_t mesh = test_patches_2d_detail::build_p1_quad_mesh(1);

    EXPECT_EQ(mesh.num_patches, 4u);
    EXPECT_EQ(mesh.num_bdy_patches, 4u);

    auto counts = test_patches_2d_detail::verify_all_patch_connectivity(mesh);

    EXPECT_EQ(counts.basics_errors, 0u);
    EXPECT_EQ(counts.surface_membership_errors, 0u);
    EXPECT_EQ(counts.internal_continuity_errors, 0u);
    EXPECT_EQ(counts.boundary_errors, 0u);
}

TEST(Patches2D, MultiElementMeshHasInternalAndBoundaryPatches)
{
    // A 2x2 mesh has both internal (shared) and boundary patches, exercising
    // the internal-patch-continuity check that the single-element mesh above
    // cannot reach.
    swage::Mesh_t mesh = test_patches_2d_detail::build_p1_quad_mesh(2);

    EXPECT_GT(mesh.num_patches, mesh.num_bdy_patches);

    auto counts = test_patches_2d_detail::verify_all_patch_connectivity(mesh);

    EXPECT_EQ(counts.basics_errors, 0u);
    EXPECT_EQ(counts.surface_membership_errors, 0u);
    EXPECT_EQ(counts.internal_continuity_errors, 0u);
    EXPECT_EQ(counts.boundary_errors, 0u);
}
