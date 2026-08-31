#include <gtest/gtest.h>
#include <matar.h>
#include "swage/unstructured_mesh.h"

// Ports the connectivity checks from
// examples/reference_element/src/ref_plus_mesh_test.cpp on the same
// 2x2x2, order-3 hex mesh (expected counts: 36 surfaces, 24 boundary
// surfaces) so the hardcoded face-node list below stays valid.

namespace test_mesh_connectivity_detail {

using namespace mtr;

struct SurfMappingCounts {
    size_t bdy_surfs_found;
    size_t invalid_elem;
    size_t invalid_face;
    size_t reverse_mismatch;
};

swage::Mesh_t build_2x2x2_p3_mesh()
{
    swage::Mesh_t mesh;

    const size_t elem_dims    = 3;
    const size_t elem_order   = 3; // cubic element: 4 nodes in 1D
    const size_t num_elems_1D = 2;
    const size_t num_gauss_1D = 2 * elem_order;

    const size_t num_elems    = num_elems_1D * num_elems_1D * num_elems_1D;
    const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
    const size_t num_nodes    = num_nodes_1D * num_nodes_1D * num_nodes_1D;

    mesh.initialize_dims(elem_dims);
    mesh.initialize_elems_Pn(num_elems, elem_order, num_gauss_1D);
    mesh.initialize_nodes(num_nodes);

    FOR_ALL(i, 0, num_elems_1D,
            j, 0, num_elems_1D,
            k, 0, num_elems_1D, {
        size_t elem_gid = i + (j + k * num_elems_1D) * num_elems_1D;
        size_t node_lid = 0;

        for (size_t kc = k * elem_order; kc <= k * elem_order + elem_order; kc++)
        for (size_t jc = j * elem_order; jc <= j * elem_order + elem_order; jc++)
        for (size_t ic = i * elem_order; ic <= i * elem_order + elem_order; ic++) {
            size_t node_gid = ic + (jc + kc * num_nodes_1D) * num_nodes_1D;
            mesh.nodes_in_elem(elem_gid, node_lid) = node_gid;
            node_lid++;
        }
    });
    Kokkos::fence();
    mesh.nodes_in_elem.update_host();

    mesh.build_corner_connectivity();
    mesh.build_elem_elem_connectivity();
    mesh.build_surf_connectivity();

    return mesh;
}

// counts, for the given elem/face, how many of the 16 target node gids
// appear (without duplicates) across all patches on that face
size_t count_face_match_nodes(swage::Mesh_t& mesh, size_t target_elem_gid,
                               size_t target_face_lid, const size_t target_nodes[16])
{
    DCArrayKokkos<size_t> target_dual(16);
    for (size_t i = 0; i < 16; i++) {
        target_dual.host(i) = target_nodes[i];
    }
    target_dual.update_device();

    DCArrayKokkos<size_t> counter(1);
    counter.set_values(0);

    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        if (elem_gid != target_elem_gid) {
            return;
        }

        for (size_t face_lid = 0; face_lid < mesh.num_surfs_in_elem; face_lid++) {
            if (face_lid != target_face_lid) {
                continue;
            }

            bool found[16];
            for (size_t i = 0; i < 16; i++) {
                found[i] = false;
            }

            size_t local_count = 0;
            for (size_t patch_surf_lid = 0; patch_surf_lid < mesh.num_patches_in_surf; patch_surf_lid++) {
                size_t patch_lid = face_lid * mesh.num_patches_in_surf + patch_surf_lid;
                size_t patch_gid = mesh.patches_in_elem(elem_gid, patch_lid);

                for (size_t node_lid = 0; node_lid < 4; node_lid++) {
                    size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);
                    for (size_t ti = 0; ti < 16; ti++) {
                        if (target_dual(ti) == node_gid && !found[ti]) {
                            found[ti] = true;
                            local_count++;
                            break;
                        }
                    }
                }
            }
            Kokkos::atomic_add(&counter(0), local_count);
        }
    });
    Kokkos::fence();
    counter.update_host();

    return counter.host(0);
}

SurfMappingCounts check_surf_elem_bidirectional_mapping(swage::Mesh_t& mesh)
{
    DCArrayKokkos<size_t> bdy_counter(1);
    bdy_counter.set_values(0);
    DCArrayKokkos<size_t> err_elem(1);
    err_elem.set_values(0);
    DCArrayKokkos<size_t> err_face(1);
    err_face.set_values(0);
    DCArrayKokkos<size_t> err_reverse(1);
    err_reverse.set_values(0);

    FOR_ALL(surf_gid, 0, mesh.num_surfs, {
        const size_t num_elems_in_surf = mesh.num_elems_in_surf(surf_gid);
        if (num_elems_in_surf == 1) {
            Kokkos::atomic_add(&bdy_counter(0), 1);
        }

        for (size_t side_lid = 0; side_lid < num_elems_in_surf; side_lid++) {
            const int elem_gid = mesh.elems_in_surf(surf_gid, side_lid);
            const int face_lid = mesh.faces_in_surf(surf_gid, side_lid);

            bool elem_ok = (elem_gid >= 0) && (elem_gid < (int)mesh.num_elems);
            bool face_ok = (face_lid >= 0) && (face_lid < (int)mesh.num_surfs_in_elem);

            if (!elem_ok) {
                Kokkos::atomic_add(&err_elem(0), 1);
            }
            if (!face_ok) {
                Kokkos::atomic_add(&err_face(0), 1);
            }

            if (elem_ok && face_ok) {
                const size_t reverse_surf_gid = mesh.surfs_in_elem((size_t)elem_gid, (size_t)face_lid);
                if (reverse_surf_gid != surf_gid) {
                    Kokkos::atomic_add(&err_reverse(0), 1);
                }
            }
        }
    });
    Kokkos::fence();

    bdy_counter.update_host();
    err_elem.update_host();
    err_face.update_host();
    err_reverse.update_host();

    return { bdy_counter.host(0), err_elem.host(0), err_face.host(0), err_reverse.host(0) };
}

} // namespace test_mesh_connectivity_detail

TEST(MeshConnectivity, SurfAndBoundaryCountsFor2x2x2P3Mesh)
{
    swage::Mesh_t mesh = test_mesh_connectivity_detail::build_2x2x2_p3_mesh();

    // 3 planes of 2x2 internal surfaces in each of x,y,z = 12+12+12 = 36 total
    EXPECT_EQ(mesh.num_surfs, 36u);
    EXPECT_EQ(mesh.num_bdy_surfs, 24u);
    EXPECT_EQ(mesh.num_bdy_patches, mesh.num_bdy_surfs * mesh.num_patches_in_surf);
}

TEST(MeshConnectivity, FaceNodesMatchBetweenAdjacentElements)
{
    swage::Mesh_t mesh = test_mesh_connectivity_detail::build_2x2x2_p3_mesh();

    // Face 1 of element 0 matches face 0 of element 1
    const size_t elem_nodes_in_face_match[16] = {
        3, 10, 17, 24, 52, 59, 66, 73, 101, 108, 115, 122, 150, 157, 164, 171
    };

    size_t count_elem0 = test_mesh_connectivity_detail::count_face_match_nodes(
        mesh, 0, 1, elem_nodes_in_face_match);
    size_t count_elem1 = test_mesh_connectivity_detail::count_face_match_nodes(
        mesh, 1, 0, elem_nodes_in_face_match);

    EXPECT_EQ(count_elem0, 16u);
    EXPECT_EQ(count_elem1, 16u);
}

TEST(MeshConnectivity, SurfElemMappingIsBidirectionallyConsistent)
{
    swage::Mesh_t mesh = test_mesh_connectivity_detail::build_2x2x2_p3_mesh();

    auto result = test_mesh_connectivity_detail::check_surf_elem_bidirectional_mapping(mesh);

    EXPECT_EQ(result.invalid_elem, 0u);
    EXPECT_EQ(result.invalid_face, 0u);
    EXPECT_EQ(result.reverse_mismatch, 0u);
    EXPECT_EQ(result.bdy_surfs_found, mesh.num_bdy_surfs);
}
