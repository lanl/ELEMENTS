#include <gtest/gtest.h>
#include <matar.h>
#include <cmath>
#include "swage/unstructured_mesh.h"

// Ports the closed-form structured-mesh topology-count formulas from
// examples/reference_element/src/remap_fv_test.cpp (3D) and
// remap_2D_fv_test.cpp (2D), verified against small meshes built directly
// here (3x3x3 3D, 4x4 2D) rather than the examples' large production sizes.

namespace test_mesh_counts_detail {

using namespace mtr;

swage::Mesh_t build_structured_3d_mesh(size_t elem_order, size_t num_elems_1D, size_t num_gauss_1D)
{
    swage::Mesh_t mesh;

    const size_t elem_dims = 3;
    const size_t num_elems = num_elems_1D * num_elems_1D * num_elems_1D;
    const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
    const size_t num_nodes = num_nodes_1D * num_nodes_1D * num_nodes_1D;

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

swage::Mesh_t build_structured_2d_mesh(size_t elem_order, size_t num_elems_1D, size_t num_gauss_1D)
{
    swage::Mesh_t mesh;

    const size_t elem_dims = 2;
    const size_t num_elems = num_elems_1D * num_elems_1D;
    const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
    const size_t num_nodes = num_nodes_1D * num_nodes_1D;

    mesh.initialize_dims(elem_dims);
    mesh.initialize_elems_Pn(num_elems, elem_order, num_gauss_1D);
    mesh.initialize_nodes(num_nodes);

    FOR_ALL(i, 0, num_elems_1D,
            j, 0, num_elems_1D, {
        size_t elem_gid = i + j * num_elems_1D;
        size_t node_lid = 0;

        for (size_t jc = j * elem_order; jc <= j * elem_order + elem_order; jc++)
        for (size_t ic = i * elem_order; ic <= i * elem_order + elem_order; ic++) {
            size_t node_gid = ic + jc * num_nodes_1D;
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

} // namespace test_mesh_counts_detail

TEST(MeshCounts, Structured3x3x3MeshMatchesClosedFormCounts)
{
    const size_t elem_order   = 2;
    const size_t num_elems_1D = 3;
    const size_t num_gauss_1D = 3;

    swage::Mesh_t mesh = test_mesh_counts_detail::build_structured_3d_mesh(
        elem_order, num_elems_1D, num_gauss_1D);

    const size_t NI = num_elems_1D;
    const size_t NJ = num_elems_1D;
    const size_t NK = num_elems_1D;

    const size_t exact_num_internal_surfs = (NI - 1) * NJ * NK + NI * (NJ - 1) * NK + NI * NJ * (NK - 1);
    const size_t exact_num_bdy_surfs = 2 * (NI * NJ + NI * NK + NJ * NK);
    const size_t exact_num_internal_patches = exact_num_internal_surfs * elem_order * elem_order;
    const size_t exact_num_bdy_patches = exact_num_bdy_surfs * elem_order * elem_order;

    EXPECT_EQ(mesh.num_bdy_surfs, exact_num_bdy_surfs);
    EXPECT_EQ(mesh.num_surfs, exact_num_bdy_surfs + exact_num_internal_surfs);
    EXPECT_EQ(mesh.num_bdy_patches, exact_num_bdy_patches);
    EXPECT_EQ(mesh.num_patches, exact_num_bdy_patches + exact_num_internal_patches);

    const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
    const size_t exact_num_nodes = num_nodes_1D * num_nodes_1D * num_nodes_1D;
    EXPECT_EQ(mesh.num_nodes, exact_num_nodes);

    const size_t exact_num_gauss_in_elem = (size_t)std::pow(num_gauss_1D, 3);
    EXPECT_EQ(mesh.num_gauss_in_elem, exact_num_gauss_in_elem);
}

TEST(MeshCounts, Structured4x4MeshMatchesClosedFormCounts2D)
{
    const size_t elem_order   = 2;
    const size_t num_elems_1D = 4;
    const size_t num_gauss_1D = 3;

    swage::Mesh_t mesh = test_mesh_counts_detail::build_structured_2d_mesh(
        elem_order, num_elems_1D, num_gauss_1D);

    const size_t NI = num_elems_1D;
    const size_t NJ = num_elems_1D;

    const size_t exact_num_internal_surfs = (NI - 1) * NJ + NI * (NJ - 1);
    const size_t exact_num_bdy_surfs = 2 * (NI + NJ);
    const size_t exact_num_internal_patches = exact_num_internal_surfs * elem_order;
    const size_t exact_num_bdy_patches = exact_num_bdy_surfs * elem_order;

    EXPECT_EQ(mesh.num_bdy_surfs, exact_num_bdy_surfs);
    EXPECT_EQ(mesh.num_surfs, exact_num_bdy_surfs + exact_num_internal_surfs);
    EXPECT_EQ(mesh.num_bdy_patches, exact_num_bdy_patches);
    EXPECT_EQ(mesh.num_patches, exact_num_bdy_patches + exact_num_internal_patches);

    const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
    const size_t exact_num_nodes = num_nodes_1D * num_nodes_1D;
    EXPECT_EQ(mesh.num_nodes, exact_num_nodes);

    const size_t exact_num_gauss_in_elem = (size_t)std::pow(num_gauss_1D, 2);
    EXPECT_EQ(mesh.num_gauss_in_elem, exact_num_gauss_in_elem);
}
