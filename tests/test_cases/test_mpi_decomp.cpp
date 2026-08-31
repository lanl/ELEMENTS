#include <gtest/gtest.h>
#include <mpi.h>

#include <cstddef>
#include <set>
#include <vector>

#include "../../examples/decomp_example/include/mesh_io.h"

namespace test_mpi_decomp_detail {

// Builds a 4x4x4 box mesh on rank 0 and partitions it across MPI_COMM_WORLD via PT-Scotch,
// mirroring examples/decomp_example/src/mesh_decomp_example.cpp but shrunk for a fast CI test.
void build_and_partition(
    swage::Mesh_t& final_mesh,
    MPICArrayKokkos<double>& final_node_coords,
    CommunicationPlan& element_communication_plan,
    CommunicationPlan& node_communication_plan,
    int world_size,
    int rank)
{
    const int num_dims = 3;
    // Pn = 1 (linear elements) keeps the global node count a simple (n*Pn+1)^3 to check against.
    const int Pn_order = 1;

    double origin[3] = {0.0, 0.0, 0.0};
    double length[3] = {1.0, 1.0, 1.0};
    int num_elems_dim[3] = {4, 4, 4};

    swage::Mesh_t initial_mesh;
    initial_mesh.num_dims = num_dims;
    MPICArrayKokkos<double> initial_node_coords;

    if (rank == 0) {
        initial_mesh.Pn = Pn_order;
        build_3d_box(initial_mesh, initial_node_coords, origin, length, num_elems_dim, Pn_order);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    element_communication_plan.initialize(MPI_COMM_WORLD);
    node_communication_plan.initialize(MPI_COMM_WORLD);

    final_mesh.num_dims = num_dims;

    if (world_size != 1) {
        elements::partition_mesh(initial_mesh, final_mesh, initial_node_coords, final_node_coords,
                                  element_communication_plan, node_communication_plan, world_size, rank);
    } else {
        final_mesh = initial_mesh;
        final_mesh.num_owned_elems = initial_mesh.num_elems;
        final_mesh.num_owned_nodes = initial_mesh.num_nodes;
        final_node_coords = initial_node_coords;
        final_mesh.num_dims = num_dims;
        final_mesh.Pn = Pn_order;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// Counts, on this rank, the owned nodes where shared_tally_owned_nodes is true. Nodes on
// partition boundaries can be locally "owned" by more than one rank; this mask marks true
// only on the single rank that should count it toward a global tally.
int count_tally_owned_nodes(swage::Mesh_t& mesh)
{
    int local_accum = 0;
    int result = 0;
    FOR_REDUCE_SUM(node_lid, 0, mesh.num_owned_nodes, local_accum, {
        if (mesh.shared_tally_owned_nodes(node_lid)) {
            local_accum++;
        }
    }, result);
    MATAR_FENCE();
    return result;
}

} // namespace test_mpi_decomp_detail

TEST(MpiDecomp, OwnershipPartitionOfUnity)
{
    int rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    swage::Mesh_t final_mesh;
    MPICArrayKokkos<double> final_node_coords;
    CommunicationPlan element_communication_plan;
    CommunicationPlan node_communication_plan;
    test_mpi_decomp_detail::build_and_partition(final_mesh, final_node_coords,
        element_communication_plan, node_communication_plan, world_size, rank);

    const long long global_num_elems = 4 * 4 * 4;
    const long long global_num_nodes = 5 * 5 * 5; // Pn = 1: num_points_i = 4*Pn+1 = 5 per axis

    long long local_owned_elems = static_cast<long long>(final_mesh.num_owned_elems);
    long long total_owned_elems = 0;
    MPI_Allreduce(&local_owned_elems, &total_owned_elems, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(total_owned_elems, global_num_elems);

    long long local_tally_nodes =
        static_cast<long long>(test_mpi_decomp_detail::count_tally_owned_nodes(final_mesh));
    long long total_tally_nodes = 0;
    MPI_Allreduce(&local_tally_nodes, &total_tally_nodes, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(total_tally_nodes, global_num_nodes);
}

TEST(MpiDecomp, LocalToGlobalMappingConsistency)
{
    int rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    swage::Mesh_t final_mesh;
    MPICArrayKokkos<double> final_node_coords;
    CommunicationPlan element_communication_plan;
    CommunicationPlan node_communication_plan;
    test_mpi_decomp_detail::build_and_partition(final_mesh, final_node_coords,
        element_communication_plan, node_communication_plan, world_size, rank);

    const size_t global_num_elems = 4 * 4 * 4;
    const size_t global_num_nodes = 5 * 5 * 5;

    final_mesh.local_to_global_elem_mapping.update_host();
    final_mesh.local_to_global_node_mapping.update_host();

    bool elem_ids_valid = true;
    bool elem_ids_injective = true;
    std::set<size_t> seen_elem_gids;
    for (size_t i = 0; i < final_mesh.num_elems; i++) {
        size_t gid = final_mesh.local_to_global_elem_mapping.host(i);
        if (gid >= global_num_elems) {
            elem_ids_valid = false;
        }
        if (!seen_elem_gids.insert(gid).second) {
            elem_ids_injective = false;
        }
    }

    bool node_ids_valid = true;
    bool node_ids_injective = true;
    std::set<size_t> seen_node_gids;
    for (size_t i = 0; i < final_mesh.num_nodes; i++) {
        size_t gid = final_mesh.local_to_global_node_mapping.host(i);
        if (gid >= global_num_nodes) {
            node_ids_valid = false;
        }
        if (!seen_node_gids.insert(gid).second) {
            node_ids_injective = false;
        }
    }

    int local_flags[4] = {
        elem_ids_valid ? 1 : 0,
        elem_ids_injective ? 1 : 0,
        node_ids_valid ? 1 : 0,
        node_ids_injective ? 1 : 0
    };
    int global_flags[4] = {0, 0, 0, 0};
    MPI_Allreduce(local_flags, global_flags, 4, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    EXPECT_EQ(global_flags[0], 1) << "local_to_global_elem_mapping produced an out-of-range global id on some rank";
    EXPECT_EQ(global_flags[1], 1) << "local_to_global_elem_mapping was not injective within some rank";
    EXPECT_EQ(global_flags[2], 1) << "local_to_global_node_mapping produced an out-of-range global id on some rank";
    EXPECT_EQ(global_flags[3], 1) << "local_to_global_node_mapping was not injective within some rank";
}

TEST(MpiDecomp, HaloCommunicationCorrectness)
{
    int rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    swage::Mesh_t final_mesh;
    MPICArrayKokkos<double> final_node_coords;
    CommunicationPlan element_communication_plan;
    CommunicationPlan node_communication_plan;
    test_mpi_decomp_detail::build_and_partition(final_mesh, final_node_coords,
        element_communication_plan, node_communication_plan, world_size, rank);

    final_mesh.local_to_global_elem_mapping.update_host();
    final_mesh.local_to_global_node_mapping.update_host();

    // Deterministic phony fields: value is a known function of the entity's global id, so a
    // ghost entry's value after communicate() can be checked against the same function applied
    // to its global id.
    auto elem_value = [](size_t gid) { return 2.0 * static_cast<double>(gid) + 1.0; };
    auto node_value  = [](size_t gid) { return 3.0 * static_cast<double>(gid) + 5.0; };

    std::vector<gauss_pt_state> gauss_states = {gauss_pt_state::fields};
    GaussPoint_t gauss_point;
    gauss_point.initialize(final_mesh.num_elems, final_mesh.num_dims, gauss_states, element_communication_plan);

    for (size_t i = 0; i < final_mesh.num_owned_elems; i++) {
        size_t gid = final_mesh.local_to_global_elem_mapping.host(i);
        gauss_point.fields.host(i) = elem_value(gid);
    }
    for (size_t i = final_mesh.num_owned_elems; i < final_mesh.num_elems; i++) {
        gauss_point.fields.host(i) = -1.0; // sentinel, must be overwritten by communicate()
    }
    gauss_point.fields.update_device();

    std::vector<node_state> node_states = {node_state::scalar_field};
    node_t final_node;
    final_node.initialize(final_mesh.num_nodes, final_mesh.num_dims, node_states, node_communication_plan);

    for (size_t i = 0; i < final_mesh.num_owned_nodes; i++) {
        size_t gid = final_mesh.local_to_global_node_mapping.host(i);
        final_node.scalar_field.host(i) = node_value(gid);
    }
    for (size_t i = final_mesh.num_owned_nodes; i < final_mesh.num_nodes; i++) {
        final_node.scalar_field.host(i) = -1.0; // sentinel, must be overwritten by communicate()
    }
    final_node.scalar_field.update_device();

    MPI_Barrier(MPI_COMM_WORLD);
    gauss_point.fields.communicate();
    final_node.scalar_field.communicate();
    MPI_Barrier(MPI_COMM_WORLD);

    gauss_point.fields.update_host();
    final_node.scalar_field.update_host();

    for (size_t i = final_mesh.num_owned_elems; i < final_mesh.num_elems; i++) {
        size_t gid = final_mesh.local_to_global_elem_mapping.host(i);
        EXPECT_DOUBLE_EQ(gauss_point.fields.host(i), elem_value(gid));
    }

    for (size_t i = final_mesh.num_owned_nodes; i < final_mesh.num_nodes; i++) {
        size_t gid = final_mesh.local_to_global_node_mapping.host(i);
        EXPECT_DOUBLE_EQ(final_node.scalar_field.host(i), node_value(gid));
    }

    // Sanity check that this decomposition actually produced inter-rank neighbors to exchange.
    long long local_ghost_elems =
        static_cast<long long>(final_mesh.num_elems - final_mesh.num_owned_elems);
    long long total_ghost_elems = 0;
    MPI_Allreduce(&local_ghost_elems, &total_ghost_elems, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_GT(total_ghost_elems, 0);
}
