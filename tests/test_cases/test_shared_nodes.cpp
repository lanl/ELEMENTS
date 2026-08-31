#include <gtest/gtest.h>
#include <algorithm>
#include <vector>
#include "swage/point_cloud.h"

using namespace mtr;
using namespace swage;

namespace test_shared_nodes_detail
{
    const size_t num_points = 16;

    // Simulates two hexahedral elements sharing a face.
    //
    //   Element 0 occupies [0,1] x [0,1] x [0,1]
    //   Element 1 occupies [1,2] x [0,1] x [0,1]
    //
    //   Each element has 8 corner points.  The two elements share 4 points
    //   on the face x=1, so 16 corner points map to 12 unique nodes:
    //
    //       Element 0 corners (points 0-7):
    //         0: (0,0,0)   1: (1,0,0)   2: (1,1,0)   3: (0,1,0)
    //         4: (0,0,1)   5: (1,0,1)   6: (1,1,1)   7: (0,1,1)
    //
    //       Element 1 corners (points 8-15):
    //         8: (1,0,0)   9: (2,0,0)  10: (2,1,0)  11: (1,1,0)
    //        12: (1,0,1)  13: (2,0,1)  14: (2,1,1)  15: (1,1,1)
    //
    //   Coincident pairs (same (x,y,z), should share a node gid):
    //     point  1 == point  8  at (1,0,0)
    //     point  2 == point 11  at (1,1,0)
    //     point  5 == point 12  at (1,0,1)
    //     point  6 == point 15  at (1,1,1)
    CArrayKokkos<double> build_corner_points()
    {
        CArrayKokkos<double> corner_pts(num_points, 3, "corner_pts");

        RUN({
            // element 0 -- unit cube [0,1]^3
            corner_pts(0,  0) = 0.0;  corner_pts(0,  1) = 0.0;  corner_pts(0,  2) = 0.0;
            corner_pts(1,  0) = 1.0;  corner_pts(1,  1) = 0.0;  corner_pts(1,  2) = 0.0;
            corner_pts(2,  0) = 1.0;  corner_pts(2,  1) = 1.0;  corner_pts(2,  2) = 0.0;
            corner_pts(3,  0) = 0.0;  corner_pts(3,  1) = 1.0;  corner_pts(3,  2) = 0.0;
            corner_pts(4,  0) = 0.0;  corner_pts(4,  1) = 0.0;  corner_pts(4,  2) = 1.0;
            corner_pts(5,  0) = 1.0;  corner_pts(5,  1) = 0.0;  corner_pts(5,  2) = 1.0;
            corner_pts(6,  0) = 1.0;  corner_pts(6,  1) = 1.0;  corner_pts(6,  2) = 1.0;
            corner_pts(7,  0) = 0.0;  corner_pts(7,  1) = 1.0;  corner_pts(7,  2) = 1.0;

            // element 1 -- cube [1,2] x [0,1] x [0,1]
            // points 8,11,12,15 are coincident with points 1,2,5,6 from element 0
            corner_pts(8,  0) = 1.0;  corner_pts(8,  1) = 0.0;  corner_pts(8,  2) = 0.0;
            corner_pts(9,  0) = 2.0;  corner_pts(9,  1) = 0.0;  corner_pts(9,  2) = 0.0;
            corner_pts(10, 0) = 2.0;  corner_pts(10, 1) = 1.0;  corner_pts(10, 2) = 0.0;
            corner_pts(11, 0) = 1.0;  corner_pts(11, 1) = 1.0;  corner_pts(11, 2) = 0.0;
            corner_pts(12, 0) = 1.0;  corner_pts(12, 1) = 0.0;  corner_pts(12, 2) = 1.0;
            corner_pts(13, 0) = 2.0;  corner_pts(13, 1) = 0.0;  corner_pts(13, 2) = 1.0;
            corner_pts(14, 0) = 2.0;  corner_pts(14, 1) = 1.0;  corner_pts(14, 2) = 1.0;
            corner_pts(15, 0) = 1.0;  corner_pts(15, 1) = 1.0;  corner_pts(15, 2) = 1.0;
        });
        Kokkos::fence();

        return corner_pts;
    }

    struct SharedNodeResult {
        CArrayKokkos<size_t> node_in_corner_point;
        size_t num_nodes;
    };

    SharedNodeResult run_shared_node_connectivity()
    {
        CArrayKokkos<double> corner_pts = build_corner_points();

        PointCloud_t pc;
        pc.initialize_shared_node_vars(1.0e-16);

        // build the bin mesh with a small padding around the domain [0,2]^3
        pc.build_bin_mesh(-0.1, -0.1, -0.1,
                            2.1,  2.1,  2.1,
                            10, 10, 10);

        SharedNodeResult result;
        pc.build_shared_node_connectivity(corner_pts, result.node_in_corner_point, result.num_nodes);
        return result;
    }

    // copies a device-only CArrayKokkos<size_t> into a host-readable DCArrayKokkos<size_t>
    std::vector<size_t> to_host_vector(const CArrayKokkos<size_t>& dev_arr, size_t num)
    {
        DCArrayKokkos<size_t> bridge(num);
        FOR_ALL(i, 0, num, {
            bridge(i) = dev_arr(i);
        });
        Kokkos::fence();
        bridge.update_host();

        std::vector<size_t> out(num);
        for (size_t i = 0; i < num; i++) {
            out[i] = bridge.host(i);
        }
        return out;
    }
} // namespace test_shared_nodes_detail

TEST(SharedNodes, NumberOfUniqueNodes)
{
    using namespace test_shared_nodes_detail;
    SharedNodeResult result = run_shared_node_connectivity();

    EXPECT_EQ(result.num_nodes, 12u);
}

TEST(SharedNodes, NodeGidsInRange)
{
    using namespace test_shared_nodes_detail;
    SharedNodeResult result = run_shared_node_connectivity();
    std::vector<size_t> node_gid = to_host_vector(result.node_in_corner_point, num_points);

    for (size_t pt = 0; pt < num_points; pt++) {
        EXPECT_LT(node_gid[pt], result.num_nodes) << "point " << pt;
    }
}

TEST(SharedNodes, CoincidentPointPairsShareGid)
{
    using namespace test_shared_nodes_detail;
    SharedNodeResult result = run_shared_node_connectivity();
    std::vector<size_t> node_gid = to_host_vector(result.node_in_corner_point, num_points);

    const size_t coincident_pairs[4][2] = {{1, 8}, {2, 11}, {5, 12}, {6, 15}};
    for (auto& pair : coincident_pairs) {
        EXPECT_EQ(node_gid[pair[0]], node_gid[pair[1]])
            << "points " << pair[0] << " and " << pair[1] << " should share a node gid";
    }
}

TEST(SharedNodes, ExactlyTwelveDistinctGids)
{
    using namespace test_shared_nodes_detail;
    SharedNodeResult result = run_shared_node_connectivity();
    std::vector<size_t> node_gid = to_host_vector(result.node_in_corner_point, num_points);

    std::vector<size_t> unique_nodes;
    for (size_t pt = 0; pt < num_points; pt++) {
        if (std::find(unique_nodes.begin(), unique_nodes.end(), node_gid[pt]) == unique_nodes.end()) {
            unique_nodes.push_back(node_gid[pt]);
        }
    }

    EXPECT_EQ(unique_nodes.size(), 12u);
}
