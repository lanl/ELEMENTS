#include <gtest/gtest.h>
#include <algorithm>
#include <vector>
#include "swage/point_cloud.h"

using namespace mtr;
using namespace swage;

namespace test_point_cloud_detail
{
    // number of points in each direction of the structured 3x3x3 cloud
    const size_t N = 3;
    const size_t num_points = N * N * N;

    // Builds a structured N x N x N point cloud spanning [0, 2]^3 and its
    // point-to-point connectivity. All output arrays inside pc are dual
    // (DCArrayKokkos / DRaggedRightArrayKokkos) and are already
    // update_host()'d by build_point_cloud_connectivity, so callers may read
    // them with .host(...) directly.
    PointCloud_t build_point_cloud()
    {
        const double h = 2.0 / ((double)N - 1);

        DCArrayKokkos<double> point_positions(num_points, 3, "point_positions");

        size_t point_gid = 0;
        for (size_t k = 0; k < N; k++)
        for (size_t j = 0; j < N; j++)
        for (size_t i = 0; i < N; i++) {
            point_positions.host(point_gid, 0) = (double)i * h;
            point_positions.host(point_gid, 1) = (double)j * h;
            point_positions.host(point_gid, 2) = (double)k * h;
            point_gid++;
        }
        point_positions.update_device();

        PointCloud_t pc;
        pc.build_bin_mesh(-0.5, -0.5, -0.5, 2.5, 2.5, 2.5, 6, 6, 6);

        const size_t min_num_points_fit = 1;
        pc.initialize_point_cloud_vars(1.5 * h, min_num_points_fit);

        pc.build_point_cloud_connectivity(point_positions);

        return pc;
    }
} // namespace test_point_cloud_detail

TEST(PointCloud, NeighborListsAreSymmetric)
{
    using namespace test_point_cloud_detail;
    PointCloud_t pc = build_point_cloud();

    for (size_t i = 0; i < num_points; i++) {
        for (size_t lid = 0; lid < pc.points_num_neighbors.host(i); lid++) {
            size_t j = pc.points_in_point.host(i, lid);

            bool found = false;
            for (size_t jlid = 0; jlid < pc.points_num_neighbors.host(j); jlid++) {
                if (pc.points_in_point.host(j, jlid) == i) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << "point " << j << " does not list point " << i << " as a neighbor";
        }
    }
}

TEST(PointCloud, ReverseNeighborLidRoundTrips)
{
    using namespace test_point_cloud_detail;
    PointCloud_t pc = build_point_cloud();

    for (size_t i = 0; i < num_points; i++) {
        for (size_t lid = 0; lid < pc.points_num_neighbors.host(i); lid++) {
            size_t j    = pc.points_in_point.host(i, lid);
            size_t jlid = pc.reverse_neighbor_lid.host(i, lid);
            EXPECT_EQ(pc.points_in_point.host(j, jlid), i)
                << "reverse map for point " << i << " lid " << lid << " does not point back";
        }
    }
}

TEST(PointCloud, NoSelfOrDuplicateNeighbors)
{
    using namespace test_point_cloud_detail;
    PointCloud_t pc = build_point_cloud();

    for (size_t i = 0; i < num_points; i++) {
        std::vector<size_t> seen;

        for (size_t lid = 0; lid < pc.points_num_neighbors.host(i); lid++) {
            size_t j = pc.points_in_point.host(i, lid);

            EXPECT_NE(j, i) << "point " << i << " lists itself as a neighbor";

            bool duplicate = (std::find(seen.begin(), seen.end(), j) != seen.end());
            EXPECT_FALSE(duplicate) << "point " << i << " has duplicate neighbor " << j;

            seen.push_back(j);
        }
    }
}

TEST(PointCloud, CenterAndCornerNeighborCounts)
{
    using namespace test_point_cloud_detail;
    PointCloud_t pc = build_point_cloud();

    // center point of the 3x3x3 grid is index 13
    EXPECT_EQ(pc.points_num_neighbors.host(13), 26u);

    // corner point (index 0) should have 2x2x2-1 = 7 neighbors
    EXPECT_EQ(pc.points_num_neighbors.host(0), 7u);
}
