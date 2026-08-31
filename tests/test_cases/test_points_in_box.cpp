#include <gtest/gtest.h>
#include "swage/point_cloud.h"

using namespace mtr;
using namespace swage;

namespace test_points_in_box_detail
{
    // number of points in each direction of the structured 4x4x4 cloud
    const size_t N = 4;
    const size_t num_points = N * N * N;
    const size_t num_boxes = 2;

    DCArrayKokkos<double> build_structured_cloud()
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

        return point_positions;
    }

    struct Bounds {
        double xmin, ymin, zmin;
        double xmax, ymax, zmax;
    };

    Bounds get_bounds(PointCloud_t& pc, const DCArrayKokkos<double>& point_positions)
    {
        Bounds b;
        pc.get_bounds_point_cloud(b.xmin, b.ymin, b.zmin,
                                   b.xmax, b.ymax, b.zmax,
                                   point_positions);
        return b;
    }

    // box 0 corners = (0,0,0) & (1,1,1); box 1 corners = (0.5,0.5,0.5) & (1.5,1.5,1.5)
    CArrayKokkos<double> build_bounding_boxes()
    {
        CArrayKokkos<double> bounding_box(num_boxes, 2, 3);

        RUN({
            bounding_box(0,0,0) = 0.0; // xmin
            bounding_box(0,0,1) = 0.0; // ymin
            bounding_box(0,0,2) = 0.0; // zmin
            bounding_box(0,1,0) = 1.0; // xmax
            bounding_box(0,1,1) = 1.0; // ymax
            bounding_box(0,1,2) = 1.0; // zmax

            bounding_box(1,0,0) = 0.5; // xmin
            bounding_box(1,0,1) = 0.5; // ymin
            bounding_box(1,0,2) = 0.5; // zmin
            bounding_box(1,1,0) = 1.5; // xmax
            bounding_box(1,1,1) = 1.5; // ymax
            bounding_box(1,1,2) = 1.5; // zmax
        });
        Kokkos::fence();

        return bounding_box;
    }

    struct PointsInBoxResult {
        RaggedRightArrayKokkos<size_t> points_in_box;
        DCArrayKokkos<size_t> num_points_in_box;
    };

    PointsInBoxResult run_points_in_box()
    {
        DCArrayKokkos<double> point_positions = build_structured_cloud();

        PointCloud_t pc;
        Bounds b = get_bounds(pc, point_positions);
        pc.build_bin_mesh(b.xmin - 0.5, b.ymin - 0.5, b.zmin - 0.5,
                           b.xmax + 0.5, b.ymax + 0.5, b.zmax + 0.5,
                           6, 6, 6);

        PointsInBoxResult result;
        result.num_points_in_box = DCArrayKokkos<size_t>(num_boxes);
        CArrayKokkos<double> bounding_box = build_bounding_boxes();

        pc.get_points_in_box(point_positions, result.points_in_box, result.num_points_in_box, bounding_box);
        result.num_points_in_box.update_host();

        return result;
    }

    // points_in_box is a device-only RaggedRightArrayKokkos (no host mirror),
    // so read back the four entries the tests care about via a device kernel.
    struct FourEntries {
        size_t box0_first;
        size_t box0_last;
        size_t box1_first;
        size_t box1_last;
    };

    FourEntries read_first_and_last_entries(const RaggedRightArrayKokkos<size_t>& points_in_box)
    {
        DCArrayKokkos<size_t> bridge(4);

        RUN({
            bridge(0) = points_in_box(0, 0);
            bridge(1) = points_in_box(0, 7);
            bridge(2) = points_in_box(1, 0);
            bridge(3) = points_in_box(1, 7);
        });
        Kokkos::fence();
        bridge.update_host();

        FourEntries entries;
        entries.box0_first = bridge.host(0);
        entries.box0_last  = bridge.host(1);
        entries.box1_first = bridge.host(2);
        entries.box1_last  = bridge.host(3);
        return entries;
    }
} // namespace test_points_in_box_detail

TEST(PointsInBox, DomainBoundsPaddedByHalfMicron)
{
    using namespace test_points_in_box_detail;

    DCArrayKokkos<double> point_positions = build_structured_cloud();
    PointCloud_t pc;
    Bounds b = get_bounds(pc, point_positions);

    // cloud spans exactly [0,2]^3; get_bounds_point_cloud pads each side by 1e-6
    EXPECT_NEAR(b.xmin, -1.0e-6, 1e-9);
    EXPECT_NEAR(b.ymin, -1.0e-6, 1e-9);
    EXPECT_NEAR(b.zmin, -1.0e-6, 1e-9);

    EXPECT_NEAR(b.xmax, 2.0 + 1.0e-6, 1e-9);
    EXPECT_NEAR(b.ymax, 2.0 + 1.0e-6, 1e-9);
    EXPECT_NEAR(b.zmax, 2.0 + 1.0e-6, 1e-9);
}

TEST(PointsInBox, EightPointsFoundInEachBox)
{
    using namespace test_points_in_box_detail;

    PointsInBoxResult result = run_points_in_box();

    EXPECT_EQ(result.num_points_in_box.host(0), 8u);
    EXPECT_EQ(result.num_points_in_box.host(1), 8u);
}

TEST(PointsInBox, PointIndicesMatchExpectedCorners)
{
    using namespace test_points_in_box_detail;

    PointsInBoxResult result = run_points_in_box();
    FourEntries entries = read_first_and_last_entries(result.points_in_box);

    EXPECT_EQ(entries.box0_first, 0u);
    EXPECT_EQ(entries.box0_last, 21u);
    EXPECT_EQ(entries.box1_first, 21u);
    EXPECT_EQ(entries.box1_last, 42u);
}
