#include <gtest/gtest.h>
#include "swage/point_cloud.h"

using namespace mtr;
using namespace swage;

TEST(BinMesh, BinSpacing)
{
    PointCloud_t pc;

    pc.build_bin_mesh(0.0, 0.0, 0.0,   // xmin, ymin, zmin
                       1.0, 1.0, 1.0,   // xmax, ymax, zmax
                       10, 10, 10);     // num bins per direction

    EXPECT_NEAR(pc.bin_dx, 0.1, 1.e-14);
    EXPECT_NEAR(pc.bin_dy, 0.1, 1.e-14);
    EXPECT_NEAR(pc.bin_dz, 0.1, 1.e-14);
}

TEST(BinMesh, GetBinKeys)
{
    PointCloud_t pc;

    pc.build_bin_mesh(0.0, 0.0, 0.0,
                       1.0, 1.0, 1.0,
                       10, 10, 10);

    // a point at (0.15, 0.25, 0.35) should land in bin (1, 2, 3)
    bin_keys_t bk = pc.get_bin_keys(0.15, 0.25, 0.35);
    EXPECT_EQ(bk.i, 1);
    EXPECT_EQ(bk.j, 2);
    EXPECT_EQ(bk.k, 3);
}

TEST(BinMesh, GetBinGidMatchesGetIdOfIjk)
{
    PointCloud_t pc;

    pc.build_bin_mesh(0.0, 0.0, 0.0,
                       1.0, 1.0, 1.0,
                       10, 10, 10);

    size_t gid = pc.get_bin_gid(0.15, 0.25, 0.35);
    EXPECT_EQ(gid, get_id_of_ijk(1, 2, 3, 10, 10));
}

TEST(BinMesh, GetBinKeysClampsBelowDomain)
{
    PointCloud_t pc;

    pc.build_bin_mesh(0.0, 0.0, 0.0,
                       1.0, 1.0, 1.0,
                       10, 10, 10);

    bin_keys_t bk = pc.get_bin_keys(-1.0, -1.0, -1.0); // below domain
    EXPECT_EQ(bk.i, 0);
    EXPECT_EQ(bk.j, 0);
    EXPECT_EQ(bk.k, 0);
}

TEST(BinMesh, GetBinKeysClampsAboveDomain)
{
    PointCloud_t pc;

    pc.build_bin_mesh(0.0, 0.0, 0.0,
                       1.0, 1.0, 1.0,
                       10, 10, 10);

    bin_keys_t bk = pc.get_bin_keys(2.0, 2.0, 2.0); // above domain
    EXPECT_EQ(bk.i, 9);
    EXPECT_EQ(bk.j, 9);
    EXPECT_EQ(bk.k, 9);
}
