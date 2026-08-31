#include <gtest/gtest.h>
#include <vector>
#include "swage/indexing_utils.h"
#include "elements/ref_elem.h"

namespace test_indexing_utils_detail {

// checks that rid(i,j) for i,j in [0,n) is a bijection onto [0, n*n)
bool qpt_rid_2d_is_bijective(size_t n)
{
    std::vector<char> seen(n * n, 0);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) {
            size_t rid = elements::get_qpt_rid(i, j, n);
            if (rid >= n * n) {
                return false;
            }
            if (seen[rid]) {
                return false;
            }
            seen[rid] = 1;
        }
    }
    return true;
}

// checks that rid(i,j,k) for i,j,k in [0,n) is a bijection onto [0, n*n*n)
bool qpt_rid_3d_is_bijective(size_t n)
{
    std::vector<char> seen(n * n * n, 0);
    for (size_t k = 0; k < n; k++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t i = 0; i < n; i++) {
                size_t rid = elements::get_qpt_rid(i, j, k, n);
                if (rid >= n * n * n) {
                    return false;
                }
                if (seen[rid]) {
                    return false;
                }
                seen[rid] = 1;
            }
        }
    }
    return true;
}

bool dof_rid_2d_is_bijective(size_t n)
{
    std::vector<char> seen(n * n, 0);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) {
            size_t rid = elements::get_dof_rid(i, j, n);
            if (rid >= n * n) {
                return false;
            }
            if (seen[rid]) {
                return false;
            }
            seen[rid] = 1;
        }
    }
    return true;
}

bool dof_rid_3d_is_bijective(size_t n)
{
    std::vector<char> seen(n * n * n, 0);
    for (size_t k = 0; k < n; k++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t i = 0; i < n; i++) {
                size_t rid = elements::get_dof_rid(i, j, k, n);
                if (rid >= n * n * n) {
                    return false;
                }
                if (seen[rid]) {
                    return false;
                }
                seen[rid] = 1;
            }
        }
    }
    return true;
}

} // namespace test_indexing_utils_detail

TEST(IndexingUtils, BubbleSortSortsAscending)
{
    size_t arr[6] = {5, 3, 8, 1, 9, 2};
    swage::bubble_sort(arr, 6);

    for (size_t i = 0; i + 1 < 6; i++) {
        EXPECT_LE(arr[i], arr[i + 1]);
    }
    EXPECT_EQ(arr[0], 1u);
    EXPECT_EQ(arr[5], 9u);
}

TEST(IndexingUtils, BubbleSortAlreadySorted)
{
    size_t arr[4] = {1, 2, 3, 4};
    swage::bubble_sort(arr, 4);

    EXPECT_EQ(arr[0], 1u);
    EXPECT_EQ(arr[1], 2u);
    EXPECT_EQ(arr[2], 3u);
    EXPECT_EQ(arr[3], 4u);
}

TEST(IndexingUtils, GetQptRid2DRoundTrip)
{
    EXPECT_TRUE(test_indexing_utils_detail::qpt_rid_2d_is_bijective(4));
    // spot check row-major ordering
    EXPECT_EQ(elements::get_qpt_rid((size_t)0, (size_t)0, (size_t)4), 0u);
    EXPECT_EQ(elements::get_qpt_rid((size_t)1, (size_t)0, (size_t)4), 1u);
    EXPECT_EQ(elements::get_qpt_rid((size_t)0, (size_t)1, (size_t)4), 4u);
}

TEST(IndexingUtils, GetQptRid3DRoundTrip)
{
    EXPECT_TRUE(test_indexing_utils_detail::qpt_rid_3d_is_bijective(3));
    EXPECT_EQ(elements::get_qpt_rid((size_t)0, (size_t)0, (size_t)0, (size_t)3), 0u);
    EXPECT_EQ(elements::get_qpt_rid((size_t)1, (size_t)0, (size_t)0, (size_t)3), 1u);
    EXPECT_EQ(elements::get_qpt_rid((size_t)0, (size_t)1, (size_t)0, (size_t)3), 3u);
    EXPECT_EQ(elements::get_qpt_rid((size_t)0, (size_t)0, (size_t)1, (size_t)3), 9u);
}

TEST(IndexingUtils, GetDofRid2DRoundTrip)
{
    EXPECT_TRUE(test_indexing_utils_detail::dof_rid_2d_is_bijective(5));
    EXPECT_EQ(elements::get_dof_rid((size_t)2, (size_t)3, (size_t)5), 2u + 3u * 5u);
}

TEST(IndexingUtils, GetDofRid3DRoundTrip)
{
    EXPECT_TRUE(test_indexing_utils_detail::dof_rid_3d_is_bijective(3));
    EXPECT_EQ(elements::get_dof_rid((size_t)1, (size_t)2, (size_t)1, (size_t)3), 1u + (2u + 1u * 3u) * 3u);
}

TEST(IndexingUtils, ZonesInElemFunctorHostIndexing)
{
    const size_t num_zones_in_elem = 8;
    swage::zones_in_elem_t zones_in_elem(num_zones_in_elem);

    EXPECT_EQ(zones_in_elem.host(0, 0), 0u);
    EXPECT_EQ(zones_in_elem.host(0, 7), 7u);
    EXPECT_EQ(zones_in_elem.host(2, 3), 2u * num_zones_in_elem + 3u);
}

TEST(IndexingUtils, GaussInElemFunctorHostIndexing)
{
    const size_t num_gauss_in_elem = 27;
    swage::gauss_in_elem_t gauss_in_elem(num_gauss_in_elem);

    EXPECT_EQ(gauss_in_elem.host(0, 0), 0u);
    EXPECT_EQ(gauss_in_elem.host(1, 5), 1u * num_gauss_in_elem + 5u);
    EXPECT_EQ(gauss_in_elem.host(3, 26), 3u * num_gauss_in_elem + 26u);
}

TEST(IndexingUtils, CornersInElemFunctorHostIndexing)
{
    const size_t num_corners_in_elem = 8;
    swage::corners_in_elem_t corners_in_elem(num_corners_in_elem);

    EXPECT_EQ(corners_in_elem.host(0, 0), 0u);
    EXPECT_EQ(corners_in_elem.host(2, 5), 2u * num_corners_in_elem + 5u);
}

TEST(IndexingUtils, PatchesInSurfFunctorHostIndexing)
{
    const size_t num_patches_in_surf = 4;
    swage::patches_in_surf_t patches_in_surf(num_patches_in_surf);

    EXPECT_EQ(patches_in_surf.host(0, 0), 0u);
    EXPECT_EQ(patches_in_surf.host(3, 2), 3u * num_patches_in_surf + 2u);
}
