#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "matar.h"

using namespace mtr;

#include "elements/ref_elem.h"

namespace test_partition_unity_detail
{
    // Builds Quadrature_t/ReferenceElement_t for the given combination and
    // returns the number of qpt/dim slots where the basis sum != 1 or the
    // grad-basis sum != 0, accumulated via atomics on device.
    size_t count_partition_unity_violations(reference_space::QuadratureType quad_type,
                                             reference_space::BasisType basis_type,
                                             size_t num_qpts_1d,
                                             size_t elem_dims,
                                             size_t p_order)
    {
        elements::Quadrature_t Quad;
        Quad.initialize_quadrature(quad_type, num_qpts_1d, elem_dims);

        elements::ReferenceElement_t RefElem;
        RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement, basis_type, Quad, p_order);

        DCArrayKokkos<size_t> errors(1);
        errors.set_values(0);

        FOR_ALL(qpt, 0, Quad.num_qpts_in_elem, {
            double sum = 0.0;
            for (size_t dof = 0; dof < RefElem.num_dofs_in_elem; dof++) {
                sum += RefElem.qpt_basis(qpt, dof);
            }
            if (fabs(sum - 1.0) > 1e-12) {
                Kokkos::atomic_add(&errors(0), (size_t)1);
            }

            for (size_t dim = 0; dim < RefElem.elem_dims; dim++) {
                double grad_sum = 0.0;
                for (size_t dof = 0; dof < RefElem.num_dofs_in_elem; dof++) {
                    grad_sum += RefElem.qpt_grad_basis(qpt, dof, dim);
                }
                if (fabs(grad_sum) > 1e-12) {
                    Kokkos::atomic_add(&errors(0), (size_t)1);
                }
            }
        });
        Kokkos::fence();

        errors.update_host();
        return errors.host(0);
    }
} // namespace test_partition_unity_detail

TEST(PartitionUnity, LegendreQuadratureLegendreBasis)
{
    const std::vector<size_t> ns = {2, 3, 5};
    const std::vector<size_t> ps = {1, 3, 5};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t n : ns) {
            for (size_t p : ps) {
                size_t violations = test_partition_unity_detail::count_partition_unity_violations(
                    reference_space::GaussLegendre, reference_space::LagrangeLegendre, n, elem_dims, p);
                EXPECT_EQ(violations, 0u) << "dims=" << elem_dims << " n=" << n << " p=" << p;
            }
        }
    }
}

TEST(PartitionUnity, LegendreQuadratureLobattoBasis)
{
    const std::vector<size_t> ns = {2, 3, 5};
    const std::vector<size_t> ps = {1, 3, 5};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t n : ns) {
            for (size_t p : ps) {
                size_t violations = test_partition_unity_detail::count_partition_unity_violations(
                    reference_space::GaussLegendre, reference_space::LagrangeLobatto, n, elem_dims, p);
                EXPECT_EQ(violations, 0u) << "dims=" << elem_dims << " n=" << n << " p=" << p;
            }
        }
    }
}

TEST(PartitionUnity, LobattoQuadratureLegendreBasis)
{
    const std::vector<size_t> ns = {2, 3, 5};
    const std::vector<size_t> ps = {1, 3, 5};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t n : ns) {
            for (size_t p : ps) {
                size_t violations = test_partition_unity_detail::count_partition_unity_violations(
                    reference_space::GaussLobatto, reference_space::LagrangeLegendre, n, elem_dims, p);
                EXPECT_EQ(violations, 0u) << "dims=" << elem_dims << " n=" << n << " p=" << p;
            }
        }
    }
}

TEST(PartitionUnity, LobattoQuadratureLobattoBasis)
{
    const std::vector<size_t> ns = {2, 3, 5};
    const std::vector<size_t> ps = {1, 3, 5};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t n : ns) {
            for (size_t p : ps) {
                size_t violations = test_partition_unity_detail::count_partition_unity_violations(
                    reference_space::GaussLobatto, reference_space::LagrangeLobatto, n, elem_dims, p);
                EXPECT_EQ(violations, 0u) << "dims=" << elem_dims << " n=" << n << " p=" << p;
            }
        }
    }
}
