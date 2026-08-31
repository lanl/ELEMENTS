#include <gtest/gtest.h>
#include <cmath>
#include "matar.h"

using namespace mtr;

#include "elements/ref_elem.h"

namespace test_kronecker_delta_detail
{
    // Collocated Lobatto quadrature/basis: num_qpts_1d == num_dofs_1d, both
    // Lobatto, so quadrature points coincide with DOF positions. Returns the
    // number of (basis, qpt) pairs where qpt_basis violates the Kronecker
    // delta property.
    size_t count_kronecker_delta_violations(size_t num_qpts_1d, size_t elem_dims)
    {
        elements::Quadrature_t Quad;
        Quad.initialize_quadrature(reference_space::GaussLobatto, num_qpts_1d, elem_dims);

        const size_t p_order = num_qpts_1d - 1;

        elements::ReferenceElement_t RefElem;
        RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                     reference_space::LagrangeLobatto, Quad, p_order);

        DCArrayKokkos<size_t> errors(1);
        errors.set_values(0);

        FOR_ALL(basis, 0, RefElem.num_dofs_in_elem, {
            for (size_t dof_pt = 0; dof_pt < RefElem.num_dofs_in_elem; dof_pt++) {
                double expected = (basis == dof_pt) ? 1.0 : 0.0;
                if (fabs(RefElem.qpt_basis(dof_pt, basis) - expected) > 1e-12) {
                    Kokkos::atomic_add(&errors(0), (size_t)1);
                }
            }
        });
        Kokkos::fence();

        errors.update_host();
        return errors.host(0);
    }
} // namespace test_kronecker_delta_detail

TEST(KroneckerDelta, CollocatedLobatto1D)
{
    for (size_t n : {2, 4, 6}) {
        size_t violations = test_kronecker_delta_detail::count_kronecker_delta_violations(n, 1);
        EXPECT_EQ(violations, 0u) << "n=" << n;
    }
}

TEST(KroneckerDelta, CollocatedLobatto2D)
{
    for (size_t n : {2, 4, 6}) {
        size_t violations = test_kronecker_delta_detail::count_kronecker_delta_violations(n, 2);
        EXPECT_EQ(violations, 0u) << "n=" << n;
    }
}

TEST(KroneckerDelta, CollocatedLobatto3D)
{
    for (size_t n : {2, 3, 4}) {
        size_t violations = test_kronecker_delta_detail::count_kronecker_delta_violations(n, 3);
        EXPECT_EQ(violations, 0u) << "n=" << n;
    }
}
