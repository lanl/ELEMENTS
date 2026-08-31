#include <gtest/gtest.h>
#include <cmath>
#include "matar.h"

using namespace mtr;

#include "elements/ref_elem.h"

namespace test_lagrange_1d_detail
{
    // Builds `num` Lobatto dof positions on [-1,1], then checks lagrange_basis_1D
    // and lagrange_derivative_1D directly (no Quadrature_t/ReferenceElement_t
    // involved):
    //   results(0) = max |lagrange_basis_1D(x_j)_i - delta_ij| over dof nodes x_j
    //   results(1) = max |sum_i lagrange_basis_1D(test_pt)_i - 1| over test points
    //   results(2) = max |sum_i lagrange_derivative_1D(test_pt)_i| over test points
    // Returned as a host-readable DCArrayKokkos<double> of length 3.
    DCArrayKokkos<double> evaluate_lagrange_1d_checks(size_t num)
    {
        CArrayKokkos<double> dof_positions_1d(num);
        RUN({
            elements::get_lobatto_nodes_1D(dof_positions_1d, num);
        });
        Kokkos::fence();

        CArrayKokkos<double> basis(num);
        CArrayKokkos<double> deriv(num);
        DCArrayKokkos<double> results(3);

        RUN({
            double max_kron_error = 0.0;
            for (size_t j = 0; j < num; j++) {
                elements::lagrange_basis_1D(basis, dof_positions_1d, dof_positions_1d(j));
                for (size_t i = 0; i < num; i++) {
                    double expected = (i == j) ? 1.0 : 0.0;
                    double err = fabs(basis(i) - expected);
                    if (err > max_kron_error) {
                        max_kron_error = err;
                    }
                }
            }

            double test_pts[4];
            test_pts[0] = 0.0;
            test_pts[1] = 0.3;
            test_pts[2] = -0.55;
            test_pts[3] = 0.8;

            double max_unity_error = 0.0;
            double max_deriv_sum_error = 0.0;
            for (size_t t = 0; t < 4; t++) {
                elements::lagrange_basis_1D(basis, dof_positions_1d, test_pts[t]);
                double sum = 0.0;
                for (size_t i = 0; i < num; i++) {
                    sum += basis(i);
                }
                double unity_err = fabs(sum - 1.0);
                if (unity_err > max_unity_error) {
                    max_unity_error = unity_err;
                }

                elements::lagrange_derivative_1D(deriv, dof_positions_1d, test_pts[t]);
                double dsum = 0.0;
                for (size_t i = 0; i < num; i++) {
                    dsum += deriv(i);
                }
                double dsum_err = fabs(dsum);
                if (dsum_err > max_deriv_sum_error) {
                    max_deriv_sum_error = dsum_err;
                }
            }

            results(0) = max_kron_error;
            results(1) = max_unity_error;
            results(2) = max_deriv_sum_error;
        });
        Kokkos::fence();

        results.update_host();
        return results;
    }
} // namespace test_lagrange_1d_detail

TEST(Lagrange1D, KroneckerDeltaAtDofNodes)
{
    for (size_t num : {3, 4, 5, 6}) {
        DCArrayKokkos<double> results = test_lagrange_1d_detail::evaluate_lagrange_1d_checks(num);
        EXPECT_NEAR(results.host(0), 0.0, 1e-12) << "num=" << num;
    }
}

TEST(Lagrange1D, PartitionOfUnityAtInteriorPoints)
{
    for (size_t num : {3, 4, 5, 6}) {
        DCArrayKokkos<double> results = test_lagrange_1d_detail::evaluate_lagrange_1d_checks(num);
        EXPECT_NEAR(results.host(1), 0.0, 1e-12) << "num=" << num;
    }
}

TEST(Lagrange1D, DerivativeSumIsZeroAtInteriorPoints)
{
    for (size_t num : {3, 4, 5, 6}) {
        DCArrayKokkos<double> results = test_lagrange_1d_detail::evaluate_lagrange_1d_checks(num);
        EXPECT_NEAR(results.host(2), 0.0, 1e-10) << "num=" << num;
    }
}
