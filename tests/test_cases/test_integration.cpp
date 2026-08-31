#include <gtest/gtest.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "matar.h"

using namespace mtr;

#include "elements/ref_elem.h"

namespace test_integration_detail
{
    // 1D tensor-product polynomial with all coefficients set to a fixed value.
    KOKKOS_INLINE_FUNCTION
    double polynomial(const CArrayKokkos<double>& coeff, const double x, const size_t p_order)
    {
        double result = 0.0;
        for (int i = 0; i <= (int)p_order; ++i) {
            result += coeff(i) * pow(x, (double)i);
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double integrate_polynomial_1D(const CArrayKokkos<double>& coeff, const size_t p_order)
    {
        double result = 0.0;
        for (int i = 0; i <= (int)p_order; ++i) {
            if (i % 2 == 0) {
                result += coeff(i) * 2.0 / (double)(i + 1);
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double polynomial(const CArrayKokkos<double>& coeff, const double x, const double y, const size_t p_order)
    {
        double result = 0.0;
        for (int j = 0; j <= (int)p_order; ++j) {
            for (int i = 0; i <= (int)p_order; ++i) {
                result += coeff(i, j) * pow(x, (double)i) * pow(y, (double)j);
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double integrate_polynomial_2D(const CArrayKokkos<double>& coeff, const size_t p_order)
    {
        double result = 0.0;
        for (int j = 0; j <= (int)p_order; ++j) {
            for (int i = 0; i <= (int)p_order; ++i) {
                double xi_integral = (i % 2 == 0) ? 2.0 / (double)(i + 1) : 0.0;
                double eta_integral = (j % 2 == 0) ? 2.0 / (double)(j + 1) : 0.0;
                result += coeff(i, j) * xi_integral * eta_integral;
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double polynomial(const CArrayKokkos<double>& coeff, const double x, const double y, const double z, const size_t p_order)
    {
        double result = 0.0;
        for (int k = 0; k <= (int)p_order; ++k) {
            for (int j = 0; j <= (int)p_order; ++j) {
                for (int i = 0; i <= (int)p_order; ++i) {
                    result += coeff(i, j, k) * pow(x, (double)i) * pow(y, (double)j) * pow(z, (double)k);
                }
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double integrate_polynomial_3D(const CArrayKokkos<double>& coeff, const size_t p_order)
    {
        double result = 0.0;
        for (int k = 0; k <= (int)p_order; ++k) {
            for (int j = 0; j <= (int)p_order; ++j) {
                for (int i = 0; i <= (int)p_order; ++i) {
                    double xi_integral = (i % 2 == 0) ? 2.0 / (double)(i + 1) : 0.0;
                    double eta_integral = (j % 2 == 0) ? 2.0 / (double)(j + 1) : 0.0;
                    double mu_integral = (k % 2 == 0) ? 2.0 / (double)(k + 1) : 0.0;
                    result += coeff(i, j, k) * xi_integral * eta_integral * mu_integral;
                }
            }
        }
        return result;
    }

    // Runs the numerical quadrature-based integral and the closed-form integral
    // of a fixed-coefficient tensor-product polynomial of degree p_order, and
    // returns the absolute error as a host double.
    double integration_error(reference_space::QuadratureType quad_type,
                              size_t num_qpts_1d,
                              size_t elem_dims,
                              size_t p_order)
    {
        elements::Quadrature_t Quad;
        Quad.initialize_quadrature(quad_type, num_qpts_1d, elem_dims);

        CArrayKokkos<double> coeff;
        if (elem_dims == 1) {
            coeff = CArrayKokkos<double>(p_order + 1);
        } else if (elem_dims == 2) {
            coeff = CArrayKokkos<double>(p_order + 1, p_order + 1);
        } else {
            coeff = CArrayKokkos<double>(p_order + 1, p_order + 1, p_order + 1);
        }
        coeff.set_values(0.78914567);

        double numerical_integral = 0.0;
        double sum_lcl = 0.0;

        FOR_REDUCE_SUM(qpt, 0, Quad.num_qpts_in_elem, sum_lcl, {
            double value = 0.0;
            if (elem_dims == 1) {
                value = polynomial(coeff, Quad.qpt_positions(qpt, 0), p_order);
            } else if (elem_dims == 2) {
                value = polynomial(coeff, Quad.qpt_positions(qpt, 0), Quad.qpt_positions(qpt, 1), p_order);
            } else {
                value = polynomial(coeff, Quad.qpt_positions(qpt, 0), Quad.qpt_positions(qpt, 1), Quad.qpt_positions(qpt, 2), p_order);
            }
            sum_lcl += value * Quad.qpt_weights(qpt);
        }, numerical_integral);
        Kokkos::fence();

        DCArrayKokkos<double> exact(1);
        RUN({
            double exact_integral = 0.0;
            if (elem_dims == 1) {
                exact_integral = integrate_polynomial_1D(coeff, p_order);
            } else if (elem_dims == 2) {
                exact_integral = integrate_polynomial_2D(coeff, p_order);
            } else {
                exact_integral = integrate_polynomial_3D(coeff, p_order);
            }
            exact(0) = exact_integral;
        });
        Kokkos::fence();
        exact.update_host();

        return fabs(numerical_integral - exact.host(0));
    }
} // namespace test_integration_detail

TEST(Integration, LegendreQuadratureExactForPolynomialsUpToBound)
{
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        size_t max_n = (elem_dims == 3) ? 4 : 8;
        for (size_t n = 1; n <= max_n; n++) {
            size_t p_bound = 2 * n - 1; // Legendre exact to 2n-1
            size_t p_order = std::min(p_bound, elem_dims == 3 ? (size_t)3 : p_bound);
            double error = test_integration_detail::integration_error(
                reference_space::GaussLegendre, n, elem_dims, p_order);
            double tol = std::max(1e-10, 1e-10 * (double)p_order);
            EXPECT_LE(error, tol) << "dims=" << elem_dims << " n=" << n << " p=" << p_order;
        }
    }
}

TEST(Integration, LobattoQuadratureExactForPolynomialsUpToBound)
{
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        size_t max_n = (elem_dims == 3) ? 5 : 8;
        for (size_t n = 2; n <= max_n; n++) {
            long p_bound_signed = 2 * (long)n - 3; // Lobatto exact to 2n-3
            size_t p_bound = (p_bound_signed < 0) ? 0 : (size_t)p_bound_signed;
            size_t p_order = std::min(p_bound, elem_dims == 3 ? (size_t)3 : p_bound);
            double error = test_integration_detail::integration_error(
                reference_space::GaussLobatto, n, elem_dims, p_order);
            double tol = std::max(1e-10, 1e-10 * (double)p_order);
            EXPECT_LE(error, tol) << "dims=" << elem_dims << " n=" << n << " p=" << p_order;
        }
    }
}
