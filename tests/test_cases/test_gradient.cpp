#include <gtest/gtest.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "matar.h"

using namespace mtr;

#include "elements/ref_elem.h"

namespace test_gradient_detail
{
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
    double polynomial(const CArrayKokkos<double>& coeff, const double x, const double y, const size_t p_order)
    {
        double result = 0.0;
        for (int j = 0; j <= (int)p_order; ++j) {
            for (int i = 0; i <= (int)p_order - j; ++i) {
                result += coeff(i, j) * pow(x, (double)i) * pow(y, (double)j);
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double polynomial(const CArrayKokkos<double>& coeff, const double x, const double y, const double z, const size_t p_order)
    {
        double result = 0.0;
        for (int k = 0; k <= (int)p_order; ++k) {
            for (int j = 0; j <= (int)p_order - k; ++j) {
                for (int i = 0; i <= (int)p_order - j - k; ++i) {
                    result += coeff(i, j, k) * pow(x, (double)i) * pow(y, (double)j) * pow(z, (double)k);
                }
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double grad_polynomial_x(const CArrayKokkos<double>& coeff, const double x, const size_t p_order)
    {
        double result = 0.0;
        if (p_order == 0) return 0.0;
        for (int i = 1; i <= (int)p_order; ++i) {
            result += (double)i * coeff(i) * pow(x, (double)i - 1);
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double grad_polynomial_x(const CArrayKokkos<double>& coeff, const double x, const double y, const size_t p_order)
    {
        double result = 0.0;
        if (p_order == 0) return 0.0;
        for (int j = 0; j <= (int)p_order; ++j) {
            for (int i = 1; i <= (int)p_order - j; ++i) {
                result += (double)i * coeff(i, j) * pow(x, (double)i - 1) * pow(y, (double)j);
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double grad_polynomial_y(const CArrayKokkos<double>& coeff, const double x, const double y, const size_t p_order)
    {
        double result = 0.0;
        if (p_order == 0) return 0.0;
        for (int j = 1; j <= (int)p_order; ++j) {
            for (int i = 0; i <= (int)p_order - j; ++i) {
                result += (double)j * coeff(i, j) * pow(x, (double)i) * pow(y, (double)j - 1);
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double grad_polynomial_x(const CArrayKokkos<double>& coeff, const double x, const double y, const double z, const size_t p_order)
    {
        double result = 0.0;
        if (p_order == 0) return 0.0;
        for (int k = 0; k <= (int)p_order; ++k) {
            for (int j = 0; j <= (int)p_order - k; ++j) {
                for (int i = 1; i <= (int)p_order - j - k; ++i) {
                    result += (double)i * coeff(i, j, k) * pow(x, (double)i - 1.0) * pow(y, (double)j) * pow(z, (double)k);
                }
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double grad_polynomial_y(const CArrayKokkos<double>& coeff, const double x, const double y, const double z, const size_t p_order)
    {
        double result = 0.0;
        if (p_order == 0) return 0.0;
        for (int k = 0; k <= (int)p_order; ++k) {
            for (int j = 1; j <= (int)p_order - k; ++j) {
                for (int i = 0; i <= (int)p_order - j - k; ++i) {
                    result += (double)j * coeff(i, j, k) * pow(x, (double)i) * pow(y, (double)j - 1.0) * pow(z, (double)k);
                }
            }
        }
        return result;
    }

    KOKKOS_INLINE_FUNCTION
    double grad_polynomial_z(const CArrayKokkos<double>& coeff, const double x, const double y, const double z, const size_t p_order)
    {
        double result = 0.0;
        if (p_order == 0) return 0.0;
        for (int k = 1; k <= (int)p_order; ++k) {
            for (int j = 0; j <= (int)p_order - k; ++j) {
                for (int i = 0; i <= (int)p_order - j - k; ++i) {
                    result += (double)k * coeff(i, j, k) * pow(x, (double)i) * pow(y, (double)j) * pow(z, (double)k - 1.0);
                }
            }
        }
        return result;
    }

    // Differentiates a fixed-coefficient polynomial of degree p_order via
    // qpt_grad_basis and returns the max absolute error (across all qpts and
    // dims) against the analytic partial derivatives, as a host double.
    double max_gradient_error(reference_space::QuadratureType quad_type,
                               reference_space::BasisType basis_type,
                               size_t num_qpts_1d,
                               size_t elem_dims,
                               size_t p_order)
    {
        elements::Quadrature_t Quad;
        Quad.initialize_quadrature(quad_type, num_qpts_1d, elem_dims);

        elements::ReferenceElement_t RefElem;
        RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement, basis_type, Quad, p_order);

        CArrayKokkos<double> coeff;
        if (elem_dims == 1) {
            coeff = CArrayKokkos<double>(p_order + 1);
        } else if (elem_dims == 2) {
            coeff = CArrayKokkos<double>(p_order + 1, p_order + 1);
        } else {
            coeff = CArrayKokkos<double>(p_order + 1, p_order + 1, p_order + 1);
        }
        coeff.set_values(0.78914567);

        CArrayKokkos<double> dof_values(RefElem.num_dofs_in_elem);

        FOR_ALL(dof, 0, RefElem.num_dofs_in_elem, {
            double value = 0.0;
            if (elem_dims == 1) {
                value = polynomial(coeff, RefElem.dof_positions(dof, 0), p_order);
            } else if (elem_dims == 2) {
                value = polynomial(coeff, RefElem.dof_positions(dof, 0), RefElem.dof_positions(dof, 1), p_order);
            } else {
                value = polynomial(coeff, RefElem.dof_positions(dof, 0), RefElem.dof_positions(dof, 1), RefElem.dof_positions(dof, 2), p_order);
            }
            dof_values(dof) = value;
        });
        Kokkos::fence();

        double max_error = 0.0;
        double max_error_lcl = 0.0;

        FOR_REDUCE_MAX(qpt, 0, Quad.num_qpts_in_elem, max_error_lcl, {
            double sum[3];
            sum[0] = 0.0;
            sum[1] = 0.0;
            sum[2] = 0.0;

            for (size_t dof = 0; dof < RefElem.num_dofs_in_elem; dof++) {
                for (size_t dim = 0; dim < RefElem.elem_dims; dim++) {
                    sum[dim] += RefElem.qpt_grad_basis(qpt, dof, dim) * dof_values(dof);
                }
            }

            double exact_value[3];
            exact_value[0] = 0.0;
            exact_value[1] = 0.0;
            exact_value[2] = 0.0;

            if (elem_dims == 1) {
                exact_value[0] = grad_polynomial_x(coeff, Quad.qpt_positions(qpt, 0), p_order);
            } else if (elem_dims == 2) {
                const double xi = Quad.qpt_positions(qpt, 0);
                const double eta = Quad.qpt_positions(qpt, 1);
                exact_value[0] = grad_polynomial_x(coeff, xi, eta, p_order);
                exact_value[1] = grad_polynomial_y(coeff, xi, eta, p_order);
            } else {
                const double xi = Quad.qpt_positions(qpt, 0);
                const double eta = Quad.qpt_positions(qpt, 1);
                const double mu = Quad.qpt_positions(qpt, 2);
                exact_value[0] = grad_polynomial_x(coeff, xi, eta, mu, p_order);
                exact_value[1] = grad_polynomial_y(coeff, xi, eta, mu, p_order);
                exact_value[2] = grad_polynomial_z(coeff, xi, eta, mu, p_order);
            }

            for (size_t dim = 0; dim < RefElem.elem_dims; dim++) {
                double err = fabs(sum[dim] - exact_value[dim]);
                if (err > max_error_lcl) {
                    max_error_lcl = err;
                }
            }
        }, max_error);
        Kokkos::fence();

        return max_error;
    }
} // namespace test_gradient_detail

TEST(Gradient, LegendreQuadratureLegendreBasis)
{
    const std::vector<size_t> ps = {1, 3, 5, 7};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t p : ps) {
            if (elem_dims == 3 && p > 4) continue;
            size_t n = p + 2;
            double err = test_gradient_detail::max_gradient_error(
                reference_space::GaussLegendre, reference_space::LagrangeLegendre, n, elem_dims, p);
            double tol = std::max(1e-10, 1e-10 * (double)p);
            EXPECT_LE(err, tol) << "dims=" << elem_dims << " p=" << p;
        }
    }
}

TEST(Gradient, LobattoQuadratureLobattoBasis)
{
    const std::vector<size_t> ps = {1, 3, 5, 7};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t p : ps) {
            if (elem_dims == 3 && p > 4) continue;
            size_t n = p + 3;
            double err = test_gradient_detail::max_gradient_error(
                reference_space::GaussLobatto, reference_space::LagrangeLobatto, n, elem_dims, p);
            double tol = std::max(1e-10, 1e-10 * (double)p);
            EXPECT_LE(err, tol) << "dims=" << elem_dims << " p=" << p;
        }
    }
}
