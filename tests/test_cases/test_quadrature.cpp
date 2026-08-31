#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "matar.h"

using namespace mtr;

#include "elements/ref_elem.h"

namespace test_quadrature_detail
{
    // Builds a Quadrature_t of the requested type/dims/num_qpts_1d and returns
    // the host-side sum of all quadrature weights.
    double sum_of_weights(reference_space::QuadratureType type,
                           size_t num_qpts_1d,
                           size_t elem_dims)
    {
        elements::Quadrature_t Quad;
        Quad.initialize_quadrature(type, num_qpts_1d, elem_dims);

        double weight_sum = 0.0;
        double weight_sum_lcl = 0.0;

        FOR_REDUCE_SUM(qpt, 0, Quad.num_qpts_in_elem, weight_sum_lcl, {
            weight_sum_lcl += Quad.qpt_weights(qpt);
        }, weight_sum);
        Kokkos::fence();

        return weight_sum;
    }

    // Reads back the 1D Lobatto nodes/weights for a given num into host arrays.
    void get_lobatto_1d_host(size_t num, CArrayKokkos<double>& nodes, CArrayKokkos<double>& weights)
    {
        nodes = CArrayKokkos<double>(num);
        weights = CArrayKokkos<double>(num);
        RUN({
            elements::get_lobatto_nodes_1D(nodes, num);
            elements::get_lobatto_weights_1D(weights, num);
        });
        Kokkos::fence();
    }

    void get_legendre_1d_host(size_t num, CArrayKokkos<double>& nodes, CArrayKokkos<double>& weights)
    {
        nodes = CArrayKokkos<double>(num);
        weights = CArrayKokkos<double>(num);
        RUN({
            elements::get_legendre_nodes_1D(nodes, num);
            elements::get_legendre_weights_1D(weights, num);
        });
        Kokkos::fence();
    }

    // Copies a small CArrayKokkos<double>(num) to a std::vector on host via a
    // DCArrayKokkos bridge (CArrayKokkos itself has no host mirror API here).
    std::vector<double> to_host_vector(const CArrayKokkos<double>& dev_arr, size_t num)
    {
        DCArrayKokkos<double> bridge(num);
        FOR_ALL(i, 0, num, {
            bridge(i) = dev_arr(i);
        });
        Kokkos::fence();
        bridge.update_host();

        std::vector<double> out(num);
        for (size_t i = 0; i < num; i++) {
            out[i] = bridge.host(i);
        }
        return out;
    }
} // namespace test_quadrature_detail

TEST(Quadrature, WeightsSumToVolumeGaussLegendre)
{
    const std::vector<size_t> ns = {2, 3, 5, 8, 12};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t n : ns) {
            double sum = test_quadrature_detail::sum_of_weights(
                reference_space::GaussLegendre, n, elem_dims);
            EXPECT_NEAR(sum, std::pow(2.0, (double)elem_dims), 1e-12)
                << "GaussLegendre dims=" << elem_dims << " n=" << n;
        }
    }
}

TEST(Quadrature, WeightsSumToVolumeGaussLobatto)
{
    const std::vector<size_t> ns = {2, 3, 5, 8, 12};
    for (size_t elem_dims = 1; elem_dims <= 3; elem_dims++) {
        for (size_t n : ns) {
            double sum = test_quadrature_detail::sum_of_weights(
                reference_space::GaussLobatto, n, elem_dims);
            EXPECT_NEAR(sum, std::pow(2.0, (double)elem_dims), 1e-12)
                << "GaussLobatto dims=" << elem_dims << " n=" << n;
        }
    }
}

TEST(Quadrature, LobattoTable1DSymmetryAndEndpoints)
{
    for (size_t n : {4, 7, 10}) {
        CArrayKokkos<double> dev_nodes, dev_weights;
        test_quadrature_detail::get_lobatto_1d_host(n, dev_nodes, dev_weights);
        std::vector<double> nodes = test_quadrature_detail::to_host_vector(dev_nodes, n);
        std::vector<double> weights = test_quadrature_detail::to_host_vector(dev_weights, n);

        for (size_t i = 0; i < n; i++) {
            EXPECT_NEAR(nodes[i], -nodes[n - 1 - i], 1e-12) << "n=" << n << " i=" << i;
            EXPECT_NEAR(weights[i], weights[n - 1 - i], 1e-12) << "n=" << n << " i=" << i;
        }
        EXPECT_DOUBLE_EQ(nodes[0], -1.0) << "n=" << n;
        EXPECT_DOUBLE_EQ(nodes[n - 1], 1.0) << "n=" << n;
    }
}

TEST(Quadrature, LegendreTable1DSymmetry)
{
    for (size_t n : {4, 7, 10}) {
        CArrayKokkos<double> dev_nodes, dev_weights;
        test_quadrature_detail::get_legendre_1d_host(n, dev_nodes, dev_weights);
        std::vector<double> nodes = test_quadrature_detail::to_host_vector(dev_nodes, n);
        std::vector<double> weights = test_quadrature_detail::to_host_vector(dev_weights, n);

        for (size_t i = 0; i < n; i++) {
            EXPECT_NEAR(nodes[i], -nodes[n - 1 - i], 1e-12) << "n=" << n << " i=" << i;
            EXPECT_NEAR(weights[i], weights[n - 1 - i], 1e-12) << "n=" << n << " i=" << i;
        }
    }
}
