#include <gtest/gtest.h>
#include <cmath>
#include "geometry/geometry.h"
#include "cramers_rule.hpp"

namespace test_jacobian_detail
{
    using namespace mtr;
    using namespace elements;

    // Ensight FE hex8 node ordering -> IJK (lexicographic) node ordering.
    const size_t convert_ensight_to_ijk[8] = {0, 1, 3, 2, 4, 5, 7, 6};

    struct JacobianResult
    {
        double J[3][3];
        double det;
    };

    // Builds a single order-1 hex8 element from Ensight-ordered node coords,
    // evaluates the (constant, single-quadrature-point) Jacobian and its
    // determinant on device, and returns both on host.
    JacobianResult compute_hex_jacobian(const double node_coords_ensight[8][3])
    {
        DCArrayKokkos<double> node_coords_dual(8, 3);
        DCArrayKokkos<double> node_in_elem_dual(8);

        for (size_t n = 0; n < 8; n++) {
            const size_t node_lid = convert_ensight_to_ijk[n];
            for (size_t i = 0; i < 3; i++) {
                node_coords_dual.host(node_lid, i) = node_coords_ensight[n][i];
            }
            node_in_elem_dual.host(node_lid) = node_lid;
        }
        node_coords_dual.update_device();
        node_in_elem_dual.update_device();

        Quadrature_t Quad;
        ReferenceElement_t RefElem;
        Quad.initialize_quadrature(reference_space::GaussLegendre, 1, 3);
        RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                     reference_space::LagrangeLobatto,
                                     Quad,
                                     1);

        DCArrayKokkos<double> jac(3, 3);
        DCArrayKokkos<double> det_arr(1);
        RUN({
            ViewCArrayKokkos<double> a_grad_basis(&RefElem.qpt_grad_basis(0, 0, 0),
                                                   RefElem.num_dofs_in_elem, 3);

            jacobian(jac, node_coords_dual, node_in_elem_dual, a_grad_basis);
            det_arr(0) = det_3x3(jac);
        });
        Kokkos::fence();
        jac.update_host();
        det_arr.update_host();

        JacobianResult result;
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                result.J[i][j] = jac.host(i, j);
            }
        }
        result.det = det_arr.host(0);
        return result;
    }

    struct ManufacturedResult
    {
        double vol_J[3][3][3];   // [point][i][j]
        double vol_det[3];
        double surf_J[18][3][3]; // [point][i][j]
        double surf_det[18];
    };

    // Builds the elem_order=1 hex8 element (IJK-ordered node coords), a
    // quadrature/reference element of degree 3, and its reference surface,
    // then evaluates the Jacobian at 3 interior test quadrature points and
    // 18 surface test quadrature points (3 per face, 6 faces), returning
    // everything on host.
    ManufacturedResult compute_manufactured_jacobians(const double node_coords_ijk[8][3])
    {
        DCArrayKokkos<double> node_coords_dual(8, 3, "node coords dual");
        DCArrayKokkos<double> node_in_elem_dual(8, "nodes in elem dual");

        for (size_t n = 0; n < 8; n++) {
            node_in_elem_dual.host(n) = n;
            for (size_t i = 0; i < 3; i++) {
                node_coords_dual.host(n, i) = node_coords_ijk[n][i];
            }
        }
        node_coords_dual.update_device();
        node_in_elem_dual.update_device();

        Quadrature_t Quad;
        ReferenceElement_t RefElem;
        Quad.initialize_quadrature(reference_space::GaussLegendre, 3, 3);
        RefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                     reference_space::LagrangeLobatto,
                                     Quad,
                                     1);

        SurfaceQuadrature_t SurfQuad;
        ReferenceSurface_t RefSurf;
        SurfQuad.initialize_quadrature(reference_space::GaussLegendre, 3, 3);
        RefSurf.initialize_ref_surf(SurfQuad, RefElem);

        DCArrayKokkos<double> vol_jac(3, 3, 3, "vol_jac");
        DCArrayKokkos<double> vol_det(3, "vol_det");

        RUN({
            // test qpt ids are 0, 13, 26 (uniform stride of 13)
            for (size_t p = 0; p < 3; p++) {
                const size_t qpt_id = 13 * p;
                ViewCArrayKokkos<double> a_grad_basis(&RefElem.qpt_grad_basis(qpt_id, 0, 0),
                                                       RefElem.num_dofs_in_elem, 3);
                ViewCArrayKokkos<double> jac(&vol_jac(p, 0, 0), 3, 3);
                jacobian(jac, node_coords_dual, node_in_elem_dual, a_grad_basis);
                vol_det(p) = det_3x3(jac);
            }
        });
        Kokkos::fence();
        vol_jac.update_host();
        vol_det.update_host();

        DCArrayKokkos<double> surf_jac(18, 3, 3, "surf_jac");
        DCArrayKokkos<double> surf_det(18, "surf_det");

        RUN({
            for (size_t side = 0; side < 6; side++) {
                size_t qpt_id = 0;
                for (size_t lid = 0; lid < 3; lid++) {
                    const size_t p = lid + side * 3;
                    ViewCArrayKokkos<double> a_grad_basis(&RefSurf.qpt_grad_basis(side, qpt_id, 0, 0),
                                                           RefElem.num_dofs_in_elem, 3);
                    ViewCArrayKokkos<double> jac(&surf_jac(p, 0, 0), 3, 3);
                    jacobian(jac, node_coords_dual, node_in_elem_dual, a_grad_basis);
                    surf_det(p) = det_3x3(jac);

                    qpt_id += 4;
                    if (qpt_id > 8) qpt_id = 0;
                }
            }
        });
        Kokkos::fence();
        surf_jac.update_host();
        surf_det.update_host();

        ManufacturedResult result;
        for (size_t p = 0; p < 3; p++) {
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++) {
                    result.vol_J[p][i][j] = vol_jac.host(p, i, j);
                }
            }
            result.vol_det[p] = vol_det.host(p);
        }
        for (size_t p = 0; p < 18; p++) {
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++) {
                    result.surf_J[p][i][j] = surf_jac.host(p, i, j);
                }
            }
            result.surf_det[p] = surf_det.host(p);
        }
        return result;
    }

    // Hex8 trilinear basis gradients, IJK reference-node ordering.
    void compute_hex8_basis_gradients_ijk(double xi, double eta, double zeta, double grad_basis[8][3])
    {
        double ref_coords[8][3] = {
            {-1, -1, -1}, { 1, -1, -1}, {-1,  1, -1}, { 1,  1, -1},
            {-1, -1,  1}, { 1, -1,  1}, {-1,  1,  1}, { 1,  1,  1}
        };

        for (int i = 0; i < 8; i++) {
            const double xi_i = ref_coords[i][0];
            const double eta_i = ref_coords[i][1];
            const double zeta_i = ref_coords[i][2];

            grad_basis[i][0] = 0.125 * xi_i * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta);
            grad_basis[i][1] = 0.125 * (1.0 + xi_i * xi) * eta_i * (1.0 + zeta_i * zeta);
            grad_basis[i][2] = 0.125 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * zeta_i;
        }
    }

    void compute_analytical_jacobian(double xi, double eta, double zeta,
                                      const double node_coords_ijk[8][3], double J[3][3])
    {
        double grad_basis[8][3];
        compute_hex8_basis_gradients_ijk(xi, eta, zeta, grad_basis);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                J[i][j] = 0.0;
                for (int k = 0; k < 8; k++) {
                    J[i][j] += node_coords_ijk[k][i] * grad_basis[k][j];
                }
            }
        }
    }

    double analytical_det_3x3(const double J[3][3])
    {
        return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
             - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
             + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    }
} // namespace test_jacobian_detail

TEST(Jacobian, AffineTransformation)
{
    // x = A*xi + b : Jacobian is constant everywhere and equals A.
    const double A[3][3] = {
        {0.5, 0.3, 0.2},
        {0.2, 0.6, 0.1},
        {0.1, 0.2, 0.4}
    };
    const double b[3] = {1.0, 2.0, 3.0};

    const double ref_corners[8][3] = {
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
        {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1}
    };

    double node_coords[8][3];
    for (size_t n = 0; n < 8; n++) {
        for (size_t i = 0; i < 3; i++) {
            node_coords[n][i] = b[i];
            for (size_t j = 0; j < 3; j++) {
                node_coords[n][i] += A[i][j] * ref_corners[n][j];
            }
        }
    }

    const double expected_det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                               - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                               + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    test_jacobian_detail::JacobianResult result =
        test_jacobian_detail::compute_hex_jacobian(node_coords);

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            EXPECT_NEAR(result.J[i][j], A[i][j], 1e-12) << "i=" << i << " j=" << j;
        }
    }
    EXPECT_NEAR(result.det, expected_det, 1e-12);
}

TEST(Jacobian, SkewedHexahedron)
{
    const double node_coords[8][3] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.2, 0.1},
        {1.1, 1.0, 0.15},
        {0.1, 0.9, 0.05},
        {0.1, 0.1, 0.9},
        {1.2, 0.1, 1.0},
        {1.3, 1.1, 1.1},
        {0.0, 1.0, 0.95}
    };

    // Trilinear-basis central-difference approximation of the Jacobian at
    // the element center (xi=eta=zeta=0).
    double J_approx[3][3];
    for (int i = 0; i < 3; i++) {
        J_approx[i][0] = 0.125 * ((node_coords[1][i] + node_coords[2][i] +
                                    node_coords[5][i] + node_coords[6][i]) -
                                   (node_coords[0][i] + node_coords[3][i] +
                                    node_coords[4][i] + node_coords[7][i]));

        J_approx[i][1] = 0.125 * ((node_coords[2][i] + node_coords[3][i] +
                                    node_coords[6][i] + node_coords[7][i]) -
                                   (node_coords[0][i] + node_coords[1][i] +
                                    node_coords[4][i] + node_coords[5][i]));

        J_approx[i][2] = 0.125 * ((node_coords[4][i] + node_coords[5][i] +
                                    node_coords[6][i] + node_coords[7][i]) -
                                   (node_coords[0][i] + node_coords[1][i] +
                                    node_coords[2][i] + node_coords[3][i]));
    }

    const double det_approx = J_approx[0][0] * (J_approx[1][1] * J_approx[2][2] - J_approx[1][2] * J_approx[2][1])
                             - J_approx[0][1] * (J_approx[1][0] * J_approx[2][2] - J_approx[1][2] * J_approx[2][0])
                             + J_approx[0][2] * (J_approx[1][0] * J_approx[2][1] - J_approx[1][1] * J_approx[2][0]);

    test_jacobian_detail::JacobianResult result =
        test_jacobian_detail::compute_hex_jacobian(node_coords);

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            EXPECT_NEAR(result.J[i][j], J_approx[i][j], 1e-10) << "i=" << i << " j=" << j;
        }
    }
    EXPECT_NEAR(result.det, det_approx, 1e-10);
}

TEST(Jacobian, RotatedScaledCube)
{
    const double theta = M_PI / 6.0;
    const double phi = M_PI / 4.0;
    const double sx = 0.5, sy = 0.8, sz = 1.2;

    const double cos_t = cos(theta), sin_t = sin(theta);
    const double cos_p = cos(phi), sin_p = sin(phi);

    const double T[3][3] = {
        {sx * cos_t * cos_p, sx * (-sin_p), sx * sin_t * cos_p},
        {sy * cos_t * sin_p, sy * cos_p,    sy * sin_t * sin_p},
        {sz * (-sin_t),      sz * 0.0,      sz * cos_t}
    };

    const double ref_verts[8][3] = {
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
        {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1}
    };

    double node_coords[8][3];
    for (int n = 0; n < 8; n++) {
        for (int i = 0; i < 3; i++) {
            node_coords[n][i] = 0.0;
            for (int j = 0; j < 3; j++) {
                node_coords[n][i] += T[i][j] * ref_verts[n][j];
            }
        }
    }

    const double expected_det = sx * sy * sz;

    test_jacobian_detail::JacobianResult result =
        test_jacobian_detail::compute_hex_jacobian(node_coords);

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            EXPECT_NEAR(result.J[i][j], T[i][j], 1e-10) << "i=" << i << " j=" << j;
        }
    }
    EXPECT_NEAR(result.det, expected_det, 1e-10);
}

TEST(Jacobian, ManufacturedSolutionVolume)
{
    const double node_coords_ensight[8][3] = {
        {0.0,  0.0,  0.0},
        {1.0,  0.1,  0.05},
        {0.9,  1.0,  0.1},
        {0.05, 0.95, 0.05},
        {0.1,  0.05, 1.0},
        {1.05, 0.15, 0.95},
        {1.0,  1.05, 1.05},
        {0.1,  1.0,  0.95}
    };

    double node_coords_ijk[8][3];
    for (size_t n = 0; n < 8; n++) {
        const size_t ijk_lid = test_jacobian_detail::convert_ensight_to_ijk[n];
        for (size_t i = 0; i < 3; i++) {
            node_coords_ijk[ijk_lid][i] = node_coords_ensight[n][i];
        }
    }

    test_jacobian_detail::ManufacturedResult result =
        test_jacobian_detail::compute_manufactured_jacobians(node_coords_ijk);

    const double qpt = 0.774596669241483377035853079956; // sqrt(3/5)
    const double test_points[3][3] = {
        {-qpt, -qpt, -qpt},
        {0.0, 0.0, 0.0},
        {qpt, qpt, qpt}
    };

    const double tol = 1e-10;
    for (int p = 0; p < 3; p++) {
        double J_analytical[3][3];
        test_jacobian_detail::compute_analytical_jacobian(
            test_points[p][0], test_points[p][1], test_points[p][2],
            node_coords_ijk, J_analytical);
        const double det_analytical = test_jacobian_detail::analytical_det_3x3(J_analytical);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                EXPECT_NEAR(result.vol_J[p][i][j], J_analytical[i][j], tol)
                    << "pt=" << p << " i=" << i << " j=" << j;
            }
        }
        EXPECT_NEAR(result.vol_det[p], det_analytical, tol) << "pt=" << p;
    }
}

TEST(Jacobian, ManufacturedSolutionSurface)
{
    const double node_coords_ensight[8][3] = {
        {0.0,  0.0,  0.0},
        {1.0,  0.1,  0.05},
        {0.9,  1.0,  0.1},
        {0.05, 0.95, 0.05},
        {0.1,  0.05, 1.0},
        {1.05, 0.15, 0.95},
        {1.0,  1.05, 1.05},
        {0.1,  1.0,  0.95}
    };

    double node_coords_ijk[8][3];
    for (size_t n = 0; n < 8; n++) {
        const size_t ijk_lid = test_jacobian_detail::convert_ensight_to_ijk[n];
        for (size_t i = 0; i < 3; i++) {
            node_coords_ijk[ijk_lid][i] = node_coords_ensight[n][i];
        }
    }

    test_jacobian_detail::ManufacturedResult result =
        test_jacobian_detail::compute_manufactured_jacobians(node_coords_ijk);

    const double qpt = 0.774596669241483377035853079956; // sqrt(3/5)
    const double surf_test_points[18][3] = {
        // side 0 (xi=-1)
        {-1, -qpt, -qpt}, {-1, 0.0, 0.0}, {-1, qpt, qpt},
        // side 1 (xi=1)
        {1, -qpt, -qpt},  {1, 0.0, 0.0},  {1, qpt, qpt},
        // side 2 (eta=-1)
        {-qpt, -1, -qpt}, {0.0, -1, 0.0}, {qpt, -1, qpt},
        // side 3 (eta=1)
        {-qpt, 1, -qpt},  {0.0, 1, 0.0},  {qpt, 1, qpt},
        // side 4 (mu=-1)
        {-qpt, -qpt, -1}, {0.0, 0.0, -1}, {qpt, qpt, -1},
        // side 5 (mu=1)
        {-qpt, -qpt, 1},  {0.0, 0.0, 1},  {qpt, qpt, 1}
    };

    const double tol = 1e-10;
    for (int p = 0; p < 18; p++) {
        double J_analytical[3][3];
        test_jacobian_detail::compute_analytical_jacobian(
            surf_test_points[p][0], surf_test_points[p][1], surf_test_points[p][2],
            node_coords_ijk, J_analytical);
        const double det_analytical = test_jacobian_detail::analytical_det_3x3(J_analytical);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                EXPECT_NEAR(result.surf_J[p][i][j], J_analytical[i][j], tol)
                    << "pt=" << p << " i=" << i << " j=" << j;
            }
        }
        EXPECT_NEAR(result.surf_det[p], det_analytical, tol) << "pt=" << p;
    }
}
