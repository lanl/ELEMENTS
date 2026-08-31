#include <gtest/gtest.h>
#include <cmath>
#include "geometry/geometry.h"
#include "cramers_rule.hpp"

namespace test_element_geometry_detail
{
    using namespace mtr;
    using namespace swage;
    using namespace elements;

    struct ElementGeometryMetrics
    {
        double max_vol_error;          // max_elem |sum_qpt det(J)*w - h^3|
        double max_normal_tally_error; // max_elem,dim |sum_face,qpt area_normal|
        double max_gradx_identity_error; // max_elem,i,j |(1/V) sum n(x)x - I|
        double max_div_error;          // max_elem |(1/V) sum n.x - num_dims|
    };

    // Builds a num_elems_1D^3 linear hex mesh on the unit cube, computes
    // per-element quantities on device (element volume via Gauss quadrature,
    // and, via Nanson's-formula face-quadrature area normals, the tally of
    // normals, the outer-product with position, and the divergence sum),
    // reduces the 4 geometric-invariant error metrics to host doubles, and
    // returns them.
    ElementGeometryMetrics compute_element_geometry_metrics(size_t num_elems_1D)
    {
        Quadrature_t Quad;
        ReferenceElement_t FERefElem;
        SurfaceQuadrature_t SurfQuad;
        ReferenceSurface_t RefSurf;
        Mesh_t Mesh;

        const size_t elem_dims = 3;
        const size_t elem_order = 1;
        const size_t num_qpts_1d = 2 * elem_order;

        Quad.initialize_quadrature(reference_space::GaussLegendre, num_qpts_1d, elem_dims);
        FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                       reference_space::LagrangeLobatto,
                                       Quad,
                                       elem_order);

        SurfQuad.initialize_quadrature(reference_space::GaussLegendre, num_qpts_1d, elem_dims);
        RefSurf.initialize_ref_surf(SurfQuad, FERefElem);

        const size_t num_elems = num_elems_1D * num_elems_1D * num_elems_1D;
        const size_t num_nodes_1D = elem_order * num_elems_1D + 1;
        const size_t num_nodes = num_nodes_1D * num_nodes_1D * num_nodes_1D;

        Mesh.initialize_dims(elem_dims);
        Mesh.initialize_elems_Pn(num_elems, elem_order, Quad.num_qpts_1d);
        Mesh.initialize_nodes(num_nodes);

        DCArrayKokkos<double> node_coords(Mesh.num_nodes, Mesh.num_dims);
        const double h = 1.0 / ((double)num_nodes_1D);

        FOR_ALL(kc, 0, num_nodes_1D,
                jc, 0, num_nodes_1D,
                ic, 0, num_nodes_1D, {
            const size_t node_gid = ic + (jc + kc * num_nodes_1D) * num_nodes_1D;
            node_coords(node_gid, 0) = (double)ic * h;
            node_coords(node_gid, 1) = (double)jc * h;
            node_coords(node_gid, 2) = (double)kc * h;
        });

        FOR_ALL(i, 0, num_elems_1D,
                j, 0, num_elems_1D,
                k, 0, num_elems_1D, {
            const size_t elem_gid = i + (j + k * num_elems_1D) * num_elems_1D;
            size_t node_lid = 0;
            for (size_t kc = k * elem_order; kc <= k * elem_order + elem_order; kc++)
            for (size_t jc = j * elem_order; jc <= j * elem_order + elem_order; jc++)
            for (size_t ic = i * elem_order; ic <= i * elem_order + elem_order; ic++) {
                const size_t node_gid = ic + (jc + kc * num_nodes_1D) * num_nodes_1D;
                Mesh.nodes_in_elem(elem_gid, node_lid) = node_gid;
                node_lid++;
            }
        });

        Mesh.build_corner_connectivity();
        Mesh.build_elem_elem_connectivity();
        Mesh.build_surf_connectivity();

        const size_t num_qpts_in_elem = Quad.num_qpts_in_elem;
        const size_t num_qpts_in_surf = SurfQuad.num_qpts_in_surf;
        const size_t num_nodes_in_elem = Mesh.num_nodes_in_elem;
        const size_t num_surfs_in_elem = Mesh.num_surfs_in_elem;

        DCArrayKokkos<double> elem_jac(num_elems, num_qpts_in_elem, elem_dims, elem_dims, "elem_jacobian");
        DCArrayKokkos<double> elem_vol(num_elems, "elem_vol");

        DCArrayKokkos<double> vol_error(num_elems, "vol_error");
        DCArrayKokkos<double> normal_error(num_elems, "normal_error");
        DCArrayKokkos<double> gradx_error(num_elems, "gradx_error");
        DCArrayKokkos<double> div_error(num_elems, "div_error");

        FOR_ALL(elem_gid, 0, num_elems, {

            elem_vol(elem_gid) = 0.0;

            ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

            for (size_t qpt_lid = 0; qpt_lid < num_qpts_in_elem; qpt_lid++) {
                ViewCArrayKokkos<double> a_grad_basis(&FERefElem.qpt_grad_basis(qpt_lid, 0, 0),
                                                       num_nodes_in_elem, 3);
                ViewCArrayKokkos<double> jac(&elem_jac(elem_gid, qpt_lid, 0, 0), 3, 3);

                jacobian(jac, node_coords, nodes_in_elem, a_grad_basis);

                const double det_jac = det_3x3(jac);
                elem_vol(elem_gid) += det_jac * Quad.qpt_weights(qpt_lid);
            } // end for qpt

            const double vol_q = elem_vol(elem_gid);
            double vol_h = 1.0;
            for (size_t dim = 0; dim < elem_dims; dim++) vol_h *= h;
            vol_error(elem_gid) = fabs(vol_q - vol_h);

            double tally_div = 0.0;
            double tally_flux[3][3];
            double tally_normal[3];
            double eye[3][3];

            for (size_t j = 0; j < elem_dims; j++) {
                for (size_t i = 0; i < elem_dims; i++) {
                    tally_flux[i][j] = 0.0;
                    eye[i][j] = 0.0;
                }
                tally_normal[j] = 0.0;
            }
            for (size_t j = 0; j < elem_dims; j++) eye[j][j] = 1.0;

            // Per-thread scratch: each element computes its own face Jacobian
            // from its own local basis, never shared with the element on the
            // other side of the face, so these must NOT be indexed by
            // surf_gid (surf_gid is shared by up to 2 elements, and two
            // FOR_ALL iterations racing on the same surf_gid slot is a data
            // race under threaded backends).
            double jac_buf[9];
            double inv_jac_buf[9];

            for (size_t face_lid = 0; face_lid < num_surfs_in_elem; face_lid++) {

                for (size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++) {

                    ViewCArrayKokkos<double> a_grad_basis(&RefSurf.qpt_grad_basis(face_lid, qpt_lid, 0, 0),
                                                           num_nodes_in_elem, 3);
                    ViewCArrayKokkos<double> a_basis(&RefSurf.qpt_basis(face_lid, qpt_lid, 0),
                                                      num_nodes_in_elem);

                    ViewCArrayKokkos<double> jac(jac_buf, 3, 3);
                    ViewCArrayKokkos<double> inv_jac(inv_jac_buf, 3, 3);

                    jacobian(jac, node_coords, nodes_in_elem, a_grad_basis);

                    const double det_jac = det_3x3(jac);
                    invert_3x3(jac, inv_jac, det_jac);

                    // Nanson's formula: s*J^-1*j*f*w
                    double area_normal[3];
                    area_normal[0] = 0.0;
                    area_normal[1] = 0.0;
                    area_normal[2] = 0.0;
                    for (size_t j = 0; j < elem_dims; j++) {
                        for (size_t i = 0; i < elem_dims; i++) {
                            area_normal[j] += RefSurf.outward_normal(face_lid, i) * inv_jac(i, j);
                        }
                        area_normal[j] *= det_jac * SurfQuad.qpt_weights(face_lid, qpt_lid);
                    }

                    for (size_t dim = 0; dim < elem_dims; dim++) {
                        tally_normal[dim] += area_normal[dim];
                    }

                    double qpt_coords[3];
                    for (size_t dim = 0; dim < elem_dims; dim++) qpt_coords[dim] = 0.0;

                    for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
                    for (size_t dim = 0; dim < elem_dims; dim++) {
                        const size_t node_gid = nodes_in_elem(node_lid);
                        qpt_coords[dim] += a_basis(node_lid) * node_coords(node_gid, dim);
                    }

                    for (size_t dim = 0; dim < elem_dims; dim++) {
                        tally_div += area_normal[dim] * qpt_coords[dim];
                    }

                    for (size_t i = 0; i < elem_dims; i++)
                    for (size_t j = 0; j < elem_dims; j++) {
                        tally_flux[i][j] += area_normal[j] * qpt_coords[i];
                    }

                } // end for qpt_lid
            } // end for face_lid

            double max_normal = 0.0;
            for (size_t dim = 0; dim < elem_dims; dim++) {
                max_normal = fmax(max_normal, fabs(tally_normal[dim]));
            }
            normal_error(elem_gid) = max_normal;

            double max_gradx = 0.0;
            for (size_t i = 0; i < elem_dims; i++)
            for (size_t j = 0; j < elem_dims; j++) {
                max_gradx = fmax(max_gradx, fabs((tally_flux[i][j] / elem_vol(elem_gid)) - eye[i][j]));
            }
            gradx_error(elem_gid) = max_gradx;

            div_error(elem_gid) = fabs((tally_div / elem_vol(elem_gid)) - 3.0);

        }); // end parallel for elems
        Kokkos::fence();

        double max_vol_error = 0.0;
        double max_normal_tally_error = 0.0;
        double max_gradx_identity_error = 0.0;
        double max_div_error = 0.0;

        FOR_REDUCE_MAX(e, 0, num_elems, max_vol_error, {
            max_vol_error = fmax(max_vol_error, vol_error(e));
        }, max_vol_error);

        FOR_REDUCE_MAX(e, 0, num_elems, max_normal_tally_error, {
            max_normal_tally_error = fmax(max_normal_tally_error, normal_error(e));
        }, max_normal_tally_error);

        FOR_REDUCE_MAX(e, 0, num_elems, max_gradx_identity_error, {
            max_gradx_identity_error = fmax(max_gradx_identity_error, gradx_error(e));
        }, max_gradx_identity_error);

        FOR_REDUCE_MAX(e, 0, num_elems, max_div_error, {
            max_div_error = fmax(max_div_error, div_error(e));
        }, max_div_error);
        Kokkos::fence();

        ElementGeometryMetrics result;
        result.max_vol_error = max_vol_error;
        result.max_normal_tally_error = max_normal_tally_error;
        result.max_gradx_identity_error = max_gradx_identity_error;
        result.max_div_error = max_div_error;
        return result;
    }
} // namespace test_element_geometry_detail

class ElementGeometryTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        const size_t num_elems_1D = 4;
        metrics_ = test_element_geometry_detail::compute_element_geometry_metrics(num_elems_1D);
    }

    test_element_geometry_detail::ElementGeometryMetrics metrics_;
};

TEST_F(ElementGeometryTest, ElementVolumeMatchesGaussQuadrature)
{
    EXPECT_NEAR(metrics_.max_vol_error, 0.0, 1e-12);
}

TEST_F(ElementGeometryTest, FaceAreaNormalsSumToZero)
{
    EXPECT_NEAR(metrics_.max_normal_tally_error, 0.0, 1e-12);
}

TEST_F(ElementGeometryTest, OuterProductOfNormalAndPositionIsIdentity)
{
    EXPECT_NEAR(metrics_.max_gradx_identity_error, 0.0, 1e-12);
}

TEST_F(ElementGeometryTest, DivergenceTheoremHolds)
{
    EXPECT_NEAR(metrics_.max_div_error, 0.0, 1e-10);
}
