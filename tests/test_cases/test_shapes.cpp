#include <gtest/gtest.h>
#include "matar.h"

using namespace mtr;

#include "geometry/shapes.h"

TEST(Shapes, PlanePointOnPlane)
{
    double pos[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 1.0};
    ViewCArrayKokkos<double> position(pos, 3);
    ViewCArrayKokkos<double> plane_normal(normal, 3);
    geometry::Plane plane(position, plane_normal);

    double on_pt[3] = {1.0, 2.0, 0.0};
    ViewCArrayKokkos<double> point_on(on_pt, 3);
    EXPECT_TRUE(plane.is_point_on_plane(point_on));

    double off_pt[3] = {1.0, 2.0, 5.0};
    ViewCArrayKokkos<double> point_off(off_pt, 3);
    EXPECT_FALSE(plane.is_point_on_plane(point_off));
}

TEST(Shapes, PlaneNormalizesNormal)
{
    double pos[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 5.0};
    ViewCArrayKokkos<double> position(pos, 3);
    ViewCArrayKokkos<double> plane_normal(normal, 3);
    geometry::Plane plane(position, plane_normal);

    EXPECT_NEAR(plane.normal(2), 1.0, 1e-12);
}

TEST(Shapes, CirclePointOnCircle)
{
    double pos[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 1.0};
    ViewCArrayKokkos<double> position(pos, 3);
    ViewCArrayKokkos<double> circle_normal(normal, 3);
    geometry::Circle circle(position, circle_normal, 2.0);

    // in-plane, within radius
    double on_pt[3] = {1.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_on(on_pt, 3);
    EXPECT_TRUE(circle.is_point_on_circle(point_on));

    // in-plane, but outside radius
    double far_pt[3] = {3.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_far(far_pt, 3);
    EXPECT_FALSE(circle.is_point_on_circle(point_far));

    // within radius in-plane distance, but off-plane
    double off_pt[3] = {1.0, 0.0, 1.0};
    ViewCArrayKokkos<double> point_off(off_pt, 3);
    EXPECT_FALSE(circle.is_point_on_circle(point_off));
}

TEST(Shapes, SphereIsPointOnSphere)
{
    double pos[3] = {0.0, 0.0, 0.0};
    ViewCArrayKokkos<double> position(pos, 3);
    geometry::Sphere sphere(position, 2.0);

    double on_pt[3] = {2.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_on(on_pt, 3);
    EXPECT_TRUE(sphere.is_point_on_sphere(point_on));

    double in_pt[3] = {1.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_in(in_pt, 3);
    EXPECT_FALSE(sphere.is_point_on_sphere(point_in));
}

TEST(Shapes, SphereIsPointInsideSphere)
{
    double pos[3] = {0.0, 0.0, 0.0};
    ViewCArrayKokkos<double> position(pos, 3);
    geometry::Sphere sphere(position, 2.0);

    double in_pt[3] = {1.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_in(in_pt, 3);
    EXPECT_TRUE(sphere.is_point_inside_sphere(point_in));

    double on_pt[3] = {2.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_on(on_pt, 3);
    EXPECT_FALSE(sphere.is_point_inside_sphere(point_on));

    double out_pt[3] = {3.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_out(out_pt, 3);
    EXPECT_FALSE(sphere.is_point_inside_sphere(point_out));
}

TEST(Shapes, SphereSignedDistance)
{
    double pos[3] = {0.0, 0.0, 0.0};
    ViewCArrayKokkos<double> position(pos, 3);
    geometry::Sphere sphere(position, 2.0);

    double on_pt[3] = {2.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_on(on_pt, 3);
    EXPECT_NEAR(sphere.signed_distance_to_sphere(point_on), 0.0, 1e-12);

    double in_pt[3] = {1.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_in(in_pt, 3);
    EXPECT_NEAR(sphere.signed_distance_to_sphere(point_in), -1.0, 1e-12);

    double out_pt[3] = {5.0, 0.0, 0.0};
    ViewCArrayKokkos<double> point_out(out_pt, 3);
    EXPECT_NEAR(sphere.signed_distance_to_sphere(point_out), 3.0, 1e-12);
}
