#ifndef SHAPES_H
#define SHAPES_H

#include "matar.h"
#include "utils/precision_utils.h"
#include <cmath>


namespace geometry{

///
/// \class Plane
///
/// \brief Represents an infinite plane in space
///
class Plane
{
public:

    int num_dims = 0; ///< Number of dimensions of the plane
    real_t epsilon = real_tol(1e-7); ///< Epsilon for floating point comparisons

    ViewCArrayKokkos<real_t> position; ///< Position of the plane
    ViewCArrayKokkos<real_t> normal; ///< Normal vector of the plane, normalize internally

    Plane(const ViewCArrayKokkos<real_t>& position_inp, const ViewCArrayKokkos<real_t>& normal_inp){

        assert(position_inp.order() == 1 && "Plane constructor requires a 1D ViewCArrayKokkos<real_t> for position");
        assert(normal_inp.order() == 1 && "Plane constructor requires a 1D ViewCArrayKokkos<real_t> for normal");
        assert(position_inp.extent() == normal_inp.extent() && "Plane constructor requires the position and normal vectors to have the same number of dimensions");
        assert(position_inp.extent() >= 2 && "Plane constructor requires at least 2 dimensions");
        
        num_dims = position_inp.extent();
        position = position_inp;
        normal = normal_inp;

        // Normalize the normal vector (just to be safe)
        real_t norm = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            norm += normal_inp(i) * normal_inp(i);
        }
        norm = std::sqrt(norm);
        if (norm > epsilon) {
            for (int i = 0; i < num_dims; ++i) {
                normal(i) = normal_inp(i) / norm;
            }
        } else {
            throw std::runtime_error("Plane constructor requires a non-zero normal vector");
        }
    }

    KOKKOS_INLINE_FUNCTION
    bool is_point_on_plane(const ViewCArrayKokkos<real_t>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<real_t> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");

        real_t dot_product = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            dot_product += (point(i) - position(i)) * normal(i);
        }
        return std::abs(dot_product) < epsilon;
    }

    ~Plane() = default;
};


///
/// \class Circle
///
/// \brief Represents a circle in space
///
class Circle
{
public:

    int num_dims = 0; ///< Number of dimensions of the circle
    real_t radius = 0.0; ///< Radius of the circle
    real_t epsilon = real_tol(1e-7); ///< Epsilon for floating point comparisons

    ViewCArrayKokkos<real_t> position; ///< Position of the circle
    ViewCArrayKokkos<real_t> normal; ///< Normal vector of the circle, normalize internally

    Circle(const ViewCArrayKokkos<real_t>& position_inp, const ViewCArrayKokkos<real_t>& normal_inp, const real_t radius_inp){

        assert(position_inp.order() == 1 && "Circle constructor requires a 1D ViewCArrayKokkos<real_t> for position");
        assert(normal_inp.order() == 1 && "Circle constructor requires a 1D ViewCArrayKokkos<real_t> for normal");
        assert(position_inp.extent() == normal_inp.extent() && "Circle constructor requires the position and normal vectors to have the same number of dimensions");
        assert(radius_inp > 0.0 && "Circle constructor requires a positive radius");
        assert(position_inp.extent() >= 2 && "Circle constructor requires at least 2 dimensions");
            

        radius = radius_inp;
        position = position_inp;
        normal = normal_inp;
        num_dims = position_inp.extent();

        // Normalize the normal vector (just to be safe)
        real_t norm = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            norm += normal_inp(i) * normal_inp(i);
        }
        norm = std::sqrt(norm);
        if (norm > epsilon) {
            for (int i = 0; i < num_dims; ++i) {
                normal(i) = normal_inp(i) / norm;
            }
        } else {
            normal(0) = 1.0;
            normal(1) = 0.0;
            normal(2) = 0.0;
        }
    }

    KOKKOS_INLINE_FUNCTION
    bool is_point_on_circle(const ViewCArrayKokkos<real_t>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<real_t> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");

        real_t dot_product = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            dot_product += (point(i) - position(i)) * normal(i);
        }

        real_t distance = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            distance += (point(i) - position(i)) * (point(i) - position(i));
        }
        distance = std::sqrt(distance);
        return std::abs(dot_product) < epsilon && distance < radius;
    }

    ~Circle() = default;
};


///
/// \class Sphere
///
/// \brief Represents a sphere in 3D space
///
class Sphere
{
public:

    int num_dims = 3; ///< Number of dimensions 
    real_t radius = 0.0; ///< Radius of the sphere
    real_t epsilon = real_tol(1e-7); ///< Epsilon for floating point comparisons

    ViewCArrayKokkos<real_t> position; ///< Position of the sphere

    Sphere(const ViewCArrayKokkos<real_t>& position_inp, const real_t radius_inp){

        assert(position_inp.order() == 1 && "Sphere constructor requires a 1D ViewCArrayKokkos<real_t> for position");
        assert(radius_inp > 0.0 && "Sphere constructor requires a positive radius");
        assert(position_inp.extent() == 3 && "Sphere constructor requires a 3D ViewCArrayKokkos<real_t> for position");
            
        radius = radius_inp;
        position = position_inp;
    

    }

    KOKKOS_INLINE_FUNCTION
    bool is_point_on_sphere(const ViewCArrayKokkos<real_t>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<real_t> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");


        real_t distance = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            distance += (point(i) - position(i)) * (point(i) - position(i));
        }
        distance = std::sqrt(distance);

        return std::abs(distance - radius) < epsilon;
    }

    KOKKOS_INLINE_FUNCTION
    bool is_point_inside_sphere(const ViewCArrayKokkos<real_t>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<real_t> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");


        real_t distance = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            distance += (point(i) - position(i)) * (point(i) - position(i));
        }
        distance = std::sqrt(distance);

        return distance < radius;
    }

    KOKKOS_INLINE_FUNCTION
    real_t signed_distance_to_sphere(const ViewCArrayKokkos<real_t>& point) const
    {
        assert(point.order() == 1 && "signed_distance_to_sphere requires a 1D ViewCArrayKokkos<real_t> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");

        real_t distance = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            distance += (point(i) - position(i)) * (point(i) - position(i));
        }
        distance = std::sqrt(distance);
        return distance - radius;
    }

    ~Sphere() = default;
};


} // end namespace geometry



#endif // SHAPES_H