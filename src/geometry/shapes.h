#ifndef SHAPES_H
#define SHAPES_H

#include "matar.h"
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
    double epsilon = 1e-7; ///< Epsilon for floating point comparisons

    ViewCArrayKokkos<double> position; ///< Position of the plane
    ViewCArrayKokkos<double> normal; ///< Normal vector of the plane, normalize internally

    Plane(const ViewCArrayKokkos<double>& position_inp, const ViewCArrayKokkos<double>& normal_inp){

        assert(position_inp.order() == 1 && "Plane constructor requires a 1D ViewCArrayKokkos<double> for position");
        assert(normal_inp.order() == 1 && "Plane constructor requires a 1D ViewCArrayKokkos<double> for normal");
        assert(position_inp.extent() == normal_inp.extent() && "Plane constructor requires the position and normal vectors to have the same number of dimensions");
        assert(position_inp.extent() >= 2 && "Plane constructor requires at least 2 dimensions");
        
        num_dims = position_inp.extent();
        position = position_inp;
        normal = normal_inp;

        // Normalize the normal vector (just to be safe)
        double norm = 0.0;
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
    bool is_point_on_plane(const ViewCArrayKokkos<double>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<double> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");

        double dot_product = 0.0;
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
    double radius = 0.0; ///< Radius of the circle
    double epsilon = 1e-7; ///< Epsilon for floating point comparisons

    ViewCArrayKokkos<double> position; ///< Position of the circle
    ViewCArrayKokkos<double> normal; ///< Normal vector of the circle, normalize internally

    Circle(const ViewCArrayKokkos<double>& position_inp, const ViewCArrayKokkos<double>& normal_inp, const double radius_inp){

        assert(position_inp.order() == 1 && "Circle constructor requires a 1D ViewCArrayKokkos<double> for position");
        assert(normal_inp.order() == 1 && "Circle constructor requires a 1D ViewCArrayKokkos<double> for normal");
        assert(position_inp.extent() == normal_inp.extent() && "Circle constructor requires the position and normal vectors to have the same number of dimensions");
        assert(radius_inp > 0.0 && "Circle constructor requires a positive radius");
        assert(position_inp.extent() >= 2 && "Circle constructor requires at least 2 dimensions");
            

        radius = radius_inp;
        position = position_inp;
        normal = normal_inp;
        num_dims = position_inp.extent();

        // Normalize the normal vector (just to be safe)
        double norm = 0.0;
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
    bool is_point_on_circle(const ViewCArrayKokkos<double>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<double> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");

        double dot_product = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            dot_product += (point(i) - position(i)) * normal(i);
        }

        double distance = 0.0;
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
    double radius = 0.0; ///< Radius of the sphere
    double epsilon = 1e-7; ///< Epsilon for floating point comparisons

    ViewCArrayKokkos<double> position; ///< Position of the sphere

    Sphere(const ViewCArrayKokkos<double>& position_inp, const double radius_inp){

        assert(position_inp.order() == 1 && "Sphere constructor requires a 1D ViewCArrayKokkos<double> for position");
        assert(radius_inp > 0.0 && "Sphere constructor requires a positive radius");
        assert(position_inp.extent() == 3 && "Sphere constructor requires a 3D ViewCArrayKokkos<double> for position");
            
        radius = radius_inp;
        position = position_inp;
    

    }

    KOKKOS_INLINE_FUNCTION
    bool is_point_on_sphere(const ViewCArrayKokkos<double>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<double> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");


        double distance = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            distance += (point(i) - position(i)) * (point(i) - position(i));
        }
        distance = std::sqrt(distance);

        return std::abs(distance - radius) < epsilon;
    }

    KOKKOS_INLINE_FUNCTION
    bool is_point_inside_sphere(const ViewCArrayKokkos<double>& point) const
    {
        assert(point.order() == 1 && "is_point_on_plane requires a 1D ViewCArrayKokkos<double> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");


        double distance = 0.0;
        for (int i = 0; i < num_dims; ++i) {
            distance += (point(i) - position(i)) * (point(i) - position(i));
        }
        distance = std::sqrt(distance);

        return distance < radius;
    }

    KOKKOS_INLINE_FUNCTION
    double signed_distance_to_sphere(const ViewCArrayKokkos<double>& point) const
    {
        assert(point.order() == 1 && "signed_distance_to_sphere requires a 1D ViewCArrayKokkos<double> for point");
        assert(point.extent() == num_dims && "Point must have same num_dims as the plane");

        double distance = 0.0;
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