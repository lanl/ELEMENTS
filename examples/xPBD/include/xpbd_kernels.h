#ifndef XPBD_KERNELS_H
#define XPBD_KERNELS_H

#include "matar.h"

using namespace mtr;

namespace xpbd {
\

KOKKOS_INLINE_FUNCTION
void invert_3x3_matrix(
    const double* input_matrix,
    double* inverse_matrix,
    double& matrix_determinant)
{
    const double a00 = input_matrix[0];
    const double a01 = input_matrix[1];
    const double a02 = input_matrix[2];
    const double a10 = input_matrix[3];
    const double a11 = input_matrix[4];
    const double a12 = input_matrix[5];
    const double a20 = input_matrix[6];
    const double a21 = input_matrix[7];
    const double a22 = input_matrix[8];

    const double cofactor_00 = a11 * a22 - a12 * a21;
    const double cofactor_01 = a12 * a20 - a10 * a22;
    const double cofactor_02 = a10 * a21 - a11 * a20;
    const double cofactor_10 = a02 * a21 - a01 * a22;
    const double cofactor_11 = a00 * a22 - a02 * a20;
    const double cofactor_12 = a01 * a20 - a00 * a21;
    const double cofactor_20 = a01 * a12 - a02 * a11;
    const double cofactor_21 = a02 * a10 - a00 * a12;
    const double cofactor_22 = a00 * a11 - a01 * a10;

    matrix_determinant = a00 * cofactor_00 + a01 * cofactor_10 + a02 * cofactor_20;

    if (std::fabs(matrix_determinant) < 1.0e-14) {
        matrix_determinant = 0.0;
        for (int index = 0; index < 9; ++index) {
            inverse_matrix[index] = 0.0;
        }
        return;
    }

    const double inverse_det = 1.0 / matrix_determinant;
    inverse_matrix[0] = cofactor_00 * inverse_det;
    inverse_matrix[1] = cofactor_10 * inverse_det;
    inverse_matrix[2] = cofactor_20 * inverse_det;
    inverse_matrix[3] = cofactor_01 * inverse_det;
    inverse_matrix[4] = cofactor_11 * inverse_det;
    inverse_matrix[5] = cofactor_21 * inverse_det;
    inverse_matrix[6] = cofactor_02 * inverse_det;
    inverse_matrix[7] = cofactor_12 * inverse_det;
    inverse_matrix[8] = cofactor_22 * inverse_det;
}

template<typename NodePositionsArray>
KOKKOS_INLINE_FUNCTION
void gather_element_node_positions(
    const NodePositionsArray& node_positions,
    const DCArrayKokkos<size_t>& nodes_in_element,
    const size_t num_nodes_in_element,
    double* element_node_positions)
{
    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        const size_t node_gid = nodes_in_element(node_lid);
        element_node_positions[3 * node_lid + 0] = node_positions(node_gid, 0);
        element_node_positions[3 * node_lid + 1] = node_positions(node_gid, 1);
        element_node_positions[3 * node_lid + 2] = node_positions(node_gid, 2);
    }
}

template<typename NodePositionsArray>
KOKKOS_INLINE_FUNCTION
void gather_element_reference_positions(
    const NodePositionsArray& reference_node_positions,
    const DCArrayKokkos<size_t>& nodes_in_element,
    const size_t num_nodes_in_element,
    double* element_reference_positions)
{
    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        const size_t node_gid = nodes_in_element(node_lid);
        element_reference_positions[3 * node_lid + 0] = reference_node_positions(node_gid, 0);
        element_reference_positions[3 * node_lid + 1] = reference_node_positions(node_gid, 1);
        element_reference_positions[3 * node_lid + 2] = reference_node_positions(node_gid, 2);
    }
}

KOKKOS_INLINE_FUNCTION
void compute_reference_to_physical_jacobian(
    const double* element_node_positions,
    const size_t num_nodes_in_element,
    const size_t quadrature_point_lid,
    const CArrayKokkos<double>& lobatto_point_grad_basis,
    const size_t num_basis,
    double* reference_to_physical_jacobian,
    double& jacobian_determinant)
{
    for (int row = 0; row < 9; ++row) {
        reference_to_physical_jacobian[row] = 0.0;
    }

    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            for (int reference_dim = 0; reference_dim < 3; ++reference_dim) {
                const double basis_gradient = lobatto_point_grad_basis(
                    quadrature_point_lid, node_lid, reference_dim);
                reference_to_physical_jacobian[3 * spatial_dim + reference_dim] +=
                    element_node_positions[3 * node_lid + spatial_dim] * basis_gradient;
            }
        }
    }

    const double a00 = reference_to_physical_jacobian[0];
    const double a01 = reference_to_physical_jacobian[1];
    const double a02 = reference_to_physical_jacobian[2];
    const double a10 = reference_to_physical_jacobian[3];
    const double a11 = reference_to_physical_jacobian[4];
    const double a12 = reference_to_physical_jacobian[5];
    const double a20 = reference_to_physical_jacobian[6];
    const double a21 = reference_to_physical_jacobian[7];
    const double a22 = reference_to_physical_jacobian[8];

    jacobian_determinant = a00 * (a11 * a22 - a12 * a21)
        - a01 * (a10 * a22 - a12 * a20)
        + a02 * (a10 * a21 - a11 * a20);
}

KOKKOS_INLINE_FUNCTION
void compute_physical_basis_gradients(
    const size_t num_nodes_in_element,
    const size_t quadrature_point_lid,
    const CArrayKokkos<double>& lobatto_point_grad_basis,
    const double* reference_to_physical_jacobian,
    const double jacobian_determinant,
    double* physical_basis_gradient)
{
    double inverse_jacobian[9];
    double ignored_det = jacobian_determinant;
    invert_3x3_matrix(reference_to_physical_jacobian, inverse_jacobian, ignored_det);

    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            physical_basis_gradient[3 * node_lid + spatial_dim] = 0.0;
            for (int reference_dim = 0; reference_dim < 3; ++reference_dim) {
                const double reference_gradient = lobatto_point_grad_basis(
                    quadrature_point_lid, node_lid, reference_dim);
                physical_basis_gradient[3 * node_lid + spatial_dim] +=
                    inverse_jacobian[3 * reference_dim + spatial_dim] * reference_gradient;
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION
void compute_deformation_gradient(
    const double* element_node_positions,
    const double* element_reference_positions,
    const size_t num_nodes_in_element,
    const double* physical_basis_gradient,
    double* deformation_gradient)
{
    for (int row = 0; row < 9; ++row) {
        deformation_gradient[row] = 0.0;
    }

    deformation_gradient[0] = 1.0;
    deformation_gradient[4] = 1.0;
    deformation_gradient[8] = 1.0;

    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            const double displacement = element_node_positions[3 * node_lid + spatial_dim]
                - element_reference_positions[3 * node_lid + spatial_dim];
            for (int reference_dim = 0; reference_dim < 3; ++reference_dim) {
                deformation_gradient[3 * spatial_dim + reference_dim] +=
                    displacement * physical_basis_gradient[3 * node_lid + reference_dim];
            }
        }
    }
}

KOKKOS_INLINE_FUNCTION
void compute_stvk_strain_energy_density(
    const double* deformation_gradient,
    const double lame_lambda,
    const double lame_mu,
    double& strain_energy_density)
{
    const double f00 = deformation_gradient[0];
    const double f01 = deformation_gradient[1];
    const double f02 = deformation_gradient[2];
    const double f10 = deformation_gradient[3];
    const double f11 = deformation_gradient[4];
    const double f12 = deformation_gradient[5];
    const double f20 = deformation_gradient[6];
    const double f21 = deformation_gradient[7];
    const double f22 = deformation_gradient[8];

    const double c00 = f00 * f00 + f10 * f10 + f20 * f20;
    const double c11 = f01 * f01 + f11 * f11 + f21 * f21;
    const double c22 = f02 * f02 + f12 * f12 + f22 * f22;
    const double c01 = f00 * f01 + f10 * f11 + f20 * f21;
    const double c02 = f00 * f02 + f10 * f12 + f20 * f22;
    const double c12 = f01 * f02 + f11 * f12 + f21 * f22;

    const double e00 = 0.5 * (c00 - 1.0);
    const double e11 = 0.5 * (c11 - 1.0);
    const double e22 = 0.5 * (c22 - 1.0);
    const double e01 = 0.5 * c01;
    const double e02 = 0.5 * c02;
    const double e12 = 0.5 * c12;

    const double trace_strain = e00 + e11 + e22;
    const double trace_strain_squared = trace_strain * trace_strain;
    const double strain_double_contraction =
        e00 * e00 + e11 * e11 + e22 * e22 + 2.0 * (e01 * e01 + e02 * e02 + e12 * e12);

    strain_energy_density =
        lame_mu * trace_strain_squared
        + lame_lambda * trace_strain_squared
        + lame_mu * strain_double_contraction;
}

KOKKOS_INLINE_FUNCTION
void compute_stvk_first_piola_stress(
    const double* deformation_gradient,
    const double lame_lambda,
    const double lame_mu,
    double* first_piola_kirchhoff_stress)
{
    const double f00 = deformation_gradient[0];
    const double f01 = deformation_gradient[1];
    const double f02 = deformation_gradient[2];
    const double f10 = deformation_gradient[3];
    const double f11 = deformation_gradient[4];
    const double f12 = deformation_gradient[5];
    const double f20 = deformation_gradient[6];
    const double f21 = deformation_gradient[7];
    const double f22 = deformation_gradient[8];

    const double c00 = f00 * f00 + f10 * f10 + f20 * f20;
    const double c11 = f01 * f01 + f11 * f11 + f21 * f21;
    const double c22 = f02 * f02 + f12 * f12 + f22 * f22;
    const double c01 = f00 * f01 + f10 * f11 + f20 * f21;
    const double c02 = f00 * f02 + f10 * f12 + f20 * f22;
    const double c12 = f01 * f02 + f11 * f12 + f21 * f22;

    const double e00 = 0.5 * (c00 - 1.0);
    const double e11 = 0.5 * (c11 - 1.0);
    const double e22 = 0.5 * (c22 - 1.0);
    const double e01 = 0.5 * c01;
    const double e02 = 0.5 * c02;
    const double e12 = 0.5 * c12;
    const double trace_strain = e00 + e11 + e22;

    const double s00 = lame_lambda * trace_strain + 2.0 * lame_mu * e00;
    const double s11 = lame_lambda * trace_strain + 2.0 * lame_mu * e11;
    const double s22 = lame_lambda * trace_strain + 2.0 * lame_mu * e22;
    const double s01 = 2.0 * lame_mu * e01;
    const double s02 = 2.0 * lame_mu * e02;
    const double s12 = 2.0 * lame_mu * e12;

    first_piola_kirchhoff_stress[0] = f00 * s00 + f01 * s01 + f02 * s02;
    first_piola_kirchhoff_stress[1] = f00 * s01 + f01 * s11 + f02 * s12;
    first_piola_kirchhoff_stress[2] = f00 * s02 + f01 * s12 + f02 * s22;
    first_piola_kirchhoff_stress[3] = f10 * s00 + f11 * s01 + f12 * s02;
    first_piola_kirchhoff_stress[4] = f10 * s01 + f11 * s11 + f12 * s12;
    first_piola_kirchhoff_stress[5] = f10 * s02 + f11 * s12 + f12 * s22;
    first_piola_kirchhoff_stress[6] = f20 * s00 + f21 * s01 + f22 * s02;
    first_piola_kirchhoff_stress[7] = f20 * s01 + f21 * s11 + f22 * s12;
    first_piola_kirchhoff_stress[8] = f20 * s02 + f21 * s12 + f22 * s22;
}

KOKKOS_INLINE_FUNCTION
void compute_constraint_at_quadrature_point(
    const double quadrature_weight_volume,
    const double strain_energy_density,
    double& constraint_value)
{
    constraint_value = quadrature_weight_volume * strain_energy_density;
}

KOKKOS_INLINE_FUNCTION
void compute_constraint_gradient_nodal(
    const double quadrature_weight_volume,
    const double* first_piola_kirchhoff_stress,
    const double* physical_basis_gradient,
    const size_t num_nodes_in_element,
    double* constraint_gradient_nodal)
{
    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            double gradient_component = 0.0;
            for (int reference_dim = 0; reference_dim < 3; ++reference_dim) {
                gradient_component += first_piola_kirchhoff_stress[3 * spatial_dim + reference_dim]
                    * physical_basis_gradient[3 * node_lid + reference_dim];
            }
            constraint_gradient_nodal[3 * node_lid + spatial_dim] =
                quadrature_weight_volume * gradient_component;
        }
    }
}

KOKKOS_INLINE_FUNCTION
double compute_constraint_denominator(
    const double* constraint_gradient_nodal,
    const DCArrayKokkos<double>& inverse_lumped_mass,
    const DCArrayKokkos<size_t>& nodes_in_element,
    const size_t num_nodes_in_element)
{
    double constraint_denominator = 0.0;
    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        const size_t node_gid = nodes_in_element(node_lid);
        const double inverse_mass = inverse_lumped_mass(node_gid);
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            const double gradient_component = constraint_gradient_nodal[3 * node_lid + spatial_dim];
            constraint_denominator += inverse_mass * gradient_component * gradient_component;
        }
    }
    return constraint_denominator;
}

template<typename NodePositionsArray>
KOKKOS_INLINE_FUNCTION
double compute_constraint_damping_numerator(
    const double* constraint_gradient_nodal,
    const NodePositionsArray& node_positions,
    const NodePositionsArray& node_positions_previous,
    const DCArrayKokkos<size_t>& nodes_in_element,
    const size_t num_nodes_in_element)
{
    double damping_numerator = 0.0;
    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        const size_t node_gid = nodes_in_element(node_lid);
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            const double gradient_component = constraint_gradient_nodal[3 * node_lid + spatial_dim];
            const double position_change = node_positions(node_gid, spatial_dim)
                - node_positions_previous(node_gid, spatial_dim);
            damping_numerator += gradient_component * position_change;
        }
    }
    return damping_numerator;
}

template<typename NodePositionsArray>
KOKKOS_INLINE_FUNCTION
void apply_constraint_position_correction(
    const double constraint_multiplier_delta,
    const double* constraint_gradient_nodal,
    const DCArrayKokkos<double>& inverse_lumped_mass,
    const DCArrayKokkos<size_t>& nodes_in_element,
    const size_t num_nodes_in_element,
    const size_t num_owned_nodes,
    NodePositionsArray& node_positions)
{
    for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
        const size_t node_gid = nodes_in_element(node_lid);
        if (node_gid >= num_owned_nodes) {
            continue;
        }
        const double inverse_mass = inverse_lumped_mass(node_gid);
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            node_positions(node_gid, spatial_dim) +=
                inverse_mass * constraint_gradient_nodal[3 * node_lid + spatial_dim]
                * constraint_multiplier_delta;
        }
    }
}

KOKKOS_INLINE_FUNCTION
void evaluate_constraint_at_quadrature_point(
    const double* element_node_positions,
    const double* element_reference_positions,
    const size_t num_nodes_in_element,
    const size_t quadrature_point_lid,
    const CArrayKokkos<double>& lobatto_point_grad_basis,
    const CArrayKokkos<double>& lobatto_point_weights,
    const double lame_lambda,
    const double lame_mu,
    double& constraint_value,
    double* constraint_gradient_nodal)
{
    double reference_to_physical_jacobian[9];
    double jacobian_determinant = 0.0;
    compute_reference_to_physical_jacobian(
        element_node_positions,
        num_nodes_in_element,
        quadrature_point_lid,
        lobatto_point_grad_basis,
        num_nodes_in_element,
        reference_to_physical_jacobian,
        jacobian_determinant);

    double physical_basis_gradient[XPBD_MAX_NODE_BUFFER];
    compute_physical_basis_gradients(
        num_nodes_in_element,
        quadrature_point_lid,
        lobatto_point_grad_basis,
        reference_to_physical_jacobian,
        jacobian_determinant,
        physical_basis_gradient);

    double deformation_gradient[9];
    compute_deformation_gradient(
        element_node_positions,
        element_reference_positions,
        num_nodes_in_element,
        physical_basis_gradient,
        deformation_gradient);

    double strain_energy_density = 0.0;
    compute_stvk_strain_energy_density(
        deformation_gradient, lame_lambda, lame_mu, strain_energy_density);

    const double quadrature_weight_volume =
        lobatto_point_weights(quadrature_point_lid) * jacobian_determinant;

    compute_constraint_at_quadrature_point(
        quadrature_weight_volume, strain_energy_density, constraint_value);

    double first_piola_kirchhoff_stress[9];
    compute_stvk_first_piola_stress(
        deformation_gradient, lame_lambda, lame_mu, first_piola_kirchhoff_stress);

    compute_constraint_gradient_nodal(
        quadrature_weight_volume,
        first_piola_kirchhoff_stress,
        physical_basis_gradient,
        num_nodes_in_element,
        constraint_gradient_nodal);
}

} // namespace xpbd

#endif
