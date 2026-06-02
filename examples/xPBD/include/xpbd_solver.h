#ifndef XPBD_SOLVER_H
#define XPBD_SOLVER_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "ELEMENTS.h"
#include "mesh_io.h"
#include "state.h"
#include "xpbd_kernels.h"

using namespace mtr;

struct XpbdMaterialProperties {
    double density = 1.0;
    double youngs_modulus = 1.0e6;
    double poissons_ratio = 0.3;
    double lame_lambda = 0.0;
    double lame_mu = 0.0;
    double compliance = 1.0e-9;
    double constraint_damping = 0.0;
};

struct XpbdSolverConfig {
    int poly_order = 1;
    int num_timesteps = 100;
    int num_solver_iterations = 20;
    int num_force_ramp_steps = 10;
    int num_outputs = 10;
    int output_interval = 5;
    double timestep_dt = 1.0e-3;
    double velocity_damping_factor = 0.01;
    double dirichlet_plane_tolerance = 1.0e-8;
    double dirichlet_plane_origin[3] = {0.0, 0.0, 0.0};
    double dirichlet_plane_normal[3] = {0.0, 0.0, 1.0};
    double gravity_vector[3] = {0.0, 0.0, -9.81};
    double surface_traction[3] = {0.0, 0.0, 0.0};
    double surface_traction_plane_normal[3] = {0.0, 0.0, 1.0};
    XpbdMaterialProperties material;
};

struct XpbdSolverState {
    elements::fe_ref_elem_t ref_elem;
    XpbdSolverConfig config;

    size_t num_quadrature_points = 0;
    size_t num_quadrature_points_1d = 0;
    size_t num_nodes_per_element = 0;
    size_t num_colors = 0;

    MPICArrayKokkos<double> node_positions;
    MPICArrayKokkos<double> reference_node_positions;
    MPICArrayKokkos<double> predicted_position;
    MPICArrayKokkos<double> nodal_velocity;
    MPICArrayKokkos<double> node_positions_previous;
    MPICArrayKokkos<double> nodal_force_external;

    DCArrayKokkos<double> lumped_mass;
    DCArrayKokkos<double> inverse_lumped_mass;
    DCArrayKokkos<double> constraint_multipliers;

    CArrayKokkos<size_t> color_element_offsets;
    CArrayKokkos<size_t> color_element_list;

    node_t node;
    GaussPoint_t gauss_point;
};

inline void compute_lame_parameters(XpbdMaterialProperties& material)
{
    const double poissons_ratio = material.poissons_ratio;
    material.lame_mu = material.youngs_modulus / (2.0 * (1.0 + poissons_ratio));
    material.lame_lambda = material.youngs_modulus * poissons_ratio
        / ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
}

inline void initialize_reference_element_tables(
    elements::fe_ref_elem_t& ref_elem,
    const int poly_order)
{
    ref_elem.init(poly_order, 3);
}

inline void build_element_graph_coloring(
    swage::Mesh& mesh,
    XpbdSolverState& solver_state)
{
    const size_t num_owned_elements = mesh.num_owned_elems;
    std::vector<int> element_color(num_owned_elements, -1);
    int max_color = -1;

    mesh.nodes_in_elem.update_host();

    for (size_t elem_gid = 0; elem_gid < num_owned_elements; ++elem_gid) {
        std::vector<int> neighbor_colors;
        const size_t num_neighbors = mesh.num_elems_in_elem(elem_gid);
        for (size_t neighbor_lid = 0; neighbor_lid < num_neighbors; ++neighbor_lid) {
            const size_t neighbor_elem_gid = mesh.elems_in_elem(elem_gid, neighbor_lid);
            if (neighbor_elem_gid < num_owned_elements && element_color[neighbor_elem_gid] >= 0) {
                neighbor_colors.push_back(element_color[neighbor_elem_gid]);
            }
        }

        std::sort(neighbor_colors.begin(), neighbor_colors.end());
        neighbor_colors.erase(std::unique(neighbor_colors.begin(), neighbor_colors.end()),
                              neighbor_colors.end());

        int chosen_color = 0;
        for (int neighbor_color : neighbor_colors) {
            if (neighbor_color == chosen_color) {
                ++chosen_color;
            }
        }
        element_color[elem_gid] = chosen_color;
        max_color = std::max(max_color, chosen_color);
    }

    solver_state.num_colors = static_cast<size_t>(max_color + 1);
    std::vector<size_t> color_counts(solver_state.num_colors, 0);
    for (size_t elem_gid = 0; elem_gid < num_owned_elements; ++elem_gid) {
        ++color_counts[element_color[elem_gid]];
    }

    std::vector<size_t> color_offsets(solver_state.num_colors + 1, 0);
    for (size_t color_index = 0; color_index < solver_state.num_colors; ++color_index) {
        color_offsets[color_index + 1] = color_offsets[color_index] + color_counts[color_index];
    }

    std::vector<size_t> color_element_list(num_owned_elements, 0);
    std::vector<size_t> color_fill(solver_state.num_colors, 0);
    for (size_t elem_gid = 0; elem_gid < num_owned_elements; ++elem_gid) {
        const size_t color_index = static_cast<size_t>(element_color[elem_gid]);
        const size_t write_index = color_offsets[color_index] + color_fill[color_index];
        color_element_list[write_index] = elem_gid;
        ++color_fill[color_index];
    }

    solver_state.color_element_offsets = CArrayKokkos<size_t>(
        solver_state.num_colors + 1, "xpbd_color_element_offsets");
    solver_state.color_element_list = CArrayKokkos<size_t>(
        num_owned_elements, "xpbd_color_element_list");

    for (size_t color_index = 0; color_index <= solver_state.num_colors; ++color_index) {
        solver_state.color_element_offsets(color_index) = color_offsets[color_index];
    }
    for (size_t elem_index = 0; elem_index < num_owned_elements; ++elem_index) {
        solver_state.color_element_list(elem_index) = color_element_list[elem_index];
    }
}

inline void assemble_lumped_mass(
    swage::Mesh& mesh,
    XpbdSolverState& solver_state)
{
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_owned_elements = mesh.num_owned_elems;
    const size_t num_nodes_in_element = mesh.num_nodes_in_elem;
    const double density = solver_state.config.material.density;

    solver_state.lumped_mass = DCArrayKokkos<double>(num_nodes, "xpbd_lumped_mass");
    solver_state.lumped_mass.set_values(0.0);

    const elements::fe_ref_elem_t& ref_elem = solver_state.ref_elem;

    FOR_ALL(elem_gid, 0, num_owned_elements, {
        double element_node_positions[xpbd::XPBD_MAX_NODE_BUFFER];
        xpbd::gather_element_node_positions(
            solver_state.node_positions,
            mesh.nodes_in_elem,
            num_nodes_in_element,
            element_node_positions);

        for (size_t node_lid = 0; node_lid < num_nodes_in_element; ++node_lid) {
            const size_t lobatto_point_lid = ref_elem.dof_lobatto_map(node_lid);

            double reference_to_physical_jacobian[9];
            double jacobian_determinant = 0.0;
            xpbd::compute_reference_to_physical_jacobian(
                element_node_positions,
                num_nodes_in_element,
                lobatto_point_lid,
                ref_elem.lobotto_point_grad_basis,
                num_nodes_in_element,
                reference_to_physical_jacobian,
                jacobian_determinant);

            const double quadrature_weight = ref_elem.lobotto_point_weights(lobatto_point_lid);
            const size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            const double mass_contribution = density * quadrature_weight * jacobian_determinant;
            Kokkos::atomic_add(&solver_state.lumped_mass(node_gid), mass_contribution);
        }
    });
    MATAR_FENCE();
}

inline void invert_lumped_mass(
    const DCArrayKokkos<double>& lumped_mass,
    DCArrayKokkos<double>& inverse_lumped_mass)
{
    const size_t num_nodes = lumped_mass.dims(0);
    inverse_lumped_mass = DCArrayKokkos<double>(num_nodes, "xpbd_inverse_lumped_mass");

    FOR_ALL(node_gid, 0, num_nodes, {
        const double node_mass = lumped_mass(node_gid);
        if (node_mass > 1.0e-14) {
            inverse_lumped_mass(node_gid) = 1.0 / node_mass;
        } else {
            inverse_lumped_mass(node_gid) = 0.0;
        }
    });
    MATAR_FENCE();
}

inline void apply_dirichlet_inverse_mass_mask(
    const swage::Mesh& mesh,
    const XpbdSolverConfig& config,
    const MPICArrayKokkos<double>& node_positions,
    DCArrayKokkos<double>& inverse_lumped_mass)
{
    const double tolerance = config.dirichlet_plane_tolerance;
    const double origin_x = config.dirichlet_plane_origin[0];
    const double origin_y = config.dirichlet_plane_origin[1];
    const double origin_z = config.dirichlet_plane_origin[2];
    const double normal_x = config.dirichlet_plane_normal[0];
    const double normal_y = config.dirichlet_plane_normal[1];
    const double normal_z = config.dirichlet_plane_normal[2];

    FOR_ALL(node_gid, 0, mesh.num_owned_nodes, {
        const double offset_x = node_positions(node_gid, 0) - origin_x;
        const double offset_y = node_positions(node_gid, 1) - origin_y;
        const double offset_z = node_positions(node_gid, 2) - origin_z;
        const double signed_distance = normal_x * offset_x + normal_y * offset_y + normal_z * offset_z;
        if (std::fabs(signed_distance) < tolerance) {
            inverse_lumped_mass(node_gid) = 0.0;
        }
    });
    MATAR_FENCE();
}

inline void pin_dirichlet_node_positions(
    const swage::Mesh& mesh,
    const XpbdSolverConfig& config,
    const MPICArrayKokkos<double>& reference_node_positions,
    MPICArrayKokkos<double>& node_positions)
{
    const double tolerance = config.dirichlet_plane_tolerance;
    const double origin_x = config.dirichlet_plane_origin[0];
    const double origin_y = config.dirichlet_plane_origin[1];
    const double origin_z = config.dirichlet_plane_origin[2];
    const double normal_x = config.dirichlet_plane_normal[0];
    const double normal_y = config.dirichlet_plane_normal[1];
    const double normal_z = config.dirichlet_plane_normal[2];

    FOR_ALL(node_gid, 0, mesh.num_owned_nodes, {
        const double offset_x = node_positions(node_gid, 0) - origin_x;
        const double offset_y = node_positions(node_gid, 1) - origin_y;
        const double offset_z = node_positions(node_gid, 2) - origin_z;
        const double signed_distance = normal_x * offset_x + normal_y * offset_y + normal_z * offset_z;
        if (std::fabs(signed_distance) < tolerance) {
            node_positions(node_gid, 0) = reference_node_positions(node_gid, 0);
            node_positions(node_gid, 1) = reference_node_positions(node_gid, 1);
            node_positions(node_gid, 2) = reference_node_positions(node_gid, 2);
        }
    });
    MATAR_FENCE();
}

inline void integrate_boundary_traction(
    swage::Mesh& mesh,
    const XpbdSolverConfig& config,
    XpbdSolverState& solver_state,
    const int timestep_index)
{
    const size_t num_nodes = mesh.num_nodes;
    const double ramp_factor = (config.num_force_ramp_steps > 0)
        ? std::min(1.0, static_cast<double>(timestep_index + 1) / config.num_force_ramp_steps)
        : 1.0;

    solver_state.nodal_force_external.set_values(0.0);

    FOR_ALL(node_gid, 0, mesh.num_owned_nodes, {
        const double node_mass = solver_state.lumped_mass(node_gid);
        solver_state.nodal_force_external(node_gid, 0) =
            ramp_factor * node_mass * config.gravity_vector[0];
        solver_state.nodal_force_external(node_gid, 1) =
            ramp_factor * node_mass * config.gravity_vector[1];
        solver_state.nodal_force_external(node_gid, 2) =
            ramp_factor * node_mass * config.gravity_vector[2];
    });
    MATAR_FENCE();

    const double traction_magnitude = std::sqrt(
        config.surface_traction[0] * config.surface_traction[0]
        + config.surface_traction[1] * config.surface_traction[1]
        + config.surface_traction[2] * config.surface_traction[2]);

    if (traction_magnitude < 1.0e-14) {
        solver_state.nodal_force_external.update_device();
        return;
    }

    const size_t num_quadrature_1d = solver_state.num_quadrature_points_1d;
    const double traction_x = ramp_factor * config.surface_traction[0];
    const double traction_y = ramp_factor * config.surface_traction[1];
    const double traction_z = ramp_factor * config.surface_traction[2];
    const double origin_x = config.dirichlet_plane_origin[0];
    const double origin_y = config.dirichlet_plane_origin[1];
    const double origin_z = config.dirichlet_plane_origin[2];
    const double normal_x = config.surface_traction_plane_normal[0];
    const double normal_y = config.surface_traction_plane_normal[1];
    const double normal_z = config.surface_traction_plane_normal[2];
    const double tolerance = config.dirichlet_plane_tolerance;

    solver_state.node_positions.update_host();

    for (size_t bdy_patch_gid = 0; bdy_patch_gid < mesh.num_bdy_patches; ++bdy_patch_gid) {
        const size_t patch_gid = mesh.bdy_patches(bdy_patch_gid);
        const size_t elem_gid = mesh.elems_in_patch(patch_gid, 0);
        (void)elem_gid;

        double patch_node_positions[4 * 3];
        for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; ++patch_node_lid) {
            const size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
            patch_node_positions[3 * patch_node_lid + 0] = solver_state.node_positions.host(node_gid, 0);
            patch_node_positions[3 * patch_node_lid + 1] = solver_state.node_positions.host(node_gid, 1);
            patch_node_positions[3 * patch_node_lid + 2] = solver_state.node_positions.host(node_gid, 2);
        }

        const double center_x = 0.25 * (patch_node_positions[0] + patch_node_positions[3]
            + patch_node_positions[6] + patch_node_positions[9]);
        const double center_y = 0.25 * (patch_node_positions[1] + patch_node_positions[4]
            + patch_node_positions[7] + patch_node_positions[10]);
        const double center_z = 0.25 * (patch_node_positions[2] + patch_node_positions[5]
            + patch_node_positions[8] + patch_node_positions[11]);
        const double signed_distance = normal_x * (center_x - origin_x)
            + normal_y * (center_y - origin_y)
            + normal_z * (center_z - origin_z);
        if (std::fabs(signed_distance) < tolerance) {
            continue;
        }

        for (size_t quadrature_j = 0; quadrature_j < num_quadrature_1d; ++quadrature_j) {
            for (size_t quadrature_i = 0; quadrature_i < num_quadrature_1d; ++quadrature_i) {
                const double xi = solver_state.ref_elem.lob_points_1D(quadrature_i);
                const double eta = solver_state.ref_elem.lob_points_1D(quadrature_j);
                const double weight_xi = solver_state.ref_elem.lob_weights_1D(quadrature_i);
                const double weight_eta = solver_state.ref_elem.lob_weights_1D(quadrature_j);

                const double shape_0 = 0.25 * (1.0 - xi) * (1.0 - eta);
                const double shape_1 = 0.25 * (1.0 + xi) * (1.0 - eta);
                const double shape_2 = 0.25 * (1.0 + xi) * (1.0 + eta);
                const double shape_3 = 0.25 * (1.0 - xi) * (1.0 + eta);

                double tangent_xi[3] = {0.0, 0.0, 0.0};
                double tangent_eta[3] = {0.0, 0.0, 0.0};
                const double dshape_dxi[4] = {
                    -0.25 * (1.0 - eta), 0.25 * (1.0 - eta),
                    0.25 * (1.0 + eta), -0.25 * (1.0 + eta)
                };
                const double dshape_deta[4] = {
                    -0.25 * (1.0 - xi), -0.25 * (1.0 + xi),
                    0.25 * (1.0 + xi), 0.25 * (1.0 - xi)
                };

                for (size_t patch_node_lid = 0; patch_node_lid < 4; ++patch_node_lid) {
                    for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
                        tangent_xi[spatial_dim] += dshape_dxi[patch_node_lid]
                            * patch_node_positions[3 * patch_node_lid + spatial_dim];
                        tangent_eta[spatial_dim] += dshape_deta[patch_node_lid]
                            * patch_node_positions[3 * patch_node_lid + spatial_dim];
                    }
                }

                const double normal_vector[3] = {
                    tangent_xi[1] * tangent_eta[2] - tangent_xi[2] * tangent_eta[1],
                    tangent_xi[2] * tangent_eta[0] - tangent_xi[0] * tangent_eta[2],
                    tangent_xi[0] * tangent_eta[1] - tangent_xi[1] * tangent_eta[0]
                };
                const double surface_jacobian = std::sqrt(
                    normal_vector[0] * normal_vector[0]
                    + normal_vector[1] * normal_vector[1]
                    + normal_vector[2] * normal_vector[2]);
                const double surface_measure = weight_xi * weight_eta * surface_jacobian;

                const double shape_values[4] = {shape_0, shape_1, shape_2, shape_3};
                for (size_t patch_node_lid = 0; patch_node_lid < 4; ++patch_node_lid) {
                    const size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
                    solver_state.nodal_force_external.host(node_gid, 0) +=
                        shape_values[patch_node_lid] * traction_x * surface_measure;
                    solver_state.nodal_force_external.host(node_gid, 1) +=
                        shape_values[patch_node_lid] * traction_y * surface_measure;
                    solver_state.nodal_force_external.host(node_gid, 2) +=
                        shape_values[patch_node_lid] * traction_z * surface_measure;
                }
            }
        }
    }

    solver_state.nodal_force_external.update_device();
}

inline void predict_unconstrained_positions(
    XpbdSolverState& solver_state,
    const double timestep_dt)
{
    const size_t num_nodes = solver_state.node_positions.dims(0);

    FOR_ALL(node_gid, 0, num_nodes, {
        const double inverse_mass = solver_state.inverse_lumped_mass(node_gid);
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            solver_state.predicted_position(node_gid, spatial_dim) =
                solver_state.node_positions(node_gid, spatial_dim)
                + timestep_dt * solver_state.nodal_velocity(node_gid, spatial_dim)
                + timestep_dt * timestep_dt * inverse_mass
                * solver_state.nodal_force_external(node_gid, spatial_dim);
        }
    });
    MATAR_FENCE();
}

inline void apply_velocity_damping(
    MPICArrayKokkos<double>& nodal_velocity,
    const double velocity_damping_factor)
{
    FOR_ALL(node_gid, 0, nodal_velocity.dims(0), {
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            nodal_velocity(node_gid, spatial_dim) *= (1.0 - velocity_damping_factor);
        }
    });
    MATAR_FENCE();
}

inline void solve_xpbd_constraints_for_color(
    swage::Mesh& mesh,
    XpbdSolverState& solver_state,
    const size_t color_index)
{
    const size_t color_start = solver_state.color_element_offsets(color_index);
    const size_t color_end = solver_state.color_element_offsets(color_index + 1);
    const size_t num_nodes_in_element = mesh.num_nodes_in_elem;
    const size_t num_quadrature_points = solver_state.num_quadrature_points;
    const double compliance = solver_state.config.material.compliance;
    const double constraint_damping = solver_state.config.material.constraint_damping;
    const double timestep_dt = solver_state.config.timestep_dt;
    const double timestep_compliance = compliance / (timestep_dt * timestep_dt);
    const double damping_gamma = compliance * constraint_damping / timestep_dt;
    const double lame_lambda = solver_state.config.material.lame_lambda;
    const double lame_mu = solver_state.config.material.lame_mu;
    const elements::fe_ref_elem_t& ref_elem = solver_state.ref_elem;

    FOR_ALL(color_elem_index, color_start, color_end, {
        const size_t elem_gid = solver_state.color_element_list(color_elem_index);

        double element_node_positions[xpbd::XPBD_MAX_NODE_BUFFER];
        double element_reference_positions[xpbd::XPBD_MAX_NODE_BUFFER];
        xpbd::gather_element_node_positions(
            solver_state.node_positions,
            mesh.nodes_in_elem,
            num_nodes_in_element,
            element_node_positions);
        xpbd::gather_element_reference_positions(
            solver_state.reference_node_positions,
            mesh.nodes_in_elem,
            num_nodes_in_element,
            element_reference_positions);

        for (size_t quadrature_point_lid = 0; quadrature_point_lid < num_quadrature_points; ++quadrature_point_lid) {
            double constraint_gradient_nodal[xpbd::XPBD_MAX_NODE_BUFFER];
            double constraint_value = 0.0;

            xpbd::evaluate_constraint_at_quadrature_point(
                element_node_positions,
                element_reference_positions,
                num_nodes_in_element,
                quadrature_point_lid,
                ref_elem.lobotto_point_grad_basis,
                ref_elem.lobotto_point_weights,
                lame_lambda,
                lame_mu,
                constraint_value,
                constraint_gradient_nodal);

            const double constraint_denominator = xpbd::compute_constraint_denominator(
                constraint_gradient_nodal,
                solver_state.inverse_lumped_mass,
                mesh.nodes_in_elem,
                num_nodes_in_element);

            const size_t constraint_gid = elem_gid * num_quadrature_points + quadrature_point_lid;
            const double constraint_multiplier = solver_state.constraint_multipliers(constraint_gid);

            const double damping_numerator = xpbd::compute_constraint_damping_numerator(
                constraint_gradient_nodal,
                solver_state.node_positions,
                solver_state.node_positions_previous,
                mesh.nodes_in_elem,
                num_nodes_in_element);

            const double constraint_numerator = constraint_value
                + timestep_compliance * constraint_multiplier
                - damping_gamma * damping_numerator;
            const double constraint_denominator_scaled =
                (1.0 + damping_gamma) * constraint_denominator + timestep_compliance;

            double constraint_multiplier_delta = 0.0;
            if (std::fabs(constraint_denominator_scaled) > 1.0e-20) {
                constraint_multiplier_delta = -constraint_numerator / constraint_denominator_scaled;
            }

            xpbd::apply_constraint_position_correction(
                constraint_multiplier_delta,
                constraint_gradient_nodal,
                solver_state.inverse_lumped_mass,
                mesh.nodes_in_elem,
                num_nodes_in_element,
                mesh.num_owned_nodes,
                solver_state.node_positions);

            solver_state.constraint_multipliers(constraint_gid) =
                constraint_multiplier + constraint_multiplier_delta;

            xpbd::gather_element_node_positions(
                solver_state.node_positions,
                mesh.nodes_in_elem,
                num_nodes_in_element,
                element_node_positions);
        }
    });
    MATAR_FENCE();
}

inline void update_nodal_velocities(
    XpbdSolverState& solver_state,
    const double timestep_dt)
{
    const size_t num_nodes = solver_state.node_positions.dims(0);
    const double inverse_dt = 1.0 / timestep_dt;

    FOR_ALL(node_gid, 0, num_nodes, {
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            solver_state.nodal_velocity(node_gid, spatial_dim) =
                inverse_dt * (solver_state.node_positions(node_gid, spatial_dim)
                    - solver_state.node_positions_previous(node_gid, spatial_dim));
        }
    });
    MATAR_FENCE();
}

inline void copy_mpic_array(
    const MPICArrayKokkos<double>& source,
    MPICArrayKokkos<double>& destination)
{
    const size_t num_nodes = source.dims(0);
    const int num_dims = static_cast<int>(source.dims(1));
    FOR_ALL(node_gid, 0, num_nodes, {
        for (int spatial_dim = 0; spatial_dim < num_dims; ++spatial_dim) {
            destination(node_gid, spatial_dim) = source(node_gid, spatial_dim);
        }
    });
    MATAR_FENCE();
}

inline void initialize_xpbd_solver(
    swage::Mesh& mesh,
    MPICArrayKokkos<double>& node_coords,
    CommunicationPlan& node_communication_plan,
    CommunicationPlan& element_communication_plan,
    XpbdSolverState& solver_state,
    const XpbdSolverConfig& config)
{
    solver_state.config = config;
    compute_lame_parameters(solver_state.config.material);

    initialize_reference_element_tables(solver_state.ref_elem, config.poly_order);
    solver_state.num_quadrature_points_1d = solver_state.ref_elem.num_lobotto_1d;
    solver_state.num_quadrature_points = solver_state.ref_elem.num_lobotto_in_elem;
    solver_state.num_nodes_per_element = solver_state.ref_elem.num_dofs_in_elem;

    std::vector<node_state> node_states = {
        node_state::coords,
        node_state::scalar_field,
        node_state::vector_field
    };
    solver_state.node.initialize(mesh.num_nodes, mesh.num_dims, node_states, node_communication_plan);
    solver_state.node.coords = node_coords;
    solver_state.node.scalar_field.set_values(0.0);
    solver_state.node.vector_field.set_values(0.0);

    solver_state.node_positions = node_coords;
    solver_state.reference_node_positions = MPICArrayKokkos<double>(
        mesh.num_nodes, mesh.num_dims, "xpbd_reference_node_positions");
    solver_state.reference_node_positions.initialize_comm_plan(node_communication_plan);
    copy_mpic_array(node_coords, solver_state.reference_node_positions);

    solver_state.predicted_position = MPICArrayKokkos<double>(
        mesh.num_nodes, mesh.num_dims, "xpbd_predicted_position");
    solver_state.predicted_position.initialize_comm_plan(node_communication_plan);
    copy_mpic_array(node_coords, solver_state.predicted_position);

    solver_state.nodal_velocity = MPICArrayKokkos<double>(
        mesh.num_nodes, mesh.num_dims, "xpbd_nodal_velocity");
    solver_state.nodal_velocity.initialize_comm_plan(node_communication_plan);
    solver_state.nodal_velocity.set_values(0.0);

    solver_state.node_positions_previous = MPICArrayKokkos<double>(
        mesh.num_nodes, mesh.num_dims, "xpbd_node_positions_previous");
    solver_state.node_positions_previous.initialize_comm_plan(node_communication_plan);
    copy_mpic_array(node_coords, solver_state.node_positions_previous);

    solver_state.nodal_force_external = MPICArrayKokkos<double>(
        mesh.num_nodes, mesh.num_dims, "xpbd_nodal_force_external");
    solver_state.nodal_force_external.initialize_comm_plan(node_communication_plan);
    solver_state.nodal_force_external.set_values(0.0);

    solver_state.constraint_multipliers = DCArrayKokkos<double>(
        mesh.num_owned_elems * solver_state.num_quadrature_points,
        "xpbd_constraint_multipliers");
    solver_state.constraint_multipliers.set_values(0.0);

    build_element_graph_coloring(mesh, solver_state);
    assemble_lumped_mass(mesh, solver_state);
    invert_lumped_mass(solver_state.lumped_mass, solver_state.inverse_lumped_mass);
    apply_dirichlet_inverse_mass_mask(
        mesh, config, solver_state.node_positions, solver_state.inverse_lumped_mass);

    std::vector<gauss_pt_state> gauss_pt_states = {
        gauss_pt_state::fields,
        gauss_pt_state::fields_vec
    };
    solver_state.gauss_point.initialize(
        mesh.num_elems, mesh.num_dims, gauss_pt_states, element_communication_plan);
    solver_state.gauss_point.fields.set_values(0.0);
    solver_state.gauss_point.fields_vec.set_values(0.0);
}

inline void populate_nodal_displacement_field(
    const swage::Mesh& mesh,
    XpbdSolverState& solver_state)
{
    FOR_ALL(node_gid, 0, mesh.num_owned_nodes, {
        double displacement_magnitude = 0.0;
        for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
            const double displacement = solver_state.node_positions(node_gid, spatial_dim)
                - solver_state.reference_node_positions(node_gid, spatial_dim);
            solver_state.node.vector_field(node_gid, spatial_dim) = displacement;
            displacement_magnitude += displacement * displacement;
        }
        solver_state.node.scalar_field(node_gid) = std::sqrt(displacement_magnitude);
    });
    MATAR_FENCE();
}

inline void write_xpbd_state(
    swage::Mesh& mesh,
    XpbdSolverState& solver_state,
    const int rank,
    MPI_Comm comm,
    const int graphics_id,
    const bool write_reference_configuration,
    const int timestep_label = -1)
{
    if (write_reference_configuration) {
        copy_mpic_array(solver_state.reference_node_positions, solver_state.node.coords);
        solver_state.node.vector_field.set_values(0.0);
        solver_state.node.scalar_field.set_values(0.0);
    } else {
        solver_state.node.coords = solver_state.node_positions;
        populate_nodal_displacement_field(mesh, solver_state);
    }

    solver_state.node.coords.update_device();
    solver_state.node.vector_field.update_device();
    solver_state.node.scalar_field.update_device();

    if (rank == 0) {
        if (write_reference_configuration) {
            std::cout << "Writing initial state to VTU (graphics_id = " << graphics_id << ")"
                      << std::endl;
        } else if (timestep_label >= 0) {
            std::cout << "Writing timestep " << timestep_label << " to VTU (graphics_id = "
                      << graphics_id << ")" << std::endl;
        } else {
            std::cout << "Writing state to VTU (graphics_id = " << graphics_id << ")"
                      << std::endl;
        }
    }

    write_vtu(mesh, solver_state.node, solver_state.gauss_point, rank, comm, graphics_id);
}

inline int compute_output_interval(const XpbdSolverConfig& config)
{
    if (config.num_outputs > 0) {
        return std::max(1, (config.num_timesteps + config.num_outputs - 1) / config.num_outputs);
    }
    return std::max(1, config.output_interval);
}

inline void run_xpbd_simulation(
    swage::Mesh& mesh,
    XpbdSolverState& solver_state,
    const int rank,
    MPI_Comm comm)
{
    const XpbdSolverConfig& config = solver_state.config;
    const int output_interval = compute_output_interval(config);
    int next_graphics_id = 1;
    bool wrote_at_final_timestep = false;

    if (rank == 0 && config.num_outputs > 0) {
        std::cout << "Writing " << config.num_outputs
                  << " simulation snapshots every " << output_interval
                  << " timesteps" << std::endl;
    }

    for (int timestep_index = 0; timestep_index < config.num_timesteps; ++timestep_index) {
        integrate_boundary_traction(mesh, config, solver_state, timestep_index);
        predict_unconstrained_positions(solver_state, config.timestep_dt);

        copy_mpic_array(solver_state.predicted_position, solver_state.node_positions);
        solver_state.node_positions.communicate();

        pin_dirichlet_node_positions(
            mesh, config, solver_state.reference_node_positions, solver_state.node_positions);

        solver_state.constraint_multipliers.set_values(0.0);

        for (int iteration_index = 0; iteration_index < config.num_solver_iterations; ++iteration_index) {
            for (size_t color_index = 0; color_index < solver_state.num_colors; ++color_index) {
                solve_xpbd_constraints_for_color(mesh, solver_state, color_index);
                solver_state.node_positions.communicate();
                pin_dirichlet_node_positions(
                    mesh, config, solver_state.reference_node_positions, solver_state.node_positions);
            }
        }

        update_nodal_velocities(solver_state, config.timestep_dt);
        copy_mpic_array(solver_state.node_positions, solver_state.node_positions_previous);
        apply_velocity_damping(solver_state.nodal_velocity, config.velocity_damping_factor);

        solver_state.node.coords = solver_state.node_positions;
        copy_mpic_array(solver_state.nodal_velocity, solver_state.node.vector_field);

        const bool is_output_timestep =
            config.num_outputs > 0 && ((timestep_index + 1) % output_interval == 0);
        if (is_output_timestep) {
            write_xpbd_state(
                mesh, solver_state, rank, comm, next_graphics_id, false, timestep_index + 1);
            ++next_graphics_id;
            if (timestep_index + 1 == config.num_timesteps) {
                wrote_at_final_timestep = true;
            }
        }
    }

    if (config.num_outputs > 0 && !wrote_at_final_timestep) {
        write_xpbd_state(
            mesh,
            solver_state,
            rank,
            comm,
            next_graphics_id,
            false,
            config.num_timesteps);
    }
}

inline void print_inverse_lumped_mass_stats(
    XpbdSolverState& solver_state,
    const swage::Mesh& mesh,
    const int rank)
{
    if (rank != 0) {
        return;
    }

    double min_inverse_mass = std::numeric_limits<double>::max();
    double max_inverse_mass = 0.0;
    solver_state.inverse_lumped_mass.update_host();
    for (size_t node_gid = 0; node_gid < mesh.num_owned_nodes; ++node_gid) {
        const double inverse_mass = solver_state.inverse_lumped_mass(node_gid);
        if (inverse_mass > 0.0) {
            min_inverse_mass = std::min(min_inverse_mass, inverse_mass);
            max_inverse_mass = std::max(max_inverse_mass, inverse_mass);
        }
    }
    std::cout << "Inverse lumped mass range: [" << min_inverse_mass << ", " << max_inverse_mass << "]"
              << std::endl;
    std::cout << "Graph colors: " << solver_state.num_colors
              << ", quadrature points per element: " << solver_state.num_quadrature_points
              << std::endl;
}

#endif
