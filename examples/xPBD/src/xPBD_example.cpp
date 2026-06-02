#include <cmath>
#include <cstdlib>
#include <iostream>

#include "ELEMENTS.h"
#include "mesh_io.h"
#include "state.h"
#include "xpbd_solver.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MATAR_INITIALIZE(argc, argv);
    {
        int world_size = 0;
        int rank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const int num_dims = 3;
        int poly_order = 1;
        int num_outputs = 10;
        if (argc > 1) {
            poly_order = std::max(1, std::atoi(argv[1]));
        }
        if (argc > 2) {
            num_outputs = std::max(0, std::atoi(argv[2]));
        }

        const double main_start_time = MPI_Wtime();

        double origin[3] = {0.0, 0.0, 0.0};
        double length[3] = {1.0, 1.0, 1.0};
        int num_elems_dim[3] = {4, 4, 8};

        swage::Mesh initial_mesh;
        MPICArrayKokkos<double> initial_node_coords;

        swage::Mesh final_mesh;
        MPICArrayKokkos<double> final_node_coords;

        CommunicationPlan element_communication_plan;
        element_communication_plan.initialize(MPI_COMM_WORLD);
        CommunicationPlan node_communication_plan;
        node_communication_plan.initialize(MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "World size: " << world_size << std::endl;
            std::cout << "Building initial 3D box mesh with poly_order = " << poly_order << std::endl;

            initial_mesh.num_dims = num_dims;
            initial_mesh.Pn = poly_order;
            build_3d_box(
                initial_mesh,
                initial_node_coords,
                origin,
                length,
                num_elems_dim,
                poly_order);

            std::cout << "Initial mesh: " << initial_mesh.num_elems << " elements, "
                      << initial_mesh.num_nodes << " nodes" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (world_size != 1) {
            elements::partition_mesh(
                initial_mesh,
                final_mesh,
                initial_node_coords,
                final_node_coords,
                element_communication_plan,
                node_communication_plan,
                world_size,
                rank);
        } else {
            final_mesh = initial_mesh;
            final_mesh.num_owned_elems = initial_mesh.num_elems;
            final_mesh.num_owned_nodes = initial_mesh.num_nodes;
            final_node_coords = initial_node_coords;
            final_mesh.num_dims = num_dims;
            final_mesh.Pn = poly_order;
        }

        if (world_size != 1) {
            element_communication_plan.verify_graph_communicator();
            node_communication_plan.verify_graph_communicator();
        }

        if (rank == 0) {
            std::cout << "Final mesh: " << final_mesh.num_elems << " elements, "
                      << final_mesh.num_nodes << " nodes" << std::endl;
        }

        // Mesh initialized, now implement the xPBD solver. 

        XpbdSolverConfig solver_config;
        solver_config.poly_order = poly_order;
        solver_config.num_timesteps = 200;
        solver_config.num_solver_iterations = 10;
        solver_config.num_force_ramp_steps = 10;
        solver_config.num_outputs = num_outputs;
        solver_config.timestep_dt = 5.0e-4;
        solver_config.velocity_damping_factor = 0.05;
        solver_config.dirichlet_plane_origin[0] = origin[0];
        solver_config.dirichlet_plane_origin[1] = origin[1];
        solver_config.dirichlet_plane_origin[2] = origin[2];
        solver_config.dirichlet_plane_normal[0] = 0.0;
        solver_config.dirichlet_plane_normal[1] = 0.0;
        solver_config.dirichlet_plane_normal[2] = 1.0;
        solver_config.surface_traction_plane_normal[0] = 0.0;
        solver_config.surface_traction_plane_normal[1] = 0.0;
        solver_config.surface_traction_plane_normal[2] = 1.0;
        solver_config.gravity_vector[0] = 0.0;
        solver_config.gravity_vector[1] = 1.0;
        solver_config.gravity_vector[2] = -9.81;
        solver_config.material.density = 1000.0;
        solver_config.material.youngs_modulus = 1.0e4;
        solver_config.material.poissons_ratio = 0.3;
        solver_config.material.compliance = 1.0 / solver_config.material.youngs_modulus;

        XpbdSolverState solver_state;
        initialize_xpbd_solver(
            final_mesh,
            final_node_coords,
            node_communication_plan,
            element_communication_plan,
            solver_state,
            solver_config);

        print_inverse_lumped_mass_stats(solver_state, final_mesh, rank);

        write_xpbd_state(
            final_mesh,
            solver_state,
            rank,
            MPI_COMM_WORLD,
            0,
            true);

        const double simulation_start_time = MPI_Wtime();
        run_xpbd_simulation(final_mesh, solver_state, rank, MPI_COMM_WORLD);
        const double simulation_end_time = MPI_Wtime();

        if (rank == 0) {
            std::cout << "XPBD simulation time: "
                      << (simulation_end_time - simulation_start_time) << " seconds"
                      << std::endl;
        }

        double max_displacement = 0.0;
        double local_max_displacement = 0.0;
        solver_state.node_positions.update_host();
        solver_state.reference_node_positions.update_host();
        for (size_t node_gid = 0; node_gid < final_mesh.num_owned_nodes; ++node_gid) {
            for (int spatial_dim = 0; spatial_dim < 3; ++spatial_dim) {
                const double displacement = solver_state.node_positions.host(node_gid, spatial_dim)
                    - solver_state.reference_node_positions.host(node_gid, spatial_dim);
                local_max_displacement = std::max(local_max_displacement, std::fabs(displacement));
            }
        }
        MPI_Allreduce(
            &local_max_displacement,
            &max_displacement,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "Maximum nodal displacement magnitude: " << max_displacement << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        const double main_end_time = MPI_Wtime();
        if (rank == 0) {
            std::cout << "Total execution time: " << (main_end_time - main_start_time)
                      << " seconds" << std::endl;
        }
    }
    MATAR_FINALIZE();
    MPI_Finalize();

    return 0;
}
