#include <cmath>
#include <cstdlib>
#include <iostream>

#include "ELEMENTS.h"
#include "mesh_io.h"
#include "state.h"
// #include "xpbd_solver.h"
#include "xpbd_kernels.h"









int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MATAR_INITIALIZE(argc, argv);
    {
        bool verbose = false;
        int print_elem_gid = 0;

        int world_size = 0;
        int rank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const int num_dims = 3;
        int poly_order = 2;
        
        double origin[3] = {0.0, 0.0, 0.0};
        double length[3] = {8.0, 1.0, 1.0};
        int num_elems_dim[3] = {16, 2, 2};

        int num_outputs = 10;

        double density = 1.0;

        int iterations = 10;   
        double final_time = 1.0;
        double time_step = 0.1;

        int write_id = 0;

        double youngs_modulus = 110.0e9; // copper younds modulus in Pa
        double poissons_ratio = 0.34; // copper poissons ratio
        
        
        
        double shear_modulus = youngs_modulus / (2.0 * (1.0 + poissons_ratio)); // copper shear modulus in Pa
        double lame_lambda = (youngs_modulus * poissons_ratio) / ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio)); // copper lame lambda in Pa
  

        

        const double main_start_time = MPI_Wtime();

        swage::Mesh initial_mesh;
        MPICArrayKokkos<double> initial_node_coords;

        swage::Mesh mesh;
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
                mesh,
                initial_node_coords,
                final_node_coords,
                element_communication_plan,
                node_communication_plan,
                world_size,
                rank);
        } else {
            mesh = initial_mesh;
            mesh.num_owned_elems = initial_mesh.num_elems;
            mesh.num_owned_nodes = initial_mesh.num_nodes;
            final_node_coords = initial_node_coords;
            mesh.num_dims = num_dims;
            mesh.Pn = poly_order;
        }

        if (world_size != 1) {
            element_communication_plan.verify_graph_communicator();
            node_communication_plan.verify_graph_communicator();
        }

        if (rank == 0) {
            std::cout << "Final mesh: " << mesh.num_elems << " elements, "
                      << mesh.num_nodes << " nodes" << std::endl;
        }

        // Build the reference element
        elements::fe_ref_elem_t ref_elem;
        ref_elem.init(mesh.Pn, mesh.num_dims);

        std::cout << "Reference element: " << ref_elem.num_gauss_in_elem << " quadrature points" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_dofs_in_elem << " DOFs" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_dg_dofs_in_elem << " DG DOFs" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_zones_in_elem << " zones" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_basis << " basis functions" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_dg_basis << " DG basis functions" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_lobotto_in_elem << " Lobatto points" << std::endl;
        std::cout << "Reference element: " << ref_elem.num_dual_lobotto_in_elem << " Dual Lobatto points" << std::endl;

        // Create the state required for the xPBD solver. 

        // Node State
        std::vector<node_state> node_states = 
            {   node_state::coords, 
                node_state::coords_reference, 
                node_state::coords_predicted, 
                node_state::velocity, 
                node_state::inv_mass, 
                node_state::force_external};
        Node_t node;
        node.initialize(mesh.num_nodes, num_dims, node_states);


        FOR_ALL(i, 0, mesh.num_nodes, {

            for(int dim = 0; dim < num_dims; dim++){
                node.coords(i, dim) = final_node_coords(i, dim);
                node.coords_reference(i, dim) = initial_node_coords(i, dim);
                node.coords_predicted(i, dim) = final_node_coords(i, dim);
                node.velocity(i, dim) = 0.0;
                node.inv_mass(i) = 1.0;
                node.force_external(i, dim) = 0.0;
            }
        });


        // Gauss Point State
        std::vector<gauss_pt_state> gauss_pt_states = 
            {   gauss_pt_state::J_inv_initial,
                gauss_pt_state::jacobian_determinant,
                gauss_pt_state::volume};
        GaussPoint_t gauss_point;
        size_t num_gauss_points = mesh.num_gauss_in_elem * mesh.num_elems;
        gauss_point.initialize(num_gauss_points, num_dims, gauss_pt_states, element_communication_plan);

        // Element State
        std::vector<element_state> element_states = 
            {   element_state::lambda_deviatoric,
                element_state::lambda_volumetric};
        Element_t element;
        element.initialize(mesh.num_elems, element_states, element_communication_plan);



        // Compute the Jacobian matrix at each gauss point
        DCArrayKokkos<double> jacobian_matrix(num_gauss_points, mesh.num_dims, mesh.num_dims);
        jacobian_matrix.set_values(0.0);
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
           
            for(int g_lid = 0; g_lid < mesh.num_gauss_in_elem; g_lid++){

                int gauss_gid = mesh.gauss_in_elem(elem_gid, g_lid);

                xpbd::compute_jacobian(mesh, jacobian_matrix, node.coords, ref_elem, num_dims, mesh.num_nodes_in_elem, elem_gid, gauss_gid, g_lid);

                double det_jacobian = 0.0;
                double J_inv_[9] = {0.0};
                ViewCArrayDevice<double> J_inv(&J_inv_[0], 3, 3);

                xpbd::compute_determinant_and_inverse(jacobian_matrix, J_inv, det_jacobian, num_dims, gauss_gid);

                // Save the initial J_inverse
                for(int dim_i = 0; dim_i < mesh.num_dims; dim_i++){
                    for(int dim_j = 0; dim_j < mesh.num_dims; dim_j++){
                        gauss_point.J_inv_initial(gauss_gid, dim_i, dim_j) = J_inv(dim_i, dim_j);
                    }
                }
                // Store the determinant in gauss_point state
                gauss_point.jacobian_determinant(gauss_gid) = det_jacobian;

            } // end loop over gauss points
        }); // end loop over elements
        jacobian_matrix.update_host();


        if (verbose) {  
            // Print out the jacobian matrix for all gauss points in element 0
            
            for(int g_lid = 0; g_lid < mesh.num_gauss_in_elem; g_lid++){
                int gauss_gid = mesh.gauss_in_elem(print_elem_gid, g_lid);
                
                std::cout << "Gauss point det_j = " <<  gauss_point.jacobian_determinant(gauss_gid) <<std::endl;
                
                std::cout << "Gauss point " << gauss_gid << " Jacobian matrix: " << std::endl;
                for(int dim_i = 0; dim_i < mesh.num_dims; dim_i++){
                    for(int dim_j = 0; dim_j < mesh.num_dims; dim_j++){
                        double value = jacobian_matrix(gauss_gid, dim_i, dim_j);
                        if (std::fabs(value) < 1E-14)
                            std::cout << 0.0 << " ";
                        else
                            std::cout << value << " ";
                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;
            }
        }


        // Compute the mass matrix
        DCArrayKokkos<double> mass_matrix(mesh.num_elems, mesh.num_nodes_in_elem, mesh.num_nodes_in_elem);
        mass_matrix.set_values(0.0);
        
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            
            for(int basis_m = 0; basis_m < ref_elem.num_basis; basis_m++){
                for(int basis_n = 0; basis_n < ref_elem.num_basis; basis_n++){
                    
                    // Integrate by looping over the quadrature points
                    for(int g_lid = 0; g_lid < ref_elem.num_gauss_in_elem; g_lid++){

                        int gauss_gid = mesh.gauss_in_elem(elem_gid, g_lid);
                        double det_j = gauss_point.jacobian_determinant(gauss_gid);
                        

                        mass_matrix(elem_gid, basis_m, basis_n) += 
                            1.0 //density
                            * ref_elem.gauss_point_basis(g_lid, basis_m)
                            * ref_elem.gauss_point_basis(g_lid, basis_n)
                            * det_j
                            * ref_elem.gauss_point_weights(g_lid);
                    } // end loop over gauss points
                } // end loop over basis functions
            } // end loop over basis function
        }); // end loop over elements

        // Print out the mass matrix for element zero

        // Print out the mass matrix for element zero
        if (verbose) {
            std::cout << "Mass matrix for element 0:" << std::endl;
            for(int m = 0; m < mesh.num_nodes_in_elem; m++) {
                for(int n = 0; n < mesh.num_nodes_in_elem; n++) {
                    double value = mass_matrix(print_elem_gid, m, n);
                    // print in a readable way, optionally thresholding tiny values
                    if (std::fabs(value) < 1E-14)
                        std::cout << 0.0 << " ";
                    else
                        std::cout << value << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // Compute the row sums of the mass matrix for element 0 and print them
        if (verbose) {
            std::cout << "Row sum (lumped mass) for each node in element 0:" << std::endl;
            for(int m = 0; m < mesh.num_nodes_in_elem; m++) {
                double row_sum = 0.0;
                for(int n = 0; n < mesh.num_nodes_in_elem; n++) {
                    row_sum += mass_matrix(print_elem_gid, m, n);
                }
                std::cout << "Node " << m << ": " << row_sum;
                if (row_sum < 0.0) {
                    std::cout << "  <--- WARNING: Negative lumped mass!";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        // Compute the mass at each corner
        std::cout << "Computing the corners masses via row tally:" << std::endl;
        DCArrayKokkos<double> corner_mass(mesh.num_corners);

        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for(int corner_lid = 0; corner_lid < mesh.num_nodes_in_elem; corner_lid++){
                int corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
                
                double mass_tally = 0.0;
                for(int n = 0; n < mesh.num_nodes_in_elem; n++) {
                    mass_tally += mass_matrix(elem_gid, corner_lid, n);
                }
                
                corner_mass(corner_gid) = mass_tally;
            } // end loop over corners
        }); // end loop over elements

        // Tally the corner masses to the node
        std::cout << "Tallying the corner masses to the nodes:" << std::endl;
        DCArrayKokkos<double> node_mass(mesh.num_nodes);
        node_mass.set_values(0.0);
        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            for(int corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++){
                int corner_gid = mesh.corners_in_node(node_gid, corner_lid);
                node_mass(node_gid) += corner_mass(corner_gid);
            } // end loop over corners
        });

        // Compute the inverse mass at the nodes
        std::cout << "Computing the inverse mass at the nodes:" << std::endl;
        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            node.inv_mass(node_gid) = 1.0 / node_mass(node_gid);
        });
        node.inv_mass.update_host();

        // verbose = true;
        if (verbose) {
            std::cout << "Inverse mass for all nodes in element 0:" << std::endl;
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                int node_gid = mesh.nodes_in_elem(print_elem_gid, node_lid);
                std::cout << "Node " << node_gid << ": " << node.inv_mass(node_gid) << std::endl;
            }
        }

        // Set the inverse mass to zero for the x=0 nodes
        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            if (std::fabs(node.coords(node_gid, 0)) < 1E-14) {
                node.inv_mass(node_gid) = 0.0;
            }
        });
        node.inv_mass.update_host();
   
        


        // Begin the xPBD solver
        std::cout << "Beginning the xPBD solver:" << std::endl;

        // Compute the body force due to gravity, safely handling zero inverse mass
        FOR_ALL(node_gid, 0, mesh.num_nodes, {

            if (fabs(node.inv_mass(node_gid)) > 1E-14) {

                node.force_external(node_gid, 2) = -1.0*(1.0 / node.inv_mass(node_gid)) * 9.81;

            }
        });
        node.force_external.update_host();


        // Write out the current state of the mesh to a VTU file. 
        std::cout << "Writing out the current state of the mesh to a VTU file:" << std::endl;
        write_vtu(mesh, node, gauss_point, element, ref_elem.num_gauss_in_elem, rank, MPI_COMM_WORLD, write_id);
        write_id++;


        double t = 0.0;
        while (t < final_time) {
            // 1. Explicit Euler Prediction
            FOR_ALL(node_gid, 0, mesh.num_nodes, {
                for(int dim = 0; dim < num_dims; dim++){
                    // v = v + dt * f_ext * m^-1
                    node.velocity(node_gid, dim) += time_step * node.force_external(node_gid, dim) * node.inv_mass(node_gid);
                    
                    // x_pred = x + dt * v
                    node.coords_predicted(node_gid, dim) = node.coords(node_gid, dim) + time_step * node.velocity(node_gid, dim);

                    node.coords(node_gid, dim) = node.coords_predicted(node_gid, dim);
                }
            });

            // Reset Lagrange multipliers for the new time step
            element.lambda.set_values(0.0);


            // // 2. Solver Iterations
            // double alpha = 1e-4; // Material compliance (inverse stiffness)
            // double alpha_tilde = alpha / (time_step * time_step); // Time-step scaled compliance

            // for (int iter = 0; iter < iterations; iter++) {
                
            //     FOR_ALL(elem_gid, 0, mesh.num_elems, {
            //         double C = 0.0;
            //         // You will need a local array to accumulate gradients for each node in the element
            //         double grad_C[27][3] = {0.0}; // Assuming max nodes per element is 27 (poly_order=2) or allocate dynamically

            //         // A. Loop over Gauss points to accumulate C and grad_C
            //         for(int g_lid = 0; g_lid < mesh.num_gauss_in_elem; g_lid++){
            //             int gauss_gid = mesh.gauss_in_elem(elem_gid, g_lid);
                        
            //             // 1. Compute current Deformation Gradient F = J(x_pred) * J(X)^-1
            //             // 2. Evaluate strain energy density (e.g., Neo-Hookean)
            //             // 3. Evaluate first Piola-Kirchhoff stress tensor
            //             // 4. Accumulate into C
            //             // 5. Accumulate into local grad_C array using shape function derivatives
            //         }
                    
            //         // B. Apply square root to finalize C (as per the paper's iso-parametric formulation)
            //         // C = sqrt(C);
            //         // Scale gradients by 1 / (2 * C)

            //         // C. Compute denominator: sum of (inv_mass * |grad_C|^2) for all nodes in the element
            //         double denom = 0.0;
            //         for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
            //             int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            //             double w_i = node.inv_mass(node_gid);
                        
            //             double grad_sq = 0.0;
            //             for(int dim = 0; dim < num_dims; dim++){
            //                 grad_sq += grad_C[node_lid][dim] * grad_C[node_lid][dim];
            //             }
            //             denom += w_i * grad_sq;
            //         }

            //         // D. Compute Delta Lambda
            //         double current_lambda = element.lambda(elem_gid);
            //         double delta_lambda = (-C - alpha_tilde * current_lambda) / (denom + alpha_tilde);
                    
            //         // Update Lagrange multiplier
            //         element.lambda(elem_gid) += delta_lambda;

            //         // E. Scatter position corrections safely using Atomic Adds
            //         for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
            //             int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            //             double w_i = node.inv_mass(node_gid);
                        
            //             if (w_i > 0.0) { // Skip fixed boundaries
            //                 for(int dim = 0; dim < num_dims; dim++){
            //                     double delta_x = w_i * grad_C[node_lid][dim] * delta_lambda;
            //                     Kokkos::atomic_add(&node.coords_predicted(node_gid, dim), delta_x);
            //                 }
            //             }
            //         }
            //     });
        
            //     // If running across MPI ranks, you must synchronize coords_predicted at partition boundaries here
            //     // node_communication_plan.exchange(node.coords_predicted);
            // } // end loop over iterations

            // Write out the current state of the mesh to a VTU file. 
            std::cout << "Writing out the current state of the mesh to a VTU file:" << std::endl;
            write_vtu(mesh, node, gauss_point, element, ref_elem.num_gauss_in_elem, rank, MPI_COMM_WORLD, write_id);
            write_id++;

            t += time_step;

        } // end loop over time steps
        




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
