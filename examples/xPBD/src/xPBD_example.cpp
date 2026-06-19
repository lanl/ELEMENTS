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
        int num_elems_dim[3] = {4, 2, 2};

        int num_outputs = 10;

        double density = 8960.0; // copper density kg/m^3

        int iterations = 2;        // nit: paper uses 1 Gauss-Seidel pass per substep
        int nsubsteps = 50;       // small steps per macro timestep (Saillant et al. Algorithm 1)
        double final_time = 10.0;   // run longer to approach static equilibrium under gravity
        double time_step = 0.1;    // macro timestep between VTK outputs

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
                            density
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

        // Body force: f_ext = m * g (gravity in -z for cantilever bending in x-z plane)
        const double gravity = 9.81;
        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            if (fabs(node.inv_mass(node_gid)) > 1E-14) {
                node.force_external(node_gid, 2) = -gravity / node.inv_mass(node_gid);
            }
        });
        node.force_external.update_host();


        // Write out the current state of the mesh to a VTU file. 
        std::cout << "Writing out the current state of the mesh to a VTU file:" << std::endl;
        write_vtu(mesh, node, gauss_point, element, ref_elem.num_gauss_in_elem, rank, MPI_COMM_WORLD, write_id);
        write_id++;


        double t = 0.0;
        const double sub_dt = time_step / static_cast<double>(nsubsteps);

        while (t < final_time) {

            for (int sub = 0; sub < nsubsteps; sub++) {

                // 1. Explicit Euler prediction (Algorithm 1, line 3)
                FOR_ALL(node_gid, 0, mesh.num_nodes, {
                    for(int dim = 0; dim < num_dims; dim++){
                        node.velocity(node_gid, dim) += sub_dt * node.force_external(node_gid, dim) * node.inv_mass(node_gid);
                        node.coords_predicted(node_gid, dim) = node.coords(node_gid, dim) + sub_dt * node.velocity(node_gid, dim);
                    }
                });

                // Reset Lagrange multipliers for each substep
                element.lambda_deviatoric.set_values(0.0);
                element.lambda_volumetric.set_values(0.0);
                MATAR_FENCE();

                // 2. Constraint solve (scaled compliance uses substep dt)
                const double alpha_deviatoric = 1.0 / shear_modulus;
                const double alpha_volumetric = 1.0 / lame_lambda;
                const double alpha_tilde_dev = alpha_deviatoric / (sub_dt * sub_dt);
                const double alpha_tilde_vol = alpha_volumetric / (sub_dt * sub_dt);

                for (int iter = 0; iter < iterations; iter++) {
                
                    // ====================================================================
                    // CONSTRAINT 1: DEVIATORIC 
                    // ====================================================================
                    FOR_ALL(elem_gid, 0, mesh.num_elems, {
                        
                        double C_dev = 0.0;
                        
                        // NOTE: the original used a variable-length array (VLA) here:
                        //   double grad_C_[mesh.num_nodes_in_elem*num_dims];
                        // VLAs are a GCC/Clang extension, not standard C++, and do NOT compile for
                        // CUDA/HIP device code. Use a fixed-size stack buffer sized for the element
                        // order instead. kMaxDofsPerElem must be >= num_nodes_in_elem*num_dims.
                        // 216 nodes covers tensor elements up to P5 hex ((5+1)^3); raise if you go higher.
                        constexpr int kMaxDofsPerElem = 216 * 3;
                        double grad_C_[kMaxDofsPerElem];
                        ViewCArrayDevice<double> grad_C(&grad_C_[0], mesh.num_nodes_in_elem, num_dims);
                        
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            for(int dim = 0; dim < num_dims; dim++){
                                grad_C(node_lid, dim) = 0.0;
                            }
                        }

                        for(int g_lid = 0; g_lid < mesh.num_gauss_in_elem; g_lid++){
                            int gauss_gid = mesh.gauss_in_elem(elem_gid, g_lid);
                            
                            // Volume uses absolute determinant to remain stable under inversion
                            double V_i = fabs(gauss_point.jacobian_determinant(gauss_gid)) * ref_elem.gauss_point_weights(g_lid);

                            xpbd::compute_jacobian(mesh, jacobian_matrix, node.coords_predicted, ref_elem, num_dims, mesh.num_nodes_in_elem, elem_gid, gauss_gid, g_lid);

                            double F_[9] = {0.0};
                            ViewCArrayDevice<double> F(&F_[0], 3, 3);

                            xpbd::compute_F(mesh, node.coords_predicted, gauss_point.J_inv_initial, ref_elem, jacobian_matrix, F, num_dims, elem_gid, gauss_gid, g_lid);
                            
                            double cofactor_F_[9] = {0.0};
                            ViewCArrayDevice<double> cofactor_F(&cofactor_F_[0], 3, 3);
                            double I2 = 0.0;
                            double I3 = 0.0;

                            xpbd::compute_cofactor_and_invariants(F, cofactor_F, num_dims, I2, I3);
                            
                            // Psi_A = I2 - 2 I3 - 1
                            double psi_A = I2 - 2.0 * I3 - 1.0;
                            if (psi_A < 0.0) psi_A = 0.0;
                            C_dev += psi_A * V_i;
                            
                            for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                                double dN_dX[3];
                                xpbd::compute_dN_dX(gauss_point.J_inv_initial, ref_elem, num_dims, gauss_gid, g_lid, node_lid, dN_dX);
                                
                                for(int i = 0; i < num_dims; i++){
                                    double P_A_i = 0.0;
                                    for(int j = 0; j < num_dims; j++){
                                        double P_A_ij = 2.0 * (F(i, j) - cofactor_F(i, j));
                                        P_A_i += P_A_ij * dN_dX[j];
                                    }
                                    grad_C(node_lid, i) += P_A_i * V_i;
                                }
                            }
                        }
                        
                        if (C_dev > 1e-14) {
                            C_dev = std::sqrt(C_dev);
                            double inv_2C = 1.0 / (2.0 * C_dev);
                            for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                                for(int i = 0; i < num_dims; i++){
                                    grad_C(node_lid, i) *= inv_2C;
                                }
                            }
                        } else {
                            C_dev = 0.0;
                            for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                                for(int i = 0; i < num_dims; i++){
                                    grad_C(node_lid, i) = 0.0;
                                }
                            }
                        }

                        double denom_dev = 0.0;
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                            double w_i = node.inv_mass(node_gid);
                            
                            double grad_sq = 0.0;
                            for(int dim = 0; dim < num_dims; dim++){
                                grad_sq += grad_C(node_lid, dim) * grad_C(node_lid, dim);
                            }
                            denom_dev += w_i * grad_sq;
                        }

                        double current_lambda = element.lambda_deviatoric(elem_gid);
                        double delta_lambda = (-C_dev - alpha_tilde_dev * current_lambda) / (denom_dev + alpha_tilde_dev);
                        element.lambda_deviatoric(elem_gid) += delta_lambda;

                        // ---- position scatter ----
                        // RACE CONDITION (parallel only): elements that share a node both read
                        // coords_predicted (to build F / grad_C above) and write it here. The atomic_add
                        // makes each write safe, but neighbouring elements still race on the
                        // read-modify-write as a whole, so this is a Jacobi/Gauss-Seidel hybrid whose
                        // result depends on thread scheduling. In SERIAL the FOR_ALL runs sequentially,
                        // giving a deterministic Gauss-Seidel sweep (each element sees the previous
                        // element's correction) -- which is correct. For a parallel run this loop needs
                        // graph colouring of the elements, or a Jacobi scheme that accumulates
                        // corrections into a separate buffer and applies them after a fence.
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                            double w_i = node.inv_mass(node_gid);
                            
                            if (w_i > 0.0) { 
                                for(int dim = 0; dim < num_dims; dim++){
                                    double delta_x = 0.99 * w_i * grad_C(node_lid, dim) * delta_lambda;
                                    Kokkos::atomic_add(&node.coords_predicted(node_gid, dim), delta_x);
                                }
                            }
                        }
                    });

                    // Deviatoric pass writes coords_predicted; the volumetric pass below reads it.
                    // This fence is required for correctness (device dependency boundary).
                    MATAR_FENCE();

                    // ====================================================================
                    // CONSTRAINT 2: VOLUMETRIC 
                    // ====================================================================
                    FOR_ALL(elem_gid, 0, mesh.num_elems, {
                        double C_vol = 0.0;
                        // See note in the deviatoric block: avoid VLA, use a fixed-size buffer.
                        constexpr int kMaxDofsPerElem = 216 * 3;
                        double grad_C_vol_[kMaxDofsPerElem];
                        ViewCArrayDevice<double> grad_C_vol(&grad_C_vol_[0], mesh.num_nodes_in_elem, num_dims);
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            for(int dim = 0; dim < num_dims; dim++){
                                grad_C_vol(node_lid, dim) = 0.0;
                            }
                        }

                        for(int g_lid = 0; g_lid < mesh.num_gauss_in_elem; g_lid++){
                            
                            int gauss_gid = mesh.gauss_in_elem(elem_gid, g_lid);
                            
                            double V_i = fabs(gauss_point.jacobian_determinant(gauss_gid)) * ref_elem.gauss_point_weights(g_lid);

                            xpbd::compute_jacobian(mesh, jacobian_matrix, node.coords_predicted, ref_elem, num_dims, mesh.num_nodes_in_elem, elem_gid, gauss_gid, g_lid);

                            double F_[9] = {0.0};
                            ViewCArrayDevice<double> F(&F_[0], 3, 3);

                            xpbd::compute_F(mesh, node.coords_predicted, gauss_point.J_inv_initial, ref_elem, jacobian_matrix, F, num_dims, elem_gid, gauss_gid, g_lid);
                            
                            double cofactor_F_[9] = {0.0};
                            ViewCArrayDevice<double> cofactor_F(&cofactor_F_[0], 3, 3);
                            double I2 = 0.0;
                            double I3 = 0.0;

                            xpbd::compute_cofactor_and_invariants(F, cofactor_F, num_dims, I2, I3);
                            
                            // Psi_B = (I3 - 1)^2
                            double psi_B = (I3 - 1.0) * (I3 - 1.0);
                            C_vol += psi_B * V_i;
                            
                            for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                                double dN_dX[3];
                                xpbd::compute_dN_dX(gauss_point.J_inv_initial, ref_elem, num_dims, gauss_gid, g_lid, node_lid, dN_dX);
                                
                                for(int i = 0; i < num_dims; i++){
                                    double P_B_i = 0.0;
                                    for(int j = 0; j < num_dims; j++){
                                        double P_B_ij = 2.0 * (I3 - 1.0) * cofactor_F(i, j);
                                        P_B_i += P_B_ij * dN_dX[j];
                                    }
                                    grad_C_vol(node_lid, i) += P_B_i * V_i;
                                }
                            }
                        }
                        
                        if (C_vol > 1e-12) {
                            C_vol = std::sqrt(C_vol);
                            double inv_2C = 1.0 / (2.0 * C_vol);
                            for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                                for(int i = 0; i < num_dims; i++){
                                    grad_C_vol(node_lid, i) *= inv_2C;
                                }
                            }
                        } else {
                            C_vol = 0.0;
                            for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                                for(int i = 0; i < num_dims; i++){
                                    grad_C_vol(node_lid, i) = 0.0;
                                }
                            }
                        }

                        double denom_vol = 0.0;
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                            double w_i = node.inv_mass(node_gid);
                            
                            double grad_sq = 0.0;
                            for(int dim = 0; dim < num_dims; dim++){
                                grad_sq += grad_C_vol(node_lid, dim) * grad_C_vol(node_lid, dim);
                            }
                            denom_vol += w_i * grad_sq;
                        }

                        double current_lambda = element.lambda_volumetric(elem_gid);
                        double delta_lambda = (-C_vol - alpha_tilde_vol * current_lambda) / (denom_vol + alpha_tilde_vol);
                        element.lambda_volumetric(elem_gid) += delta_lambda;

                        // ---- position scatter ----
                        // RACE CONDITION (parallel only): same hazard as the deviatoric scatter above
                        // (shared nodes read+written across elements). Correct in serial; needs element
                        // colouring or a Jacobi accumulate-then-apply scheme to be correct in parallel.
                        for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                            double w_i = node.inv_mass(node_gid);
                            
                            if (w_i > 0.0) { 
                                for(int dim = 0; dim < num_dims; dim++){
                                    double delta_x = 0.99 * w_i * grad_C_vol(node_lid, dim) * delta_lambda;
                                    Kokkos::atomic_add(&node.coords_predicted(node_gid, dim), delta_x);
                                }
                            }
                        }
                    }); // end 
        
                // If running across MPI ranks, you must synchronize coords_predicted at partition boundaries here
                } // end loop over iterations

                // 3. Kinematic update per substep (Algorithm 1, line 14)
                FOR_ALL(node_gid, 0, mesh.num_nodes, {
                    for(int dim = 0; dim < num_dims; dim++){
                        node.velocity(node_gid, dim) = (node.coords_predicted(node_gid, dim) - node.coords(node_gid, dim)) / sub_dt;
                        node.coords(node_gid, dim) = node.coords_predicted(node_gid, dim);
                    }
                });

            } // end loop over substeps

            // node.velocity.update_host();
            // node.coords.update_host();

            // if (rank == 0) {
            //     node.coords.update_host();
            //     initial_node_coords.update_host();
            //     double tip_z_disp = 0.0;
            //     double max_tip_z = 0.0;
            //     int tip_count = 0;
            //     for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
            //         if (node.coords.host(node_gid, 0) > length[0] - 1e-6) {
            //             double dz = node.coords.host(node_gid, 2) - initial_node_coords.host(node_gid, 2);
            //             tip_z_disp += dz;
            //             max_tip_z = std::max(max_tip_z, std::fabs(dz));
            //             tip_count++;
            //         }
            //     }
            //     const double tip_avg = tip_count > 0 ? tip_z_disp / tip_count : 0.0;
            //     std::cout << "t=" << t << " tip avg z-disp=" << tip_avg
            //               << " max |z-disp|=" << max_tip_z << std::endl;
            // }

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