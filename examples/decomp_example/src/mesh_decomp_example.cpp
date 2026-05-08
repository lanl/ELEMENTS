
//
// Mesh Decomposition Example
//
// This example demonstrates how to:
//   1. Create initial meshes (3D box or 2D polar) on MPI rank 0.
//   2. Partition the mesh across multiple MPI ranks using parallel mesh decomposition (e.g., PT-Scotch).
//   3. Load and distribute node coordinates and mesh connectivity to each rank.
//   4. (Optionally) visualize or post-process partitioned mesh data.
//
// The goal is to show how to use Swage and associated mesh utilities for distributed memory mesh setup,
// which is a typical workflow in scalable finite element or particle-based simulations.
//
// Usage: mpirun -n <num_ranks> ./mesh_decomp_example
//
// Dependencies:
//   - MPI
//   - MATAR (container library)
//   - Swage mesh + decomposition utilities
//   - Optionally: PT-Scotch for mesh partitioning
//


#include <cmath> // for sin


#include "ELEMENTS.h"
#include "state.h"
#include "mesh_io.h"


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MATAR_INITIALIZE(argc, argv);
    { // MATAR scope

    int world_size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_dims = 3;
    int Pn_order = 2;

    double t_main_start = MPI_Wtime();

    // Mesh size for 3D box
    double origin[3] = {0.0, 0.0, 0.0};
    double length[3] = {1.0, 1.0, 1.0};
    int num_elems_dim[3] = {20, 20, 20};
   
    // Mesh size for 2D polar
    double inner_radius = 1.0;
    double outer_radius = 2.0;
    double start_angle = 0.0;
    double end_angle = 45.0;
    int num_elems_i = 100;
    int num_elems_j = 100;

    // Initial mesh built on rank zero
    swage::Mesh initial_mesh;
    MPICArrayKokkos<double> initial_node_coords;

    // Mesh partitioned by pt-scotch, including ghost
    swage::Mesh final_mesh;
    node_t final_node;
    MPICArrayKokkos<double> final_node_coords;

    GaussPoint_t gauss_point;

// ********************************************************  
//              Build the initial mesh
// ********************************************************  

    double t_init_mesh_start = MPI_Wtime();
    if (rank == 0) {
        std::cout<<"World size: "<<world_size<<std::endl;
        std::cout<<"Rank "<<rank<<" Building initial mesh"<<std::endl;

        std::cout<<"Initializing mesh"<<std::endl;
        if(num_dims == 3) {
            build_3d_box(initial_mesh,  initial_node_coords, origin, length, num_elems_dim, Pn_order);
        } else if(num_dims == 2) {
            build_2d_polar(initial_mesh,  initial_node_coords, inner_radius, outer_radius, start_angle, end_angle, num_elems_i, num_elems_j);
        }

        // Read the mesh from a file
        // read_vtk_mesh(initial_mesh, initial_node, 3, "/full/path/to/mesh/file.vtk");


        // Morph the inital mesh to show curvature
        bool morph_mesh = false;
        if(morph_mesh && Pn_order > 1) {

            FOR_ALL(i, 0, initial_mesh.num_nodes, {
                initial_node_coords(i, 0) += 0.0; //0.02 * std::sin(10.0 * initial_node_coords(i, 0));
                initial_node_coords(i, 1) += 0.03 * ( std::sin(10* initial_node_coords(i, 0)) + std::sin(16 * initial_node_coords(i, 0))); 
                initial_node_coords(i, 2) += 0.05 * std::sin(12 * std::sqrt(initial_node_coords(i, 0)*initial_node_coords(i, 0) + initial_node_coords(i, 1)*initial_node_coords(i, 1)));
            });
            initial_node_coords.update_device();
        }

        double t_init_mesh_end = MPI_Wtime();
        std::cout << "Initial mesh build time: " << (t_init_mesh_end - t_init_mesh_start) << " seconds" << std::endl;
        std::cout << "Initial mesh has " << initial_mesh.num_elems << " elements and " << initial_mesh.num_nodes << " nodes" << std::endl;
        std::cout <<" Num_nodes_in_elem: " << initial_mesh.num_nodes_in_elem << std::endl;
        std::cout <<" Num_dims: " << initial_mesh.num_dims << std::endl;
        std::cout <<" Num_elems: " << initial_mesh.num_elems << std::endl;
        std::cout <<" Num_nodes: " << initial_mesh.num_nodes << std::endl;
        std::cout <<" Num_nodes_in_elem: " << initial_mesh.num_nodes_in_elem << std::endl;
        std::cout <<" Num_nodes_in_elem: " << initial_mesh.num_nodes_in_elem << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
// ********************************************************  
//             Partition and balance the mesh
// ********************************************************  
    double t_partition_start = MPI_Wtime();
    // Create communicaion plans
    CommunicationPlan element_communication_plan;
    element_communication_plan.initialize(MPI_COMM_WORLD);
    CommunicationPlan node_communication_plan;
    node_communication_plan.initialize(MPI_COMM_WORLD);

    if(world_size != 1) { // pass through the partitioning function if not a single rank
        elements::partition_mesh(initial_mesh, final_mesh, initial_node_coords, final_node_coords, element_communication_plan, node_communication_plan, world_size, rank);   
    } else {
        final_mesh = initial_mesh;
        final_mesh.num_owned_elems = initial_mesh.num_elems;
        final_mesh.num_owned_nodes = initial_mesh.num_nodes;
        final_node_coords = initial_node_coords;
    }
    double t_partition_end = MPI_Wtime();

    std::cout<<"Final mesh has " << final_mesh.num_elems << " elements and " << final_mesh.num_nodes << " nodes" << std::endl;
    std::cout<<"Final mesh has " << final_mesh.num_owned_elems << " owned elements and " << final_mesh.num_owned_nodes << " owned nodes" << std::endl;

    // Verify communicaiton plans
    if(world_size != 1) element_communication_plan.verify_graph_communicator();
    if(world_size != 1) node_communication_plan.verify_graph_communicator();
   
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0) {
        printf("Mesh partitioning time: %.2f seconds\n", t_partition_end - t_partition_start);
    }


// ****************************************************************************************** 
//     Test element communication using MPI_Neighbor_alltoallv
// ****************************************************************************************** 
    // Gauss points share the same communication plan as elements.
    // This test initializes gauss point fields on owned elements and exchanges them with ghost elements.

    if(world_size != 1) {
        // Test that the shared_tally_owned_nodes mask works correctly by counting all nodes across all ranks and verifying that the number of unique nodes is equal to the number of owned nodes.
        int total_num_nodes = 0;
        int total_local_nodes = 0;
        int total_global_nodes = 0;

    
        FOR_REDUCE_SUM(node_gid, 0, final_mesh.num_owned_nodes, total_local_nodes, {

            if(final_mesh.shared_tally_owned_nodes(node_gid)){
                total_local_nodes++;
            }   

        }, total_num_nodes);
        MATAR_FENCE();
    
        MPI_Allreduce(&total_num_nodes, &total_global_nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(rank == 0){
            std::cout<<"Total number of nodes: "<<total_global_nodes<<std::endl;
            std::cout<<"Error in node count = "<<total_global_nodes - initial_mesh.num_nodes<<std::endl;
        }
    }



    std::cout<<"Generating phony data for gauss points to test MPI communications on rank "<<rank<<std::endl;

    std::vector<gauss_pt_state> gauss_pt_states = {gauss_pt_state::fields, gauss_pt_state::fields_vec};
    gauss_point.initialize(final_mesh.num_elems, final_mesh.num_dims, gauss_pt_states, element_communication_plan); // , &element_communication_plan

    // Initialize the gauss point fields on each rank
    // Set owned elements to rank number, ghost elements to -1 (to verify communication)
    for (int i = 0; i < final_mesh.num_owned_elems; i++) {
        gauss_point.fields.host(i) = static_cast<double>(rank);
        for(int dim = 0; dim < final_mesh.num_dims; dim++){
            gauss_point.fields_vec.host(i, dim) = static_cast<double>(rank);
        }
    }
    for (int i = final_mesh.num_owned_elems; i < final_mesh.num_elems; i++) {
        gauss_point.fields.host(i) = -1.0;  // Ghost elements should be updated
        for(int dim = 0; dim < final_mesh.num_dims; dim++){
            gauss_point.fields_vec.host(i, dim) = -100.0;
        }
    }
    gauss_point.fields.update_device();
    gauss_point.fields_vec.update_device();

    MPI_Barrier(MPI_COMM_WORLD);
    
    gauss_point.fields.communicate();
    gauss_point.fields_vec.communicate();

    MPI_Barrier(MPI_COMM_WORLD);

    CArrayKokkos <double> tmp(final_mesh.num_elems);
    
    // Loop over all elements and average the values of elements connected to that element
    FOR_ALL(i, 0, final_mesh.num_elems, {
        double value = 0.0;
        for (int j = 0; j < final_mesh.num_elems_in_elem(i); j++) {
            value += gauss_point.fields(final_mesh.elems_in_elem(i, j));
        }
        value /= final_mesh.num_elems_in_elem(i);

        tmp(i) = value;
        

        value = 0.0;
        for (int j = 0; j < final_mesh.num_elems_in_elem(i); j++) {
            value += gauss_point.fields_vec(final_mesh.elems_in_elem(i, j), 0);
        }
        value /= final_mesh.num_elems_in_elem(i);

        for(int dim = 0; dim < final_mesh.num_dims; dim++){
            gauss_point.fields_vec(i, dim) = value;
        }
    });
    MATAR_FENCE();

    FOR_ALL(i, 0, final_mesh.num_elems, {
        gauss_point.fields(i) = tmp(i);
    });
    MATAR_FENCE();

    gauss_point.fields.update_host();
    gauss_point.fields_vec.update_host();

    std::cout<<"Generating phony data for nodes to test MPI communications on rank "<<rank<<std::endl;

    // Test node communication using MPI_Neighbor_alltoallv
    std::vector<node_state> node_states = {node_state::coords, node_state::scalar_field, node_state::vector_field};
    final_node.initialize(final_mesh.num_nodes, final_mesh.num_dims, node_states, node_communication_plan);


    // Copy the final node coordinates to the final node
    final_node.coords = final_node_coords;
    final_node.coords.update_device();
    
    for (int i = 0; i < final_mesh.num_owned_nodes; i++) {
        final_node.scalar_field.host(i) = static_cast<double>(rank);
        for(int dim = 0; dim < final_mesh.num_dims; dim++){
            final_node.vector_field.host(i, dim) = static_cast<double>(rank);
        }
    }
    for (int i = final_mesh.num_owned_nodes; i < final_mesh.num_nodes; i++) {
        final_node.scalar_field.host(i) = -100.0;
        for(int dim = 0; dim < final_mesh.num_dims; dim++){
            final_node.vector_field.host(i, dim) = -100;
        }
    }

    final_node.coords.update_device();
    final_node.scalar_field.update_device();
    final_node.vector_field.update_device();
    MATAR_FENCE();
    MPI_Barrier(MPI_COMM_WORLD);

    final_node.scalar_field.communicate();
    final_node.vector_field.communicate();
    
    MATAR_FENCE();
    MPI_Barrier(MPI_COMM_WORLD);

    DCArrayKokkos <double> tmp_too(final_mesh.num_nodes);
    for(int smooth = 0; smooth < 3; smooth++){
        FOR_ALL(i, 0, final_mesh.num_nodes, {

            double value = final_node.scalar_field(i);
            for(int j = 0; j < final_mesh.num_nodes_in_node(i); j++){
                value += final_node.scalar_field(final_mesh.nodes_in_node(i, j));
            }
            value /= final_mesh.num_nodes_in_node(i) + 1;
            tmp_too(i) = value;
        });
        MATAR_FENCE();

        FOR_ALL(i, 0, final_mesh.num_nodes, {
            final_node.scalar_field(i) = tmp_too(i);
            for(int dim = 0; dim < final_mesh.num_dims; dim++){
                final_node.vector_field(i, dim) = tmp_too(i);
            }
        });
        MATAR_FENCE();
    }

    final_node.scalar_field.update_host();

    MATAR_FENCE();


    MPI_Barrier(MPI_COMM_WORLD);
    std::cout<<"Writing VTU file for rank "<<rank<<std::endl;
    write_vtu(final_mesh, final_node, gauss_point, rank, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

    // Stop timer and get execution time
    double t_main_end = MPI_Wtime();
    
    if(rank == 0) {
        printf("Total execution time: %.2f seconds\n", t_main_end - t_main_start);
    }

    } // end MATAR scope
    MATAR_FINALIZE();
    MPI_Finalize();

    return 0;
}