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

    int num_dims = 2;
    int Pn_order = 1;

    double t_main_start = MPI_Wtime();

    // Mesh size for 3D box
    double origin[3] = {0.0, 0.0, 0.0};
    double length[3] = {1.0, 1.0, 1.0};
    int num_elems_dim[3] = {20, 20, 20};

    // Initial mesh built on rank zero
    swage::Mesh initial_mesh;
    MPICArrayKokkos<double> initial_node_coords;

    // Mesh partitioned by pt-scotch, including ghost
    swage::Mesh final_mesh;
    node_t final_node;
    MPICArrayKokkos<double> final_node_coords;

    // ********************************************************  
    //              Build the initial mesh
    // ********************************************************  

    double t_init_mesh_start = MPI_Wtime();
    if (rank == 0) {
        std::cout<<"World size: "<<world_size<<std::endl;
        std::cout<<"Rank "<<rank<<" Building initial mesh"<<std::endl;

        std::cout<<"Initializing mesh"<<std::endl;
        initial_mesh.num_dims = num_dims;
        initial_mesh.Pn = Pn_order;
        
        if(num_dims == 3) {
            build_3d_box(initial_mesh,  initial_node_coords, origin, length, num_elems_dim, Pn_order);
        } else if(num_dims == 2) {
            build_2d_polar(initial_mesh,  initial_node_coords, inner_radius, outer_radius, start_angle, end_angle, num_elems_i, num_elems_j);
        }

        // Read the mesh from a file
        // read_vtk_mesh(initial_mesh, initial_node, 3, "/full/path/to/mesh/file.vtk");

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

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0) printf("After partitioning\n");

    } else {
        final_mesh = initial_mesh;
        final_mesh.num_owned_elems = initial_mesh.num_elems;
        final_mesh.num_owned_nodes = initial_mesh.num_nodes;
        final_node_coords = initial_node_coords;
        final_mesh.num_dims = num_dims;
        final_mesh.Pn = Pn_order;
    }
    double t_partition_end = MPI_Wtime();

    std::cout<<"Final mesh has " << final_mesh.num_elems << " elements and " << final_mesh.num_nodes << " nodes" << std::endl;
    std::cout<<"Final mesh has " << final_mesh.num_owned_elems << " owned elements and " << final_mesh.num_owned_nodes << " owned nodes" << std::endl;
    std::cout<<"Final mesh num_dims = " << final_mesh.num_dims << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // Verify communicaiton plans
    if(world_size != 1) element_communication_plan.verify_graph_communicator();
    if(world_size != 1) node_communication_plan.verify_graph_communicator();

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0) {
        printf("Mesh partitioning time: %.2f seconds\n", t_partition_end - t_partition_start);
    }


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