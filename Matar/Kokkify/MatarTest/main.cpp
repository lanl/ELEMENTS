#include <iostream>
#include <chrono>         // To access timing calipers 
#include "pseudo_mesh.hpp"
#include "inherited_inits.hpp"

int main() {

    int size_i = 5, size_j = 4, size_k = 3;
    int big_i = 5000, big_j = 2000, big_k = 3000;


    // Kokkos GPU test
   
    Kokkos::initialize();
    {

    using policy2D = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;
    policy2D array_type = policy2D({0,0}, {size_i, size_j});
    policy2D matrix_type = policy2D({1,1}, {size_i+1, size_j+1});
    policy2D splice_type = policy2D({0,0}, {size_j, size_k});
    using policy3D = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
    policy3D array_type3 = policy3D({0,0,0}, {size_i, size_j, size_k});
    policy3D matrix_type3 = policy3D({1,1,1}, {size_i+1, size_j+1, size_k+1});
    

    printf("~~~~~~~~~~~~~~~~~GPU TEST~~~~~~~~~~~~~~~\n"); 


    int num_parent = 2; // number of materials
    Parent1D parent("parent", num_parent); // Initialize Kokkos View on the GPU of type material, size num_parent
    auto h_parent = HostMirror(parent); // Create a host view of the Kokkos View


    AllocateHost(h_parent, 0, BABY2_SIZE); // Function performed on Host to do raw Kokkos allocation of baby2 GPU space inside of Host data structure
    AllocateHost(h_parent, 1, BABY1_SIZE); // Function performed on Host to do raw Kokkos allocation of baby1 GPU space inside of Host data structure

    Kokkos::deep_copy(parent, h_parent); // deep copy Host data (allocated above) to the GPU Kokkos View. GPU View now has the class space allocated

    InitChildModels(parent, 0, baby2{}); // Kokkos Function to create new instances of the baby2 model on the GPU
    InitChildModels(parent, 1, baby1{1.4,1.0}); // Kokkos Function to create new instances of the baby1 models on the GPU

    // Model test, also shows a Kokkos reduction
    double value_1;
    Kokkos::parallel_reduce(
        "CheckValues",                             
        num_parent,                   
        KOKKOS_LAMBDA(const int idx, real_t &lsum)
        {
             lsum += parent(idx).child->math(2.0, 4.0);
        }
        , value_1);                                

    printf("value %f\n", value_1);



    ClearDeviceModels(parent); // Kokkos Function to call deconstructors of objects on the GPU

    FreeHost(h_parent); // Function performed on Host to free the allocated GPU classes inside of the Host mirror

    


printf("Pseudo Mesh Kokkos\n");
    pseudo_mesh pmesh;
    pmesh.init(size_i, size_j);

    Kokkos::parallel_for("PseudoMesh", 1, KOKKOS_LAMBDA(const int&) {
            //pmesh->init(size_i, size_j);
        });

    //using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    //Kokkos::parallel_for ("PseudoMesh", TeamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
        //pmesh->init(size_i, size_j);
    //    });


    printf("\nCArray\n");
    Kokkos::parallel_for("CArrayValues", array_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.carray(i, j) = (real_t) j + i * size_j;
            printf("(%d, %d) %lf\n", i, j, pmesh.carray(i, j));
        });
    Kokkos::fence();
    printf("\nCMatrix\n");
    Kokkos::parallel_for("CMatrixValues", matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.cmatrix(i, j) = (real_t) (j - 1) + (i - 1) * size_j;
            printf("(%d, %d) %lf\n", i, j, pmesh.cmatrix(i, j));
        });
    Kokkos::fence();
    printf("\nFArray\n");
    Kokkos::parallel_for("FArrayValues", array_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.farray(i, j) = (real_t) i + j * size_i;
            printf("(%d, %d) %lf\n", i, j, pmesh.farray(i, j));
        });
    Kokkos::fence();
    printf("\nFMatrix\n");
    Kokkos::parallel_for("FMatrixValues", matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.fmatrix(i, j) = (real_t) (i - 1) + (j - 1) * size_i;
            printf("(%d, %d) %lf\n", i, j, pmesh.fmatrix(i, j));
        });
    Kokkos::fence();

    printf("\nRaggedRight\n");
    Kokkos::parallel_for ("RaggedRightTeam", TeamPolicy(size_i, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int i = teamMember.league_rank();
            const int i_stride = pmesh.raggedright.stride(i);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, i_stride), [=] (const int j) {
                pmesh.raggedright(i, j) = (double) i + j * i_stride; // not the exact placement, needs start index
                printf("%d %d %lf\n", i, j, pmesh.raggedright(i, j));
            });
        });
    Kokkos::fence();
    
    printf("\nRaggedDown\n");
    Kokkos::parallel_for ("RaggedDownTeam", TeamPolicy(size_j, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int j = teamMember.league_rank();
            const int j_stride = pmesh.raggeddown.stride(j);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, j_stride), [=] (const int i) {
                pmesh.raggeddown(i, j) = (double) j + i * j_stride; // not the exact placement, needs start index
                printf("%d %d %lf\n", i, j, pmesh.raggeddown(i, j));
            });
        });
    Kokkos::fence();

    /*------------------------- Stream Benchmarks ----------------------------------------*/

    // Kokkos stream benchmark test (for FMatrixKokkos and CArrayKokkos objects)
    printf("\n Kokkos stream benchmark (FMatrixKokkos and CArrayKokkos)\n");
    
    double scalar = 3.0;
    size_t nsize = 90000000;
    pseudo_mesh bench;
    bench.init(nsize);

    printf("Stream benchmark with %d elements.\n", nsize);
    
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr1(i) = 1.0;
            bench.arr2(i) = 2.0;
            bench.arr3(i) = 0.0;
            });
    Kokkos::fence();
   
    // Perform copy benchmark
    auto begin = std::chrono::high_resolution_clock::now();
     
    Kokkos::parallel_for("Copy", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr3(i) = bench.arr1(i);
            });
    Kokkos::fence();
    

    std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;

    // Perform scaling benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Scale", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr2(i) = scalar * bench.arr3(i);
            });
    Kokkos::fence();
    
    std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;

    // Perform sum benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Sum", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr3(i) = bench.arr1(i) + bench.arr2(i);
            });
    Kokkos::fence();
    
    std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;

    // Perform triad benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
            bench.arr1(i) = bench.arr2(i) + 
                                  (scalar * bench.arr3(i));
            });
    Kokkos::fence();
    
    std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;
    
    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
