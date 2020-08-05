#include <iostream>
#include <chrono>         // To access timing calipers 
#include "pseudo_mesh.hpp"
#include "kokkos_simple.hpp"

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


    int num_mats = 2; // number of materials
    Material1D mats("mats", num_mats); // Initialize Kokkos View on the GPU of type material, size num_mats
    auto h_mats = HostMirror(mats); // Create a host view of the Kokkos View


    AllocateHost(h_mats, 0, SESAME_SIZE); // Function performed on Host to do raw Kokkos allocation of sesame GPU space inside of Host data structure
    AllocateHost(h_mats, 1, GAMMA_LAW_SIZE); // Function performed on Host to do raw Kokkos allocation of gamma_law GPU space inside of Host data structure

    Kokkos::deep_copy(mats, h_mats); // deep copy Host data (allocated above) to the GPU Kokkos View. GPU View now has the class space allocated

    InitEOSModels(mats, 0, sesame{}); // Kokkos Function to create new instances of the sesame model on the GPU
    InitEOSModels(mats, 1, gamma_law{1.4,1.0}); // Kokkos Function to create new instances of the gamma_law models on the GPU

    // Model test, also shows a Kokkos reduction
    double value_1;
    Kokkos::parallel_reduce(
        "CheckValues",                             
        num_mats,                   
        KOKKOS_LAMBDA(const int idx, real_t &lsum)
        {
             lsum += mats(idx).eos->pres_from_den_sie(2.0, 4.0);
        }
        , value_1);                                

    printf("value %f\n", value_1);



    ClearDeviceModels(mats); // Kokkos Function to call deconstructors of objects on the GPU

    FreeHost(h_mats); // Function performed on Host to free the allocated GPU classes inside of the Host mirror

    


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

    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
