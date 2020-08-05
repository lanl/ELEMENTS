#include <iostream>
#include <chrono>         // To access timing calipers 
#include "pseudo_mesh.hpp"

int main() {

    int size_i = 5, size_j = 2, size_k = 3;
    int big_i = 5000, big_j = 2000, big_k = 3000;

    printf("~~~~~~~~~~~~~~~~~GPU TEST~~~~~~~~~~~~~~~\n"); 

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
    

printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
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
    


    /*------------------------- Stream Benchmarks ----------------------------------------*/

    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
