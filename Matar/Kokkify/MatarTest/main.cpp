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


    printf("CArray\n");
    Kokkos::parallel_for("CArrayValues", array_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.carray(i, j) = (real_t) j + i * size_j;
            printf("(%d, %d) %lf\n", i, j, pmesh.carray(i, j));
        });
    Kokkos::fence();
    printf("CMatrix\n");
    Kokkos::parallel_for("CMatrixValues", matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.cmatrix(i, j) = (real_t) (j - 1) + (i - 1) * size_j;
            printf("(%d, %d) %lf\n", i, j, pmesh.cmatrix(i, j));
        });
    Kokkos::fence();
    printf("FArray\n");
    Kokkos::parallel_for("FArrayValues", array_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.farray(i, j) = (real_t) i + j * size_i;
            printf("(%d, %d) %lf\n", i, j, pmesh.farray(i, j));
        });
    Kokkos::fence();
    printf("FMatrix\n");
    Kokkos::parallel_for("FMatrixValues", matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
            pmesh.fmatrix(i, j) = (real_t) (i - 1) + (j - 1) * size_i;
            printf("(%d, %d) %lf\n", i, j, pmesh.fmatrix(i, j));
        });
    Kokkos::fence();
   /*
    Kokkos::parallel_for ("RaggedDownTeam", TeamPolicy(size_j, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int j = teamMember.league_rank();
            const int j_stride = pmesh->var.stride(j);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, j_stride), [=] (const int i) {
                pmesh->var(i, j) = (double) j + i * j_stride; // not the exact placement, needs start index
                printf("%d %d %lf\n", i, j, pmesh->var(i, j));
            });
        });
    Kokkos::fence();
    */


    /*------------------------- Stream Benchmarks ----------------------------------------*/

#ifdef PSEUDOVAR
    // Kokkos stream benchmark test (for FMatrixKokkos and CArrayKokkos objects)
    printf("\n Kokkos stream benchmark (FMatrixKokkos and CArrayKokkos)\n");
    
    double scalar = 3.0;
    size_t nsize = 10000000;
    pseudo_var* pvar = (pseudo_var*) Kokkos::kokkos_malloc<Kokkos::CudaSpace> 
                                             (sizeof(pseudo_var));

    // Initialize three 1D FMatrixKokkos and CArrayKokkos objects,
    // respectively
    Kokkos::parallel_for("StreamTriadInit", 1, KOKKOS_LAMBDA(const int&) {
            pvar->init(nsize);
            });

    printf("Stream benchmark with %d elements.\n", nsize);
    /*
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var1(i) = 1.0;
            pvar->c_arr_var2(i) = 2.0;
            pvar->c_arr_var3(i) = 0.0;
            });
    */
    // Perform copy benchmark
    auto begin = std::chrono::high_resolution_clock::now();
    /* 
    Kokkos::parallel_for("Copy", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var3(i) = pvar->c_arr_var1(i);
            });
    */

    std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;

    // Perform scaling benchmark
    begin = std::chrono::high_resolution_clock::now();
    /*
    Kokkos::parallel_for("Scale", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var2(i) = scalar * pvar->c_arr_var3(i);
            });
    */
    std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;

    // Perform sum benchmark
    begin = std::chrono::high_resolution_clock::now();
    /*
    Kokkos::parallel_for("Sum", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var3(i) = pvar->c_arr_var1(i) + pvar->c_arr_var2(i);
            });
    */
    std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;

    // Perform triad benchmark
    begin = std::chrono::high_resolution_clock::now();
    /*
    Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var1(i) = pvar->c_arr_var2(i) + 
                                  (scalar * pvar->c_arr_var3(i));
            });
    */
    std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;
    


    // Test RaggedRight 
    /*
    // First, create RaggedRight using CArray
    auto c_array_strides_array = CArray<size_t>(size_i);
    for (size_t i = 0; i < size_i; i++) {
        c_array_strides_array(i) = i;
    }
    
    auto rr_array = RaggedRightArray<double>(c_array_strides_array);

    int val = 1;
    for (size_t i = 0; i < size_i; i++) {
        for (size_t j = 0; j < rr_array.stride(i); j++) {
            rr_array(i, j) = val;
            val++;
        }
    }
    
    for (size_t i = 0; i < size_i; i++) {
        for (size_t j = 0; j < rr_array.stride(i); j++) {
            printf("rr_array(%d, %d) = %f\n", i, j, rr_array(i, j));
        }
    }
    
    // Then, create RaggedRight using ViewCArray
    auto c_array_strides_array_2 = CArray<size_t>(size_i + 2);
    for (size_t i = 0; i < (size_i + 2); i++) {
        c_array_strides_array_2(i) = 2*i;
    }

    auto view_c_array_strides_array = ViewCArray<size_t>(&c_array_strides_array_2(0),
                                                          size_i + 2);

    auto rr_array_vs = RaggedRightArray<double>(view_c_array_strides_array);

    val = 1;
    for (size_t i = 0; i < (size_i + 2); i++) {
        for (size_t j = 0; j < rr_array_vs.stride(i); j++) {
            rr_array_vs(i, j) = val;
            val++;
        }
    }

    printf("\n");

    for (size_t i = 0; i < (size_i + 2); i++) {
        for (size_t j = 0; j < rr_array_vs.stride(i); j++) {
            printf("rr_array_vs(%d, %d) = %f\n", i, j, rr_array_vs(i, j));
        }
    }

    // Create RaggedRight using regular C array
    size_t num_strides = 4;
    size_t strides[num_strides] = {4, 5, 1, 2}; 
    
    auto rr_array_reg = RaggedRightArray<double>(strides, num_strides);

    val = 1;
    for (size_t i = 0; i < num_strides; i++) {
        for (size_t j = 0; j < rr_array_reg.stride(i); j++) {
            rr_array_reg(i, j) = val;
            val++;
        }
    }

    printf("\n");

    for (size_t i = 0; i < num_strides; i++) {
        for (size_t j = 0; j < rr_array_reg.stride(i); j++) {
            printf("rr_array_reg(%d, %d) = %f\n", i, j, rr_array_reg(i, j));
        }
    }
    */
#endif
    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
