#include <iostream>
#include <chrono>         // To access timing calipers 
#include "matar.h"

int main() {

    int size_i = 5, size_j = 4, size_k = 3;


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
    
    /*------------------------- Stream Benchmarks ----------------------------------------*/

    // Kokkos stream benchmark test (for FMatrixKokkos and CArrayKokkos objects)
    printf("\n Kokkos stream benchmark CArrayKokkos\n");
    
    double scalar = 3.0;
    size_t nsize = 64 * 64 * 64 * 64;

    auto arr1 = CArrayKokkos <real_t> (nsize);
    auto arr2 = CArrayKokkos <real_t> (nsize);
    auto arr3 = CArrayKokkos <real_t> (nsize);

    printf("Stream benchmark with %d elements.\n", nsize);
    
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
            arr1(i) = 1.0;
            arr2(i) = 2.0;
            arr3(i) = 0.0;
            });
    Kokkos::fence();
   
    // Perform copy benchmark
    auto begin = std::chrono::high_resolution_clock::now();
     
    Kokkos::parallel_for("Copy", nsize, KOKKOS_LAMBDA(const int i) {
            arr3(i) = arr1(i);
            });
    Kokkos::fence();
    

    std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;

    // Perform scaling benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Scale", nsize, KOKKOS_LAMBDA(const int i) {
            arr2(i) = scalar * arr3(i);
            });
    Kokkos::fence();
    
    std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;

    // Perform sum benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Sum", nsize, KOKKOS_LAMBDA(const int i) {
            arr3(i) = arr1(i) + arr2(i);
            });
    Kokkos::fence();
    
    std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;

    // Perform triad benchmark
    begin = std::chrono::high_resolution_clock::now();
    
    Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
            arr1(i) = arr2(i) + 
                                  (scalar * arr3(i));
            });
    Kokkos::fence();
    
    std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;


    // Perform matrix matrix multiply benchmark

    int matrix_size = 64 * 64;
    auto mat1 = CArrayKokkos <real_t> (matrix_size, matrix_size);
    auto mat2 = CArrayKokkos <real_t> (matrix_size, matrix_size);
    auto mat3 = CArrayKokkos <real_t> (matrix_size, matrix_size);

    policy2D mmm_type = policy2D({0,0}, {matrix_size,matrix_size});
    Kokkos::parallel_for("MatrixInit", mmm_type, KOKKOS_LAMBDA(const int i, const int j) {
            mat1(i, j) = (real_t) (i + 1) * (j + 1);
            mat2(i, j) = (real_t) (i + 2) * (j + 2);
            //printf("Mat1 (%d, %d) %lf\n", i, j, mat1(i, j));
            //printf("Mat2 (%d, %d) %lf\n", i, j, mat2(i, j));
        });
    Kokkos::fence();

    begin = std::chrono::high_resolution_clock::now();

    Kokkos::parallel_for ("RaggedDownTeam", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int i = teamMember.league_rank();

            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, matrix_size), [=] (const int j) {
                double temp_var = 0.0;

                Kokkos::parallel_reduce (Kokkos::ThreadVectorRange (teamMember, matrix_size), [=] (const int k, double &mat_val) {
                    mat_val += mat1(i, k) * mat2(k, j);
                }, temp_var);

                mat3(i, j) = temp_var;
                //printf("Mat3 (%d, %d) %lf\n", i, j, mat3(i, j));
            });
        });
    
    Kokkos::fence();
    
    std::chrono::duration<double> mmm_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "MMM time: " << mmm_time.count() << " s." << std::endl;



    
    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
