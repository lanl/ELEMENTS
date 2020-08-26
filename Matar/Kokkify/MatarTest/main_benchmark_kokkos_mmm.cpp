#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <chrono>         // To access timing calipers
#include <cmath>
#include <limits>
#include <cstdlib>
#include "matar.h"

// For LIKWID

#ifdef LIKWID_PERFMON
#include "likwid.h"
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

// If we are compiling with LIKWID support and if we want to conditionally
// enable LIKWID markers, use the following macros
// #define PERF_INIT       if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_INIT;
// #define PERF_START(tag) if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_START(tag);
// #define PERF_STOP(tag)  if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_STOP(tag);
// #define PERF_CLOSE      if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_CLOSE;

int main(int argc, char** argv) {
    // Start LIKWID
    LIKWID_MARKER_INIT;

    const unsigned int default_max_iter = 3;
    const unsigned int repeat = (argc > 1) ? std::atoi(argv[1]) : default_max_iter;
    
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "Kokkos test" << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    // Declare timers
    auto begin = std::chrono::high_resolution_clock::now();
    auto end   = std::chrono::high_resolution_clock::now();

    // Kokkos GPU test
   
    Kokkos::initialize();
    {
    using policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    using policy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

    ////////////////////////////////////////////////////////////////////////////
    // 2D matrix-matrix multiplication (MMM) suite
    ////////////////////////////////////////////////////////////////////////////

    // Perform matrix matrix multiply benchmark
    int matrix_size = 64 * 64;
    int matrix_total_size = (matrix_size * matrix_size);
    auto cak_mat1 = CArrayKokkos <real_t> (matrix_size, matrix_size);
    auto cak_mat2 = CArrayKokkos <real_t> (matrix_size, matrix_size);
    auto cak_mat3 = CArrayKokkos <real_t> (matrix_size, matrix_size);

    RMatrix2D kv_mat1("kv_mat1", matrix_size, matrix_size); 
    RMatrix2D kv_mat2("kv_mat2", matrix_size, matrix_size);
    RMatrix2D kv_mat3("kv_mat3", matrix_size, matrix_size);

    std::streamsize ss = std::cout.precision();

    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "2D matrix-matrix multiplication (MMM) suite" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Running 2D (MMM) benchmark: " 
              << repeat << " times"
              << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << std::endl;

    std::cout << "Number of 2D array elements: " << matrix_total_size 
              << " (" << matrix_size << " elements in each dimension)"
              << std::endl;   
    std::cout << std::setprecision(1) << std::fixed
              << "2D array size: " << (matrix_total_size * sizeof(real_t) * 1.0E-6) << "MB"
              << " (" << (matrix_total_size * sizeof(real_t) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size (2D): " << (3.0 * matrix_total_size * sizeof(real_t) * 1.0E-6) << "MB"
              << " (" << (3.0 * matrix_total_size * sizeof(real_t) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout << std::endl;

    std::cout.precision(ss);

    policy2D mmm_type = policy2D({0,0}, {matrix_size,matrix_size});
    Kokkos::parallel_for("MatrixInit", mmm_type, KOKKOS_LAMBDA(const int i, const int j) {
            // Initialize 2D CArrayKokkos objects (arrays)
            cak_mat1(i, j) = (real_t) (i + 1) * (j + 1);
            cak_mat2(i, j) = (real_t) (i + 2) * (j + 2);

            // Initialize 2D Kokkos Views objects (arrays)
            kv_mat1(i, j) = (real_t) (i + 1) * (j + 1);
            kv_mat2(i, j) = (real_t) (i + 2) * (j + 2);
        });
    Kokkos::fence();

    std::vector<double> mmm_timings;
    
    //LIKWID_MARKER_START("2D_CAK_MMM");
    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for ("MMM", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
                const int i = teamMember.league_rank();

                Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, matrix_size), [=] (const int j) {
                    double temp_var = 0.0;

                    Kokkos::parallel_reduce (Kokkos::ThreadVectorRange (teamMember, matrix_size), [=] (const int k, double &mat_val) {
                        mat_val += (cak_mat1(i, k) * cak_mat2(k, j));
                    }, temp_var);

                    cak_mat3(i, j) = temp_var;
                    //printf("Mat3 (%d, %d) %lf\n", i, j, mat3(i, j));
                });
            });
        Kokkos::fence();
    
        end = std::chrono::high_resolution_clock::now();

        mmm_timings.push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }
    //LIKWID_MARKER_STOP("2D_CAK_MMM");

    // Print tables
   
    // Total size of the three matrices used in the MMM benchmark 
    size_t sizes_MMM = (3 * sizeof(real_t) * matrix_total_size);

    // Calculate minimum and maximum times taken on matrix-matrix
    // multiplication
    // (ignore the first result)
    auto minmax_mmm = std::minmax_element(mmm_timings.begin() + 1,
                                          mmm_timings.end());

    // Calculate average time taken on matrix matrix multiplication
    // (ignore the first result)
    double average_mmm = std::accumulate(mmm_timings.begin() + 1,
                                         mmm_timings.end(),
                                         0.0) / (double) (repeat - 1);
    
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "2D CArrayKokkos MMM benchmark results" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    std::cout << std::endl;

    std::cout << std::left << std::setw(12) << "MMM" 
              << std::left << std::setw(12) << std::setprecision(3) <<
              ((1.0E-6 * sizes_MMM)/ (*minmax_mmm.first))
              << std::left << std::setw(12) << std::setprecision(5) << *minmax_mmm.first
              << std::left << std::setw(12) << std::setprecision(5) << *minmax_mmm.second
              << std::left << std::setw(12) << std::setprecision(5) << average_mmm
              << std::endl
              << std::endl;

    // Kokkos View matrix-matrix multiplication benchmark
    std::vector<double> mmm_kv_timings;

    //LIKWID_MARKER_START("2D_KV_MMM");
    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for ("MMM (KV)", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
                const int i = teamMember.league_rank();

                Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, matrix_size), [=] (const int j) {
                    double temp_var = 0.0;

                    Kokkos::parallel_reduce (Kokkos::ThreadVectorRange (teamMember, matrix_size), [=] (const int k, double &mat_val) {
                        mat_val += (kv_mat1(i, k) * kv_mat2(k, j));
                    }, temp_var);

                    kv_mat3(i, j) = temp_var;
                    //printf("Mat3 (%d, %d) %lf\n", i, j, mat3(i, j));
                });
            });
        Kokkos::fence();
    
        end = std::chrono::high_resolution_clock::now();

        mmm_kv_timings.push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }
    //LIKWID_MARKER_STOP("2D_KV_MMM");

    // Calculate minimum and maximum times taken on matrix-matrix
    // multiplication
    // (ignore the first result)
    auto minmax_mmm_kv = std::minmax_element(mmm_kv_timings.begin() + 1,
                                             mmm_kv_timings.end());

    // Calculate average time taken on matrix matrix multiplication
    // (ignore the first result)
    double average_mmm_kv = std::accumulate(mmm_kv_timings.begin() + 1,
                                            mmm_kv_timings.end(),
                                            0.0) / (double) (repeat - 1);

    // std::cout << "MMM (KV) time: " << mmm_time_kv.count() << " s." << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "2D Kokkos View MMM benchmark results" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    std::cout << std::endl;

    std::cout << std::left << std::setw(12) << "MMM" 
              << std::left << std::setw(12) << std::setprecision(3) <<
              ((1.0E-6 * sizes_MMM)/ (*minmax_mmm_kv.first))
              << std::left << std::setw(12) << std::setprecision(5) << *minmax_mmm_kv.first
              << std::left << std::setw(12) << std::setprecision(5) << *minmax_mmm_kv.second
              << std::left << std::setw(12) << std::setprecision(5) << average_mmm_kv
              << std::endl
              << std::endl;
    }

    Kokkos::finalize();


    printf("--- finished ---\n");
    
    // Stop LIKWID
    LIKWID_MARKER_CLOSE;

    return 0;
}