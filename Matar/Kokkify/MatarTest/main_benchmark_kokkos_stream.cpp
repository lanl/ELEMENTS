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

    int size_i = 5, size_j = 4, size_k = 3;
    const int size3 = 256;
    const int size1 = size3 * size3 * size3;

    const unsigned int default_max_iter = 3;
    const unsigned int repeat = (argc > 1) ? std::atoi(argv[1]) : default_max_iter;

    ///////////////////////////////////////////////////////////////////////////
    // The following code runs the STREAM benchmark suite on the following
    // data types:
    //      1. "Conventionally" allocated multi-dimensional arrays
    //      2. Kokkos Views
    //      3. MATAR's Kokkos-specific classes
    //         (in particular, the CArrayKokkos classes)
    //
    // In (1)-(3), the STREAM benchmark suite is run on one-dimensional (1D) 
    // and three-dimensional (3D) instances of the types, i.e., the STREAM
    // benchmark suite is run on:
    //      1. 1D and 3D "conventionally" allocated multi-dimensional arrays
    //      2. 1D and 3D Kokkos Views
    //      3. 1D and 3D CArrayKokkos objects
    //
    // In addition, the code also runs the dot product kernel, adapted from the
    // BabelStream project (https://github.com/UoB-HPC/BabelStream)
    ///////////////////////////////////////////////////////////////////////////

    // std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~GPU  TESTS~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "Kokkos test" << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    // Number of entries for 1D types (by default, it is 2e25)
    // size_t ARRAY_SIZE_1D = 33554432;
    size_t nsize = 64 * 64 * 64 * 64;
    size_t nsize_3D = 512;
    size_t ARRAY_SIZE_3D = (nsize_3D * nsize_3D * nsize_3D); 

    // Create vector of vectors that will store timings for kernel runs
    // Rows correspond to kernels:
    //      1. Copy kernel
    //      2. Scale kernel
    //      3. Sum kernel
    //      4. Triad kernel
    //      5. Dot product kernel
    
    int num_kernels = 5;

    // std::vector<std::vector<double>> reg_arr_timings(num_kernels);

    // std::vector<std::vector<double>> reg_arr_timings_3D(num_kernels);

    // Declare timers
    auto begin = std::chrono::high_resolution_clock::now();
    auto end   = std::chrono::high_resolution_clock::now();

    std::string labels[num_kernels] = {"Copy", "Mul", "Add", "Triad", "Dot"};

    size_t sizes[num_kernels] = {
        2 * sizeof(real_t) * nsize,
        2 * sizeof(real_t) * nsize,
        3 * sizeof(real_t) * nsize,
        3 * sizeof(real_t) * nsize,
        2 * sizeof(real_t) * nsize
    };

    size_t sizes_3D[num_kernels] = {
        2 * sizeof(real_t) * ARRAY_SIZE_3D,
        2 * sizeof(real_t) * ARRAY_SIZE_3D,
        3 * sizeof(real_t) * ARRAY_SIZE_3D,
        3 * sizeof(real_t) * ARRAY_SIZE_3D,
        2 * sizeof(real_t) * ARRAY_SIZE_3D
    };

    // Kokkos GPU test
   
    Kokkos::initialize();
    {
    using policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    using policy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
    
    /***************************************************************************
     * 1D array STREAM benchmark suite
     **************************************************************************/

    // There are three vectors involved in the STREAM benchmark: a, b, and c
    // THe following are their respective initial values
    real_t arr1_init_val = 0.1;
    real_t arr2_init_val = 0.2;
    real_t arr3_init_val = 0.0;

    real_t scalar = 0.4;

    // The following variables contain their respective final values after
    // all iterations of the STREAM benchmark
    real_t arr1_fin_val = arr1_init_val;
    real_t arr2_fin_val = arr2_init_val;
    real_t arr3_fin_val = arr3_init_val;

    real_t dot_1D_fin_val = 0.0;

    for (int iter = 0; iter < repeat; iter++) {
        arr3_fin_val = arr1_fin_val;
        arr2_fin_val = scalar * arr3_fin_val;
        arr3_fin_val = arr1_fin_val + arr2_fin_val;
        arr1_fin_val = arr2_fin_val + (scalar * arr3_fin_val);
    }

    dot_1D_fin_val = arr1_fin_val * arr2_fin_val * nsize;

    // Create "regularly" allocated 1D C++ arrays
    //real_t* reg_arr1 = new real_t[nsize];
    //real_t* reg_arr2 = new real_t[nsize];
    //real_t* reg_arr3 = new real_t[nsize];
   
    // Create 1D FArrayKokkos objects
    auto fak_arr1 = FArrayKokkos <real_t> (nsize);
    auto fak_arr2 = FArrayKokkos <real_t> (nsize);
    auto fak_arr3 = FArrayKokkos <real_t> (nsize);

    // Create 1D ViewFArrayKokkos objects
    auto vfak_arr1 = ViewFArrayKokkos <real_t> (fak_arr1.pointer(), nsize);
    auto vfak_arr2 = ViewFArrayKokkos <real_t> (fak_arr2.pointer(), nsize);
    auto vfak_arr3 = ViewFArrayKokkos <real_t> (fak_arr3.pointer(), nsize);

    // Create 1D CArrayKokkos objects
    auto cak_arr1 = CArrayKokkos <real_t> (nsize);
    auto cak_arr2 = CArrayKokkos <real_t> (nsize);
    auto cak_arr3 = CArrayKokkos <real_t> (nsize);

    // Create 1D ViewCArrayKokkos objects
    auto vcak_arr1 = ViewCArrayKokkos <real_t> (cak_arr1.pointer(), nsize);
    auto vcak_arr2 = ViewCArrayKokkos <real_t> (cak_arr2.pointer(), nsize);
    auto vcak_arr3 = ViewCArrayKokkos <real_t> (cak_arr3.pointer(), nsize);

    // Create 1D Kokkos View objects
    RMatrix1D kv_arr1("kv_arr1", nsize);
    RMatrix1D kv_arr2("kv_arr2", nsize);
    RMatrix1D kv_arr3("kv_arr3", nsize);

    // printf("Stream benchmark with %d elements.\n", nsize);

    // Initialize 1D FArrayKokkos, 1D CArrayKokkos, and 1D Kokkos View objects
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
                // Initialize "regularly" allocated 1D C++ arrays
                //reg_arr1[i] = 1.0;
                //reg_arr2[i] = 2.0;
                //reg_arr3[i] = 0.0;

                // Initialize 1D FArrayKokkos objects
                fak_arr1(i) = arr1_init_val;
                fak_arr2(i) = arr2_init_val;
                fak_arr3(i) = arr3_init_val;

                // Initialize 1D CArrayKokkos objects
                cak_arr1(i) = arr1_init_val;
                cak_arr2(i) = arr2_init_val;
                cak_arr3(i) = arr3_init_val;

                // Initialize 1D Kokkos View objects
                kv_arr1(i) = arr1_init_val;
                kv_arr2(i) = arr2_init_val;
                kv_arr3(i) = arr3_init_val;
                });
    Kokkos::fence();

    std::streamsize ss = std::cout.precision();

    std::cout << "-------------------------------------------------------------"
              << std::endl;
    std::cout << "1D array STREAM benchmark suite" << std::endl;
    std::cout << "-------------------------------------------------------------"
              << std::endl
              << std::endl;

    std::cout << "Running kernels " << repeat << " times" << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << "Number of 1D array elements: " << nsize << std::endl;
    std::cout << std::endl;

    // Conversion factor between megabytes (MB) and bytes: 1 byte = 10^(-6) MB
    // Conversion factor between bytes and gigabytes (GB): 1 byte = 10^(-9) GB
    // std::cout << "1D information" << std::endl;
    // std::cout << "--------------" << std::endl;
    std::cout << std::setprecision(1) << std::fixed
              << "1D array size: " << (nsize * sizeof(real_t) * 1.0E-6) << "MB"
              << " (" << (nsize * sizeof(real_t) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size (1D): " << (3.0 * nsize * sizeof(real_t) * 1.0E-6) << "MB"
              << " (" << (3.0 * nsize * sizeof(real_t) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout << std::endl;

    std::cout.precision(ss);

    ////////////////////////////////////////////////////////////////////////////
    // 1D "conventionally" allocated C++ array STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // Free allocated memory for "regularly" allocated C++ arrays
    //delete[] reg_arr1;
    //delete[] reg_arr2;
    //delete[] reg_arr3;

    ////////////////////////////////////////////////////////////////////////////
    // 1D FArrayKokkos STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // Vector that stores the times taken by the various kernel calls on the
    // 1D FArrayKokkos objects
    std::vector<std::vector<double>> fak_1D_timings(num_kernels);

    // Variable the stores the dot product of vectors a and b
    real_t fak_dot_1D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 1D FArrayKokkos copy kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (1D FAK)", nsize, KOKKOS_LAMBDA(const int i) {
                fak_arr3(i) = fak_arr1(i);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record copy kernel timing
        fak_1D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 1D FArrayKokkos scale kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale (1D FAK)", nsize, KOKKOS_LAMBDA(const int i) {
                fak_arr2(i) = (scalar * fak_arr3(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record scale kernel timing
        fak_1D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // PERF_START(tag) // tag is a string, e.g., "copy"
        // 1D FArrayKokkos sum kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Sum (1D FAK)", nsize, KOKKOS_LAMBDA(const int i) {
                fak_arr3(i) = (fak_arr1(i) + fak_arr2(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record sum kernel timing
        fak_1D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count()); 

        // 1D FArrayKokkos triad kernel

        //LIKWID_MARKER_START("1D_FAK_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad (1D FAK)", nsize, KOKKOS_LAMBDA(const int i) {
                fak_arr1(i) = (fak_arr2(i) + (scalar * fak_arr3(i)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_FAK_TRIAD");

        // Record triad kernel timing
        fak_1D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
       
        // 1D FArrayKokkos dot product kernel

        //LIKWID_MARKER_START("1D_FAK_DOT");
        fak_dot_1D_fin_val = 0.0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (1D FAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += (fak_arr1(i) * fak_arr2(i));
        }, fak_dot_1D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_FAK_DOT");

        // Record dot product kernel timing
        fak_1D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 1D FArrayKokkos STREAM benchmark results
    real_t fak_arr1_1D_err = 0;
    real_t fak_arr2_1D_err = 0;
    real_t fak_arr3_1D_err = 0;
    real_t fak_dot_1D_err = std::fabs(dot_1D_fin_val - fak_dot_1D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (1D FAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (fak_arr1(i) - arr1_fin_val) >= 0 ? (fak_arr1(i) - arr1_fin_val) : (arr1_fin_val - fak_arr1(i));
    }, fak_arr1_1D_err);
    Kokkos::fence();

    fak_arr1_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D FAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (fak_arr2(i) - arr2_fin_val) >= 0 ? (fak_arr2(i) - arr2_fin_val) : (arr2_fin_val - fak_arr2(i));
    }, fak_arr2_1D_err);
    Kokkos::fence();

    fak_arr2_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D FAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (fak_arr3(i) - arr3_fin_val) >= 0 ? (fak_arr3(i) - arr3_fin_val) : (arr3_fin_val - fak_arr3(i));
    }, fak_arr3_1D_err);
    Kokkos::fence();

    fak_arr3_1D_err /= nsize;

    real_t epsi = (std::numeric_limits<real_t>::epsilon() * 100.0);

    if (fak_arr1_1D_err > epsi) {
        std::cout << "Validation failed on fak_arr1. Average error "
                  << fak_arr1_1D_err << std::endl << std::endl;
    }

    if (fak_arr2_1D_err > epsi) {
        std::cout << "Validation failed on fak_arr2. Average error "
                  << fak_arr2_1D_err << std::endl << std::endl;
    }

    if (fak_arr3_1D_err > epsi) {
        std::cout << "Validation failed on fak_arr3. Average error "
                  << fak_arr3_1D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (fak_dot_1D_err > 1.0E-8) {
        std::cout << "Validation failed on 1D FAK dot product kernel. Error is "
                  << fak_dot_1D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << fak_dot_1D_fin_val 
                  << " but should be "  << dot_1D_fin_val
                  << std::endl << std::endl;
    }


    // Print kernel computation memory bandwidth table header
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "1D FArrayKokkos STREAM benchmark results" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(fak_1D_timings[ker].begin() + 1,
                                          fak_1D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(fak_1D_timings[ker].begin() + 1,
                                         fak_1D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // 1D ViewFArrayKokkos STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////
    
    // Initialize the 1D ViewFArrayKokkos objects
    Kokkos::parallel_for("Initialize (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i) {
            vfak_arr1(i) = arr1_init_val;
            vfak_arr2(i) = arr2_init_val;
            vfak_arr3(i) = arr3_init_val;
            });
    Kokkos::fence();

    // Vector that stores the times taken by the various kernel calls on the
    // 1D ViewFArrayKokkos objects
    std::vector<std::vector<double>> vfak_1D_timings(num_kernels);

    // Variable the stores the dot product of vectors a and b
    real_t vfak_dot_1D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 1D ViewFArrayKokkos copy kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vfak_arr3(i) = vfak_arr1(i);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record copy kernel timing
        vfak_1D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 1D ViewFArrayKokkos scale kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vfak_arr2(i) = (scalar * vfak_arr3(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record scale kernel timing
        vfak_1D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // PERF_START(tag) // tag is a string, e.g., "copy"
        // 1D ViewFArrayKokkos sum kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Sum (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vfak_arr3(i) = (vfak_arr1(i) + vfak_arr2(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record sum kernel timing
        vfak_1D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count()); 

        // 1D ViewFArrayKokkos triad kernel

        //LIKWID_MARKER_START("1D_VFAK_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vfak_arr1(i) = (vfak_arr2(i) + (scalar * vfak_arr3(i)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_VFAK_TRIAD");

        // Record triad kernel timing
        vfak_1D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
       
        // 1D ViewFArrayKokkos dot product kernel

        //LIKWID_MARKER_START("1D_VFAK_DOT");
        vfak_dot_1D_fin_val = 0.0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += (vfak_arr1(i) * vfak_arr2(i));
        }, vfak_dot_1D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_VFAK_DOT");

        // Record dot product kernel timing
        vfak_1D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 1D ViewFArrayKokkos STREAM benchmark results
    real_t vfak_arr1_1D_err = 0;
    real_t vfak_arr2_1D_err = 0;
    real_t vfak_arr3_1D_err = 0;
    real_t vfak_dot_1D_err = std::fabs(dot_1D_fin_val - vfak_dot_1D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (vfak_arr1(i) - arr1_fin_val) >= 0 ? (vfak_arr1(i) - arr1_fin_val) : (arr1_fin_val - vfak_arr1(i));
    }, vfak_arr1_1D_err);
    Kokkos::fence();

    vfak_arr1_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (vfak_arr2(i) - arr2_fin_val) >= 0 ? (vfak_arr2(i) - arr2_fin_val) : (arr2_fin_val - vfak_arr2(i));
    }, vfak_arr2_1D_err);
    Kokkos::fence();

    vfak_arr2_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D VFAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (vfak_arr3(i) - arr3_fin_val) >= 0 ? (vfak_arr3(i) - arr3_fin_val) : (arr3_fin_val - vfak_arr3(i));
    }, vfak_arr3_1D_err);
    Kokkos::fence();

    vfak_arr3_1D_err /= nsize;

    if (vfak_arr1_1D_err > epsi) {
        std::cout << "Validation failed on vfak_arr1. Average error "
                  << vfak_arr1_1D_err << std::endl << std::endl;
    }

    if (vfak_arr2_1D_err > epsi) {
        std::cout << "Validation failed on vfak_arr2. Average error "
                  << vfak_arr2_1D_err << std::endl << std::endl;
    }

    if (vfak_arr3_1D_err > epsi) {
        std::cout << "Validation failed on vfak_arr3. Average error "
                  << vfak_arr3_1D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (vfak_dot_1D_err > 1.0E-8) {
        std::cout << "Validation failed on 1D VFAK dot product kernel. Error is "
                  << vfak_dot_1D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << vfak_dot_1D_fin_val 
                  << " but should be "  << dot_1D_fin_val
                  << std::endl << std::endl;
    }

    // Print kernel computation memory bandwidth table header
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "1D ViewFArrayKokkos STREAM benchmark results" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(vfak_1D_timings[ker].begin() + 1,
                                          vfak_1D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(vfak_1D_timings[ker].begin() + 1,
                                         vfak_1D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;
    
    ////////////////////////////////////////////////////////////////////////////
    // 1D CArrayKokkos STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // Vector that stores the times taken by the various kernel calls on the
    // 1D CArrayKokkos objects
    std::vector<std::vector<double>> cak_1D_timings(num_kernels);

    real_t cak_dot_1D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 1D CArrayKokkos copy kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (1D CAK)", nsize, KOKKOS_LAMBDA(const int i) {
                cak_arr3(i) = cak_arr1(i);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record copy kernel timing
        cak_1D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    
        // 1D CArrayKokkos scale kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale (1D CAK)", nsize, KOKKOS_LAMBDA(const int i) {
                cak_arr2(i) = (scalar * cak_arr3(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record scale kernel timing
        cak_1D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    
        // PERF_START(tag) // tag is a string, e.g., "copy"
        // 1D CArrayKokkos sum kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Sum (1D CAK)", nsize, KOKKOS_LAMBDA(const int i) {
                cak_arr3(i) = (cak_arr1(i) + cak_arr2(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record sum kernel timing
        cak_1D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count()); 

        // 1D CArrayKokkos triad kernel

        //LIKWID_MARKER_START("1D_CAK_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad (1D CAK)", nsize, KOKKOS_LAMBDA(const int i) {
                cak_arr1(i) = (cak_arr2(i) + (scalar * cak_arr3(i)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_CAK_TRIAD");

        // Record triad kernel timing
        cak_1D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
   
        // 1D CArrayKokkos dot product kernel

        //LIKWID_MARKER_START("1D_CAK_DOT");
        cak_dot_1D_fin_val = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (1D CAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += (cak_arr1(i) * cak_arr2(i));
        }, cak_dot_1D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_CAK_DOT");

        // Record dot product kernel timing
        cak_1D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 1D CArrayKokkos STREAM benchmark results
    real_t cak_arr1_1D_err = 0;
    real_t cak_arr2_1D_err = 0;
    real_t cak_arr3_1D_err = 0;
    real_t cak_dot_1D_err = std::fabs(dot_1D_fin_val - cak_dot_1D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (1D CAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (cak_arr1(i) - arr1_fin_val) >= 0 ? (cak_arr1(i) - arr1_fin_val) : (arr1_fin_val - cak_arr1(i));
    }, cak_arr1_1D_err);
    Kokkos::fence();

    cak_arr1_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D CAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (cak_arr2(i) - arr2_fin_val) >= 0 ? (cak_arr2(i) - arr2_fin_val) : (arr2_fin_val - cak_arr2(i));
    }, cak_arr2_1D_err);
    Kokkos::fence();

    cak_arr2_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D CAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (cak_arr3(i) - arr3_fin_val) >= 0 ? (cak_arr3(i) - arr3_fin_val) : (arr3_fin_val - cak_arr3(i));
    }, cak_arr3_1D_err);
    Kokkos::fence();

    cak_arr3_1D_err /= nsize;

    if (cak_arr1_1D_err > epsi) {
        std::cout << "Validation failed on cak_arr1. Average error "
                  << cak_arr1_1D_err << std::endl << std::endl;
    }

    if (cak_arr2_1D_err > epsi) {
        std::cout << "Validation failed on cak_arr2. Average error "
                  << cak_arr2_1D_err << std::endl << std::endl;
    }

    if (cak_arr3_1D_err > epsi) {
        std::cout << "Validation failed on cak_arr3. Average error "
                  << cak_arr3_1D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (cak_dot_1D_err > 1.0E-8) {
        std::cout << "Validation failed on 1D CAK dot product kernel. Error is "
                  << cak_dot_1D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << cak_dot_1D_fin_val 
                  << " but should be "  << dot_1D_fin_val
                  << std::endl << std::endl;
    }

    // Print kernel computation memory bandwidth table header
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "1D CArrayKokkos STREAM benchmark results" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(cak_1D_timings[ker].begin() + 1,
                                          cak_1D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(cak_1D_timings[ker].begin() + 1,
                                         cak_1D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // 1D ViewCArrayKokkos STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////
    
    // Initialize the 1D ViewCArrayKokkos objects
    Kokkos::parallel_for("Initialize (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i) {
            vcak_arr1(i) = arr1_init_val;
            vcak_arr2(i) = arr2_init_val;
            vcak_arr3(i) = arr3_init_val;
            });
    Kokkos::fence();

    // Vector that stores the times taken by the various kernel calls on the
    // 1D ViewCArrayKokkos objects
    std::vector<std::vector<double>> vcak_1D_timings(num_kernels);

    // Variable the stores the dot product of vectors a and b
    real_t vcak_dot_1D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 1D ViewCArrayKokkos copy kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vcak_arr3(i) = vcak_arr1(i);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record copy kernel timing
        vcak_1D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 1D ViewCArrayKokkos scale kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vcak_arr2(i) = (scalar * vcak_arr3(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record scale kernel timing
        vcak_1D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // PERF_START(tag) // tag is a string, e.g., "copy"
        // 1D ViewCArrayKokkos sum kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Sum (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vcak_arr3(i) = (vcak_arr1(i) + vcak_arr2(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record sum kernel timing
        vcak_1D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count()); 

        // 1D ViewFArrayKokkos triad kernel

        //LIKWID_MARKER_START("1D_VCAK_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i) {
                vcak_arr1(i) = (vcak_arr2(i) + (scalar * vcak_arr3(i)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_VCAK_TRIAD");

        // Record triad kernel timing
        vcak_1D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
       
        // 1D ViewFArrayKokkos dot product kernel

        //LIKWID_MARKER_START("1D_VCAK_DOT");
        vcak_dot_1D_fin_val = 0.0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += (vcak_arr1(i) * vcak_arr2(i));
        }, vcak_dot_1D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_VCAK_DOT");

        // Record dot product kernel timing
        vcak_1D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 1D ViewCArrayKokkos STREAM benchmark results
    real_t vcak_arr1_1D_err = 0;
    real_t vcak_arr2_1D_err = 0;
    real_t vcak_arr3_1D_err = 0;
    real_t vcak_dot_1D_err = std::fabs(dot_1D_fin_val - vcak_dot_1D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (vcak_arr1(i) - arr1_fin_val) >= 0 ? (vcak_arr1(i) - arr1_fin_val) : (arr1_fin_val - vcak_arr1(i));
    }, vcak_arr1_1D_err);
    Kokkos::fence();

    vcak_arr1_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (vcak_arr2(i) - arr2_fin_val) >= 0 ? (vcak_arr2(i) - arr2_fin_val) : (arr2_fin_val - vcak_arr2(i));
    }, vcak_arr2_1D_err);
    Kokkos::fence();

    vcak_arr2_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D VCAK)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (vcak_arr3(i) - arr3_fin_val) >= 0 ? (vcak_arr3(i) - arr3_fin_val) : (arr3_fin_val - vcak_arr3(i));
    }, vcak_arr3_1D_err);
    Kokkos::fence();

    vcak_arr3_1D_err /= nsize;

    if (vcak_arr1_1D_err > epsi) {
        std::cout << "Validation failed on vcak_arr1. Average error "
                  << vcak_arr1_1D_err << std::endl << std::endl;
    }

    if (vcak_arr2_1D_err > epsi) {
        std::cout << "Validation failed on vcak_arr2. Average error "
                  << vcak_arr2_1D_err << std::endl << std::endl;
    }

    if (vcak_arr3_1D_err > epsi) {
        std::cout << "Validation failed on vcak_arr3. Average error "
                  << vcak_arr3_1D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (vcak_dot_1D_err > 1.0E-8) {
        std::cout << "Validation failed on 1D VCAK dot product kernel. Error is "
                  << vcak_dot_1D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << vcak_dot_1D_fin_val 
                  << " but should be "  << dot_1D_fin_val
                  << std::endl << std::endl;
    }

    // Print kernel computation memory bandwidth table header
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "1D ViewCArrayKokkos STREAM benchmark results" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(vcak_1D_timings[ker].begin() + 1,
                                          vcak_1D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(vcak_1D_timings[ker].begin() + 1,
                                         vcak_1D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // 1D Kokkos View STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // Vector that stores the times taken by the various kernel calls on the
    // 1D Kokkos Views objects
    std::vector<std::vector<double>> kv_1D_timings(num_kernels);

    real_t kv_dot_1D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 1D Kokkos View copy kernel 

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr3(i) = kv_arr1(i);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record copy kernel timing
        kv_1D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
   
        // 1D Kokkos View scale kernel

        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr2(i) = (scalar * kv_arr3(i));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        
        // Record scale kernel timing
        kv_1D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 1D Kokkos View sum kernel

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Sum (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr3(i) = kv_arr1(i) + kv_arr2(i);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record 1D Kokkos View sum kernel timing
        kv_1D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());   

        // 1D Kokkos View triad kernel

        //LIKWID_MARKER_START("1D_KV_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr1(i) = (kv_arr2(i) + (scalar * kv_arr3(i)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_KV_TRIAD");

        // Record triad kernel timing
        kv_1D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 1D Kokkos View dot product kernel

        //LIKWID_MARKER_START("1D_KV_DOT");
        kv_dot_1D_fin_val = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (1D KV)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += (kv_arr1(i) * kv_arr2(i));
        }, kv_dot_1D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("1D_KV_DOT");

        // Record dot product kernel timing
        kv_1D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 1D Kokkos View STREAM benchmark results
    real_t kv_arr1_1D_err = 0;
    real_t kv_arr2_1D_err = 0;
    real_t kv_arr3_1D_err = 0;
    real_t kv_dot_1D_err = std::fabs(dot_1D_fin_val - kv_dot_1D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (1D KV)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (kv_arr1(i) - arr1_fin_val) >= 0 ? (kv_arr1(i) - arr1_fin_val) : (arr1_fin_val - kv_arr1(i));
    }, kv_arr1_1D_err);
    Kokkos::fence();

    kv_arr1_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D KV)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (kv_arr2(i) - arr2_fin_val) >= 0 ? (kv_arr2(i) - arr2_fin_val) : (arr2_fin_val - kv_arr2(i));
    }, kv_arr2_1D_err);
    Kokkos::fence();

    kv_arr2_1D_err /= nsize;

    Kokkos::parallel_reduce("arr1 Error (1D KV)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
            tmp += (kv_arr3(i) - arr3_fin_val) >= 0 ? (kv_arr3(i) - arr3_fin_val) : (arr3_fin_val - kv_arr3(i));
    }, kv_arr3_1D_err);
    Kokkos::fence();

    kv_arr3_1D_err /= nsize;

    if (kv_arr1_1D_err > epsi) {
        std::cout << "Validation failed on kv_arr1. Average error "
                  << kv_arr1_1D_err << std::endl << std::endl;
    }

    if (kv_arr2_1D_err > epsi) {
        std::cout << "Validation failed on kv_arr2. Average error "
                  << kv_arr2_1D_err << std::endl << std::endl;
    }

    if (kv_arr3_1D_err > epsi) {
        std::cout << "Validation failed on kv_arr3. Average error "
                  << kv_arr3_1D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (kv_dot_1D_err > 1.0E-8) {
        std::cout << "Validation failed on 1D KV dot product kernel. Error is "
                  << kv_dot_1D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << kv_dot_1D_fin_val 
                  << " but should be "  << dot_1D_fin_val
                  << std::endl << std::endl;
    }

    // Print kernel computation memory bandwidth table header
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "1D Kokkos View STREAM benchmark results" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(kv_1D_timings[ker].begin() + 1,
                                          kv_1D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(kv_1D_timings[ker].begin() + 1,
                                         kv_1D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    /***************************************************************************
     * 3D array STREAM benchmark suite
     **************************************************************************/

    real_t dot_3D_fin_val = arr1_fin_val * arr2_fin_val * ARRAY_SIZE_3D;

    std::cout << "-------------------------------------------------------------"
              << std::endl;
    std::cout << "3D array STREAM benchmark suite" << std::endl;
    std::cout << "-------------------------------------------------------------"
              << std::endl
              << std::endl;

    std::cout << "Running kernels " << repeat << " times" << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << "Number of 3D array elements: " << ARRAY_SIZE_3D 
              << " (" << nsize_3D << " elements in each dimension)" 
              << std::endl;
    std::cout << std::endl;

    // std::cout << "3D information" << std::endl;
    // std::cout << "--------------" << std::endl;
    std::cout << std::setprecision(1) << std::fixed
              << "3D array size: " << (ARRAY_SIZE_3D * sizeof(real_t) * 1.0E-6) << "MB"
              << " (" << (ARRAY_SIZE_3D * sizeof(real_t) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size (3D): " << (3.0 * ARRAY_SIZE_3D * sizeof(real_t) * 1.0E-6) << "MB"
              << " (" << (3.0 * ARRAY_SIZE_3D * sizeof(real_t) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout << std::endl;

    std::cout.precision(ss);

    // Create "regularly" allocated 3D C++ arrays
    //real_t*** reg_arr1_3D = new real_t**[nsize_3D];
    //real_t*** reg_arr2_3D = new real_t**[nsize_3D];
    //real_t*** reg_arr3_3D = new real_t**[nsize_3D];

    /*
    for (int i = 0; i < nsize_3D; i++) {
        reg_arr1_3D[i] = new real_t*[nsize_3D];
        reg_arr2_3D[i] = new real_t*[nsize_3D];
        reg_arr3_3D[i] = new real_t*[nsize_3D];

        for (int j = 0; j < nsize_3D; j++) {
            reg_arr1_3D[i][j] = new real_t[nsize_3D];
            reg_arr2_3D[i][j] = new real_t[nsize_3D];
            reg_arr3_3D[i][j] = new real_t[nsize_3D];
        }
    }
    */

    // Create 3D FArrayKokkos objects
    auto fak_arr1_3D = FArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto fak_arr2_3D = FArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto fak_arr3_3D = FArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);

    // Create 3D ViewFArrayKokkos objects
    auto vfak_arr1_3D = ViewFArrayKokkos <real_t> (fak_arr1_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);
    auto vfak_arr2_3D = ViewFArrayKokkos <real_t> (fak_arr2_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);
    auto vfak_arr3_3D = ViewFArrayKokkos <real_t> (fak_arr3_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);

    // Create 3D CArrayKokkos objects
    auto cak_arr1_3D = CArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto cak_arr2_3D = CArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto cak_arr3_3D = CArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);

    // Create 3D ViewCArrayKokkos objects
    auto vcak_arr1_3D = ViewCArrayKokkos <real_t> (cak_arr1_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);
    auto vcak_arr2_3D = ViewCArrayKokkos <real_t> (cak_arr2_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);
    auto vcak_arr3_3D = ViewCArrayKokkos <real_t> (cak_arr3_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);

    // Create 3D Kokkos View objects
    RMatrix3D kv_arr1_3D("kv_arr1_3D", nsize_3D, nsize_3D, nsize_3D);
    RMatrix3D kv_arr2_3D("kv_arr2_3D", nsize_3D, nsize_3D, nsize_3D);
    RMatrix3D kv_arr3_3D("kv_arr3_3D", nsize_3D, nsize_3D, nsize_3D);

    policy3D array_type_STREAM = policy3D({0, 0, 0},
                                          {nsize_3D, nsize_3D, nsize_3D});

    // Initialize 3D FArrayKokkos objects
    Kokkos::parallel_for("Initialize (3D FAK)", array_type_STREAM,
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
            // Initialize 3D FArrayKokkos objects
            fak_arr1_3D(i, j, k) = arr1_init_val;
            fak_arr2_3D(i, j, k) = arr2_init_val;
            fak_arr3_3D(i, j, k) = arr3_init_val;
            });
    Kokkos::fence();

    // Initialize 3D CArrayKokkos and 3D Kokkos View objects
    Kokkos::parallel_for("Initialize (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
            // Initialize "regularly" allocated 3D C++ arrays
            //reg_arr1_3D[i][j][k] = 1.0;
            //reg_arr2_3D[i][j][k] = 2.0;
            //reg_arr3_3D[i][j][k] = 0.0;

            // Initialize 3D CArrayKokkos objects
            cak_arr1_3D(i, j, k) = arr1_init_val;
            cak_arr2_3D(i, j, k) = arr2_init_val;
            cak_arr3_3D(i, j, k) = arr3_init_val;
    
            // Initialize 3D Kokkos View objects
            kv_arr1_3D(i, j, k) = arr1_init_val;
            kv_arr2_3D(i, j, k) = arr2_init_val;
            kv_arr3_3D(i, j, k) = arr3_init_val;
            });
    Kokkos::fence();

    ////////////////////////////////////////////////////////////////////////////
    // 3D "conventionally" allocated C++ array STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // De-allocate memory used for "regularly" allocated 3D C++ arrays
    /*
    for (int i = 0; i < nsize_3D; i++) {
        for (int j = 0; j < nsize_3D; j++) {
            delete[] reg_arr1_3D[i][j];
            delete[] reg_arr2_3D[i][j];
            delete[] reg_arr3_3D[i][j];
        }
        delete[] reg_arr1_3D[i];
        delete[] reg_arr2_3D[i];
        delete[] reg_arr3_3D[i];
    }

    delete[] reg_arr1_3D;
    delete[] reg_arr2_3D;
    delete[] reg_arr3_3D;
    */

    ////////////////////////////////////////////////////////////////////////////
    // 3D FArrayKokkos STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////    

    // Vector that stores the times taken by the various kernel calls on the
    // 3D FArrayKokkos objects
    std::vector<std::vector<double>> fak_3D_timings(num_kernels);

    real_t fak_dot_3D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 3D FArrayKokkos copy kernel  

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (3D FAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                fak_arr3_3D(i, j, k) = fak_arr1_3D(i, j, k);
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();

        // Record 3D FArrayKokkos copy kernel timing
        fak_3D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D FArrayKokkos scale kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Scale (3D FAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                fak_arr2_3D(i, j, k) = (scalar * fak_arr3_3D(i, j, k));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();

        // Record 3D FArrayKokkos scale kernel timing
        fak_3D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D FArrayKokkos sum kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Sum (3D FAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                fak_arr3_3D(i, j, k) = (fak_arr1_3D(i, j, k) + fak_arr2_3D(i, j, k));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();

        // Record 3D FArrayKokkos sum kernel timing
        fak_3D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D FArrayKokkos triad kernel

        //LIKWID_MARKER_START("3D_FAK_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Triad (3D FAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                fak_arr1_3D(i, j, k) = (fak_arr2_3D(i, j, k) + (scalar * fak_arr3_3D(i, j, k)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("3D_FAK_TRIAD");

        // Record 3D FArrayKokkos triad kernel timing
        fak_3D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D FArrayKokkos dot product kernel

        //LIKWID_MARKER_START("3D_FAK_DOT");
        fak_dot_3D_fin_val = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (3D FAK)", array_type_STREAM, 
                                KOKKOS_LAMBDA(const int i, const int j, 
                                              const int k, real_t& tmp) {
                tmp += (fak_arr1_3D(i, j, k) * fak_arr2_3D(i, j, k));
        }, fak_dot_3D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("3D_FAK_DOT");

        // Record dot product kernel timing
        fak_3D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 3D FArrayKokkos STREAM benchmark results
    real_t fak_arr1_3D_err = 0;
    real_t fak_arr2_3D_err = 0;
    real_t fak_arr3_3D_err = 0;
    real_t fak_dot_3D_err = std::fabs(dot_3D_fin_val - fak_dot_3D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (3D FAK)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int i, const int j, 
                                          const int k, real_t& tmp) {
            tmp += (fak_arr1_3D(i, j, k) - arr1_fin_val) >= 0
                   ? (fak_arr1_3D(i, j, k) - arr1_fin_val)
                   : (arr1_fin_val - fak_arr1_3D(i, j, k));
    }, fak_arr1_3D_err);
    Kokkos::fence();

    fak_arr1_3D_err /= ARRAY_SIZE_3D;

    Kokkos::parallel_reduce("arr2 Error (3D FAK)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int i, const int j, 
                                          const int k, real_t& tmp) {
            tmp += (fak_arr2_3D(i, j, k) - arr2_fin_val) >= 0
                   ? (fak_arr2_3D(i, j, k) - arr2_fin_val)
                   : (arr2_fin_val - fak_arr2_3D(i, j, k));
    }, fak_arr2_3D_err);
    Kokkos::fence();

    fak_arr2_3D_err /= ARRAY_SIZE_3D;

    Kokkos::parallel_reduce("arr3 Error (3D FAK)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int i, const int j, 
                                          const int k, real_t& tmp) {
            tmp += (fak_arr3_3D(i, j, k) - arr3_fin_val) >= 0
                   ? (fak_arr3_3D(i, j, k) - arr3_fin_val)
                   : (arr3_fin_val - fak_arr3_3D(i, j, k));
    }, fak_arr3_3D_err);
    Kokkos::fence();

    fak_arr3_3D_err /= ARRAY_SIZE_3D;

    if (fak_arr1_3D_err > epsi) {
        std::cout << "Validation failed on fak_arr1_3D. Average error "
                  << fak_arr1_3D_err << std::endl << std::endl;
    }

    if (fak_arr2_3D_err > epsi) {
        std::cout << "Validation failed on fak_arr2_3D. Average error "
                  << fak_arr2_3D_err << std::endl << std::endl;
    }

    if (fak_arr3_3D_err > epsi) {
        std::cout << "Validation failed on fak_arr3_3D. Average error "
                  << fak_arr3_3D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (fak_dot_3D_err > 1.0E-8) {
        std::cout << "Validation failed on 3D FAK dot product kernel. Error is "
                  << fak_dot_3D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << fak_dot_3D_fin_val 
                  << " but should be "  << dot_3D_fin_val
                  << std::endl << std::endl;
    }
    
    // Print kernel computation memory bandwidth table header
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "3D FArrayKokkos STREAM benchmark results" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(fak_3D_timings[ker].begin() + 1,
                                          fak_3D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(fak_3D_timings[ker].begin() + 1,
                                         fak_3D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes_3D[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // 3D CArrayKokkos STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // Vector that stores the times taken by the various kernel calls on the
    // 3D CArrayKokkos objects
    std::vector<std::vector<double>> cak_3D_timings(num_kernels);   

    real_t cak_dot_3D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 3D CArrayKokkos copy kernel  

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (3D CAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int k, const int j, const int i) {
                cak_arr3_3D(i, j, k) = cak_arr1_3D(i, j, k);
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();

        // Record 3D CArrayKokkos copy kernel timing
        cak_3D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D CArrayKokkos scale kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Scale (3D CAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int k, const int j, const int i) {
                cak_arr2_3D(i, j, k) = (scalar * cak_arr3_3D(i, j, k));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();

        // Record 3D CArrayKokkos scale kernel timing
        cak_3D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D CArrayKokkos sum kernel

        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Sum (3D CAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int k, const int j, const int i) {
                cak_arr3_3D(i, j, k) = (cak_arr1_3D(i, j, k) + cak_arr2_3D(i, j, k));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();

        // Record 3D CArrayKokkos sum kernel timing
        cak_3D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D CArrayKokkos triad kernel

        //LIKWID_MARKER_START("3D_CAK_TRIAD");
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Triad (3D CAK)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int k, const int j, const int i) {
                cak_arr1_3D(i, j, k) = (cak_arr2_3D(i, j, k) + (scalar * cak_arr3_3D(i, j, k)));
                });
        Kokkos::fence();
        
        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("3D_CAK_TRIAD");

        // Record 3D CArrayKokkos triad kernel timing
        cak_3D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D CArrayKokkos dot product kernel

        //LIKWID_MARKER_START("3D_CAK_DOT");
        cak_dot_3D_fin_val = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (3D CAK)", array_type_STREAM, 
                                KOKKOS_LAMBDA(const int k, const int j, 
                                              const int i, real_t& tmp) {
                tmp += (cak_arr1_3D(i, j, k) * cak_arr2_3D(i, j, k));
        }, cak_dot_3D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("3D_CAK_DOT");

        // Record dot product kernel timing
        cak_3D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 3D CArrayKokkos STREAM benchmark results
    real_t cak_arr1_3D_err = 0;
    real_t cak_arr2_3D_err = 0;
    real_t cak_arr3_3D_err = 0;
    real_t cak_dot_3D_err = std::fabs(dot_3D_fin_val - cak_dot_3D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (3D CAK)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int k, const int j, 
                                          const int i, real_t& tmp) {
            tmp += (cak_arr1_3D(i, j, k) - arr1_fin_val) >= 0
                   ? (cak_arr1_3D(i, j, k) - arr1_fin_val)
                   : (arr1_fin_val - cak_arr1_3D(i, j, k));
    }, cak_arr1_3D_err);
    Kokkos::fence();

    cak_arr1_3D_err /= ARRAY_SIZE_3D;

    Kokkos::parallel_reduce("arr2 Error (3D CAK)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int k, const int j, 
                                          const int i, real_t& tmp) {
            tmp += (cak_arr2_3D(i, j, k) - arr2_fin_val) >= 0
                   ? (cak_arr2_3D(i, j, k) - arr2_fin_val)
                   : (arr2_fin_val - cak_arr2_3D(i, j, k));
    }, cak_arr2_3D_err);
    Kokkos::fence();

    cak_arr2_3D_err /= ARRAY_SIZE_3D;

    Kokkos::parallel_reduce("arr3 Error (3D CAK)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int k, const int j, 
                                          const int i, real_t& tmp) {
            tmp += (cak_arr3_3D(i, j, k) - arr3_fin_val) >= 0
                   ? (cak_arr3_3D(i, j, k) - arr3_fin_val)
                   : (arr3_fin_val - cak_arr3_3D(i, j, k));
    }, cak_arr3_3D_err);
    Kokkos::fence();

    cak_arr3_3D_err /= ARRAY_SIZE_3D;

    if (cak_arr1_3D_err > epsi) {
        std::cout << "Validation failed on cak_arr1_3D. Average error "
                  << cak_arr1_3D_err << std::endl << std::endl;
    }

    if (cak_arr2_3D_err > epsi) {
        std::cout << "Validation failed on cak_arr2_3D. Average error "
                  << cak_arr2_3D_err << std::endl << std::endl;
    }

    if (cak_arr3_3D_err > epsi) {
        std::cout << "Validation failed on cak_arr3_3D. Average error "
                  << cak_arr3_3D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (cak_dot_3D_err > 1.0E-8) {
        std::cout << "Validation failed on 3D CAK dot product kernel. Error is "
                  << cak_dot_3D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << cak_dot_3D_fin_val 
                  << " but should be "  << dot_3D_fin_val
                  << std::endl << std::endl;
    }
    
    // Print kernel computation memory bandwidth table header
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "3D CArrayKokkos STREAM benchmark results" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(cak_3D_timings[ker].begin() + 1,
                                          cak_3D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(cak_3D_timings[ker].begin() + 1,
                                         cak_3D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes_3D[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // 3D Kokkos View STREAM benchmark suite 
    ////////////////////////////////////////////////////////////////////////////

    // Vector that stores the times taken by the various kernel calls on the
    // 3D Kokkos Views objects
    std::vector<std::vector<double>> kv_3D_timings(num_kernels);

    real_t kv_dot_3D_fin_val;

    for (int iter = 0; iter < repeat; iter++) {
        // 3D Kokkos View copy kernel

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Copy (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr3_3D(i, j, k) = kv_arr1_3D(i, j, k);
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record 3D Kokkos View copy kernel timing
        kv_3D_timings[0].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D Kokkos View scale kernel

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Scale (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr2_3D(i, j, k) = (scalar * kv_arr3_3D(i, j, k));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record 3D Kokkos View scale kernel timing
        kv_3D_timings[1].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D Kokkos View sum kernel

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Sum (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr3_3D(i, j, k) = (kv_arr1_3D(i, j, k) + kv_arr2_3D(i, j, k));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();

        // Record 3D Kokkos View sum kernel timing
        kv_3D_timings[2].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D Kokkos View triad kernel

        //LIKWID_MARKER_START("3D_KV_TRIAD");
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Triad (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr1_3D(i, j, k) = (kv_arr2_3D(i, j, k) + (scalar * kv_arr3_3D(i, j, k)));
                });
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("3D_KV_TRIAD");

        // Record 3D Kokkos View triad kernel timing
        kv_3D_timings[3].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());

        // 3D Kokkos View dot product kernel

        //LIKWID_MARKER_START("3D_KV_DOT");
        kv_dot_3D_fin_val = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (3D KV)", array_type_STREAM, 
                               KOKKOS_LAMBDA(const int i, const int j,
                                             const int k, real_t& tmp) {
                tmp += (kv_arr1_3D(i, j, k) * kv_arr2_3D(i, j, k));
        }, kv_dot_3D_fin_val);
        Kokkos::fence();

        end = std::chrono::high_resolution_clock::now();
        //LIKWID_MARKER_STOP("3D_KV_DOT");

        kv_3D_timings[4].push_back(std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count());
    }

    // Verify 3D CArrayKokkos STREAM benchmark results
    real_t kv_arr1_3D_err = 0;
    real_t kv_arr2_3D_err = 0;
    real_t kv_arr3_3D_err = 0;
    real_t kv_dot_3D_err = std::fabs(dot_3D_fin_val - kv_dot_3D_fin_val);

    Kokkos::parallel_reduce("arr1 Error (3D KV)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int i, const int j, 
                                          const int k, real_t& tmp) {
            tmp += (kv_arr1_3D(i, j, k) - arr1_fin_val) >= 0
                   ? (kv_arr1_3D(i, j, k) - arr1_fin_val)
                   : (arr1_fin_val - kv_arr1_3D(i, j, k));
    }, kv_arr1_3D_err);
    Kokkos::fence();

    kv_arr1_3D_err /= ARRAY_SIZE_3D;

    Kokkos::parallel_reduce("arr2 Error (3D KV)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int i, const int j, 
                                          const int k, real_t& tmp) {
            tmp += (kv_arr2_3D(i, j, k) - arr2_fin_val) >= 0
                   ? (kv_arr2_3D(i, j, k) - arr2_fin_val)
                   : (arr2_fin_val - kv_arr2_3D(i, j, k));
    }, kv_arr2_3D_err);
    Kokkos::fence();

    kv_arr2_3D_err /= ARRAY_SIZE_3D;

    Kokkos::parallel_reduce("arr3 Error (3D KV)", array_type_STREAM, 
                            KOKKOS_LAMBDA(const int i, const int j, 
                                          const int k, real_t& tmp) {
            tmp += (kv_arr3_3D(i, j, k) - arr3_fin_val) >= 0
                   ? (kv_arr3_3D(i, j, k) - arr3_fin_val)
                   : (arr3_fin_val - kv_arr3_3D(i, j, k));
    }, kv_arr3_3D_err);
    Kokkos::fence();

    kv_arr3_3D_err /= ARRAY_SIZE_3D;

    if (kv_arr1_3D_err > epsi) {
        std::cout << "Validation failed on kv_arr1_3D. Average error "
                  << kv_arr1_3D_err << std::endl << std::endl;
    }

    if (kv_arr2_3D_err > epsi) {
        std::cout << "Validation failed on kv_arr2_3D. Average error "
                  << kv_arr2_3D_err << std::endl << std::endl;
    }

    if (kv_arr3_3D_err > epsi) {
        std::cout << "Validation failed on kv_arr3_3D. Average error "
                  << kv_arr3_3D_err << std::endl << std::endl;
    }

    // Check the dot product error up to 8 decimal places
    if (kv_dot_3D_err > 1.0E-8) {
        std::cout << "Validation failed on 3D KV dot product kernel. Error is "
                  << kv_dot_3D_err << std::endl << std::setprecision(15)
                  << "Dot product was " << kv_dot_3D_fin_val 
                  << " but should be "  << dot_3D_fin_val
                  << std::endl << std::endl;
    }

    // Print kernel computation memory bandwidth table header
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "3D Kokkos View STREAM benchmark results" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec)"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    std::cout << std::endl;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(kv_3D_timings[ker].begin() + 1,
                                          kv_3D_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(kv_3D_timings[ker].begin() + 1,
                                         kv_3D_timings[ker].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[ker]
                  << std::left << std::setw(12) << std::setprecision(3) <<
                  ((1.0E-6 * sizes_3D[ker]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.first
                  << std::left << std::setw(12) << std::setprecision(5) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5) << average
                  << std::endl;
    }
    
    std::cout << std::endl;

    }

    Kokkos::finalize();


    printf("--- finished ---\n");
    
    // Stop LIKWID
    LIKWID_MARKER_CLOSE;

    return 0;
}