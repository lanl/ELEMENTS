#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <chrono>         // To access timing calipers 
#include "matar.h"

// For LIKWID
// #ifdef LIKWID_PERFMON
// #include "likwid.h"
// #else
// #define LIKWID_PERFMON 1
// #define LIKWID_MARKER_INIT
// #define LIKWID_MARKER_THREADINIT
// #define LIKWID_MARKER_REGISTER(regionTag)
// #define LIKWID_MARKER_START(regionTag)
// #define LIKWID_MARKER_STOP(regionTag)
// #define LIKWID_MARKER_CLOSE
// #endif

// If we are compiling with LIKWID support and if we want to conditionally
// enable LIKWID markers, use the following macros
// #define PERF_INIT       if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_INIT;
// #define PERF_START(tag) if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_START(tag);
// #define PERF_STOP(tag)  if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_STOP(tag);
// #define PERF_CLOSE      if (LIKWID_PERFMON && LIKWID_ENABLED) LIKWID_MARKER_CLOSE;

int main() {
    int size_i = 5, size_j = 4, size_k = 3;
	const int size3 = 256;
	const int size1 = size3 * size3 * size3;
	const int repeat = 10;

	std::cout<<"Size of 1D problem: "<<size1<<"\n";
	std::cout<<"Size of 3D problem (each dimension): "<<size3<<"\n";
//#ifdef CPUBENCH
	//===========================================================
	//	cpu stream benchmarks
	//===========================================================	

	//1. create carrays

	//1D
	auto c_arr1 = CArray<double>(size1);
	auto c_arr2 = CArray<double>(size1);
	auto c_arr3 = CArray<double>(size1);

	//3D
	auto c_arr1_3d = CArray<double>(size3, size3, size3);
	auto c_arr2_3d = CArray<double>(size3, size3, size3);
	auto c_arr3_3d = CArray<double>(size3, size3, size3);

	//1D
	double* reg_arr1_1d = new double[size1];
	double* reg_arr2_1d = new double[size1];
	double* reg_arr3_1d = new double[size1];

	//3D
	double*** reg_arr1_3d = new double**[size3];
	double*** reg_arr2_3d = new double**[size3];
	double*** reg_arr3_3d = new double**[size3];
	
	for(size_t i = 0; i < size3; i++) {
	 reg_arr1_3d[i] = new double*[size3];	
	 reg_arr2_3d[i] = new double*[size3];	
	 reg_arr3_3d[i] = new double*[size3];	
	  for(size_t j = 0; j < size3; j++) {
		reg_arr1_3d[i][j] = new double[size3];
		reg_arr2_3d[i][j] = new double[size3];
		reg_arr3_3d[i][j] = new double[size3];
	   }
	  }

	//=========================================================================
	//	initialize
	//=========================================================================

	//1D
	std::cout<<"initializing 1D arrays\n";
#pragma omp simd
	for(size_t i = 0; i < size1; i++) {
		c_arr1(i) = 1.0;
		reg_arr1_1d[i] = 1.0;
	  }

	//3D
	std::cout<<"initializing 3D arrays\n";
#pragma omp simd collapse(3)
	for(size_t i = 0; i < size3; i++) {
	 for(size_t j = 0; j < size3; j++) {
	  for(size_t k = 0; k < size3; k++) {
		c_arr1_3d(i,j,k) = 1.0;
		reg_arr1_3d[i][j][k] = 1.0;
		}
	   }
	  }

	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~COPY TESTS~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";
	
	//=========================================================================
	//	COPY TEST
	// array2 = array1
	//=========================================================================
	
	//==============================
	//	CArrays
	//==============================	

	//1D
	double carr_copy_times[repeat];	//array to hold times for carray copy
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_c_1d = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
		c_arr2(i) = c_arr1(i);
		}
		auto end_copy_c_1d = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_copy_c_1d = end_copy_c_1d - start_copy_c_1d;
		carr_copy_times[t] = tot_copy_c_1d.count();
	}
	
	//calculate average
	double avg_copy_c, tot_copy_c;
	for(size_t i = 0; i < repeat; i++) {
		tot_copy_c += carr_copy_times[i];
	}	
	avg_copy_c = tot_copy_c / (double) repeat;
	std::cout<<"Average elapsed time for 1D Carray copy = "<<avg_copy_c<<"s\n";

	//3D
	double carr_copy_3dtimes[repeat]; 
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_c_3d = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			c_arr2_3d(i,j,k) = c_arr1_3d(i,j,k);
		   }
		  }
		 }
	   auto end_copy_c_3d = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_copy_c_3d = end_copy_c_3d - start_copy_c_3d;
	   carr_copy_3dtimes[t] = tot_copy_c_3d.count();
	}

	//calculate average
	double avg_copy_3dc, tot_copy_3dc;
	for(size_t i = 0; i < repeat; i++) {
		tot_copy_3dc += carr_copy_3dtimes[i];
	}
	avg_copy_3dc = tot_copy_3dc / (double) repeat;
	std::cout<<"Average elapsed time for 3D Carray copy = "<<avg_copy_3dc<<"s\n";

	std::cout<<"============================================================\n";
	
	//==============================
	//	Traditional Arrays
	//==============================	

	//1D
	double reg_copy_times[repeat];	//array to hold times for carray copy
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_r_1d = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
		reg_arr2_1d[i] = reg_arr1_1d[i];
		}
		auto end_copy_r_1d = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_copy_r_1d = end_copy_r_1d - start_copy_r_1d;
		reg_copy_times[t] = tot_copy_r_1d.count();
	}
	
	//calculate average
	double avg_copy_r, tot_copy_r;
	for(size_t i = 0; i < repeat; i++) {
		tot_copy_r += reg_copy_times[i];
	}	
	avg_copy_r = tot_copy_r / (double) repeat;
	std::cout<<"Average elapsed time for 1D trad copy = "<<avg_copy_r<<"s\n";

	//3D
	double reg_copy_3dtimes[repeat]; 
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_r_3d = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			reg_arr2_3d[i][j][k] = reg_arr1_3d[i][j][k];
		   }
		  }
		 }
	   auto end_copy_r_3d = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_copy_r_3d = end_copy_r_3d - start_copy_r_3d;
	   reg_copy_3dtimes[t] = tot_copy_r_3d.count();
	}

	//calculate average
	double avg_copy_3dr, tot_copy_3dr;
	for(size_t i = 0; i < repeat; i++) {
		tot_copy_3dr += reg_copy_3dtimes[i];
	}
	avg_copy_3dr = tot_copy_3dr / (double) repeat;
	std::cout<<"Average elapsed time for 3D trad copy = "<<avg_copy_3dr<<"s\n";


	std::cout<<"============================================================\n";

	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~SCALE TESTS~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	//=========================================================================
	//	SCALE TEST
	//	array2 = scale*array2
	//=========================================================================
	
	//1. CArrays
	
	//1D.
	double c_scale1d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale1dc = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
		 c_arr2(i) = 2.0*c_arr2(i);
		}
	   auto end_scale1dc = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_scale1dc = end_scale1dc - start_scale1dc;
	   c_scale1d_times[t] = tot_scale1dc.count();
	}

	//compute average
	double avg_scale1dc, tot_scale_1dc;
	for(size_t i = 0; i < repeat; i++) {
		tot_scale_1dc += c_scale1d_times[i];
	}
	avg_scale1dc = tot_scale_1dc / (double) repeat;
	std::cout<<"Averaged elapsed time for 1D Carray scale = "<<avg_scale1dc<<"s\n";

	//3D
	double c_scale3d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale3dc = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			c_arr2_3d(i,j,k) = 2.0*c_arr2_3d(i,j,k);
		   }
		  }
		 }
	   auto end_scale3dc = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_scale3dc = end_scale3dc - start_scale3dc;
	   c_scale3d_times[t] = tot_scale3dc.count();
	}

	//compute average
	double avg_scale3dc, tot_scale_3dc;
	for(size_t i = 0; i < repeat; i++) {
		tot_scale_3dc += c_scale3d_times[i];
	}
	avg_scale3dc = tot_scale_3dc / (double) repeat;
	std::cout<<"Average elapsed time for 3D Carray scale = "<<avg_scale3dc<<"s\n";

	std::cout<<"============================================================\n";

	//2. regular c++ arrays

	//1D.
	double reg_scale1d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale1dr = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
		 reg_arr2_1d[i] = 2.0*reg_arr2_1d[i];
		}
	   auto end_scale1dr = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_scale1dr = end_scale1dr - start_scale1dr;
	   reg_scale1d_times[t] = tot_scale1dr.count();
	}

	//compute average
	double avg_scale1dr, tot_scale_1dr;
	for(size_t i = 0; i < repeat; i++) {
		tot_scale_1dr += reg_scale1d_times[i];
	}
	avg_scale1dr = tot_scale_1dr / (double) repeat;
	std::cout<<"Averaged elapsed time for 1D trad scale = "<<avg_scale1dr<<"s\n";

	//3D
	double reg_scale3d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale3dr = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			reg_arr2_3d[i][j][k] = 2.0*reg_arr2_3d[i][j][k];
		   }
		  }
		 }
	   auto end_scale3dr = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_scale3dr = end_scale3dr - start_scale3dr;
	   reg_scale3d_times[t] = tot_scale3dr.count();
	}

	//compute average
	double avg_scale3dr, tot_scale_3dr;
	for(size_t i = 0; i < repeat; i++) {
		tot_scale_3dr += reg_scale3d_times[i];
	}
	avg_scale3dr = tot_scale_3dr / (double) repeat;
	std::cout<<"Average elapsed time for 3D trad scale = "<<avg_scale3dr<<"s\n";

	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~SUM  TESTS~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";
	//=========================================================================
	//	STREAM SUM  MATAR
	// arr3 = arr2 + arr1
	//=========================================================================

	//1. CARRAY
	
	//1D
	double sum_1dc[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_sum1dc = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
			c_arr3(i) = c_arr1(i) + c_arr2(i);
		}
	   auto end_sum1dc = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_1dc_sum = end_sum1dc - start_sum1dc;
	   sum_1dc[t] = tot_1dc_sum.count();
	}

	//compute average
	double avg_sum1dc, tot_sum_1dc;
	for(size_t i = 0; i < repeat; i++) {
		tot_sum_1dc += sum_1dc[i];
	}
	avg_sum1dc = tot_sum_1dc / (double) repeat;
	std::cout<<"Averate elapsed time for 1D CArray sum = " <<avg_sum1dc<<"s\n";

	//3d	
	double sum_3dc[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_sum3dc = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			c_arr3_3d(i,j,k) = c_arr2_3d(i,j,k) + c_arr1_3d(i,j,k);
		  }
		 }
		}
	   auto end_sum3dc = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_3dc_sum = end_sum3dc - start_sum3dc;
	   sum_3dc[t] = tot_3dc_sum.count();
	}

	//compute average
	double avg_sum3dc, tot_sum_3dc;
	for(size_t i = 0; i < repeat; i++) {
		tot_sum_3dc += sum_3dc[i];
	}
	avg_sum3dc = tot_sum_3dc / (double) repeat;
	std::cout<<"Averape elapsed time for 3D CArray sum = "<<avg_sum3dc<<"s\n";

	std::cout<<"=============================================================\n";

//#endif







	//=========================================================================
	//	STREAM TRIAD MATAR
	//=========================================================================

	//==============================
	//	CArrays
	//=============================	

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
    double scalar = 3.0;
    size_t nsize = 64 * 64 * 64 * 64;
    size_t nsize_3D = 512;
    size_t ARRAY_SIZE_3D = (nsize_3D * nsize_3D * nsize_3D); 

    std::streamsize ss = std::cout.precision();

    std::cout << "Running kernels " << repeat << " times" << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << "Number of 1D array elements: " << nsize << std::endl;
    std::cout << std::endl;

    // Conversion factor between megabytes (MB) and bytes: 1 byte = 10^(-6) MB
    // Conversion factor between bytes and gigabytes (GB): 1 byte = 10^(-9) GB
    std::cout << "1D information" << std::endl;
    std::cout << "--------------" << std::endl;
    std::cout << std::setprecision(1) << std::fixed
              << "1D array size: " << (nsize * sizeof(double) * 1.0E-6) << "MB"
              << " (" << (nsize * sizeof(double) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size (1D): " << (3.0 * nsize * sizeof(double) * 1.0E-6) << "MB"
              << " (" << (3.0 * nsize * sizeof(double) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout << std::endl;

    std::cout.precision(ss);

    // Create vector of vectors that will store timings for kernel runs
    // Rows correspond to kernels:
    //      1. Copy kernel
    //      2. Scale kernel
    //      3. Sum kernel
    //      4. Triad kernel
    //      5. Dot product kernel
    
    int num_kernels = 5;

    std::vector<std::vector<double>> matar_kokkos_timings(num_kernels);
    std::vector<std::vector<double>> kokkos_views_timings(num_kernels);

    std::vector<std::vector<double>> matar_kokkos_timings_3D(num_kernels);
    std::vector<std::vector<double>> kokkos_views_timings_3D(num_kernels);

    // Declare timers
    // std::chrono::high_resolution_clock::time_point begin;
    auto begin = std::chrono::high_resolution_clock::now();

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

    using policy2D = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;
    policy2D array_type = policy2D({0,0}, {size_i, size_j});
    policy2D matrix_type = policy2D({1,1}, {size_i+1, size_j+1});
    policy2D splice_type = policy2D({0,0}, {size_j, size_k});
    using policy3D = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
    policy3D array_type3 = policy3D({0,0,0}, {size_i, size_j, size_k});
    policy3D matrix_type3 = policy3D({1,1,1}, {size_i+1, size_j+1, size_k+1});
    
    /*------------------------- Stream Benchmarks ----------------------------------------*/

    // Kokkos stream benchmark test (for FMatrixKokkos and CArrayKokkos objects)
    // printf("\n Kokkos stream benchmark CArrayKokkos\n");

    // Create 1D FArrayKokkos objects
    auto arr1 = FArrayKokkos <real_t> (nsize);
    auto arr2 = FArrayKokkos <real_t> (nsize);
    auto arr3 = FArrayKokkos <real_t> (nsize);

    // Create 1D ViewFArrayKokkos objects
    auto arr1_v = ViewFArrayKokkos <real_t> (arr1.pointer(), nsize);
    auto arr2_v = ViewFArrayKokkos <real_t> (arr2.pointer(), nsize);
    auto arr3_v = ViewFArrayKokkos <real_t> (arr3.pointer(), nsize);

    // Create 1D Kokkos View objects
    RMatrix1D kv_arr1("kv_arr1", nsize);
    RMatrix1D kv_arr2("kv_arr2", nsize);
    RMatrix1D kv_arr3("kv_arr3", nsize);

    // Create 3D FArrayKokkos objects
    auto arr1_3D = FArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto arr2_3D = FArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);
    auto arr3_3D = FArrayKokkos <real_t> (nsize_3D, nsize_3D, nsize_3D);

    // Create 3D ViewFArrayKokkos objects
    auto arr1_3D_v = ViewFArrayKokkos <real_t> (arr1_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);
    auto arr2_3D_v = ViewFArrayKokkos <real_t> (arr2_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);
    auto arr3_3D_v = ViewFArrayKokkos <real_t> (arr3_3D.pointer(), nsize_3D, nsize_3D, nsize_3D);

    // Create 3D Kokkos View objects
    RMatrix3D kv_arr1_3D("kv_arr1_3D", nsize_3D, nsize_3D, nsize_3D);
    RMatrix3D kv_arr2_3D("kv_arr2_3D", nsize_3D, nsize_3D, nsize_3D);
    RMatrix3D kv_arr3_3D("kv_arr3_3D", nsize_3D, nsize_3D, nsize_3D);

    // printf("Stream benchmark with %d elements.\n", nsize);

    // Initialize 1D CArrayKokkos and 1D Kokkos View objects
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
                arr1(i) = 1.0;
                arr2(i) = 2.0;
                arr3(i) = 0.0;

                // Initialize 1D Kokkos View objects
                kv_arr1(i) = 1.0;
                kv_arr2(i) = 2.0;
                kv_arr3(i) = 0.0;
                });
    Kokkos::fence();

    // Initialize 3D CArrayKokkos and 3D Kokkos View objects
    policy3D array_type_STREAM = policy3D({0, 0, 0},
                                          {nsize_3D, nsize_3D, nsize_3D});

    Kokkos::parallel_for("Initialize (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr1_3D(i, j, k) = 1.0;
                arr2_3D(i, j, k) = 2.0;
                arr3_3D(i, j, k) = 0.0;

                kv_arr1_3D(i, j, k) = 1.0;
                kv_arr2_3D(i, j, k) = 2.0;
                kv_arr3_3D(i, j, k) = 0.0;
                });
    Kokkos::fence();

    // Perform 1D copy kernel 
    for (int iter = 0; iter < repeat; iter++) {
        // 1D CArrayKokkos copy kernel
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (1D)", nsize, KOKKOS_LAMBDA(const int i) {
                arr3(i) = arr1(i);
                });
        std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();
        

        // Record copy kernel timing
        matar_kokkos_timings[0].push_back(copy_time.count());

        // std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;
       
        // 1D Kokkos View copy kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Copy (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr3(i) = kv_arr1(i);
                });
        std::chrono::duration<double> copy_time_kv_1D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();


        // Record 1D Kokkos View copy kernel timing
        kokkos_views_timings[0].push_back(copy_time_kv_1D.count());
    }

    // Perform 3D copy kernel 
    for (int iter = 0; iter < repeat; iter++) {
        // 3D CArrayKokkos copy kernel
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr3_3D(i, j, k) = arr1_3D(i, j, k);
                });
        std::chrono::duration<double> copy_time_3D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();
        

        // Record 3D CArrayKokkos copy kernel timing
        matar_kokkos_timings_3D[0].push_back(copy_time_3D.count());

        // 3D Kokkos View copy kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Copy (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr3_3D(i, j, k) = kv_arr1_3D(i, j, k);
                });
        std::chrono::duration<double> copy_time_kv_3D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();


        // Record 3D Kokkos View copy kernel timing
        kokkos_views_timings_3D[0].push_back(copy_time_kv_3D.count());
    }

    // Perform 1D scaling kernel
    for (int iter = 0; iter < repeat; iter++) {
        // 1D CArrayKokkos scale kernel
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale (1D)", nsize, KOKKOS_LAMBDA(const int i) {
                arr2(i) = scalar * arr3(i);
                });
        std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();
        

        // Record scale kernel timing
        matar_kokkos_timings[1].push_back(scale_time.count());

        // std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;
       
        // 1D Kokkos View scale kernel 
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Scale (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr2(i) = scalar * kv_arr3(i);
                });
        std::chrono::duration<double> scale_time_kv_1D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();


        // Record 1D Kokkos View scale kernel timing
        kokkos_views_timings[1].push_back(scale_time_kv_1D.count());   
    }

    // Perform 3D scale kernel 
    for (int iter = 0; iter < repeat; iter++) {
        // 3D CArrayKokkos scale kernel
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Scale (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr2_3D(i, j, k) = scalar * arr3_3D(i, j, k);
                });
        std::chrono::duration<double> scale_time_3D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();
        

        // Record 3D CArrayKokkos scale kernel timing
        matar_kokkos_timings_3D[1].push_back(scale_time_3D.count());

        // 3D Kokkos View scale kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Scale (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr2_3D(i, j, k) = scalar * kv_arr3_3D(i, j, k);
                });
        std::chrono::duration<double> scale_time_kv_3D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();


        // Record 3D Kokkos View scale kernel timing
        kokkos_views_timings_3D[1].push_back(scale_time_kv_3D.count());
    }

    // Perform 1D sum kernel
    // PERF_START(tag) // tag is a string, e.g., "copy"
    for (int iter = 0; iter < repeat; iter++) {
        // 1D CArrayKokkos sum kernel
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Sum (1D)", nsize, KOKKOS_LAMBDA(const int i) {
                arr3(i) = arr1(i) + arr2(i);
                });
        std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();
        

        // Record sum kernel timing
        // matar_kokkos_timings(2, iter) = sum_time.count();
        matar_kokkos_timings[2].push_back(sum_time.count()); 

        // std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;

        // 1D Kokkos View sum kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Sum (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr3(i) = kv_arr1(i) + kv_arr2(i);
                });
        std::chrono::duration<double> sum_time_kv_1D = std::chrono::high_resolution_clock::now() - begin;
        Kokkos::fence();


        // Record 1D Kokkos View sum kernel timing
        kokkos_views_timings[2].push_back(sum_time_kv_1D.count());   
    }

    // Perform 3D sum kernel 
    for (int iter = 0; iter < repeat; iter++) {
        // 3D CArrayKokkos sum kernel
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Sum (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr3_3D(i, j, k) = arr1_3D(i, j, k) + arr2_3D(i, j, k);
                });
        Kokkos::fence();
        
        std::chrono::duration<double> sum_time_3D = std::chrono::high_resolution_clock::now() - begin;

        // Record 3D CArrayKokkos sum kernel timing
        matar_kokkos_timings_3D[2].push_back(sum_time_3D.count());

        // 3D Kokkos View sum kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Sum (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr3_3D(i, j, k) = kv_arr1_3D(i, j, k) + kv_arr2_3D(i, j, k);
                });
        Kokkos::fence();

        std::chrono::duration<double> sum_time_kv_3D = std::chrono::high_resolution_clock::now() - begin;

        // Record 3D Kokkos View sum kernel timing
        kokkos_views_timings_3D[2].push_back(sum_time_kv_3D.count());
    }

    // Perform triad kernel
    for (int iter = 0; iter < repeat; iter++) {
        // 1D CArrayKokkos triad kernel
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
                arr1(i) = arr2(i) + 
                                      (scalar * arr3(i));
                });
        Kokkos::fence();
        
        std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

        // Record triad kernel timing
        matar_kokkos_timings[3].push_back(triad_time.count());

        // std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;
        
        // 1D Kokkos View triad kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Triad (1D KV)", nsize, KOKKOS_LAMBDA(const int i) {
                kv_arr1(i) = kv_arr2(i) + (scalar * kv_arr3(i));
                });
        Kokkos::fence();

        std::chrono::duration<double> triad_time_kv_1D = std::chrono::high_resolution_clock::now() - begin;

        // Record 1D Kokkos View sum kernel timing
        kokkos_views_timings[3].push_back(triad_time_kv_1D.count());
    }

    // Perform 3D triad kernel 
    for (int iter = 0; iter < repeat; iter++) {
        // 3D CArrayKokkos triad kernel
        begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Triad (3D)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                arr1_3D(i, j, k) = arr2_3D(i, j, k) + (scalar * arr3_3D(i, j, k));
                });
        Kokkos::fence();
        
        std::chrono::duration<double> triad_time_3D = std::chrono::high_resolution_clock::now() - begin;

        // Record 3D CArrayKokkos triad kernel timing
        matar_kokkos_timings_3D[3].push_back(triad_time_3D.count());

        // 3D Kokkos View triad kernel
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for("Triad (3D KV)", array_type_STREAM, 
                         KOKKOS_LAMBDA(const int i, const int j, const int k) {
                kv_arr1_3D(i, j, k) = kv_arr2_3D(i, j, k) + (scalar * kv_arr3_3D(i, j, k));
                });
        Kokkos::fence();

        std::chrono::duration<double> triad_time_kv_3D = std::chrono::high_resolution_clock::now() - begin;

        // Record 3D Kokkos View triad kernel timing
        kokkos_views_timings_3D[3].push_back(triad_time_kv_3D.count());
    }

    // Perform 1D dot product kernel
    for (int iter = 0; iter < repeat; iter++) {
        // 1D CArrayKokkos dot product kernel
        real_t sum = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += arr3(i) * arr3(i);
        }, sum);
        Kokkos::fence();

        std::chrono::duration<double> dot_time = std::chrono::high_resolution_clock::now() - begin;

        // Record dot product kernel timing
        matar_kokkos_timings[4].push_back(dot_time.count());

        // 1D Kokkos View dot product kernel
        sum = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (1D KV)", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += kv_arr3(i) * kv_arr3(i);
        }, sum);
        Kokkos::fence();

        std::chrono::duration<double> dot_time_kv_1D = std::chrono::high_resolution_clock::now() - begin;

        kokkos_views_timings[4].push_back(dot_time_kv_1D.count());
    }

    // Perform 3D dot product kernel
    for (int iter = 0; iter < repeat; iter++) {
        // 3D CArrayKokkos dot product kernel
        real_t sum = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (3D)", array_type_STREAM, 
                                KOKKOS_LAMBDA(const int i, const int j, 
                                              const int k, real_t& tmp) {
                tmp += (arr3_3D(i, j, k) * arr3_3D(i, j, k));
        }, sum);
        Kokkos::fence();

        std::chrono::duration<double> dot_time_3D = std::chrono::high_resolution_clock::now() - begin;

        // Record dot product kernel timing
        matar_kokkos_timings_3D[4].push_back(dot_time_3D.count());

        // 3D Kokkos View dot product kernel
        sum = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product (3D KV)", array_type_STREAM, 
                               KOKKOS_LAMBDA(const int i, const int j,
                                             const int k, real_t& tmp) {
                tmp += (kv_arr3_3D(i, j, k) * kv_arr3_3D(i, j, k));
        }, sum);
        Kokkos::fence();

        std::chrono::duration<double> dot_time_kv_3D = std::chrono::high_resolution_clock::now() - begin;

        kokkos_views_timings_3D[4].push_back(dot_time_kv_3D.count());
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
        auto minmax = std::minmax_element(matar_kokkos_timings[ker].begin() + 1,
                                          matar_kokkos_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(matar_kokkos_timings[ker].begin() + 1,
                                         matar_kokkos_timings[ker].end(),
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
        auto minmax = std::minmax_element(kokkos_views_timings[ker].begin() + 1,
                                          kokkos_views_timings[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(kokkos_views_timings[ker].begin() + 1,
                                         kokkos_views_timings[ker].end(),
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

    std::cout << "Running kernels " << repeat << " times" << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << "Number of 3D array elements: " << ARRAY_SIZE_3D 
              << " (" << nsize_3D << " elements in each dimension)" 
              << std::endl;
    std::cout << std::endl;

    std::cout << "3D information" << std::endl;
    std::cout << "--------------" << std::endl;
    std::cout << std::setprecision(1) << std::fixed
              << "3D array size: " << (ARRAY_SIZE_3D * sizeof(double) * 1.0E-6) << "MB"
              << " (" << (ARRAY_SIZE_3D * sizeof(double) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size (3D): " << (3.0 * ARRAY_SIZE_3D * sizeof(double) * 1.0E-6) << "MB"
              << " (" << (3.0 * ARRAY_SIZE_3D * sizeof(double) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout << std::endl;

    std::cout.precision(ss);

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
        auto minmax = std::minmax_element(matar_kokkos_timings_3D[ker].begin() + 1,
                                          matar_kokkos_timings_3D[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(matar_kokkos_timings_3D[ker].begin() + 1,
                                         matar_kokkos_timings_3D[ker].end(),
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
        auto minmax = std::minmax_element(kokkos_views_timings_3D[ker].begin() + 1,
                                          kokkos_views_timings_3D[ker].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(kokkos_views_timings_3D[ker].begin() + 1,
                                         kokkos_views_timings_3D[ker].end(),
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

//#ifdef MMM
    // Perform matrix matrix multiply benchmark
    int matrix_size = 64 * 64;
    int matrix_total_size = (matrix_size * matrix_size);
    auto mat1 = CArrayKokkos <real_t> (matrix_size, matrix_size);
    auto mat2 = CArrayKokkos <real_t> (matrix_size, matrix_size);
    auto mat3 = CArrayKokkos <real_t> (matrix_size, matrix_size);

    RMatrix2D kv_mat1("kv_mat1", matrix_size, matrix_size); 
    RMatrix2D kv_mat2("kv_mat2", matrix_size, matrix_size);
    RMatrix2D kv_mat3("kv_mat3", matrix_size, matrix_size);

    std::cout << std::endl; 
    std::cout << "Running 2D matrix-matrix multiplication (MMM) benchmark: " 
              << repeat << " times"
              << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << "Number of 2D array elements: " << matrix_total_size 
              << " (" << matrix_size << " elements in each dimension)"
              << std::endl;
    std::cout << std::endl;

    std::cout << "2D information" << std::endl;
    std::cout << "--------------" << std::endl;
    std::cout << std::setprecision(1) << std::fixed
              << "2D array size: " << (matrix_total_size* sizeof(double) * 1.0E-6) << "MB"
              << " (" << (matrix_total_size * sizeof(double) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size (2D): " << (3.0 * matrix_total_size * sizeof(double) * 1.0E-6) << "MB"
              << " (" << (3.0 * matrix_total_size * sizeof(double) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout << std::endl;

    std::cout.precision(ss);

    policy2D mmm_type = policy2D({0,0}, {matrix_size,matrix_size});
    Kokkos::parallel_for("MatrixInit", mmm_type, KOKKOS_LAMBDA(const int i, const int j) {
            mat1(i, j) = (real_t) (i + 1) * (j + 1);
            mat2(i, j) = (real_t) (i + 2) * (j + 2);

            // Initialize matrices
            kv_mat1(i, j) = (real_t) (i + 1) * (j + 1);
            kv_mat2(i, j) = (real_t) (i + 2) * (j + 2);

            //printf("Mat1 (%d, %d) %lf\n", i, j, mat1(i, j));
            //printf("Mat2 (%d, %d) %lf\n", i, j, mat2(i, j));
        });
    Kokkos::fence();

    std::vector<double> mmm_timings;

    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for ("MMM", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
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

        mmm_timings.push_back(mmm_time.count());
    }

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

    // std::cout << "MMM time: " << mmm_time.count() << " s." << std::endl;
    
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "2D CArrayKokkos MMM benchmark results" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MMM min time: " << *minmax_mmm.first << " sec" << std::endl;
    std::cout << "MMM max time: " << *minmax_mmm.second << " sec" << std::endl;
    std::cout << "MMM avg time: " << average_mmm << " sec" << std::endl;
    std::cout << std::endl;

    // Kokkos View matrix-matrix multiplication benchmark
    std::vector<double> mmm_kv_timings;

    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_for ("MMM (KV)", TeamPolicy(matrix_size, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
                const int i = teamMember.league_rank();

                Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, matrix_size), [=] (const int j) {
                    double temp_var = 0.0;

                    Kokkos::parallel_reduce (Kokkos::ThreadVectorRange (teamMember, matrix_size), [=] (const int k, double &mat_val) {
                        mat_val += kv_mat1(i, k) * kv_mat2(k, j);
                    }, temp_var);

                    kv_mat3(i, j) = temp_var;
                    //printf("Mat3 (%d, %d) %lf\n", i, j, mat3(i, j));
                });
            });
    
        Kokkos::fence();
    
        std::chrono::duration<double> mmm_time_kv = std::chrono::high_resolution_clock::now() - begin;

        mmm_kv_timings.push_back(mmm_time_kv.count());
    }

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
    std::cout << "MMM min time: " << *minmax_mmm_kv.first << " sec" << std::endl;
    std::cout << "MMM max time: " << *minmax_mmm_kv.second << " sec" << std::endl;
    std::cout << "MMM avg time: " << average_mmm_kv << " sec" << std::endl;
    std::cout << std::endl;
//#endif
    
    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
