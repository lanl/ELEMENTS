#include <iostream>
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
	const int size1 = 12000000;
	const int size3 = 256;
	const int repeat = 101;

	std::cout<<"Size of 1D problem: "<<size1<<"\n";
	std::cout<<"Size of 3D problem (each dimension): "<<size3<<"\n";

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

    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << "Kokkos test" << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    // Number of entries for 1D types (by default, it is 2e25)
    size_t ARRAY_SIZE_1D = 33554432;
    double scalar = 3.0;
    size_t nsize = 64 * 64 * 64 * 64;

    std::streamsize ss = std::cout.precision();

    std::cout << "Running kernels " << repeat << " times" << std::endl;
    std::cout << "Precision: double" << std::endl;
    std::cout << std::endl;

    // Conversion factor between megabytes (MB) and bytes: 1 byte = 10^(-6) MB
    // Conversion factor between bytes and gigabytes (GB): 1 byte = 10^(-9) GB
    std::cout << std::setprecision(1) << std::fixed
              << "Array size: " << (nsize * sizeof(double) * 1.0E-6) << "MB"
              << " (=" << (nsize * sizeof(double) * 1.0E-9) << " GB)"
              << std::endl;

    std::cout << "Total size: " << (3.0 * nsize * sizeof(double) * 1.0E-6) << "MB"
              << " (=" << (3.0 * nsize * sizeof(double) * 1.0E-9) << "GB)"
              << std::endl;

    std::cout.precision(ss);

    // Create array that will store timings for kernel runs
    // Rows correspond to kernels:
    //      1. Copy kernel
    //      2. Scale kernel
    //      3. Sum kernel
    //      4. Triad kernel
    //      5. Dot product kernel
    
    int num_kernels = 5;

    // auto matar_kokkos_timings = CArray<double> (5, repeat);
    std::vector<std::vector<double>> matar_kokkos_timings(num_kernels);
    std::vector<std::vector<double>> kokkos_views_timings(num_kernels);

    std::string labels[num_kernels] = {"Copy", "Mul", "Add", "Triad", "Dot"};

    size_t sizes[num_kernels] = {
        2 * sizeof(real_t) * nsize,
        2 * sizeof(real_t) * nsize,
        3 * sizeof(real_t) * nsize,
        3 * sizeof(real_t) * nsize,
        2 * sizeof(real_t) * nsize
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

    auto arr1 = CArrayKokkos <real_t> (nsize);
    auto arr2 = CArrayKokkos <real_t> (nsize);
    auto arr3 = CArrayKokkos <real_t> (nsize);

    // printf("Stream benchmark with %d elements.\n", nsize);

    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
                arr1(i) = 1.0;
                arr2(i) = 2.0;
                arr3(i) = 0.0;
                });
    Kokkos::fence();

    // Perform copy kernel 
    for (int iter = 0; iter < repeat; iter++) {
        auto begin = std::chrono::high_resolution_clock::now();
         
        Kokkos::parallel_for("Copy", nsize, KOKKOS_LAMBDA(const int i) {
                arr3(i) = arr1(i);
                });
        Kokkos::fence();
        

        std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;

        // Record copy kernel timing
        // matar_kokkos_timings(0, iter) = copy_time.count();
        matar_kokkos_timings[0].push_back(copy_time.count());

        // std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;
    }

    // Perform scaling kernel
    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Scale", nsize, KOKKOS_LAMBDA(const int i) {
                arr2(i) = scalar * arr3(i);
                });
        Kokkos::fence();
        
        std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;

        // Record scale kernel timing
        // matar_kokkos_timings(1, iter) = scale_time.count();
        matar_kokkos_timings[1].push_back(scale_time.count());

        // std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;
    }

    // Perform sum kernel
    PERF_START(tag) // tag is a string, e.g., "copy"
    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Sum", nsize, KOKKOS_LAMBDA(const int i) {
                arr3(i) = arr1(i) + arr2(i);
                });
        Kokkos::fence();
        
        std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin;

        // Record sum kernel timing
        // matar_kokkos_timings(2, iter) = sum_time.count();
        matar_kokkos_timings[2].push_back(sum_time.count()); 

        // std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;
    }

    // Perform triad kernel
    for (int iter = 0; iter < repeat; iter++) {
        begin = std::chrono::high_resolution_clock::now();
        
        Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
                arr1(i) = arr2(i) + 
                                      (scalar * arr3(i));
                });
        Kokkos::fence();
        
        std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

        // Record triad kernel timing
        // matar_kokkos_timings(3, iter) = triad_time.count();
        matar_kokkos_timings[3].push_back(triad_time.count());

        // std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;
    }

    // Perform dot product kernel
    for (int iter = 0; iter < repeat; iter++) {
        real_t sum = 0;

        begin = std::chrono::high_resolution_clock::now();

        Kokkos::parallel_reduce("Dot product", nsize, KOKKOS_LAMBDA(const int i, real_t& tmp) {
                tmp += arr3(i) * arr3(i);
        }, sum);
        Kokkos::fence();

        std::chrono::duration<double> dot_time = std::chrono::high_resolution_clock::now() - begin;

        // Record dot product kernel timing
        // matar_kokkos_timings(4, iter) = dot_time.count();
        matar_kokkos_timings[4].push_back(dot_time.count());
    }

    // Print kernel computation memory bandwidth table header
    std::cout << std::left << std::setw(12) << "Kernel"
              << std::left << std::setw(12) << "MBytes/sec"
              << std::left << std::setw(12) << "Min (sec)"
              << std::left << std::setw(12) << "Max (sec"
              << std::left << std::setw(12) << "Average (sec)"
              << std::endl
              << std::fixed;

    // Calculate kernel memory bandwidths
    for (int ker = 0; ker < num_kernels; ker++) {
        // Get min/max times taken on kernel computation
        // (ignore the first result)
        auto minmax = std::minmax_element(matar_kokkos_timings[i].begin() + 1,
                                          matar_kokkos_timings[i].end());

        // Calculate average time taken on kernel computation
        // (ignore the first result)
        double average = std::accumulate(matar_kokkos_timings[i].begin() + 1,
                                         matar_kokkos_timings[i].end(),
                                         0.0) / (double) (repeat - 1);

        // Print kernel computation memory bandwidth statistics
        std::cout << std::left << std::setw(12) << labels[i]
                  << std::left << ((1.0E-6 * sizes[i]) / (*minmax.first))
                  << std::left << std::setw(12) << std::setprecision(5) 
                                                << *minmax.first
                  << std::left << std::setw(12) << *minmax.second
                  << std::left << std::setw(12) << std::setprecision(5)
                                                << average
                  << std::endl;
    }
    
    std::cout << std::endl;

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
