#include <iostream>
#include <chrono>         // To access timing calipers 
#include "matar.h"

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

	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~GPU  TESTS~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

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
