#include <iostream>
#include <chrono>         // To access timing calipers 
#include <bits/stdc++.h>

#include "matar.h"

//-----likwid preprocessors----
#ifdef LIKWID_PERFMON
#include "likwid.h"
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#endif

int main() {

//	LIKWID_MARKER_INIT;

    int size_i = 5, size_j = 4, size_k = 3;
	const int size1 = 12000000;
	const int size3 = 256;
	const int repeat = 101;

	std::cout<<"Size of 1D problem: "<<size1<<"\n";
	std::cout<<"Size of 3D problem (each dimension): "<<size3<<"\n";

	//=========================================================================
	//	cpu stream benchmarks
	//=========================================================================	

	//1. create carrays

	//1D
	auto c_arr1 = CArray<double>(size1);
	auto c_arr2 = CArray<double>(size1);
	auto c_arr3 = CArray<double>(size1);
	auto c_arr4 = CArray<double>(size1);

	//3D
	auto c_arr1_3d = CArray<double>(size3, size3, size3);
	auto c_arr2_3d = CArray<double>(size3, size3, size3);
	auto c_arr3_3d = CArray<double>(size3, size3, size3);
	auto c_arr4_3d = CArray<double>(size3, size3, size3);

	//1D
	double* reg_arr1_1d = new double[size1];
	double* reg_arr2_1d = new double[size1];
	double* reg_arr3_1d = new double[size1];
	double* reg_arr4_1d = new double[size1];

	//3D
	double*** reg_arr1_3d = new double**[size3];
	double*** reg_arr2_3d = new double**[size3];
	double*** reg_arr3_3d = new double**[size3];
	double*** reg_arr4_3d = new double**[size3];

	for(size_t i = 0; i < size3; i++) {
	 reg_arr1_3d[i] = new double*[size3];	
	 reg_arr2_3d[i] = new double*[size3];	
	 reg_arr3_3d[i] = new double*[size3];	
	 reg_arr4_3d[i] = new double*[size3];	
	  for(size_t j = 0; j < size3; j++) {
		reg_arr1_3d[i][j] = new double[size3];
		reg_arr2_3d[i][j] = new double[size3];
		reg_arr3_3d[i][j] = new double[size3];
		reg_arr4_3d[i][j] = new double[size3];
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

	
	//=========================================================================
	//	1D CArray Tests
	//=========================================================================

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~~~~1D CArray TESTS~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";


	//=======================
	//	COPY
	//=======================	
	
	double carr_copy_times[repeat];	//array to hold times for carray copy
//	LIKWID_MARKER_START("1D_CARRAY_COPY");
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_c_1d = std::chrono::steady_clock::now();
#pragma omp simd	   
	   for(size_t i = 0; i < size1; i++) {
		 c_arr2(i) = c_arr1(i);
		}
		auto end_copy_c_1d = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_copy_c_1d = end_copy_c_1d - start_copy_c_1d;
		carr_copy_times[t] = tot_copy_c_1d.count();
	}
//	LIKWID_MARKER_START("1D_CARRAY_COPY");
	/* N will be used in the function to get min and max elements in the array*/	
	int N = sizeof(carr_copy_times) / sizeof(carr_copy_times[0]);


	//calculate average
	double avg_copy_c, tot_copy_c;
	for(size_t i = 0; i < repeat; i++) {
		tot_copy_c += carr_copy_times[i];
	}	
	avg_copy_c = tot_copy_c / (double) repeat;
	double min1dc_copy;
	double max1dc_copy;
	min1dc_copy = *std::min_element(carr_copy_times, carr_copy_times + N);
	max1dc_copy = *std::max_element(carr_copy_times, carr_copy_times + N);
	std::cout<<"Average elapsed time for 1D Carray copy = "<<avg_copy_c<<"s\n";
	std::cout<<"Min 1D CArray copy time = "<<min1dc_copy<<"s\n";
	std::cout<<"Max 1D CArray copy time = "<<max1dc_copy<<"s\n";

	//=======================
	//	SCALE
	//=======================	
	
	double c_scale1d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale1dc = std::chrono::steady_clock::now();
#pragma omp simd
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
	double min1dc_scale, max1dc_scale;
	min1dc_scale = *std::min_element(c_scale1d_times, c_scale1d_times + N);
	max1dc_scale = *std::max_element(c_scale1d_times, c_scale1d_times + N);
	std::cout<<"Averaged elapsed time for 1D Carray scale = "<<avg_scale1dc<<"s\n";
	std::cout<<"Min 1D CArray scale time = " <<min1dc_scale<<"\n";
	std::cout<<"Max 1D CArray scale time = " <<max1dc_scale<<"\n";

	//=======================
	//	SUM
	//=======================	

	double sum_1dc[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_sum1dc = std::chrono::steady_clock::now();
#pragma omp simd
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
	double min1dc_sum, max1dc_sum;
	min1dc_sum = *std::min_element(sum_1dc, sum_1dc + N);
	max1dc_sum = *std::max_element(sum_1dc, sum_1dc + N);
	std::cout<<"Averate elapsed time for 1D CArray sum = " <<avg_sum1dc<<"s\n";
	std::cout<<"Min 1D Carray sum time = "<<min1dc_sum<<"s\n";
	std::cout<<"Max 1D Carray sum time = "<<max1dc_sum<<"s\n";

	//=======================
	//	TRIAD
	//=======================	
	
	double triad_1dc[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_triad1dc = std::chrono::steady_clock::now();
#pragma omp simd
	   for(size_t i = 0; i < size1; i++) {
		c_arr4(i) = 2.0*c_arr1(i) + c_arr2(i);
		}
	   auto end_triad1dc = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_triad1dc = end_triad1dc - start_triad1dc;
	   triad_1dc[t] = tot_triad1dc.count();
	}

	//compute average
	double avg_triad_1dc, tot_1dctriad;
	for(size_t i = 0; i < repeat; i++) {
		tot_1dctriad += triad_1dc[i];
	}
	avg_triad_1dc = tot_1dctriad / (double) repeat;
	std::cout<<"Average elapsed time for 1D Carray triad = "<<avg_triad_1dc<<"s\n";
	double min1dc_triad, max1dc_triad;
	min1dc_triad = *std::min_element(triad_1dc, triad_1dc + N);
	max1dc_triad = *std::min_element(triad_1dc, triad_1dc + N);
	std::cout<<"Min 1D CArray triad time = "<<min1dc_triad<<"s\n";
	std::cout<<"Max 1D CArray triad time = "<<max1dc_triad<<"s\n";

	//=======================
	//	DOT PRODUCT
	//=======================	

	double dot_1dc[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   double dot_sum = 0.0;
	   auto start_1dcdot = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
		dot_sum += c_arr1(i)*c_arr1(i);
	   }
	   auto end_1dcdot = std::chrono::steady_clock::now();
	   std::chrono::duration<double>tot_1dcdot = end_1dcdot - start_1dcdot;
	   dot_1dc[t] = tot_1dcdot.count();
	}

	//compute averages, min, max
	double avg_1dcdot, tot_dot1dc;
	for(size_t i = 0; i < repeat; i++) {
		tot_dot1dc += dot_1dc[i];
	}
	avg_1dcdot = tot_dot1dc / (double) repeat;
	std::cout<<"Average elapsed time for 1D Carray dot product = "<<avg_1dcdot<<"s\n";
	double min1dc_dot, max1dc_dot;
	min1dc_dot = *std::min_element(dot_1dc, dot_1dc + N);
	max1dc_dot = *std::max_element(dot_1dc, dot_1dc + N);
	std::cout<<"Min 1D Carray dot product time = "<<min1dc_dot<<"s\n";
	std::cout<<"Max 1D Carray dot product time = "<<max1dc_dot<<"s\n";


	std::cout<<"=============================================================\n";


	//=========================================================================
	//	1D Regular Array Tests
	//=========================================================================

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~~~~1D Regular Array TESTS~~~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";


	//=======================
	//	COPY
	//=======================	

	double reg_copy_times[repeat];	//array to hold times for carray copy
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_r_1d = std::chrono::steady_clock::now();
#pragma omp simd
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
	double min1dr_copy, max1dr_copy;
	min1dr_copy = *std::min_element(reg_copy_times, reg_copy_times + N);
	max1dr_copy = *std::max_element(reg_copy_times, reg_copy_times + N);
	std::cout<<"Average elapsed time for 1D trad copy = "<<avg_copy_r<<"s\n";
	std::cout<<"Min 1D trad copy time = "<<min1dr_copy<<"\n";
	std::cout<<"Max 1D trad copy time = "<<max1dr_copy<<"\n";

	//=======================
	//	SCALE
	//=======================	

	double reg_scale1d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale1dr = std::chrono::steady_clock::now();
#pragma omp simd
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
	double min1dr_scale, max1dr_scale;
	min1dr_scale = *std::min_element(reg_scale1d_times, reg_scale1d_times + N);
	max1dr_scale = *std::max_element(reg_scale1d_times, reg_scale1d_times + N);
	std::cout<<"Averaged elapsed time for 1D trad scale = "<<avg_scale1dr<<"s\n";
	std::cout<<"Min 1D trad scale time = "<<min1dr_scale<<"s\n";
	std::cout<<"Max 1D trad scale time = "<<max1dr_scale<<"s\n";

	//=======================
	//	SUM
	//=======================	

	double sum_1dr[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_sum1dr = std::chrono::steady_clock::now();
#pragma omp simd
	   for(size_t i = 0; i < size1; i++) {
			reg_arr3_1d[i] = reg_arr1_1d[i] + reg_arr2_1d[i];
		}
	   auto end_sum1dr = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_1dr_sum = end_sum1dr - start_sum1dr;
	   sum_1dr[t] = tot_1dr_sum.count();
	}

	//compute average
	double avg_sum1dr, tot_sum_1dr;
	for(size_t i = 0; i < repeat; i++) {
		tot_sum_1dr += sum_1dr[i];
	}
	avg_sum1dr = tot_sum_1dr / (double) repeat;
	double min1dr_sum, max1dr_sum;
	min1dr_sum = *std::min_element(sum_1dr, sum_1dr + N);
	max1dr_sum = *std::min_element(sum_1dr, sum_1dr + N);
	std::cout<<"Averate elapsed time for 1D trad sum = " <<avg_sum1dr<<"s\n";
	std::cout<<"Min 1D trad sum time = "<<min1dr_sum<<"s\n";
	std::cout<<"Max 1D trad sum time = "<<max1dr_sum<<"s\n";

	//=======================
	//	TRIAD
	//=======================	

	double triad_1dr[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_triad1dr = std::chrono::steady_clock::now();
#pragma omp simd
	   for(size_t i = 0; i < size1; i++) {
		reg_arr4_1d[i] = 2.0*reg_arr1_1d[i] + reg_arr2_1d[i];
		}
	   auto end_triad1dr = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_triad1dr = end_triad1dr - start_triad1dr;
	   triad_1dr[t] = tot_triad1dr.count();
	}

	//compute average
	double avg_triad_1dr, tot_1drtriad;
	for(size_t i = 0; i < repeat; i++) {
		tot_1drtriad += triad_1dr[i];
	}
	avg_triad_1dr = tot_1drtriad / (double) repeat;
	std::cout<<"Average elapsed time for 1D trad triad = "<<avg_triad_1dr<<"s\n";
	double min1dr_triad, max1dr_triad;
	min1dr_triad = *std::min_element(triad_1dr, triad_1dr + N);
	max1dr_triad = *std::min_element(triad_1dr, triad_1dr + N);
	std::cout<<"Min 1D trad triad time = "<<min1dr_triad<<"s\n";
	std::cout<<"Max 1D trad triad time = "<<max1dr_triad<<"s\n";

	//=======================
	//	DOT PRODUCT
	//=======================	

	double dot_1dr[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   double dot_sum = 0.0;
	   auto start_1drdot = std::chrono::steady_clock::now();
	   for(size_t i = 0; i < size1; i++) {
		dot_sum += reg_arr1_1d[i]*reg_arr1_1d[i];
	   }
	   auto end_1drdot = std::chrono::steady_clock::now();
	   std::chrono::duration<double>tot_1drdot = end_1drdot - start_1drdot;
	   dot_1dr[t] = tot_1drdot.count();
	}

	//compute averages, min, max
	double avg_1drdot, tot_dot1dr;
	for(size_t i = 0; i < repeat; i++) {
		tot_dot1dr += dot_1dr[i];
	}
	avg_1drdot = tot_dot1dr / (double) repeat;
	std::cout<<"Average elapsed time for 1D trad dot product = "<<avg_1drdot<<"s\n";
	double min1dr_dot, max1dr_dot;
	min1dr_dot = *std::min_element(dot_1dr, dot_1dr + N);
	max1dr_dot = *std::max_element(dot_1dr, dot_1dr + N);
	std::cout<<"Min 1D trad dot product time = "<<min1dr_dot<<"s\n";
	std::cout<<"Max 1D trad dot product time = "<<max1dr_dot<<"s\n";


	std::cout<<"=============================================================\n";


	//=========================================================================
	//	3D CArray Tests
	//=========================================================================

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~3D CArray TESTS~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	//=============================
	//	COPY
	//=============================	

	double carr_copy_3dtimes[repeat]; 
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_c_3d = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
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
	double min3dc_copy, max3dc_copy;
	min3dc_copy = *std::min_element(carr_copy_3dtimes, carr_copy_3dtimes + N);
	max3dc_copy = *std::max_element(carr_copy_3dtimes, carr_copy_3dtimes + N);
	std::cout<<"Average elapsed time for 3D Carray copy = "<<avg_copy_3dc<<"s\n";
	std::cout<<"Min 3D CArray copy time = "<<min3dc_copy<<"\n";
	std::cout<<"Max 3D Carray copy time = "<<max3dc_copy<<"\n";

	//=============================
	//	SCALE
	//=============================	

	double c_scale3d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale3dc = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
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
	double min3dc_scale, max3dc_scale;
	min3dc_scale = *std::min_element(c_scale3d_times, c_scale3d_times + N);
	max3dc_scale = *std::max_element(c_scale3d_times, c_scale3d_times + N);
	std::cout<<"Average elapsed time for 3D Carray scale = "<<avg_scale3dc<<"s\n";
	std::cout<<"Min 3D Carray scale time = "<<min3dc_scale<<"\n";
	std::cout<<"Max 3D CArray scale time = "<<max3dc_scale<<"\n";
	
	//=============================
	//	SUM
	//=============================	

	double sum_3dc[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_sum3dc = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
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
	double min3dc_sum, max3dc_sum;
	min3dc_sum = *std::min_element(sum_3dc, sum_3dc + N);
	max3dc_sum = *std::max_element(sum_3dc, sum_3dc + N);
	std::cout<<"Averape elapsed time for 3D CArray sum = "<<avg_sum3dc<<"s\n";
	std::cout<<"Min 3D Carray sum time = "<<min3dc_sum<<"s\n";
	std::cout<<"Max 3D Carray sum time = "<<max3dc_sum<<"s\n";

	//=============================
	//	TRIAD
	//=============================	

	double triad_3dc[repeat];
	for(size_t t = 0; t < repeat; t++ ) {
	   auto start_triad3dc = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			c_arr4_3d(i,j,k) = 2.0*c_arr1_3d(i,j,k) + c_arr2_3d(i,j,k);
		   }
		  }
		 }
		auto end_triad3dc = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_triad3dc = end_triad3dc - start_triad3dc;
		triad_3dc[t] = tot_triad3dc.count();
	}
	
	//compute average
	double avg_triad_3dc, tot_3dctriad;
	for(size_t i = 0; i < repeat; i++) {
		tot_3dctriad += triad_3dc[i];
	}
	avg_triad_3dc = tot_3dctriad / (double) repeat;
	std::cout<<"Average elapsed time for 3D Carray triad = "<<avg_triad_3dc<<"s\n";
	double min3dc_triad, max3dc_triad;
	min3dc_triad = *std::min_element(triad_3dc, triad_3dc + N);
	max3dc_triad = *std::min_element(triad_3dc, triad_3dc + N);
	std::cout<<"Min 3D CArray triad time = "<<min3dc_triad<<"s\n";
	std::cout<<"Max 3D CArray triad time = "<<max3dc_triad<<"s\n";

	//=============================
	//	DOT PRODUCT
	//=============================	


	std::cout<<"=============================================================\n";


	//=========================================================================
	//	3D TRAD Tests
	//=========================================================================

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~~3D TRADITIONAL ARRAY TESTS~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	//=============================
	//	COPY
	//=============================	

	double reg_copy_3dtimes[repeat]; 
	for(size_t t = 0; t < repeat; t++) {
	   auto start_copy_r_3d = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
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
	double min3dr_copy, max3dr_copy;
	min3dr_copy = *std::min_element(reg_copy_3dtimes, reg_copy_3dtimes + N);
	max3dr_copy = *std::max_element(reg_copy_3dtimes, reg_copy_3dtimes + N);
	std::cout<<"Average elapsed time for 3D trad copy = "<<avg_copy_3dr<<"s\n";
	std::cout<<"Min 3D trad copy time = "<<min3dr_copy<<"\n";
	std::cout<<"Max 3D trad copy time = "<<max3dr_copy<<"\n";


	//=============================
	//	SCALE
	//=============================	

	double reg_scale3d_times[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_scale3dr = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
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
	double min3dr_scale, max3dr_scale;
	min3dr_scale = *std::min_element(reg_scale3d_times, reg_scale3d_times + N);
	max3dr_scale = *std::max_element(reg_scale3d_times, reg_scale3d_times + N);
	std::cout<<"Average elapsed time for 3D trad scale = "<<avg_scale3dr<<"s\n";
	std::cout<<"Min 3D trad scale time = "<<min3dr_scale<<"s\n";
	std::cout<<"Max 3D trad scale time = "<<max3dr_scale<<"s\n";


	//=============================
	//	SUM
	//=============================	

	double sum_3dr[repeat];
	for(size_t t = 0; t < repeat; t++) {
	   auto start_sum3dr = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			reg_arr3_3d[i][j][k] = reg_arr2_3d[i][j][k] + reg_arr1_3d[i][j][k];
		  }
		 }
		}
	   auto end_sum3dr = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_3dr_sum = end_sum3dr - start_sum3dr;
	   sum_3dr[t] = tot_3dr_sum.count();
	}

	//compute average
	double avg_sum3dr, tot_sum_3dr;
	for(size_t i = 0; i < repeat; i++) {
		tot_sum_3dr += sum_3dr[i];
	}
	avg_sum3dr = tot_sum_3dr / (double) repeat;
	std::cout<<"Averape elapsed time for 3D trad sum = "<<avg_sum3dr<<"s\n";
	double min3dr_sum, max3dr_sum;
	min3dr_sum = *std::min_element(sum_3dr, sum_3dr + N);
	max3dr_sum = *std::min_element(sum_3dr, sum_3dr + N);
	std::cout<<"Min 3D trad sum time = "<<min3dr_sum<<"s\n";
	std::cout<<"Max 3D trad sum time = "<<max3dr_sum<<"s\n";


	//=============================
	//	TRIAD
	//=============================	

	double triad_3dr[repeat];
	for(size_t t = 0; t < repeat; t++ ) {
	   auto start_triad3dr = std::chrono::steady_clock::now();
#pragma omp simd collapse(3)
	   for(size_t i = 0; i < size3; i++) {
		for(size_t j = 0; j < size3; j++) {
		 for(size_t k = 0; k < size3; k++) {
			reg_arr4_3d[i][j][k] = 2.0*reg_arr1_3d[i][j][k] + reg_arr2_3d[i][j][k];
		   }
		  }
		 }
		auto end_triad3dr = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_triad3dr = end_triad3dr - start_triad3dr;
		triad_3dr[t] = tot_triad3dr.count();
	}
	
	//compute average
	double avg_triad_3dr, tot_3drtriad;
	for(size_t i = 0; i < repeat; i++) {
		tot_3drtriad += triad_3dr[i];
	}
	avg_triad_3dr = tot_3drtriad / (double) repeat;
	std::cout<<"Average elapsed time for 3D trad triad = "<<avg_triad_3dr<<"s\n";
	double min3dr_triad, max3dr_triad;
	min3dr_triad = *std::min_element(triad_3dr, triad_3dr + N);
	max3dr_triad = *std::min_element(triad_3dr, triad_3dr + N);
	std::cout<<"Min 3D trad triad time = "<<min3dr_triad<<"s\n";
	std::cout<<"Max 3D trad triad time = "<<max3dr_triad<<"s\n";


	std::cout<<"=============================================================\n";
















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

//	LIKWID_MARKER_CLOSE;

    printf("--- finished ---\n");

    return 0;
}
