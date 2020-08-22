
#include <bits/stdc++.h>
#include <stdio.h> //for rand

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

int main() 
{

	std::cout<<"=============================================================\n";

	LIKWID_MARKER_INIT;

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~CREATE AND INIT ARRAYS~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	//included different repeat
	 int repeat = 101;	

	const int size2 = 3200;
	const int N = repeat;
//	const int nn = repeat1;

	//create matrices
	double** A = new double*[size2];
	double** B = new double*[size2];
	double** C = new double*[size2];

	//create vector for matrix-vector product
	double* x = new double[size2];
	double* b = new double[size2];

	for(size_t i = 0; i < size2; i++) {
		A[i] = new double[size2];
		B[i] = new double[size2];
		C[i] = new double[size2];
	}

	// CArrays
	auto AC = CArray<double>(size2, size2);
	auto BC = CArray<double>(size2, size2);
	auto CC = CArray<double>(size2, size2);
	auto DC = CArray<double>(size2, size2);

	auto xc = CArray<double>(size2);	//vector for Ax = b
	auto bc = CArray<double>(size2);	//holds solution for Ax = b

	//FArray
	auto BF = FArray<double>(size2, size2);
	auto xf = FArray<double>(size2);	//vector for Ax = b
	auto bf = FArray<double>(size2);	//holds solution for Ax = b for the fcase

	//initialize 
	std::cout<<"Initialize matrices!\n";
#pragma omp simd collapse(2)
	for(size_t i = 0; i < size2; i++ ) {
	 for(size_t j = 0; j < size2; j++) {
		A[i][j] = 1.0;
		B[i][j] = 1.0;
		AC(i,j) = 1.0;
		BC(i,j) = 1.0;
		BF(i,j) = 1.0;
	  }
	 }
	
	std::cout<<"Initialize vectors!\n";
#pragma omp simd
	for(size_t i = 0; i < size2; i++) {
		xc(i) = 1.0;
		xf(i) = 1.0;
		x[i] = 1.0;
	}

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~MAT-MAT MULT REGULAR TESTS~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	//1. MATRIX mult 
	std::cout<<"Starting regular matmult\n";
	double reg_mm[repeat];
	LIKWID_MARKER_START("MATMULT_REG");
	for(size_t t = 0; t < repeat; t++) {
	 std::cout<<"Iteration = "<<t<<"\n";
	 auto start_regmm = std::chrono::steady_clock::now();
	  for(size_t i = 0; i < size2; i++ ){
	   for(size_t j = 0; j < size2; j++) {
	    for(size_t k = 0; k < size2; k++) {
			C[i][j] += A[i][k] * B[k][j];
		  }
		 }
		}
	   auto end_regmm = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_regmm = end_regmm - start_regmm;
	   reg_mm[t] = tot_regmm.count();
	}
	LIKWID_MARKER_STOP("MATMULT_REG");
	std::cout<<"Finished regular matmult\n";
	//compute averages, min, max
	double rand_elem_r = C[1][1];
	double rand_elem_r2 = C[50][101];
	std::cout<<"rand_elem_r = "<<rand_elem_r<<"\n";
	std::cout<<"rand_elem_r2 = "<<rand_elem_r2<<"\n";

	double avg_rmm, min_rmm, max_rmm, tot_rmm;
	for(size_t i = 1; i < repeat; i++) {
		tot_rmm += reg_mm[i];
	}
	avg_rmm = tot_rmm / (double) (repeat - 1);
	min_rmm = *std::min_element(reg_mm + 1, reg_mm + repeat);
	max_rmm = *std::max_element(reg_mm + 1, reg_mm + repeat);

	std::cout<<"Average run time for 2D regular mat-mat mult = "<<avg_rmm<<"s\n";
	std::cout<<"Min run time for 2D regular mat-mat mult = "<<min_rmm<<"s\n";
	std::cout<<"Max run time for 2D regular mat-mat mult = "<<max_rmm<<"s\n";

/*

	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~MAT-MAT MULT CArray TESTS~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	std::cout<<"Starting CArray matmult\n";
	double c_mm[repeat];
	LIKWID_MARKER_START("MATMULT_CARR");
	for(size_t t = 0; t < repeat; t++) {
	 auto start_cmm = std::chrono::steady_clock::now();
	  for(size_t i = 0; i < size2; i++ ){
	   for(size_t j = 0; j < size2; j++) {
	    for(size_t k = 0; k < size2; k++) {
			CC(i,j) += AC(i,k) * BC(k,j);
		  }
		 }
		}
	   auto end_cmm = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_cmm = end_cmm - start_cmm;
	   c_mm[t] = tot_cmm.count();
	}
	LIKWID_MARKER_STOP("MATMULT_CARR");
	std::cout<<"FInished carray matmult\n";
	//compute averages, min, max
	double rand_elem_c = CC(1,1);
	std::cout<<"rand_elem_c = "<<rand_elem_c<<"\n";

	double avg_cmm, min_cmm, max_cmm, tot_cmm;
	for(size_t i = 1; i < repeat; i++) {
		tot_cmm += c_mm[i];
	}
	avg_cmm = tot_cmm / (double) (repeat - 1);
	min_cmm = *std::min_element(c_mm + 1, c_mm + repeat);
	max_cmm = *std::max_element(c_mm + 1, c_mm + repeat);

	std::cout<<"Average run time for 2D CArray mat-mat mult = "<<avg_cmm<<"s\n";
	std::cout<<"Min run time for 2D CArray mat-mat mult = "<<min_cmm<<"s\n";
	std::cout<<"Max run time for 2D CArray mat-mat mult = "<<max_cmm<<"s\n";


	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~MAT-MAT MULT C&F Array TESTS~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	std::cout<<"Starting c&f matmult\n";
	double cf_mm[repeat];
	LIKWID_MARKER_START("MATMULT_C_F");
	for(size_t t = 0; t < repeat; t++) {
	 std::cout<<"Iteration "<<t<<"\n";
	 auto start_cfmm = std::chrono::steady_clock::now();
	  for(size_t i = 0; i < size2; i++ ){
	   for(size_t j = 0; j < size2; j++) {
	    for(size_t k = 0; k < size2; k++) {
			DC(i,j) += AC(i,k) * BF(k,j);
		  }
		 }
		}
	   auto end_cfmm = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_cfmm = end_cfmm - start_cfmm;
	   cf_mm[t] = tot_cfmm.count();
	}
	LIKWID_MARKER_STOP("MATMULT_C_F");
	std::cout<<"Finished c&f matmult\n";
	//compute averages, min, max
	double rand_elem_cf = DC(1,1);
	std::cout<<"rand_elem_cf = "<<rand_elem_cf<<"\n";

	double avg_cfmm, min_cfmm, max_cfmm, tot_cfmm;
	for(size_t i = 1; i < repeat; i++) {
		tot_cfmm += cf_mm[i];
	}
	avg_cfmm = tot_cfmm / (double) (repeat - 1);
	min_cfmm = *std::min_element(cf_mm + 1, cf_mm + N);
	max_cfmm = *std::max_element(cf_mm + 1, cf_mm + N);

	std::cout<<"Average run time for 2D C&F mat-mat mult = "<<avg_cfmm<<"s\n";
	std::cout<<"Min run time for 2D C&F mat-mat mult = "<<min_cfmm<<"s\n";
	std::cout<<"Max run time for 2D C&F mat-mat mult = "<<max_cfmm<<"s\n";

*/

	//percent difference table 

	std::cout<<"=============================================================\n";
//	std::cout<<"Test	CArray		C&F\n";
	std::cout<<"=============================================================\n";

	// matrix-vector product
	// Ax = b
	// where 
	// (1) A and x are traditional c++ arrays
	// (2) A and x are Matar CArrays
	// (3) A is a Carray and x is an farray
	// using the same A matrix as above (all entries are 1 for simplicity)

 /*
	
	std::cout<<"=============================================================\n";
	std::cout<<"~~~~~~~~~~~~~~~~~Matrix-vector TESTS~~~~~~~~~~~~~~~~~\n";
	std::cout<<"=============================================================\n";

	std::cout<<"begin regular mv test\n";
	//regular arrays
	double reg_mv[repeat];
	LIKWID_MARKER_START("REG_MATVEC");
	for(size_t t = 0; t < repeat; t++) {
	  auto start_reg_mv = std::chrono::steady_clock::now();
	  for(size_t i = 0; i < size2; i++) {
	   for(size_t j = 0; j < size2; j++) {
		b[i] += A[i][j]*x[j];
		}
	   }
	   auto end_reg_mv = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_reg_mv = end_reg_mv - start_reg_mv;
	   reg_mv[t] = tot_reg_mv.count();
	}
	LIKWID_MARKER_STOP("REG_MATVEC");

	//compute average, min, max
	double avg_reg_mv, min_reg_mv, max_reg_mv, tot_reg_mvtime;
	for(size_t i = 1; i < repeat; i++) {
		tot_reg_mvtime += reg_mv[i];
	}
	avg_reg_mv = tot_reg_mvtime / (double)(repeat - 1);
	min_reg_mv = *std::min_element( reg_mv +1, reg_mv + repeat);
	max_reg_mv = *std::max_element( reg_mv +1, reg_mv + repeat);
	

	// carray test
	std::cout<<"Begin CArray MV test\n";
	double c_mv[repeat];
	LIKWID_MARKER_START("C_MATVEC");
	for(size_t t = 0; t < repeat; t++) {
	  auto start_c_mv = std::chrono::steady_clock::now();
	  for(size_t i = 0; i < size2; i++) {
	   for(size_t j = 0; j < size2; j++) {
		bc(i) += AC(i,j)*xc(j);
		}
	   }
	   auto end_c_mv = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_c_mv = end_c_mv - start_c_mv;
	  c_mv[t] = tot_c_mv.count();
	}
	LIKWID_MARKER_STOP("C_MATVEC");

	//compute average, min, max
	double avg_c_mv, min_c_mv, max_c_mv, tot_c_mvtime;
	for(size_t i = 1; i < repeat; i++) {
		tot_c_mvtime += c_mv[i];
	}
	avg_c_mv = tot_c_mvtime / (double)(repeat - 1);
	min_c_mv = *std::min_element( c_mv +1, c_mv + repeat);
	max_c_mv = *std::max_element( c_mv +1, c_mv + repeat);


	// cxf test
	std::cout<<"Begin C&F MV Test\n";
	double cf_mv[repeat];
	LIKWID_MARKER_START("C_F_MATVEC");
	for(size_t t = 0; t < repeat; t++) {
	  auto start_cf_mv = std::chrono::steady_clock::now();
	  for(size_t i = 0; i < size2; i++) {
	   for(size_t j = 0; j < size2; j++) {
		bf(i) += AC(i,j)*xf(i);
		}
	   }
	   auto end_cf_mv = std::chrono::steady_clock::now();
	   std::chrono::duration<double> tot_cf_mv = end_cf_mv - start_cf_mv;
	   cf_mv[t] = tot_cf_mv.count();
	}
	LIKWID_MARKER_STOP("C_F_MATVEC");

	//compute average, min, max
	double avg_cf_mv, min_cf_mv, max_cf_mv, tot_cf_mvtime;
	for(size_t i = 1; i < repeat; i++) {
		tot_cf_mvtime += cf_mv[i];
	}
	avg_cf_mv = tot_cf_mvtime / (double)(repeat - 1);
	min_cf_mv = *std::min_element( cf_mv +1, cf_mv + repeat);
	max_cf_mv = *std::max_element( cf_mv +1, cf_mv + repeat);

	std::cout<<"=============================================================\n";
	std::cout<<"			Matrix-Vector Product Times\n";
	std::cout<<"=============================================================\n";
	std::cout<<" Type		Avg.(s) 		Min.(s) 		Max.(s)\n";
	std::cout<<" Reg		"<<avg_reg_mv<<"	"<<min_reg_mv<<"	"<<max_reg_mv<<"\n";
	std::cout<<" CArray		"<<avg_c_mv<<"	"<<min_c_mv<<"	"<<max_c_mv<<"\n";
	std::cout<<" C&F array	"<<avg_cf_mv<<"	"<<min_cf_mv<<"	"<<max_cf_mv<<"\n";
	std::cout<<"=============================================================\n";

	


	std::cout<<"Some output of each vector\n";
	for(size_t i = 0; i < 5; i++) {
	  std::cout<<"b["<<i<<"] = "<<b[i]<<"\n";
	  std::cout<<"bc["<<i<<"] = "<<bc(i)<<"\n";
	  std::cout<<"bf["<<i<<"] = "<<bf(i)<<"\n";
	  std::cout<<"--------------------------------\n";
	}
 */


	LIKWID_MARKER_CLOSE;

}















