#include <iostream>
#include <chrono>
#include <stdio.h> 

#include "matar.h"

//---likwid preprocessors---
#ifdef LIKWID_PERFMON
#include "likwid.h"
#else
#define LIKWID_MARKER_INIT;
#define LIKWID_MARKER_START(regionTag);
#define LIKWID_MARKER_STOP(regionTag);
#define LIKWID_MARKER_THREADINIT;
#define LIKWID_MARKER_REGISTER(regionTag);
#define LIKWID_MARKER_CLOSE;
#endif


int main()
{

	LIKWID_MARKER_INIT;
	
	//size of the matrix
	const size_t size2 = 3200;

	//========================================================================
	//create regular and carray mat
	//========================================================================
	auto cmat1 = CArray<double>(size2, size2);
	auto cmat2 = CArray<double>(size2, size2);
	auto cmat3 = CArray<double>(size2, size2);

	double** rmat1 = new double*[size2];
	double** rmat2 = new double*[size2];
	double** rmat3 = new double*[size2];

	for(size_t i = 0; i < size2; i++) {
	 rmat1[i] = new double[size2];
	 rmat2[i] = new double[size2];
	 rmat3[i] = new double[size2];
	}

	//========================================================================
	//initialize
	//========================================================================
#pragma omp simd collapse(2)
	for(size_t i = 0; i < size2; i++) {
	 for(size_t j = 0; j < size2; j++) {
		cmat1(i,j) = 1.0;
		cmat2(i,j) = 1.0;
		rmat1[i][j] = 1.0;
		rmat2[i][j] = 1.0;
	  }
	 }

	//========================================================================
	//	mat mult
	//========================================================================	


	//================================
	//	carray
	//===============================
	auto startc = std::chrono::steady_clock::now();
	LIKWID_MARKER_START("MATMULT_CARRAY");
	for(size_t i = 0; i < size2; i++) {
	 for(size_t j = 0; j < size2; j++) {
	  for(size_t k = 0; k < size2; k++) {
		cmat3(i,j) += cmat1(i,k) * cmat2(k,j);
		}
	   }
	  }
	LIKWID_MARKER_STOP("MATMULT_CARRAY");
	auto endc = std::chrono::steady_clock::now();
	std::chrono::duration<double> totc = endc - startc;
	std::cout<<"Elapsed time for matmult with carray " <<totc.count() <<"s\n";

	double rand1 = cmat3(1,1) + cmat3(3,50);
	std::cout<<"rand1 = "<<rand1<<"\n";

	//================================
	//	regular array
	//===============================
	auto startr = std::chrono::steady_clock::now();
	LIKWID_MARKER_START("MATMULT_REGULAR");
	for(size_t i = 0; i < size2; i++) {
	 for(size_t j = 0; j < size2; j++) {
	  for(size_t k = 0; k < size2; k++) {
		rmat3[i][j] += rmat1[i][k] * rmat2[k][j];
		}
	   }
	  }
	LIKWID_MARKER_STOP("MATMULT_REGULAR");
	auto endr = std::chrono::steady_clock::now();
	std::chrono::duration<double> totr = endr - startr;
	std::cout<<"Elapsed time for matmult with regular " <<totr.count() <<"s\n";
	
	double rand2 = rmat3[1][1] + rmat3[4][56];
	std::cout<<"rand2 = "<<rand2<<"\n";

	//begin Kokkos!
	Kokkos::initialize();
	{
		//setting the indices for matrices and arrays
		using policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
		policy2D array_type = policy2D({0,0},{size2,size2});
		policy2D matrix_type = policy2D({1,1}, {size2+1, size2+1});

		auto ckmat1 = CArrayKokkos<double>(size2, size2);
		auto ckmat2 = CArrayKokkos<double>(size2, size2);
		auto ckmat3 = CArrayKokkos<double>(size2, size2);
		
		//initialize mat1 and mat2 array
		Kokkos::parallel_for("InitailizeData", array_type, KOKKOS_LAMBDA (const int i, const int j) {
			ckmat1(i,j) = 1.0;
		 	ckmat2(i,j) = 1.0;
		});


		//matmult
		auto start_matmult = std::chrono::steady_clock::now();
		LIKWID_MARKER_START("2D_MATMULT_CK");
		for(size_t i = 0; i < size2; i++) {
		  for(size_t j = 0; j < size2; j++) {
			double temp_var = 0.0;
		    Kokkos::parallel_reduce("MatMult", size2, KOKKOS_LAMBDA(const int k, double &mat_val) {
				mat_val += ckmat1(i,k) * ckmat2(k,j);
			}, temp_var);
			ckmat3(i,j) = temp_var;
		  }
		 }
		LIKWID_MARKER_STOP("2D_MATMULT_CK");
		auto end_matmult = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_matmult = end_matmult - start_matmult;
		printf("Elapsed time for matmult with CArrayKokkos = %lf\n", tot_matmult.count() );
		std::cout<<"Elapsed time for matmult with CArrayKokkos = "<<tot_matmult.count()<<"s\n";

		//assign a random element from resulting matrix so compiler doesn't 
		//ignore loop
		double rand1 = rmat3[10][1];

		//outputi just to see if ouput is correct
		Kokkos::parallel_for("printVals", 10, KOKKOS_LAMBDA(const int i) {
		  printf("mat3.cak(1,%d) = %lf\n", i, rmat3[1][i] );
		});


	} //end kokkos
	Kokkos::finalize();

	
	LIKWID_MARKER_CLOSE;



} //---end---
