#include <iostream>
#include <chrono>
#include <stdio.h> 

#include "pseudo_mesh_cpu_benchmark.hpp"

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

	//size of the matrix
	const size_t size2 = 3200;

	LIKWID_MARKER_INIT;

	//begin Kokkos!
	Kokkos::initialize();
	{
		//setting the indices for matrices and arrays
		using policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
		policy2D array_type = policy2D({0,0},{size2,size2});
		policy2D matrix_type = policy2D({1,1}, {size2+1, size2+1});

		//creating 3 instaces of the pseudo_mesh...cpu class
		pmcb mat1;
		pmcb mat2;
		pmcb mat3;
	
		//create the matar-kokkos data
		mat1.init(size2, size2);
		mat2.init(size2, size2);
		mat3.init(size2, size2);
	
		//initialize mat1 and mat2 array
		Kokkos::parallel_for("InitailizeData", array_type, KOKKOS_LAMBDA (const int i, const int j) {
			mat1.cak(i,j) = 1.0;
			mat2.cak(i,j) = 1.0;
		});


		//matmult
		auto start_matmult = std::chrono::steady_clock::now();
		LIKWID_MARKER_START("2D_MATMULT");
		for(size_t i = 0; i < size2; i++) {
		  for(size_t j = 0; j < size2; j++) {
			double temp_var = 0.0;
		    Kokkos::parallel_reduce("MatMult", size2, KOKKOS_LAMBDA(const int k, double &mat_val) {
				mat_val += mat1.cak(i,k) * mat2.cak(k,j);
			}, temp_var);
			mat3.cak(i,j) = temp_var;
		  }
		 }
		LIKWID_MARKER_STOP("2D_MATMULT");
		auto end_matmult = std::chrono::steady_clock::now();
		std::chrono::duration<double> tot_matmult = end_matmult - start_matmult;
		printf("Elapsed time for matmult with CArrayKokkos = %ld\n", tot_matmult.count() );
		std::cout<<"Elapsed time for matmult with CArrayKokkos = "<<tot_matmult.count()<<"s\n";

		//assign a random element from resulting matrix so compiler doesn't 
		//ignore loop
		double rand1 = mat3.cak(10,1);

		//outputi just to see if ouput is correct
		Kokkos::parallel_for("printVals", 10, KOKKOS_LAMBDA(const int i) {
		  printf("mat3.cak(1,%d) = %lf\n", i, mat3.cak(1,i) );
		});


	} //end kokkos
	Kokkos::finalize();

	
	LIKWID_MARKER_CLOSE;



} //---end---
