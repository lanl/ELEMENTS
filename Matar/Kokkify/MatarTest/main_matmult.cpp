#include <iostream>
#include <chrono>

#include "pseudo_mesh_cpu_benchmark.hpp"

int main()
{

	//size of the matrix
	size_t size2 = 3200;

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
	
		//initialize mat1 and mat2
		Kokkos::parallel_for("InitailizeData", size2, KOKKOS_LAMBDA(const int i, const int j) {
			mat1.cak(i,j) = 1.0;
			mat2.cak(i,j) = 1.0;
		});
		


	} //end kokkos
	Kokkos::finalize();

	




} //---end---
