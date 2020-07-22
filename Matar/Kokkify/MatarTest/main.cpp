#include <iostream>
#include <chrono>         // To access timing calipers 
#include "pseudo_mesh.hpp"

int main() {

    int size_i = 5, size_j = 2, size_k = 3;

    // CPU test
    
    printf("~~~~~~~~~~~~~~~Daniel's Example :)~~~~~~~~~~~~~~~~~\n");

    ///*
    auto og3d = CArray<int> (size_i, size_j, size_k); // 3 dimensional original

    for (int i = 0; i < size_i; i++) {
        for (int j = 0; j < size_j; j++) {
            for (int k = 0; k < size_k; k++) {
                //og3d[i][j][k] = k + j * size_k + i * size_j * size_k;
                og3d(i, j, k) = k + j * size_k + i * size_j * size_k;
            }
        }
    }

    auto view2d_0 = ViewCArray<int> (&og3d(0, 0, 0), size_j, size_k); // 2D slice at first index
    auto view2d_1 = ViewCArray<int> (&og3d(1, 0, 0), size_j, size_k); // 2D slice at second index

    // Both print statements should have OG be the same as new_
    for (int j = 0; j < size_j; j++) {
        for (int k = 0; k < size_k; k++) {
            printf("OG %d -> new_0 %d\n", og3d(0, j, k), view2d_0(j, k));
        }
    }

    for (int j = 0; j < size_j; j++) {
        for (int k = 0; k < size_k; k++) {
            printf("OG %d -> new_1 %d\n", og3d(1, j, k), view2d_1(j, k));
        }
    }
    //*/

    
    /*
    pseudo_mesh pmesh;
    pmesh.init(size_i, size_j);
    for (int i = 0; i < size_i; i++) {
        for (int j = 0; j < size_j; j++) {
            pmesh.var(i, j) = (double) j + size_j * i;
            printf("%d %d %lf\n", i, j, pmesh.var(i, j));
        }
    }
    */

/*
    //create a copy of test carray
    auto testCopy = c_array_t <double>(test);

    for(int j = 0; j < size_j; j++) {
     for(int i = 0; i < size_i; i++) {
	printf("testCopy(%d,%d) = %lf\n",i,j,testCopy(i,j));
	}
     }
*/

    /*
    printf("~~~~~~~~~~~~~~~~~Begin Test for CMatrix class!~~~~~~~~~~~~~~~~~!\n");

    //initialize the FArrays
    int size1 = 5; int size2 = 4; int size3 = 2; //just for the loops, should of initialized the arrays like this 
	
	//creating cmatrices of different dimension
	auto my2dcm = CMatrix <int> (4,4); //2d
	auto my3dcm = CMatrix <int> (2,2,2); //3d
	auto my4dcm = CMatrix <int> (2,2,2,2); //4d
	auto my5dcm = CMatrix <int> (2,2,2,2,2); //5d
	auto my6dcm = CMatrix <int> (2,2,2,2,2,2); //6d 

	//initialize 2d matrix
    	printf("Initialize and print 2D CMatrix!\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	for(int j = 1; j < 5; j++) {
		for(int i = 1; i < 5; i++) {
		  my2dcm(i,j) = (i-1)*size2 + (j-1);
		  printf("row %d col %d = %d\n",i,j,my2dcm(i,j));
		}
	   }

	printf("try to be out of bounds my2dcm(5,2) = %d\n",my2dcm(5,2));
*/
/*
	//initialize 3D matrix
    	printf("Initialize and print 3D CMatrix!\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	for(int i = 1; i < 3; i++) {
	  for(int j = 1; j < 3; j++){
	    for(int k = 1; k < 3; k++) {
		my3dcm(i,j,k) = (k-1) + (j-1)*2 + (i-1)*2*2;
		printf("i = %d, j = %d, k = %d, val = %d\n", i, j, k, my3dcm(i,j,k));
		}
	     }
	  }

	printf("Initialize and print 4D CMatrix\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	for(int i =1; i < size3+1; i++) {
	 for(int j = 1; j < size3+1; j++){
	  for(int k = 1; k < size3+1; k++) {
	    for(int l = 1; l < size3+1; l++) {
		my4dcm(i,j,k,l) = (l-1) + (k-1)*size3 + (j-1)*size3*size3 + (i-1)*size3*size3*size3;
		printf("my4dcm(%d,%d,%d,%d) = %d\n",i,j,k,l,my4dcm(i,j,k,l));
		}
	      }
	    }
	   }
	
	printf("Initialize and print 5D CMatrix\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	for(int i = 1; i < size3+1; i++) {
	 for(int j = 1; j < size3+1; j++) {
	  for(int k = 1; k < size3+1; k++) {
	   for(int l = 1; l < size3+1; l++) {
	    for(int m = 1; m < size3+1; m++) {
	     	my5dcm(i,j,k,l,m) = (m-1) + (l-1)*size3 + (k-1)*size3*size3 + (j-1)*size3*size3*size3 + (i-1)*size3*size3*size3*size3;
		printf("my5dcm(%d,%d,%d,%d,%d) = %d\n",i,j,k,l,m,my5dcm(i,j,k,l,m));
		}
	       }
	      }
	     }
	    }

	printf("Initialize and print 6D CMatrix\n");
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	for(int i = 1; i < size3+1; i++) {
	 for(int j = 1; j < size3+1;j++) {
	  for(int k = 1; k<size3+1; k++) {
	   for(int l = 1; l<size3+1; l++){
	    for(int m = 1; m<size3+1;m++){
	     for(int n=1 ;n<size3+1; n++){
	     	my5dcm(i,j,k,l,m,n) = (n-1)+ (m-1)*size3 + (l-1)*size3*size3 + (k-1)*size3*size3*size3 + (j-1)*size3*size3*size3*size3 + (i-1)*size3*size3*size3*size3*size3;
		printf("my6dcm(%d,%d,%d,%d,%d,%d) = %d\n",i,j,k,l,m,n,my5dcm(i,j,k,l,m,n));
		}
	       }
	      }
	     }
	    }
	  }

    printf("~~~~~~~~~~~~~~~~~Test FArray Class~~~~~~~~~~~~~~\n");

    auto farr1d = FArray<int>(5); //1d array
    auto farr2d = FArray<int>(4,4); //2d array
    auto farr3d = FArray<int>(2,2,2); //3d array 
    auto farr4d = FArray<int>(2,2,2,2); //4d array
    auto farr5d = FArray<int>(2,2,2,2,2); //5d array
    auto farr6d = FArray<int>(2,2,2,2,2,2); //6d array


    //1d
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Initializing and printing 1D FArray elements!\n");
    for(int i = 0; i< size1; i++) {
     farr1d(i) = i;
     printf("farr1d(%d) = %d\n",i,farr1d(i));
    }

    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Initializing and printing 2D FAray elements!\n");
    for(int i = 0; i < size2; i++) {
     for(int j = 0; j < size2; j++) {
	farr2d(i,j) = j*size2 + i;
	printf("farr2d(%d,%d) = %d\n",i,j,farr2d(i,j));
     }
    }
   
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Initializing and printing 3D FArray elements!\n");
    for(int i = 0; i < size3; i++) {
     for(int j = 0; j < size3; j++) {
      for(int k = 0; k < size3; k++) {
	farr3d(i,j,k) = i + j*size3 + k*size3*size3;
	printf("far3d(%d,%d,%d) = %d\n",i,j,k,farr3d(i,j,k));
	}
       }
      }

    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Initializing and printing 4D FArray elements!\n");
    for(int i = 0; i < size3; i++) {
     for(int j = 0; j < size3; j++) {
      for(int k = 0; k < size3; k++) {
	for(int l = 0; l < size3; l++) {
	  farr4d(i,j,k,l) = i + j*size3 + k*size3*size3 + l*size3*size3*size3;
	  printf("farr4d(%d,%d,%d,%d) = %d\n", i,j,k,l,farr4d(i,j,k,l));
	  }
	 }
	}
      }

    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Initializing and printing 5D FArray elements!\n");
    for(int i = 0; i < size3; i++) {
     for(int j = 0; j < size3; j++) {
      for(int k = 0; k < size3; k++) {
	for(int l = 0; l < size3; l++) {
	  for(int m = 0; m < size3; m++) {
	    farr5d(i,j,k,l,m) = i + j*size3 + k*size3*size3 + l*size3*size3*size3 + m*size3*size3*size3*size3;
	    printf("farr5d(%d,%d,%d,%d,%d) = %d\n", i,j,k,l,m,farr5d(i,j,k,l,m));
	   }
	  }
	 }
	}
      }

    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Initializing and printing 6D FArray elements!\n");
    for(int i = 0; i < size3; i++) {
     for(int j = 0; j < size3; j++) {
      for(int k = 0; k < size3; k++) {
	for(int l = 0; l < size3; l++) {
	  for(int m = 0; m < size3; m++) {
	    for(int n = 0; n < size3; n++) {
		  farr6d(i,j,k,l,m,n) = i + j*size3 + k*size3*size3 + l*size3*size3*size3 + m*size3*size3*size3*size3 + n*size3*size3*size3*size3*size3;
		   printf("farr6d(%d,%d,%d,%d,%d,%d) = %d\n",i,j,k,l,m,n,farr6d(i,j,k,l,m,n));
		}
	    }
	  }
	 }
	}
      }

    printf("~~~~~~~~~~~~END OF TESTING FARRAY and CMATRIX~~~~~~~~~~~~~\n");
   
*/
    
    printf("\nTesting CArray and ViewCArray classes\n");

    // Create and test 1D CArray
    printf("1D case\n");

    size_t size_i_2 = 100;

    auto c_array_1d = CArray<double> (size_i_2);

    for (size_t i = 0; i < size_i_2; i++) {
        c_array_1d(i) = i;
        // printf("c_array_1d(%d) = %f\n", i, c_array_1d(i));
    }
    
    printf("\n");

    // Create and test multidimensional ViewCArray for 1D CArray
    auto view_c_array_md = ViewCArray<double> (&c_array_1d(0), 5, 10, 2);

    for (size_t i = 0; i < 5; i++) {
        for (size_t j = 0; j < 10; j++) {
            for (size_t k = 0; k < 2; k++) {
                view_c_array_md(i, j, k) = (2 * view_c_array_md(i, j, k));
            }
        }
    }

    for (size_t i = 0; i < size_i_2; i++) {
        // printf("c_array_1d(%d) = %f\n", c_array_1d(i));
    }

    printf("\n");

    // Create and test ViewCArray of a multidimensional ViewCArray
    // For example, modify the entries in the 5th row and 10th column
    // of the two arrays in the multidimensional view
    auto view_c_array_1d = ViewCArray<double> (&view_c_array_md(4, 9, 0), 2);

    for (size_t k = 0; k < 2; k++) {
        view_c_array_1d(k) = k;
    }

    for (size_t i = 0; i < size_i_2; i++) {
        // printf("c_array_1d(%d) = %f\n", c_array_1d(i));
    }
    
    printf("\n");

    // Create and test 2D CArray
    size_t size_j_2 = 55;

    auto c_array_2d = CArray<double> (size_i_2, size_j_2);

    for (size_t i = 0; i < size_i_2; i++) {
        for (size_t j = 0; j < size_j_2; j++) {
            c_array_2d(i, j) = j + (i * size_j_2);
            // printf("c_array_2d(%d, %d) = %f\n", i, j, c_array_2d(i, j));
        }
    }

    printf("\n");

    // Create and and test a ViewCArray object that treats a 2D CArray object
    // as if it were 1D
    auto view_c_array_1d_2 = ViewCArray<double> (&c_array_2d(0, 0), 
                                                 (size_i_2 * size_j_2));

    for (size_t i = 0; i < (size_i_2 * size_j_2); i++) {
        view_c_array_1d_2(i) = (2 * view_c_array_1d_2(i));
    }

    for (size_t i = 0; i < size_i_2; i++) {
        for (size_t j = 0; j < size_j_2; j++) {
            // printf("c_array_2d(%d, %d) = %f\n", i, j, c_array_2d(i, j));
        }
    }

    printf("\n");

    // Create and test 3D CArray
    size_t size_k_2 = 20;

    auto c_array_3d = CArray<double> (size_i_2, size_j_2, size_k_2);

    for (size_t i = 0; i < size_i_2; i++) {
        for (size_t j = 0; j < size_j_2; j++) {
            for (size_t k = 0; k < size_k_2; k++) {
                c_array_3d(i, j, k) = (k + (j * size_k_2)
                                         + (i * size_k_2 * size_j_2));
             }
        }
    }

    // Create and test a ViewCArray object that treats a 3D CArray object
    // as 

    printf("\n---------------------------------------\n");

    printf("\nTesting FMatrix and ViewFMatrix classes\n");
    
    size_t size_l = 5;
    size_t size_m = 5;
    size_t size_n = 5;

    // Create and test 1D FMatrix
    auto f_matrix_1d = FMatrix<double>(size_i);

    printf("\n");
    /*
    for (size_t i = 1; i < (size_i + 1); i++) {
        f_matrix_1d(i) = i;
        std::cout << "f_matrix_1d(" << i << ") = " << f_matrix_1d(i)
                  << std::endl;
    }
    */
    printf("\n");

    // Test view for 1D FMatrix
    auto view_f_matrix_1d = ViewFMatrix<double>(&f_matrix_1d(1), size_i);
    /*
    for (size_t i = 1; i < (size_i + 1); i++) {
        printf("f_matrix_1d(%d) %f  == view_f_matrix_1d(%d) %f\n",
               i, f_matrix_1d(i), i, view_f_matrix_1d(i));
    } 
    */ 
    auto f_matrix_2d = FMatrix<double>(size_i, size_j);
    /*
    for(size_t j = 1; j < (size_j + 1); j++) {
        for (size_t i = 1; i < (size_i + 1); i++) {
            f_matrix_2d(i, j) = i + ((j - 1) * size_i);
            printf("f_matrix_2d(%d, %d) = %f\n", i, j, 
                    f_matrix_2d(i, j));
        }
    }
    */ 
    printf("\n");

    // Test view for 2D FMatrix
    // First, slice in i-direction, i.e., grab first column
    auto view_f_matrix_2d_first_col = ViewFMatrix<double>(&f_matrix_2d(1, 1), 
                                                          size_i);
    /* 
    for (size_t i = 1; i < (size_i + 1); i++) {
        printf("f_matrix_2d(%d, 1) %f == view_f_matrix_2d_first_col(%d) %f\n",
                i, f_matrix_2d(i, 1), i, view_f_matrix_2d_first_col(i));
    }
    */ 
    printf("\n");

    // Then, slice in j-direction, i.e., grab first row
    auto view_f_matrix_2d_first_row = ViewFMatrix<double>(&f_matrix_2d(1, 1),
                                                          size_j);
    /* 
    for (size_t j = 1; j < (size_j + 1); j++) {
        printf("f_matrix_2d(1, %d) %f == view_f_matrix_2d_first_row(%d) %f\n",
                j, f_matrix_2d(1, j), j, view_f_matrix_2d_first_row(j));
    }
    */ 
    printf("\n");

    auto f_matrix_3d = FMatrix<double>(size_i, size_j, size_k);
    /* 
    for (size_t k = 1; k < (size_k + 1); k++) {
        for (size_t j = 1; j < (size_j + 1); j++) {
            for (size_t i = 1; i < (size_i + 1); i++) {
                f_matrix_3d(i, j, k) = i + ((j - 1) * size_i)
                                         + ((k - 1) * size_i * size_j);
                printf("f_matrix_3d(%d, %d, %d) = %f\n",
                        i, j, k, f_matrix_3d(i, j, k));
            }
        }
    } 
    */
    printf("\n");

    // Test view for 3D FMatrix
    // First, slice in the i-direction, i.e., grab the first column of every
    // matrix in the "stack"
    auto view_f_matrix_3d_first_col_mat = ViewFMatrix<double>(&f_matrix_3d(1, 1, 1),
                                                              size_j, size_k);
    /*
    for (size_t k = 1; k < (size_k + 1); k++) {
        for (size_t i = 1; i < (size_i + 1); i++) {
            printf("f_matrix_3d(%d, 1, %d) %f == view_f_matrix_3d_first_col_mat(%d, 1, %d) %f\n",
                    i, k, f_matrix_3d(i, 1, k), i, k, 
                    view_f_matrix_3d_first_col_mat(i, 1, k));
        }
    }
    */

    auto f_matrix_4d = FMatrix<double>(size_i, size_j, size_k, size_l);
    /*
    for (size_t l = 1; l < (size_l + 1); l++) {
        for (size_t k = 1; k < (size_k + 1); k++) {
            for (size_t j = 1; j < (size_j + 1); j++) {
                for (size_t i = 1; i < (size_i + 1); i++) {
                    f_matrix_4d(i, j, k, l) = i + ((j - 1) * size_i)
                                                + ((k - 1) * size_i * size_j)
                                                + ((l - 1) * size_i _ size_j * size_k);
                    std::cout << "f_matrix_4d(" << i << ", " << j << ", " 
                              << k << ", " << l << ") = " 
                              << f_matrix_4d(i, j, k, l) 
                              << std::endl;
                }
            }
        }
    }
    */

    std::cout << std::endl;

    auto f_matrix_5d = FMatrix<double>(size_i, size_j, size_k, size_l,
                                       size_m);
    auto f_matrix_6d = FMatrix<double>(size_i, size_j, size_k, size_l,
                                       size_m, size_n);


    printf("~~~~~~~~~~~~~~~~~GPU TEST~~~~~~~~~~~~~~~\n"); 

    // Kokkos GPU test
   
    Kokkos::initialize();
    {

    using policy2D = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;
    policy2D array_type = policy2D({0,0}, {size_i, size_j});
    policy2D matrix_type = policy2D({1,1}, {size_i+1, size_j+1});
    policy2D view_type = policy2D({0,0}, {size_j, size_k});
    using policy3D = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
    policy3D array_type3 = policy3D({0,0,0}, {size_i, size_j, size_k});
    policy3D matrix_type3 = policy3D({1,1,1}, {size_i+1, size_j+1, size_k+1});
    

printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Daniel's ragged example\n"); 
    pseudo_mesh * pmesh = (pseudo_mesh *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(sizeof(pseudo_mesh));
    // Template of stream benchmark code
    // size_t STREAM_ARRAY_SIZE = 100; // Dummy value; could be whatever we want
    // pseudo_mesh_variant* pmesh_variant = ... // Create CArrayKokkos a, b, c

    Kokkos::parallel_for("PseudoMesh", 1, KOKKOS_LAMBDA(const int&) {
        pmesh->init(size_i, size_j);
        });

    //using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;
    Kokkos::parallel_for ("PseudoMesh", TeamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
        //pmesh->init(size_i, size_j);
        });


   // Kokkos::parallel_for("Values", array_type3, KOKKOS_LAMBDA(const int i, const int j, const int k) {
   //     pmesh->var(i,j) = (double) (k + j*size_k + i*size_k*size_j);
   //     });
    Kokkos::parallel_for ("RaggedDownTeam", TeamPolicy(size_j, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamPolicy::member_type &teamMember) {
            const int j = teamMember.league_rank();
            const int j_stride = pmesh->var.stride(j);
            Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, j_stride), [=] (const int i) {
                pmesh->var(i, j) = (double) j + i * j_stride; // not the exact placement, needs start index
                printf("%d %d %lf\n", i, j, pmesh->var(i, j));
            });
        });
    Kokkos::fence();





printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Daniel's large/stream example\n"); 
    int big_i = 5000, big_j = 2000;
    k_test_t <real_t> kview;
    kview = k_test_t <real_t> (big_i, big_j);
    Kokkos::parallel_for("Values", mdrange_policy2({0,0}, {big_i,big_j}), KOKKOS_LAMBDA(const int i, const int j) {
            kview(i,j) = (real_t) j + i * big_j;
            // Only print for small examples, e.g. 5 x 2, not these. Takes a long time to print
            if (i < size_i && j < size_j)printf("%d) %lf\n", j + i * big_j, kview(i,j));
        });



    //********************This is just an example (untested) for a view_c type.******************************
    //Kokkos::parallel_for("ValueCheck", view_type, KOKKOS_LAMBDA(const int j, const int k) {
    //    auto var_view2d = view_c_array <double> (&pmesh->var(1, 0, 0), size_j, size_k);
    //    printf("%lf ====== %lf\n", pmesh->var(1,j,k), var_view2d(j,k));
    //    });
    
    /*
    Kokkos::parallel_for("ValueCheck", array_type, KOKKOS_LAMBDA(const int i, const int j) {
        auto var_view2d = ViewFArrayKokkos<double> (&pmesh->var(0, 0, 1), size_i, size_j);
        printf("%lf ====== %lf\n", pmesh->var(i,j,1), var_view2d(i,j));
        });
    */
/* 
    Kokkos::parallel_for("ValCM",matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
	pmesh->var1(i,j) = (j-1) + (i-1)*size_i;
	printf("pmeshvar1CM(%d,%d) = %d\n",i,j,pmesh->var1(i,j));
	});   

    //blocking out the copy fornow
    //create another 2D mesh 
    pseudo_mesh * my2dmesh = (pseudo_mesh *)Kokkos::kokkos_malloc<Kokkos::CudaSpace>(sizeof(pseudo_mesh));

    

    //print values from 2dmesh
    Kokkos::parallel_for("printCopyMesh",matrix_type, KOKKOS_LAMBDA(const int i, const int j) {
    printf("my2dmesh(%d,%d) = %d\n", i,j,my2dmesh->var1(i,j));
    } );

	//create a mesh to test 3D FArray
    	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    	printf("Creating 3D FArray Mesh\n");
	pseudo_mesh * my3dmesh = (pseudo_mesh *)Kokkos::kokkos_malloc<Kokkos::CudaSpace>(sizeof(pseudo_mesh));

	Kokkos::parallel_for("Init3dMesh",1,KOKKOS_LAMBDA(const int&) {
	my3dmesh->init(2,2,2);
	});

	Kokkos::parallel_for("SetVal3d", Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{2,2,2}), KOKKOS_LAMBDA(const int i, const int j, const int k) {
	  my3dmesh->var(i,j,k) = (double)(i + j*2 + k*2*2);
	  printf("my3dmesh(%d,%d,%d) = %lf\n", i,j,k,my3dmesh->var(i,j,k));
	  } ) ;
	
	Kokkos::parallel_for("SetVal3dCM",Kokkos::MDRangePolicy<Kokkos::Rank<3>>({1,1,1},{3,3,3}),KOKKOS_LAMBDA(const int i, const int j, const int k){
	  my3dmesh->var1(i,j,k) = (k-1) + (j-1)*2 + (i-1)*2*2;
	  printf("my3dmeshVar1(%d,%d,%d) = %d\n",i,j,k,my3dmesh->var1(i,j,k));
	});
  
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Creating 4D FArray Mesh\n");
    pseudo_mesh * my4dmesh = (pseudo_mesh *)Kokkos::kokkos_malloc<Kokkos::CudaSpace>(sizeof(pseudo_mesh));

    Kokkos::parallel_for("Init4DMesh",1,KOKKOS_LAMBDA(const int&) {
    my4dmesh->init(2,2,2,2);
    });

    Kokkos::parallel_for("SetVal4D",Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0,0,0,0},{2,2,2,2}),KOKKOS_LAMBDA(const int i, const int j, const int k, const int l){
    my4dmesh->var(i,j,k,l) =  (double)(i + j*2 + k*2*2 + l*2*2*2);
    printf("my4dmesh(%d,%d,%d,%d) = %lf\n",i,j,k,l,my4dmesh->var(i,j,k,l));
    });

    Kokkos::parallel_for("SetVal4dC",Kokkos::MDRangePolicy<Kokkos::Rank<4>>({1,1,1,1},{3,3,3,3}),KOKKOS_LAMBDA(const int i, const int j, const int k, const int l){
	my4dmesh->var1(i,j,k,l) = (l-1) + (k-1)*2 + (j-1)*2*2 + (i-1)*2*2*2;
	printf("my4dmeshCM(%d,%d,%d,%d) = %d\n",i,j,k,l,my4dmesh->var1(i,j,k,l));
	});

 
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Creating 5D FArray Mesh\n");
    pseudo_mesh * my5dmesh = (pseudo_mesh*)Kokkos::kokkos_malloc<Kokkos::CudaSpace>(sizeof(pseudo_mesh));

    Kokkos::parallel_for("Init5DMesh",1, KOKKOS_LAMBDA(const int&){
    my5dmesh->init(2,2,2,2,2);
    });

    Kokkos::parallel_for("Set5DVal",Kokkos::MDRangePolicy<Kokkos::Rank<5>>({0,0,0,0,0},{2,2,2,2,2}),KOKKOS_LAMBDA(const int i, const int j, const int k, const int l, const int m) {
    my5dmesh->var(i,j,k,l,m) = (double)(i + j*2 + k*2*2 + l*2*2*2 + m*2*2*2*2);
    printf("my5dmesh(%d,%d,%d,%d,%d) = %f\n", i,j,k,l,m,my5dmesh->var(i,j,k,l,m));
    });

    Kokkos::parallel_for("set5dcm",Kokkos::MDRangePolicy<Kokkos::Rank<5>>({1,1,1,1,1},{3,3,3,3,3}),KOKKOS_LAMBDA(const int i, const int j, const int k, const int l, const int m){
     my5dmesh->var1(i,j,k,l,m) = (m-1) + (l-1)*2 + (k-1)*2*2 + (j-1)*2*2*2 + (i-1)*2*2*2*2;
     printf("mymesh5dvar1(%d,%d,%d,%d,%d) = %d\n",i,j,k,l,m,my5dmesh->var1(i,j,k,l,m));
     });

    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Creating 6D FArray Mesh\n");
    pseudo_mesh *my6dmesh = (pseudo_mesh*)Kokkos::kokkos_malloc<Kokkos::CudaSpace>(sizeof(pseudo_mesh));

    Kokkos::parallel_for("Init5dMesh",1,KOKKOS_LAMBDA(const int&) {
    my6dmesh->init(2,2,2,2,2,2);
    });

    Kokkos::parallel_for("Set6dVals",Kokkos::MDRangePolicy<Kokkos::Rank<6>>({0,0,0,0,0,0},{2,2,2,2,2,2}),KOKKOS_LAMBDA(const int i, const int j, const int k, const int l, const int m, const int n) {
    my6dmesh->var(i,j,k,l,m,n) = (double)(i + j*2 + k*2*2 + l*2*2*2 + m*2*2*2*2 + n*2*2*2*2*2);
    printf("my6dmesh(%d,%d,%d,%d,%d,%d) = %lf\n",i,j,k,l,m,n,my6dmesh->var(i,j,k,l,m,n));
    });

    Kokkos::parallel_for("Set6dCM",Kokkos::MDRangePolicy<Kokkos::Rank<6>>({1,1,1,1,1,1},{3,3,3,3,3,3}),KOKKOS_LAMBDA(const int i, const int j, const int k, const int l, const int m, const int n){
    my6dmesh->var1(i,j,k,l,m,n) = (n-1) + (m-1)*2 + (l-1)*2*2 + (k-1)*2*2*2 + (j-1)*2*2*2*2 + (i-1)*2*2*2*2*2;
    printf("mymesh6dVar1(%d,%d,%d,%d,%d,%d) = %d\n",i,j,k,l,m,n,my6dmesh->var1(i,j,k,l,m,n));
     });

 */
	} //end kokkos
    Kokkos::finalize();

 	
#ifdef PSEUDOVAR
    // Kokkos stream benchmark test (for FMatrixKokkos and CArrayKokkos objects)
    printf("\n Kokkos stream benchmark (FMatrixKokkos and CArrayKokkos)\n");
    
    double scalar = 3.0;
    size_t nsize = 10000000;
    pseudo_var* pvar = (pseudo_var*) Kokkos::kokkos_malloc<Kokkos::CudaSpace> 
                                             (sizeof(pseudo_var));

    // Initialize three 1D FMatrixKokkos and CArrayKokkos objects,
    // respectively
    Kokkos::parallel_for("StreamTriadInit", 1, KOKKOS_LAMBDA(const int&) {
            pvar->init(nsize);
            });

    printf("Stream benchmark with %d elements.\n", nsize);
    /*
    Kokkos::parallel_for("Initialize", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var1(i) = 1.0;
            pvar->c_arr_var2(i) = 2.0;
            pvar->c_arr_var3(i) = 0.0;
            });
    */
    // Perform copy benchmark
    auto begin = std::chrono::high_resolution_clock::now();
    /* 
    Kokkos::parallel_for("Copy", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var3(i) = pvar->c_arr_var1(i);
            });
    */

    std::chrono::duration<double> copy_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Copy time: " << copy_time.count() << " s." << std::endl;

    // Perform scaling benchmark
    begin = std::chrono::high_resolution_clock::now();
    /*
    Kokkos::parallel_for("Scale", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var2(i) = scalar * pvar->c_arr_var3(i);
            });
    */
    std::chrono::duration<double> scale_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Scale time: " << scale_time.count() << " s." << std::endl;

    // Perform sum benchmark
    begin = std::chrono::high_resolution_clock::now();
    /*
    Kokkos::parallel_for("Sum", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var3(i) = pvar->c_arr_var1(i) + pvar->c_arr_var2(i);
            });
    */
    std::chrono::duration<double> sum_time = std::chrono::high_resolution_clock::now() - begin; 

    std::cout << "Sum time: " << sum_time.count() << " s." << std::endl;

    // Perform triad benchmark
    begin = std::chrono::high_resolution_clock::now();
    /*
    Kokkos::parallel_for("Triad", nsize, KOKKOS_LAMBDA(const int i) {
            pvar->c_arr_var1(i) = pvar->c_arr_var2(i) + 
                                  (scalar * pvar->c_arr_var3(i));
            });
    */
    std::chrono::duration<double> triad_time = std::chrono::high_resolution_clock::now() - begin;

    std::cout << "Triad time: " << triad_time.count() << " s." << std::endl;
    


    // Test RaggedRight 
    /*
    // First, create RaggedRight using CArray
    auto c_array_strides_array = CArray<size_t>(size_i);
    for (size_t i = 0; i < size_i; i++) {
        c_array_strides_array(i) = i;
    }
    
    auto rr_array = RaggedRightArray<double>(c_array_strides_array);

    int val = 1;
    for (size_t i = 0; i < size_i; i++) {
        for (size_t j = 0; j < rr_array.stride(i); j++) {
            rr_array(i, j) = val;
            val++;
        }
    }
    
    for (size_t i = 0; i < size_i; i++) {
        for (size_t j = 0; j < rr_array.stride(i); j++) {
            printf("rr_array(%d, %d) = %f\n", i, j, rr_array(i, j));
        }
    }
    
    // Then, create RaggedRight using ViewCArray
    auto c_array_strides_array_2 = CArray<size_t>(size_i + 2);
    for (size_t i = 0; i < (size_i + 2); i++) {
        c_array_strides_array_2(i) = 2*i;
    }

    auto view_c_array_strides_array = ViewCArray<size_t>(&c_array_strides_array_2(0),
                                                          size_i + 2);

    auto rr_array_vs = RaggedRightArray<double>(view_c_array_strides_array);

    val = 1;
    for (size_t i = 0; i < (size_i + 2); i++) {
        for (size_t j = 0; j < rr_array_vs.stride(i); j++) {
            rr_array_vs(i, j) = val;
            val++;
        }
    }

    printf("\n");

    for (size_t i = 0; i < (size_i + 2); i++) {
        for (size_t j = 0; j < rr_array_vs.stride(i); j++) {
            printf("rr_array_vs(%d, %d) = %f\n", i, j, rr_array_vs(i, j));
        }
    }

    // Create RaggedRight using regular C array
    size_t num_strides = 4;
    size_t strides[num_strides] = {4, 5, 1, 2}; 
    
    auto rr_array_reg = RaggedRightArray<double>(strides, num_strides);

    val = 1;
    for (size_t i = 0; i < num_strides; i++) {
        for (size_t j = 0; j < rr_array_reg.stride(i); j++) {
            rr_array_reg(i, j) = val;
            val++;
        }
    }

    printf("\n");

    for (size_t i = 0; i < num_strides; i++) {
        for (size_t j = 0; j < rr_array_reg.stride(i); j++) {
            printf("rr_array_reg(%d, %d) = %f\n", i, j, rr_array_reg(i, j));
        }
    }
    */
#endif


    // Test CSR
    

    printf("--- finished ---\n");

    return 0;
}
