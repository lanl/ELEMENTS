
/* Stream Triad Test for performance analyis on MATAR data
 *  The Stream Triad will test copy, scale and triad.
 */

#include <stdio.h> 
#include <chrono>
#include <iostream>

#include "cpu_test_class.hpp"  //includes matar data

#define size1 4000000
#define size2 2000 
#define size3 150 
#define size4 45 
#define size5 21 
#define size6 13 
#define ntimes 20

#define ndbug //stop asserts

using namespace std;
using namespace std::chrono;


int main() 
{

    int scale = 2;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1D tests~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   printf("~~~~~~~~~1D Data~~~~~~~~~~\n");

    cpu_test_class test_1;
    test_1.init(size1);

    //initialize 1D
    for(int i = 0; i < size1; i++) {
	test_1.f_arr1(i) = 1;
	test_1.c_arr1(i) = 1;
	}

    for(int i = 1; i < size1+1; i++) {
	test_1.f_mat1(i) = 1;
	test_1.c_mat1(i) = 1;
	}


  // Test for FArray 1D 
    
    //test copy
    auto start_f_1d_copy = high_resolution_clock::now();
    test_1.test_copy_farr(ntimes, 1);
    auto stop_f_1d_copy = high_resolution_clock::now();
    auto tot_f_1d_copy = duration_cast<microseconds>(stop_f_1d_copy - start_f_1d_copy);
    cout<< "Average time for an FArray copy of size " << size1<< " is " << tot_f_1d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_f_1d_scale = high_resolution_clock::now();
    test_1.test_scale_farr(ntimes, 1, scale);
    auto end_f_1d_scale = high_resolution_clock::now();
    auto tot_f_1d_scale = duration_cast<microseconds>(end_f_1d_scale - start_f_1d_scale);
    cout<<"Average time for FArray scale of size " << size1 << " is " << tot_f_1d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for farray
    auto start_f_1d_triad = high_resolution_clock::now();
    test_1.test_triad_farr(ntimes, 1, scale);
    auto end_f_1d_triad = high_resolution_clock::now();
    auto tot_f_1d_triad = duration_cast<microseconds>(end_f_1d_triad - start_f_1d_triad);
    cout<<"Average time for FArray triad of size " << size1 << " is " << tot_f_1d_triad.count()/ntimes << " microseconds!"<<endl;

    printf("\n\n");

   //CArray 1D test
    auto start_c_1d_copy = high_resolution_clock::now();
    test_1.test_copy_carr(ntimes, 1);
    auto stop_c_1d_copy = high_resolution_clock::now();
    auto tot_c_1d_copy = duration_cast<microseconds>(stop_c_1d_copy - start_c_1d_copy);
    cout<< "Average time for an CArray copy of size " << size1<< " is " << tot_c_1d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_c_1d_scale = high_resolution_clock::now();
    test_1.test_scale_carr(ntimes, 1, scale);
    auto end_c_1d_scale = high_resolution_clock::now();
    auto tot_c_1d_scale = duration_cast<microseconds>(end_c_1d_scale - start_c_1d_scale);
    cout<<"Average time for CArray scale of size " << size1 << " is " << tot_c_1d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for carray
    auto start_c_1d_triad = high_resolution_clock::now();
    test_1.test_triad_carr(ntimes, 1, scale);
    auto end_c_1d_triad = high_resolution_clock::now();
    auto tot_c_1d_triad = duration_cast<microseconds>(end_c_1d_triad - start_c_1d_triad);
    cout<<"Average time for CArray triad of size " << size1 << " is " << tot_c_1d_triad.count()/ntimes << " microseconds!"<<endl;
  
    printf("\n\n");

    //CMatrix

    //copy
    auto start_cm_1d_copy = high_resolution_clock::now();
    test_1.test_copy_cmat(ntimes, 1);
    auto stop_cm_1d_copy = high_resolution_clock::now();
    auto tot_cm_1d_copy = duration_cast<microseconds>(stop_cm_1d_copy - start_cm_1d_copy);
    cout<< "Average time for an CMatrix copy of size " << size1<< " is " << tot_cm_1d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_cm_1d_scale = high_resolution_clock::now();
    test_1.test_scale_cmat(ntimes, 1, scale);
    auto end_cm_1d_scale = high_resolution_clock::now();
    auto tot_cm_1d_scale = duration_cast<microseconds>(end_cm_1d_scale - start_cm_1d_scale);
    cout<<"Average time for CMatrix scale of size " << size1 << " is " << tot_cm_1d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for cmatrix
    auto start_cm_1d_triad = high_resolution_clock::now();
    test_1.test_triad_cmat(ntimes, 1, scale);
    auto end_cm_1d_triad = high_resolution_clock::now();
    auto tot_cm_1d_triad = duration_cast<microseconds>(end_cm_1d_triad - start_cm_1d_triad);
    cout<<"Average time for CMatrix triad of size " << size1 << " is " << tot_cm_1d_triad.count()/ntimes << " microseconds!"<<endl;
  //end CMatrix 1D test

    printf("\n\n");

    //FMatrix
    
    //copy
    auto start_fm_1d_copy = high_resolution_clock::now();
    test_1.test_copy_fmat(ntimes, 1);
    auto stop_fm_1d_copy = high_resolution_clock::now();
    auto tot_fm_1d_copy = duration_cast<microseconds>(stop_fm_1d_copy - start_fm_1d_copy);
    cout<< "Average time for an FMatrix copy of size " << size1<< " is " << tot_fm_1d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_fm_1d_scale = high_resolution_clock::now();

    //test scale
    test_1.test_scale_fmat(ntimes, 1, scale);

    auto end_fm_1d_scale = high_resolution_clock::now();
    auto tot_fm_1d_scale = duration_cast<microseconds>(end_fm_1d_scale - start_fm_1d_scale);
    cout<<"Average time for FMatrix scale of size " << size1 << " is " << tot_fm_1d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for fmatrix
    auto start_fm_1d_triad = high_resolution_clock::now();
    test_1.test_triad_fmat(ntimes, 1, scale);
    auto end_fm_1d_triad = high_resolution_clock::now();
    auto tot_fm_1d_triad = duration_cast<microseconds>(end_fm_1d_triad - start_fm_1d_triad);
    cout<<"Average time for FMatrix triad of size " << size1 << " is " << tot_fm_1d_triad.count()/ntimes << " microseconds!"<<endl;


//~~~~~~~~~~~~~~~~~~~~~end of 1D tests~~~~~~~~~~~~~~~~~~~~~



//~~~~~~~~~create objects for 2D - 6D

 /*
  int size2 = 2000; //4mil
  int size3 = 150; //3.75 mil
  int size4 = 45; 
  int size5 = 21; 
  int size6 = 13; 
 */
   
  //1. create objects
  //the number is associated with dimension
  cpu_test_class test_2;
  test_2.init(size2, size2);

  cpu_test_class test_3;
  test_3.init(size3, size3, size3);

  cpu_test_class test_4;
  test_4.init(size4, size4, size4, size4);

  cpu_test_class test_5;
  test_5.init(size5, size5, size5, size5, size5);

  cpu_test_class test_6;
  test_6.init(size6, size6, size6, size6, size6, size6);

  //1. initialize arrays and matrices

  //2D
  for(int i = 0; i < size2; i++) {
  // printf("in row %d \n", i);
   for(int j = 0; j < size2; j++) {
	test_2.f_arr1(i,j) = 1;
	test_2.c_arr1(i,j) = 1;
	}
      }

  for(int i = 1; i < size2+1; i++) {
   for(int j = 1; j < size2+1; j++) {
	test_2.f_mat1(i,j) = 1;
	test_2.c_mat1(i,j) = 1;
	}
       }

  
  //3D
  
	for(int i = 0; i < size3; i++) {
		//printf("i = %d\n", i);
	 for(int j = 0; j < size3; j++) {
	  for(int k = 0; k < size3; k++) {
		test_3.f_arr1(i,j,k) = 1;
		test_3.c_arr1(i,j,k) = 1;
		}
	    }
	   }
 
	for(int i = 1; i < size3+1; i++) {
	 for(int j = 1; j < size3+1; j++) {
	  for(int k = 1; k < size3+1; k++) {
		test_3.f_mat1(i,j,k) = 1;
		test_3.c_mat1(i,j,k) = 1;
		}
	      }
	    }

  //4D
	for(int i = 0; i < size4; i++) {
	 for(int j = 0; j < size4; j++) {
	  for(int k = 0; k < size4; k++) {
	   for(int l = 0; l < size4; l++) {
		test_4.f_arr1(i,j,k,l) = 1;
		test_4.c_arr1(i,j,k,l) = 1;
		}
	       }
	      }
	     }
  
	for(int i = 1; i < size4+1; i++) {
	 for(int j = 1; j < size4+1; j++) {
	  for(int k = 1; k < size4+1; k++) {
	   for(int l = 1; l < size4+1; l++) {
		test_4.f_mat1(i,j,k,l) = 1;
		test_4.c_mat1(i,j,k,l) = 1;
		}
	       }
	      }
	     }

  //5D
	for(int i = 0; i < size5; i++) {
	 for(int j = 0; j < size5; j++) {
	  for(int k = 0; k < size5; k++) {
	   for(int l = 0; l < size5; l++) {
	     for(int m = 0; m < size5; m++) {
		test_5.f_arr1(i,j,k,l,m) = 1;
		test_5.c_arr1(i,j,k,l,m) = 1;
		}
	       }
	      }
	     }
	    }

	for(int i = 1; i < size5+1; i++) {
	 for(int j = 1; j < size5+1; j++) {
	  for(int k = 1; k < size5+1; k++) {
	   for(int l = 1; l < size5+1; l++) {
	     for(int m = 1; m < size5+1; m++) {
		test_5.f_mat1(i,j,k,l,m) = 1;
		test_5.c_mat1(i,j,k,l,m) = 1;
		}
	       }
	      }
	     }
	    }

  //6D 
	for(int i = 0; i < size6; i++) {
	 for(int j = 0; j < size6; j++) {
	  for(int k = 0; k < size6; k++) {
	   for(int l = 0; l < size6; l++) {
	    for(int m = 0; m < size6; m++) {
	     for(int n = 0; n < size6; n++) {
		test_6.f_arr1(i,j,k,l,m,n) = 1;
		test_6.c_arr1(i,j,k,l,m,n) = 1;
		}
	       }
	      }
	     }
	    }
	   }

	for(int i = 1; i < size6+1; i++) {
	 for(int j = 1; j < size6+1; j++) {
	  for(int k = 1; k < size6+1; k++) {
	   for(int l = 1; l < size6+1; l++) {
	    for(int m = 1; m < size6+1; m++) {
	     for(int n = 1; n < size6+1; n++) {
		test_6.f_mat1(i,j,k,l,m,n) = 1;
		test_6.c_mat1(i,j,k,l,m,n) = 1;
		}
	       }
	      }
	     }
	    }
	   }

   //~~~~time tests for 2D data~~~~~~~
  printf("~~~~~~~2D DATA~~~~~~~~~~\n");
  //FArray
  //copy
  auto start_f_2d_copy = high_resolution_clock::now();
  test_2.test_copy_farr(ntimes, 2);
  auto end_f_2d_copy = high_resolution_clock::now();
  auto tot_f_2d_copy = duration_cast<microseconds>(end_f_2d_copy - start_f_2d_copy);
  cout<<"Average time for 2D FArray copy of size " << size2 << " is " << tot_f_2d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_f_2d_scale = high_resolution_clock::now();
  test_2.test_scale_farr(ntimes,2, scale);
  auto end_f_2d_scale = high_resolution_clock::now();
  auto tot_f_2d_scale = duration_cast<microseconds>(end_f_2d_scale - start_f_2d_scale);
  cout<<"Average time for 2D FArray scale of size " <<size2<< " is " <<tot_f_2d_scale.count()/ntimes << "microseconds!"<<endl;

  //triad
  auto start_f_2d_triad = high_resolution_clock::now();
  test_2.test_triad_farr(ntimes, 2, scale);
  auto end_f_2d_triad = high_resolution_clock::now();
  auto tot_f_2d_triad = duration_cast<microseconds>(end_f_2d_triad - start_f_2d_triad);
  cout << " Average time for 2D FArray triad of size " << size2 << " is " << tot_f_2d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");

  //CArray
  //copy
  auto start_c_2d_copy = high_resolution_clock::now();
  test_2.test_copy_carr(ntimes, 2);
  auto end_c_2d_copy = high_resolution_clock::now();
  auto tot_c_2d_copy = duration_cast<microseconds>(end_c_2d_copy - start_c_2d_copy);
  cout<<"Average time for 2D CArray copy of size " << size2 << " is " << tot_c_2d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_c_2d_scale = high_resolution_clock::now();
  test_2.test_scale_carr(ntimes,2, scale);
  auto end_c_2d_scale = high_resolution_clock::now();
  auto tot_c_2d_scale = duration_cast<microseconds>(end_c_2d_scale - start_c_2d_scale);
  cout<<"Average time for 2D CArray scale of size " <<size2<< " is " <<tot_c_2d_scale.count()/ntimes << "microseconds!"<<endl;

  //triad
  auto start_c_2d_triad = high_resolution_clock::now();
  test_2.test_triad_carr(ntimes, 2, scale);
  auto end_c_2d_triad = high_resolution_clock::now();
  auto tot_c_2d_triad = duration_cast<microseconds>(end_c_2d_triad - start_c_2d_triad);
  cout << " Average time for 2D CArray triad of size " << size2 << " is " << tot_c_2d_triad.count()/ntimes << " microseconds!" <<endl;

    printf("\n\n");
   
    //CMatrix
    auto start_cm_2d_copy = high_resolution_clock::now();
    test_2.test_copy_cmat(ntimes, 2);
    auto stop_cm_2d_copy = high_resolution_clock::now();
    auto tot_cm_2d_copy = duration_cast<microseconds>(stop_cm_2d_copy - start_cm_2d_copy);
    cout<< "Average time for a 2D CMatrix copy of size " << size2<< " is " << tot_cm_2d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_cm_2d_scale = high_resolution_clock::now();
    test_2.test_scale_cmat(ntimes, 2, scale);
    auto end_cm_2d_scale = high_resolution_clock::now();
    auto tot_cm_2d_scale = duration_cast<microseconds>(end_cm_2d_scale - start_cm_2d_scale);
    cout<<"Average time for a 2D CMatrix scale of size " << size2 << " is " << tot_cm_2d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad
    auto start_cm_2d_triad = high_resolution_clock::now();
    test_2.test_triad_cmat(ntimes, 2, scale);
    auto end_cm_2d_triad = high_resolution_clock::now();
    auto tot_cm_2d_triad = duration_cast<microseconds>(end_cm_2d_triad - start_cm_2d_triad);
    cout<<"Average time for 2D CMatrix triad of size " << size2 << " is " << tot_cm_2d_triad.count()/ntimes << " microseconds!"<<endl;

  printf("\n\n");

    //FMatrix
    auto start_fm_2d_copy = high_resolution_clock::now();
    test_2.test_copy_fmat(ntimes, 2);
    auto stop_fm_2d_copy = high_resolution_clock::now();
    auto tot_fm_2d_copy = duration_cast<microseconds>(stop_fm_2d_copy - start_fm_2d_copy);
    cout<< "Average time for a 2D FMatrix copy of size " << size2<< " is " << tot_fm_2d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_fm_2d_scale = high_resolution_clock::now();
    test_2.test_scale_fmat(ntimes, 2, scale);
    auto end_fm_2d_scale = high_resolution_clock::now();
    auto tot_fm_2d_scale = duration_cast<microseconds>(end_fm_2d_scale - start_fm_2d_scale);
    cout<<"Average time for 2D FMatrix scale of size " << size2 << " is " << tot_fm_2d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for farray
    auto start_fm_2d_triad = high_resolution_clock::now();
    test_2.test_triad_fmat(ntimes, 2, scale);
    auto end_fm_2d_triad = high_resolution_clock::now();
    auto tot_fm_2d_triad = duration_cast<microseconds>(end_fm_2d_triad - start_fm_2d_triad);
    cout<<"Average time for 2D FMatrix triad of size " << size2 << " is " << tot_fm_2d_triad.count()/ntimes << " microseconds!"<<endl;


  printf("\n\n");

  //~~~~~3D~~~~~~
   
  printf("~~~~~~3D~~~~~~~\n");
    
  //FArray
  //copy
  auto start_f_3d_copy = high_resolution_clock::now();
  test_3.test_copy_farr(ntimes, 3);
  auto end_f_3d_copy = high_resolution_clock::now();
  auto tot_f_3d_copy = duration_cast<microseconds>(end_f_3d_copy - start_f_3d_copy);
  cout<<"Average time for 3D FArray copy of size " << size3 << " is " << tot_f_3d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_f_3d_scale = high_resolution_clock::now();
  test_3.test_scale_farr(ntimes,3, scale);
  auto end_f_3d_scale = high_resolution_clock::now();
  auto tot_f_3d_scale = duration_cast<microseconds>(end_f_3d_scale - start_f_3d_scale);
  cout<<"Average time for 3D FArray scale of size " <<size3<< " is " <<tot_f_3d_scale.count()/ntimes << "microseconds!"<<endl;

  //triad
  auto start_f_3d_triad = high_resolution_clock::now();
  test_3.test_triad_farr(ntimes, 3, scale);
  auto end_f_3d_triad = high_resolution_clock::now();
  auto tot_f_3d_triad = duration_cast<microseconds>(end_f_3d_triad - start_f_3d_triad);
  cout << " Average time for 3D FArray triad of size " << size3 << " is " << tot_f_3d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");

  //CArray
  //copy
  auto start_c_3d_copy = high_resolution_clock::now();
  test_3.test_copy_carr(ntimes, 3);
  auto end_c_3d_copy = high_resolution_clock::now();
  auto tot_c_3d_copy = duration_cast<microseconds>(end_c_3d_copy - start_c_3d_copy);
  cout<<"Average time for 3D CArray copy of size " << size3 << " is " << tot_c_3d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_c_3d_scale = high_resolution_clock::now();
  test_3.test_scale_carr(ntimes,3, scale);
  auto end_c_3d_scale = high_resolution_clock::now();
  auto tot_c_3d_scale = duration_cast<microseconds>(end_c_3d_scale - start_c_3d_scale);
  cout<<"Average time for 3D CArray scale of size " <<size3<< " is " <<tot_c_3d_scale.count()/ntimes << "microseconds!"<<endl;

  //triad
  auto start_c_3d_triad = high_resolution_clock::now();
  test_3.test_triad_carr(ntimes, 3, scale);
  auto end_c_3d_triad = high_resolution_clock::now();
  auto tot_c_3d_triad = duration_cast<microseconds>(end_c_3d_triad - start_c_3d_triad);
  cout << " Average time for 3D CArray triad of size " << size3 << " is " << tot_c_3d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");


    //CMatrix
    auto start_cm_3d_copy = high_resolution_clock::now();
    test_3.test_copy_cmat(ntimes, 3);
    auto stop_cm_3d_copy = high_resolution_clock::now();
    auto tot_cm_3d_copy = duration_cast<microseconds>(stop_cm_3d_copy - start_cm_3d_copy);
    cout<< "Average time for a 3D CMatrix copy of size " << size3<< " is " << tot_cm_3d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_cm_3d_scale = high_resolution_clock::now();
    test_3.test_scale_cmat(ntimes, 3, scale);
    auto end_cm_3d_scale = high_resolution_clock::now();
    auto tot_cm_3d_scale = duration_cast<microseconds>(end_cm_3d_scale - start_cm_3d_scale);
    cout<<"Average time for a 3D CMatrix scale of size " << size3 << " is " << tot_cm_3d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad
    auto start_cm_3d_triad = high_resolution_clock::now();
    test_3.test_triad_cmat(ntimes, 3, scale);
    auto end_cm_3d_triad = high_resolution_clock::now();
    auto tot_cm_3d_triad = duration_cast<microseconds>(end_cm_3d_triad - start_cm_3d_triad);
    cout<<"Average time for 3D CMatrix triad of size " << size3 << " is " << tot_cm_3d_triad.count()/ntimes << " microseconds!"<<endl;

    printf("\n\n");

    //FMatrix
    auto start_fm_3d_copy = high_resolution_clock::now();
    test_3.test_copy_fmat(ntimes, 3);
    auto stop_fm_3d_copy = high_resolution_clock::now();
    auto tot_fm_3d_copy = duration_cast<microseconds>(stop_fm_3d_copy - start_fm_3d_copy);
    cout<< "Average time for a 3D FMatrix copy of size " << size3<< " is " << tot_fm_3d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_fm_3d_scale = high_resolution_clock::now();
    test_3.test_scale_fmat(ntimes, 3, scale);
    auto end_fm_3d_scale = high_resolution_clock::now();
    auto tot_fm_3d_scale = duration_cast<microseconds>(end_fm_3d_scale - start_fm_3d_scale);
    cout<<"Average time for 3D FMatrix scale of size " << size3 << " is " << tot_fm_3d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for farray
    auto start_fm_3d_triad = high_resolution_clock::now();
    test_3.test_triad_fmat(ntimes, 3, scale);
    auto end_fm_3d_triad = high_resolution_clock::now();
    auto tot_fm_3d_triad = duration_cast<microseconds>(end_fm_3d_triad - start_fm_3d_triad);
    cout<<"Average time for 3D FMatrix triad of size " << size3 << " is " << tot_fm_3d_triad.count()/ntimes << " microseconds!"<<endl;

    printf("\n\n");
  

  //~~~~~~~4D~~~~~
  printf("~~~~~~~~4D DATA~~~~~~~\n");
  //FArray
  //copy
  auto start_f_4d_copy = high_resolution_clock::now();
  test_4.test_copy_farr(ntimes, 4);
  auto end_f_4d_copy = high_resolution_clock::now();
  auto tot_f_4d_copy = duration_cast<microseconds>(end_f_4d_copy - start_f_4d_copy);
  cout<<"Average time for 4D FArray copy of size " << size4 << " is " << tot_f_4d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_f_4d_scale = high_resolution_clock::now();
  test_4.test_scale_farr(ntimes, 4, scale);
  auto end_f_4d_scale = high_resolution_clock::now();
  auto tot_f_4d_scale = duration_cast<microseconds>(end_f_4d_scale - start_f_4d_scale);
  cout<<"Average time for 4D FArray scale of size " <<size4<< " is " <<tot_f_4d_scale.count()/ntimes << " microseconds!"<<endl;

  //triad
  auto start_f_4d_triad = high_resolution_clock::now();
  test_4.test_triad_farr(ntimes, 4, scale);
  auto end_f_4d_triad = high_resolution_clock::now();
  auto tot_f_4d_triad = duration_cast<microseconds>(end_f_4d_triad - start_f_4d_triad);
  cout << " Average time for 4D FArray triad of size " << size4 << " is " << tot_f_4d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");

  //CArray
  //copy
  auto start_c_4d_copy = high_resolution_clock::now();
  test_4.test_copy_carr(ntimes, 4);
  auto end_c_4d_copy = high_resolution_clock::now();
  auto tot_c_4d_copy = duration_cast<microseconds>(end_c_4d_copy - start_c_4d_copy);
  cout<<"Average time for 4D CArray copy of size " << size4 << " is " << tot_c_4d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_c_4d_scale = high_resolution_clock::now();
  test_4.test_scale_carr(ntimes, 4, scale);
  auto end_c_4d_scale = high_resolution_clock::now();
  auto tot_c_4d_scale = duration_cast<microseconds>(end_c_4d_scale - start_c_4d_scale);
  cout<<"Average time for 4D CArray scale of size " <<size4<< " is " <<tot_c_4d_scale.count()/ntimes << "microseconds!"<<endl;

  //triad
  auto start_c_4d_triad = high_resolution_clock::now();
  test_4.test_triad_carr(ntimes, 4, scale);
  auto end_c_4d_triad = high_resolution_clock::now();
  auto tot_c_4d_triad = duration_cast<microseconds>(end_c_4d_triad - start_c_4d_triad);
  cout << " Average time for 4D CArray triad of size " << size4 << " is " << tot_c_4d_triad.count()/ntimes << " microseconds!" <<endl;

    printf("\n\n");
    
    //CMatrix
    auto start_cm_4d_copy = high_resolution_clock::now();
    test_4.test_copy_cmat(ntimes, 4);
    auto stop_cm_4d_copy = high_resolution_clock::now();
    auto tot_cm_4d_copy = duration_cast<microseconds>(stop_cm_4d_copy - start_cm_4d_copy);
    cout<< "Average time for an 4D CMatrix copy of size " << size4<< " is " << tot_cm_4d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_cm_4d_scale = high_resolution_clock::now();
    test_4.test_scale_cmat(ntimes, 4, scale);
    auto end_cm_4d_scale = high_resolution_clock::now();
    auto tot_cm_4d_scale = duration_cast<microseconds>(end_cm_4d_scale - start_cm_4d_scale);
    cout<<"Average time for 4D CMatrix scale of size " << size4 << " is " << tot_cm_4d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad
    auto start_cm_4d_triad = high_resolution_clock::now();
    test_4.test_triad_cmat(ntimes, 4, scale);
    auto end_cm_4d_triad = high_resolution_clock::now();
    auto tot_cm_4d_triad = duration_cast<microseconds>(end_cm_4d_triad - start_cm_4d_triad);
    cout<<"Average time for4D CMatrix triad of size " << size4 << " is " << tot_cm_4d_triad.count()/ntimes << " microseconds!"<<endl;

    printf("\n\n");
   
    //FMatrix
    auto start_fm_4d_copy = high_resolution_clock::now();
    test_4.test_copy_fmat(ntimes, 4);
    auto stop_fm_4d_copy = high_resolution_clock::now();
    auto tot_fm_4d_copy = duration_cast<microseconds>(stop_fm_4d_copy - start_fm_4d_copy);
    cout<< "Average time for a 4D FMatrix copy of size " << size4<< " is " << tot_fm_4d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_fm_4d_scale = high_resolution_clock::now();
    test_4.test_scale_fmat(ntimes, 4, scale);
    auto end_fm_4d_scale = high_resolution_clock::now();
    auto tot_fm_4d_scale = duration_cast<microseconds>(end_fm_4d_scale - start_fm_4d_scale);
    cout<<"Average time for 4D FMatrix scale of size " << size4 << " is " << tot_fm_4d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad
    auto start_fm_4d_triad = high_resolution_clock::now();
    test_4.test_triad_fmat(ntimes, 4, scale);
    auto end_fm_4d_triad = high_resolution_clock::now();
    auto tot_fm_4d_triad = duration_cast<microseconds>(end_fm_4d_triad - start_fm_4d_triad);
    cout<<"Average time for 4D FMatrix triad of size " << size4 << " is " << tot_fm_4d_triad.count()/ntimes << " microseconds!"<<endl;

    printf("\n\n");

  //~~~~~~~~5D~~~~~~~
  printf("~~~~~~~5D DATA~~~~~~~~~\n");
  //FArray
  //copy
  auto start_f_5d_copy = high_resolution_clock::now();
  test_5.test_copy_farr(ntimes, 5);
  auto end_f_5d_copy = high_resolution_clock::now();
  auto tot_f_5d_copy = duration_cast<microseconds>(end_f_5d_copy - start_f_5d_copy);
  cout<<"Average time for 5D FArray copy of size " << size5 << " is " << tot_f_5d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_f_5d_scale = high_resolution_clock::now();
  test_5.test_scale_farr(ntimes,5, scale);
  auto end_f_5d_scale = high_resolution_clock::now();
  auto tot_f_5d_scale = duration_cast<microseconds>(end_f_5d_scale - start_f_5d_scale);
  cout<<"Average time for 5D FArray scale of size " <<size5<< " is " <<tot_f_5d_scale.count()/ntimes << " microseconds!"<<endl;

  //triad
  auto start_f_5d_triad = high_resolution_clock::now();
  test_5.test_triad_farr(ntimes, 5, scale);
  auto end_f_5d_triad = high_resolution_clock::now();
  auto tot_f_5d_triad = duration_cast<microseconds>(end_f_5d_triad - start_f_5d_triad);
  cout << " Average time for 5D FArray triad of size " << size5 << " is " << tot_f_5d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");

  //CArray
  //copy
  auto start_c_5d_copy = high_resolution_clock::now();
  test_5.test_copy_carr(ntimes, 5);
  auto end_c_5d_copy = high_resolution_clock::now();
  auto tot_c_5d_copy = duration_cast<microseconds>(end_c_5d_copy - start_c_5d_copy);
  cout<<"Average time for 5D CArray copy of size " << size5 << " is " << tot_c_5d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_c_5d_scale = high_resolution_clock::now();
  test_5.test_scale_carr(ntimes,5, scale);
  auto end_c_5d_scale = high_resolution_clock::now();
  auto tot_c_5d_scale = duration_cast<microseconds>(end_c_5d_scale - start_c_5d_scale);
  cout<<"Average time for 5D CArray scale of size " <<size5<< " is " <<tot_c_5d_scale.count()/ntimes << " microseconds!"<<endl;

  //triad
  auto start_c_5d_triad = high_resolution_clock::now();
  test_5.test_triad_carr(ntimes, 5, scale);
  auto end_c_5d_triad = high_resolution_clock::now();
  auto tot_c_5d_triad = duration_cast<microseconds>(end_c_5d_triad - start_c_5d_triad);
  cout << " Average time for 5D CArray triad of size " << size5 << " is " << tot_c_5d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");

    //CMatrix
    auto start_cm_5d_copy = high_resolution_clock::now();
    test_5.test_copy_cmat(ntimes, 5);
    auto stop_cm_5d_copy = high_resolution_clock::now();
    auto tot_cm_5d_copy = duration_cast<microseconds>(stop_cm_5d_copy - start_cm_5d_copy);
    cout<< "Average time for a 5D CMatrix copy of size " << size5<< " is " << tot_cm_5d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_cm_5d_scale = high_resolution_clock::now();
    test_5.test_scale_cmat(ntimes, 5, scale);
    auto end_cm_5d_scale = high_resolution_clock::now();
    auto tot_cm_5d_scale = duration_cast<microseconds>(end_cm_5d_scale - start_cm_5d_scale);
    cout<<"Average time for 5D CMatrix scale of size " << size5 << " is " << tot_cm_5d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad
    auto start_cm_5d_triad = high_resolution_clock::now();
    test_5.test_triad_cmat(ntimes, 5, scale);
    auto end_cm_5d_triad = high_resolution_clock::now();
    auto tot_cm_5d_triad = duration_cast<microseconds>(end_cm_5d_triad - start_cm_5d_triad);
    cout<<"Average time for 5D CMatrix triad of size " << size5 << " is " << tot_cm_5d_triad.count()/ntimes << " microseconds!"<<endl;
  //end CMatrix 1D test

    printf("\n\n");

    //FMatrix
    auto start_fm_5d_copy = high_resolution_clock::now();
    test_5.test_copy_fmat(ntimes, 5);
    auto stop_fm_5d_copy = high_resolution_clock::now();
    auto tot_fm_5d_copy = duration_cast<microseconds>(stop_fm_5d_copy - start_fm_5d_copy);
    cout<< "Average time for an 5D FMatrix copy of size " << size5<< " is " << tot_fm_5d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_fm_5d_scale = high_resolution_clock::now();
    test_5.test_scale_fmat(ntimes, 5, scale);
    auto end_fm_5d_scale = high_resolution_clock::now();
    auto tot_fm_5d_scale = duration_cast<microseconds>(end_fm_5d_scale - start_fm_5d_scale);
    cout<<"Average time for 5D FMatrix scale of size " << size5 << " is " << tot_fm_5d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for farray
    auto start_fm_5d_triad = high_resolution_clock::now();
    test_5.test_triad_fmat(ntimes, 5, scale);
    auto end_fm_5d_triad = high_resolution_clock::now();
    auto tot_fm_5d_triad = duration_cast<microseconds>(end_fm_5d_triad - start_fm_5d_triad);
    cout<<"Average time for 5D FMatrix triad of size " << size5 << " is " << tot_fm_5d_triad.count()/ntimes << " microseconds!"<<endl;

  printf("\n\n");

  //~~~~~~6D~~~~~~~
  printf("~~~~~~6D~~~~~~~\n");
  //FArray
  //copy
  auto start_f_6d_copy = high_resolution_clock::now();
  test_6.test_copy_farr(ntimes, 6);
  auto end_f_6d_copy = high_resolution_clock::now();
  auto tot_f_6d_copy = duration_cast<microseconds>(end_f_6d_copy - start_f_6d_copy);
  cout<<"Average time for 6D FArray copy of size " << size6 << " is " << tot_f_6d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_f_6d_scale = high_resolution_clock::now();
  test_6.test_scale_farr(ntimes,6, scale);
  auto end_f_6d_scale = high_resolution_clock::now();
  auto tot_f_6d_scale = duration_cast<microseconds>(end_f_6d_scale - start_f_6d_scale);
  cout<<"Average time for 6D FArray scale of size " <<size6<< " is " <<tot_f_6d_scale.count()/ntimes << " microseconds!"<<endl;

  //triad
  auto start_f_6d_triad = high_resolution_clock::now();
  test_6.test_triad_farr(ntimes, 6, scale);
  auto end_f_6d_triad = high_resolution_clock::now();
  auto tot_f_6d_triad = duration_cast<microseconds>(end_f_6d_triad - start_f_6d_triad);
  cout << " Average time for 6D FArray triad of size " << size6 << " is " << tot_f_6d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");

  //CArray
  //copy
  auto start_c_6d_copy = high_resolution_clock::now();
  test_6.test_copy_carr(ntimes, 6);
  auto end_c_6d_copy = high_resolution_clock::now();
  auto tot_c_6d_copy = duration_cast<microseconds>(end_c_6d_copy - start_c_6d_copy);
  cout<<"Average time for 6D CArray copy of size " << size6 << " is " << tot_c_6d_copy.count()/ntimes << " microseconds!"<<endl;

  //scale
  auto start_c_6d_scale = high_resolution_clock::now();
  test_6.test_scale_carr(ntimes, 6, scale);
  auto end_c_6d_scale = high_resolution_clock::now();
  auto tot_c_6d_scale = duration_cast<microseconds>(end_c_6d_scale - start_c_6d_scale);
  cout<<"Average time for 6D CArray scale of size " <<size6<< " is " <<tot_c_6d_scale.count()/ntimes << "microseconds!"<<endl;

  //triad
  auto start_c_6d_triad = high_resolution_clock::now();
  test_6.test_triad_carr(ntimes, 6, scale);
  auto end_c_6d_triad = high_resolution_clock::now();
  auto tot_c_6d_triad = duration_cast<microseconds>(end_c_6d_triad - start_c_6d_triad);
  cout << " Average time for 6D CArray triad of size " << size6 << " is " << tot_c_6d_triad.count()/ntimes << " microseconds!" <<endl;

  printf("\n\n");
    
    //CMatrix
    auto start_cm_6d_copy = high_resolution_clock::now();
    test_6.test_copy_cmat(ntimes, 6);
    auto stop_cm_6d_copy = high_resolution_clock::now();
    auto tot_cm_6d_copy = duration_cast<microseconds>(stop_cm_6d_copy - start_cm_6d_copy);
    cout<< "Average time for a 6D CMatrix copy of size " << size6 << " is " << tot_cm_6d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_cm_6d_scale = high_resolution_clock::now();
    test_6.test_scale_cmat(ntimes, 6, scale);
    auto end_cm_6d_scale = high_resolution_clock::now();
    auto tot_cm_6d_scale = duration_cast<microseconds>(end_cm_6d_scale - start_cm_6d_scale);
    cout<<"Average time for 6D CMatrix scale of size " << size6 << " is " << tot_cm_6d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad
    auto start_cm_6d_triad = high_resolution_clock::now();
    test_6.test_triad_cmat(ntimes, 6, scale);
    auto end_cm_6d_triad = high_resolution_clock::now();
    auto tot_cm_6d_triad = duration_cast<microseconds>(end_cm_6d_triad - start_cm_6d_triad);
    cout<<"Average time for 6D CMatrix triad of size " << size6 << " is " << tot_cm_6d_triad.count()/ntimes << " microseconds!"<<endl;

    printf("\n\n");

    //FMatrix
    auto start_fm_6d_copy = high_resolution_clock::now();
    test_6.test_copy_fmat(ntimes, 6);
    auto stop_fm_6d_copy = high_resolution_clock::now();
    auto tot_fm_6d_copy = duration_cast<microseconds>(stop_fm_6d_copy - start_fm_6d_copy);
    cout<< "Average time for an 6D FMatrix copy of size " << size6<< " is " << tot_fm_6d_copy.count()/ntimes<< " microseconds!"<<endl;

    //test scale
    auto start_fm_6d_scale = high_resolution_clock::now();
    test_6.test_scale_fmat(ntimes, 6, scale);
    auto end_fm_6d_scale = high_resolution_clock::now();
    auto tot_fm_6d_scale = duration_cast<microseconds>(end_fm_6d_scale - start_fm_6d_scale);
    cout<<"Average time for 6D FMatrix scale of size " << size6 << " is " << tot_fm_6d_scale.count()/ntimes << " microseconds!"<<endl;

    //test the triad for farray
    auto start_fm_6d_triad = high_resolution_clock::now();
    test_6.test_triad_fmat(ntimes, 6, scale);
    auto end_fm_6d_triad = high_resolution_clock::now();
    auto tot_fm_6d_triad = duration_cast<microseconds>(end_fm_6d_triad - start_fm_6d_triad);
    cout<<"Average time for 6D FMatrix triad of size " << size6 << " is " << tot_fm_6d_triad.count()/ntimes << " microseconds!"<<endl;


  printf("\n\n");

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DONE!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~commending out print statements to check if the test functions worked~~~~~~
  /*
    //check if it worked
    for(int i = 0; i < 15; i++) {
	printf("test1.f_arr2(%d) = %d\n", i, test_1.f_arr2(i) );
	printf("test1.arr1[%d] = %d\n", i, test_1.arr1[i] );
	}


    //test scale function
    test_1.test_scale_farr(ntimes, 1, 2);

     //check if it worked
     for(int i = 0; i < 15; i++) {
	printf("test1.f_arr2(%d) = %d (after scaled) \n", i, test_1.f_arr2(i) );
	}

    //test triad
    test_1.test_triad_farr(ntimes, 1, 2);
    
    //check if it worked
    for(int i = 0; i < 15; i++) {
	printf("test1.f_arr3(%d) = %d\n", i, test_1.f_arr3(i) );
	}

    */


 }//done! end main















