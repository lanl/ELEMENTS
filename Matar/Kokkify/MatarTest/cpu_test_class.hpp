

#ifndef CPU_TEST_CLASS_H
#define CPU_TEST_CLASS_H

#include <chrono>
#include <iostream>

#include "matar.h"

using namespace std;
using namespace std::chrono;

class cpu_test_class {

    public:
	int size1, size2, size3, size4, size5, size6;
	int length;

	//arrays for testing


	//FArray
	FArray <int> f_arr1;
	FArray <int> f_arr2;
	FArray <int> f_arr3;

	//CArray
	CArray <int> c_arr1;
	CArray <int> c_arr2;
	CArray <int> c_arr3;

	//CMatrix
	CMatrix <int> c_mat1;
	CMatrix <int> c_mat2;
	CMatrix <int> c_mat3;

	//FMatrix
	FMatrix <int> f_mat1;
	FMatrix <int> f_mat2;
	FMatrix <int> f_mat3;
		
	//ViewFArray
	ViewFArray <int> v_farr1;
	ViewFArray <int> v_farr2;
	ViewFArray <int> v_farr3;

	//ViewCMatrix
	ViewCMatrix <int> v_cmat1;
	ViewCMatrix <int> v_cmat2;
	ViewCMatrix <int> v_cmat3;

	//1D
	void init( int pnt1);

	//2D
	void init(int pnt1, int pnt2);

	//3D
	void init(int pnt1, int pnt2, int pnt3);

	//4D
	void init(int pnt1, int pnt2, int pnt3, int pnt4);

	//5D
	void init(int pnt1, int pnt2, int pnt3, int pnt4, int pnt5);

	//6D
	void init(int pnt1, int pnt2, int pnt3, int pnt4, int pnt5, int pnt6);
	

 	//~~~test functions~~~~~~
 	
	//FArray
	void test_copy_farr(int, int);
	void test_scale_farr(int, int, int);
	void test_triad_farr( int, int, int);

	//CArray
	void test_copy_carr(int, int);
	void test_scale_carr(int, int, int);
	void test_triad_carr(int, int, int);

	//CMatrix
	void test_copy_cmat(int, int);
	void test_scale_cmat(int,int, int);
	void test_triad_cmat(int, int, int);

	//FMatrix
	void test_copy_fmat(int, int);
	void test_scale_fmat(int, int, int);
	void test_triad_fmat(int, int, int);	

	//destructor
	~cpu_test_class() {};


};

//~~~~~begin class declarations~~~~~~~~~

//1D init
void cpu_test_class::init(int pnt1)
{
	size1 = pnt1;
	length = size1;

	//FArray
	f_arr1 = FArray<int>(size1);
	f_arr2 = FArray<int>(size1);
	f_arr3 = FArray<int>(size1);

	//CMatrix
	c_mat1 = CMatrix<int>(size1);
	c_mat2 = CMatrix<int>(size1);
	c_mat3 = CMatrix<int>(size1);

	//CArray
	c_arr1 = CArray<int>(size1);
	c_arr2 = CArray<int>(size1);
	c_arr3 = CArray<int>(size1);

	//FMatrix
	f_mat1 = FMatrix<int>(size1);
	f_mat2 = FMatrix<int>(size1);
	f_mat3 = FMatrix<int>(size1);

} //end 1D

//2D
void cpu_test_class::init(int pnt1, int pnt2) {
 
  size1 = pnt1;
  size2 = pnt2;
  length = size1*size2;


  //FArray
  f_arr1 = FArray<int>(size1,size2);
  f_arr2 = FArray<int>(size1,size2);
  f_arr3 = FArray<int>(size1,size2);

  //CMatrix
  c_mat1 = CMatrix<int>(size1,size2);
  c_mat2 = CMatrix<int>(size1,size2);
  c_mat3 = CMatrix<int>(size1,size2);

  //CArray
  c_arr1 = CArray<int>(size1, size2);
  c_arr2 = CArray<int>(size1, size2);
  c_arr3 = CArray<int>(size1, size2);

  //FMatrix
  f_mat1 = FMatrix<int>(size1, size2);
  f_mat2 = FMatrix<int>(size1, size2);
  f_mat3 = FMatrix<int>(size1, size2);
}

//3D
void cpu_test_class::init(int pnt1, int pnt2, int pnt3) {

  size1 = pnt1;
  size2 = pnt2;
  size3 = pnt3;
  length = size1*size2*size3;

  //FArray
  f_arr1 = FArray<int>(size1,size2,size3);
  f_arr2 = FArray<int>(size1,size2,size3);
  f_arr3 = FArray<int>(size1,size2,size3);

  //CMatrix
  c_mat1 = CMatrix<int>(size1,size2,size3);
  c_mat2 = CMatrix<int>(size1,size2,size3);
  c_mat3 = CMatrix<int>(size1,size2,size3);

  //CArray
  c_arr1 = CArray<int>(size1, size2, size3);
  c_arr2 = CArray<int>(size1, size2, size3);
  c_arr3 = CArray<int>(size1, size2, size3);

  //FMatrix
  f_mat1 = FMatrix<int>(size1, size2, size3);
  f_mat2 = FMatrix<int>(size1, size2, size3);
  f_mat3 = FMatrix<int>(size1, size2, size3);
} //end 3D

//4D
void cpu_test_class::init(int pnt1, int pnt2, int pnt3, int pnt4)
{

  size1 = pnt1;
  size2 = pnt2;
  size3 = pnt3;
  size4 = pnt4;
  length = size1*size2*size3*size4;

  //FArray
  f_arr1 = FArray<int>(size1,size2,size3,size4);
  f_arr2 = FArray<int>(size1,size2,size3,size4);
  f_arr3 = FArray<int>(size1,size2,size3,size4);

  //CMatrix
  c_mat1 = CMatrix<int>(size1,size2,size3,size4);
  c_mat2 = CMatrix<int>(size1,size2,size3,size4);
  c_mat3 = CMatrix<int>(size1,size2,size3,size4);

  //CArray
  c_arr1 = CArray<int>(size1, size2, size3, size4);
  c_arr2 = CArray<int>(size1, size2, size3, size4);
  c_arr3 = CArray<int>(size1, size2, size3, size4);

  //FMatrix
  f_mat1 = FMatrix<int>(size1, size2, size3, size4);
  f_mat2 = FMatrix<int>(size1, size2, size3, size4);
  f_mat3 = FMatrix<int>(size1, size2, size3, size4);
} //end 4D

//5D
void cpu_test_class::init(int pnt1, int pnt2, int pnt3, \
 int pnt4, int pnt5)
{

  size1 = pnt1;
  size2 = pnt2;
  size3 = pnt3;
  size4 = pnt4;
  size5 = pnt5;
  length = size1*size2*size3*size4*size5;

  //FArray
  f_arr1 = FArray<int>(size1,size2,size3,size4,size5);
  f_arr2 = FArray<int>(size1,size2,size3,size4,size5);
  f_arr3 = FArray<int>(size1,size2,size3,size4,size5);

  //CMatrix
  c_mat1 = CMatrix<int>(size1,size2,size3,size4,size5);
  c_mat2 = CMatrix<int>(size1,size2,size3,size4,size5);
  c_mat3 = CMatrix<int>(size1,size2,size3,size4,size5);
  //CArray
  c_arr1 = CArray<int>(size1, size2, size3, size4, size5);
  c_arr2 = CArray<int>(size1, size2, size3, size4, size5);
  c_arr3 = CArray<int>(size1, size2, size3, size4, size5);

  //FMatrix
  f_mat1 = FMatrix<int>(size1, size2, size3, size4, size5);
  f_mat2 = FMatrix<int>(size1, size2, size3, size4, size5);
  f_mat3 = FMatrix<int>(size1, size2, size3, size4, size5);

} //end 5D

//6D
void cpu_test_class::init(int pnt1, int pnt2, int pnt3, \
 int pnt4, int pnt5, int pnt6)
{

  size1 = pnt1;
  size2 = pnt2;
  size3 = pnt3;
  size4 = pnt4;
  size5 = pnt5;
  size6 = pnt6;
  length = size1*size2*size3*size4*size5*size6;

  //FArray
  f_arr1 = FArray<int>(size1,size2,size3,size4,size5,size6);
  f_arr2 = FArray<int>(size1,size2,size3,size4,size5,size6);
  f_arr3 = FArray<int>(size1,size2,size3,size4,size5,size6);

  //CMatrix
  c_mat1 = CMatrix<int>(size1,size2,size3,size4,size5,size6);
  c_mat2 = CMatrix<int>(size1,size2,size3,size4,size5,size6);
  c_mat3 = CMatrix<int>(size1,size2,size3,size4,size5,size6);

  //CArray
  c_arr1 = CArray<int>(size1, size2, size3, size4, size5, size6);
  c_arr2 = CArray<int>(size1, size2, size3, size4, size5, size6);
  c_arr3 = CArray<int>(size1, size2, size3, size4, size5, size6);

  //FMatrix
  f_mat1 = FMatrix<int>(size1, size2, size3, size4, size5, size6);
  f_mat2 = FMatrix<int>(size1, size2, size3, size4, size5, size6);
  f_mat3 = FMatrix<int>(size1, size2, size3, size4, size5, size6);

} //end 6D

/* Note: for all the Stream test, the array ordering is as follows:
 * 	1. for copy: array 1 gets coppied to 2; arr2 = arr1;
 * 	2. for scale: array 2 gets scaled
 * 	3. 3 holds the sum; arr3 = const*arr2 + arr1
 * This is just to hold consistency to the numbering of the arrays/matrices. 
 */ 	


//~~~~~~~~~~stream triad functions~~~~~~

//FArray tests

//copy
void cpu_test_class::test_copy_farr(int ntimes, int dims)
{


   switch(dims) {
	
     //1D
     case(1):
       
       //auto start_time = high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 0; j < size1; j++) {
	  f_arr2(j) = f_arr1(j);
	  }
	}

	//auto stop_time = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(stop_time - start_time);
	//cout<<"Time for FArray copy an array of size "<< size1<< " is " << duration.count() <<"microseconds!"<<endl;
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
		f_arr2(i,j) = f_arr1(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int k = 0; k < size3; k++) {
	   for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
	      f_arr2(i,j,k) = f_arr1(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int l = 0; l < size4; l++) {
	  for(int k = 0; k < size3; k++) {
	   for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
	 	f_arr2(i,j,k,l) = f_arr1(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int m = 0; m < size5; m++) {
	   for(int l = 0; l < size4; l++) {
	    for(int k = 0; k < size3; k++) {
	     for(int j = 0; j < size2; j++) {
	      for(int i = 0; i < size1; i++) {
		f_arr2(i,j,k,l,m) = f_arr1(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int n = 0; n < size6; n++) {
	  for(int m = 0; m < size5; m++) {
	   for(int l = 0; l < size4; l++) {
	    for(int k = 0; k < size3; k++) {
	     for(int j = 0; j < size2; j++) {
	      for(int i = 0; i < size1; i++) {
		f_arr2(i,j,k,l,m,n) = f_arr1(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end copy test


//scale
//f_arr2(i) = const*f_arr2(i)
void cpu_test_class::test_scale_farr(int ntimes, int dims, int scalar)
{

   switch(dims) {
	
     //1D
     case(1):
       for(int i = 0; i < ntimes; i++) {
	//printf("at iteration %d\n", i);
	for(int j = 0; j < size1; j++) {
	  //printf("f_arr2(%d) = %d\n", j, f_arr2(j));
	  f_arr2(j) = scalar*f_arr2(j);
	  }
	}
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
		f_arr2(i,j) = scalar*f_arr2(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int k = 0; k < size3; k++) {
	   for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
	      f_arr2(i,j,k) = scalar*f_arr2(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int l = 0; l < size4; l++) {
	  for(int k = 0; k < size3; k++) {
	   for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
	 	f_arr2(i,j,k,l) = scalar*f_arr2(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int m = 0; m < size5; m++) {
	   for(int l = 0; l < size4; l++) {
	    for(int k = 0; k < size3; k++) {
	     for(int j = 0; j < size2; j++) {
	      for(int i = 0; i < size1; i++) {
		f_arr2(i,j,k,l,m) = scalar*f_arr2(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int n = 0; n < size6; n++) {
	  for(int m = 0; m < size5; m++) {
	   for(int l = 0; l < size4; l++) {
	    for(int k = 0; k < size3; k++) {
	     for(int j = 0; j < size2; j++) {
	      for(int i = 0; i < size1; i++) {
		f_arr2(i,j,k,l,m,n) = scalar*f_arr2(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end scale test

//triad test
//arr3 = scale*arr2 + arr1
void cpu_test_class::test_triad_farr(int ntimes, int dims, int scalar)
{

   switch(dims) {
	
     //1D
     case(1):
       for(int i = 0; i < ntimes; i++) {
	for(int j = 0; j < size1; j++) {

	  f_arr3(j) = scalar*f_arr2(j) + f_arr1(j);
	  }
	}
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
		  f_arr3(i,j) = scalar*f_arr2(i,j) + f_arr1(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int k = 0; k < size3; k++) {
	   for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
		  f_arr3(i,j,k) = scalar*f_arr2(i,j,k) + f_arr1(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int l = 0; l < size4; l++) {
	  for(int k = 0; k < size3; k++) {
	   for(int j = 0; j < size2; j++) {
	    for(int i = 0; i < size1; i++) {
		  f_arr3(i,j,k,l) = scalar*f_arr2(i,j,k,l) + f_arr1(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int m = 0; m < size5; m++) {
	   for(int l = 0; l < size4; l++) {
	    for(int k = 0; k < size3; k++) {
	     for(int j = 0; j < size2; j++) {
	      for(int i = 0; i < size1; i++) {
		  f_arr3(i,j,k,l,m) = scalar*f_arr2(i,j,k,l,m) + f_arr1(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int n = 0; n < size6; n++) {
	  for(int m = 0; m < size5; m++) {
	   for(int l = 0; l < size4; l++) {
	    for(int k = 0; k < size3; k++) {
	     for(int j = 0; j < size2; j++) {
	      for(int i = 0; i < size1; i++) {
		  f_arr3(i,j,k,l,m,n) = scalar*f_arr2(i,j,k,l,m,n) + f_arr1(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end triad test

//~~~~~~test for CArray~~~~~~~~~

//copy
void cpu_test_class::test_copy_carr(int ntimes, int dims)
{

  //Note: timer doesn't work here :(
  // auto start = std::chrono::high_resolution_clock::now();

   switch(dims) {
	
     //1D
     case(1):
      // auto start = std::chrono::high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 0; j < size1; j++) {
	  c_arr2(j) = c_arr1(j);
	  }
	}
	//std::cout<< "Elapsed time of an average of " <<ntimes<< "runs with FArray of size "<< length<< "is: "<< watch.ElapsedMilliseconds() << " ms!\n";
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	//std::cout<<"Average runtime for a copy of FArray of size "<< length " with" << ntimes << "runs is" << duration.count() <<"s\n";
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int i = 0; i < size1; i++) {
	    for(int j = 0; j < size2; j++) {
		c_arr2(i,j) = c_arr1(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int i = 0; i < size1; i++) {
	   for(int j = 0; j < size2; j++) {
	    for(int k = 0; k < size3; k++) {
	      c_arr2(i,j,k) = c_arr1(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int i = 0; i < size1; i++) {
	  for(int j = 0; j < size2; j++) {
	   for(int k = 0; k < size3; k++) {
	    for(int l = 0; l < size4; l++) {
	 	c_arr2(i,j,k,l) = c_arr1(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int i = 0; i < size1; i++) {
	   for(int j = 0; j < size2; j++) {
	    for(int k = 0; k < size3; k++) {
	     for(int l = 0; l < size4; l++) {
	      for(int m = 0; m < size5; m++) {
		c_arr2(i,j,k,l,m) = c_arr1(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int i = 0; i < size1; i++) {
	  for(int j = 0; j < size2; j++) {
	   for(int k = 0; k < size3; k++) {
	    for(int l = 0; l < size4; l++) {
	     for(int m = 0; m < size5; m++) {
	      for(int n = 0; n < size6; n++) {
		c_arr2(i,j,k,l,m,n) = c_arr1(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end copy test


//scale
void cpu_test_class::test_scale_carr(int ntimes, int dims, int scale)
{

  // auto start = std::chrono::high_resolution_clock::now();

   switch(dims) {
	
     //1D
     case(1):
      // auto start = std::chrono::high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 0; j < size1; j++) {
	  c_arr2(j) = scale*c_arr2(j);
	  }
	}
	//std::cout<< "Elapsed time of an average of " <<ntimes<< "runs with FArray of size "<< length<< "is: "<< watch.ElapsedMilliseconds() << " ms!\n";
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	//std::cout<<"Average runtime for a copy of FArray of size "<< length " with" << ntimes << "runs is" << duration.count() <<"s\n";
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int i = 0; i < size1; i++) {
	    for(int j = 0; j < size2; j++) {
		c_arr2(i,j) = scale*c_arr2(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int i = 0; i < size1; i++) {
	   for(int j = 0; j < size2; j++) {
	    for(int k = 0; k < size3; k++) {
	      c_arr2(i,j,k) = scale*c_arr2(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int i = 0; i < size1; i++) {
	  for(int j = 0; j < size2; j++) {
	   for(int k = 0; k < size3; k++) {
	    for(int l = 0; l < size4; l++) {
	 	c_arr2(i,j,k,l) = scale*c_arr2(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int i = 0; i < size1; i++) {
	   for(int j = 0; j < size2; j++) {
	    for(int k = 0; k < size3; k++) {
	     for(int l = 0; l < size4; l++) {
	      for(int m = 0; m < size5; m++) {
		c_arr2(i,j,k,l,m) = scale*c_arr2(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int i = 0; i < size1; i++) {
	  for(int j = 0; j < size2; j++) {
	   for(int k = 0; k < size3; k++) {
	    for(int l = 0; l < size4; l++) {
	     for(int m = 0; m < size5; m++) {
	      for(int n = 0; n < size6; n++) {
		c_arr2(i,j,k,l,m,n) = scale*c_arr2(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end scale test

//triad
void cpu_test_class::test_triad_carr(int ntimes, int dims, int scale)
{

  // auto start = std::chrono::high_resolution_clock::now();

   switch(dims) {
	
     //1D
     case(1):
      // auto start = std::chrono::high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 0; j < size1; j++) {
	  c_arr3(j) = scale*c_arr2(j) + c_arr1(j);
	  }
	}
	//watch.Stop();
	//std::cout<< "Elapsed time of an average of " <<ntimes<< "runs with FArray of size "<< length<< "is: "<< watch.ElapsedMilliseconds() << " ms!\n";
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	//std::cout<<"Average runtime for a copy of FArray of size "<< length " with" << ntimes << "runs is" << duration.count() <<"s\n";
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int i = 0; i < size1; i++) {
	    for(int j = 0; j < size2; j++) {
		c_arr3(i,j) = scale*c_arr2(i,j) + c_arr1(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int i = 0; i < size1; i++) {
	   for(int j = 0; j < size2; j++) {
	    for(int k = 0; k < size3; k++) {
	      c_arr3(i,j,k) = scale*c_arr2(i,j,k) + c_arr1(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int i = 0; i < size1; i++) {
	  for(int j = 0; j < size2; j++) {
	   for(int k = 0; k < size3; k++) {
	    for(int l = 0; l < size4; l++) {
	 	c_arr3(i,j,k,l) = scale*c_arr2(i,j,k,l) + c_arr1(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int i = 0; i < size1; i++) {
	   for(int j = 0; j < size2; j++) {
	    for(int k = 0; k < size3; k++) {
	     for(int l = 0; l < size4; l++) {
	      for(int m = 0; m < size5; m++) {
		c_arr3(i,j,k,l,m) = scale*c_arr2(i,j,k,l,m) + c_arr1(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int i = 0; i < size1; i++) {
	  for(int j = 0; j < size2; j++) {
	   for(int k = 0; k < size3; k++) {
	    for(int l = 0; l < size4; l++) {
	     for(int m = 0; m < size5; m++) {
	      for(int n = 0; n < size6; n++) {
		c_arr3(i,j,k,l,m,n) = scale*c_arr2(i,j,k,l,m,n) + c_arr1(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end triad test

//~~~~~~~test for CMatrix~~~~~~~~~
//need a different class because the indices are from 1,...N

//copy
void cpu_test_class::test_copy_cmat(int ntimes, int dims)
{

     switch(dims) {

       //1D
       case(1):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++ ){
	  //iterate over dimensions
	  for(int i = 1; i < size1+1; i++) {
	   c_mat2(i) = c_mat1(i);
	   }
	  }
	 break;

	//2D
	case(2):
	
	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++){
	    for(int j = 1; j < size2+1; j++) {
	     c_mat2(i,j) = c_mat1(i,j);
	     }
	    }
	   }
	  break;

	//3D
	case(3):

	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
		c_mat2(i,j,k) = c_mat1(i,j,k);
		}
	       }
	      }
	     }
	    break;

	//4D
	case(4):

	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterator over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
		c_mat2(i,j,k,l) = c_mat1(i,j,k,l);
		}
	       }
	      }
	     }
	    }
	   break;

	//5D
	case(5):
	
	  //runs ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
	       for(int m = 1; m < size5+1; m++) {
		c_mat2(i,j,k,l,m) = c_mat1(i,j,k,l,m);
		}
	       }
	      }
	     }
	    }
	   }
	  break;

	//6D
	case(6):

	  //runs ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
	       for(int m = 1; m < size5+1; m++) {
		for(int n = 1; n < size6+1; n++) {
		 c_mat2(i,j,k,l,m,n) = c_mat1(i,j,k,l,m,n);
		}
	       }
	      }
	     }
	    }
	   }
	  }
	  break;

     } //end switch cases

} //end CMatrix copy test

//~~~test scale for CMatrix~~~
void cpu_test_class::test_scale_cmat(int ntimes, int dims, int scale)
{

     switch(dims) {

       //1D
       case(1):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++ ){
	  //iterate over dimensions
	  for(int i = 1; i < size1+1; i++) {
	   c_mat2(i) = scale*c_mat2(i);
	   }
	  }
	 break;

	//2D
	case(2):
	
	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++){
	    for(int j = 1; j < size2+1; j++) {
	     c_mat2(i,j) = scale*c_mat2(i,j);
	     }
	    }
	   }
	  break;

	//3D
	case(3):

	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
		c_mat2(i,j,k) = scale*c_mat2(i,j,k);
		}
	       }
	      }
	     }
	    break;

	//4D
	case(4):

	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterator over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
		c_mat2(i,j,k,l) = scale*c_mat2(i,j,k,l);
		}
	       }
	      }
	     }
	    }
	   break;

	//5D
	case(5):
	
	  //runs ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
	       for(int m = 1; m < size5+1; m++) {
		c_mat2(i,j,k,l,m) = scale*c_mat2(i,j,k,l,m);
		}
	       }
	      }
	     }
	    }
	   }
	  break;

	//6D
	case(6):

	  //runs ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
	       for(int m = 1; m < size5+1; m++) {
		for(int n = 1; n < size6+1; n++) {
		 c_mat2(i,j,k,l,m,n) = scale*c_mat2(i,j,k,l,m,n);
		}
	       }
	      }
	     }
	    }
	   }
	  }
	  break;

     } //end switch cases

} //end CMatrix scale test

//~~~~~~CMatrix triad test~~~~~~
void cpu_test_class::test_triad_cmat(int ntimes, int dims, int scale)
{

     switch(dims) {

       //1D
       case(1):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++ ){
	  //iterate over dimensions
	  for(int i = 1; i < size1+1; i++) {
	   c_mat3(i) = scale*c_mat2(i) + c_mat1(i);
	   }
	  }
	 break;

	//2D
	case(2):
	
	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++){
	    for(int j = 1; j < size2+1; j++) {
	     c_mat3(i,j) = c_mat2(i,j)*scale + c_mat1(i,j);
	     }
	    }
	   }
	  break;

	//3D
	case(3):

	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
		c_mat3(i,j,k) = scale*c_mat2(i,j,k) +  c_mat1(i,j,k);
		}
	       }
	      }
	     }
	    break;

	//4D
	case(4):

	  //run ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterator over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
		c_mat3(i,j,k,l) = scale*c_mat2(i,j,k,l) + c_mat1(i,j,k,l);
		}
	       }
	      }
	     }
	    }
	   break;

	//5D
	case(5):
	
	  //runs ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
	       for(int m = 1; m < size5+1; m++) {
		c_mat3(i,j,k,l,m) = scale*c_mat2(i,j,k,l,m) + c_mat1(i,j,k,l,m);
		}
	       }
	      }
	     }
	    }
	   }
	  break;

	//6D
	case(6):

	  //runs ntimes
	  for(int runs = 0; runs < ntimes; runs++) {
	   //iterate over dimensions
	   for(int i = 1; i < size1+1; i++) {
	    for(int j = 1; j < size2+1; j++) {
	     for(int k = 1; k < size3+1; k++) {
	      for(int l = 1; l < size4+1; l++) {
	       for(int m = 1; m < size5+1; m++) {
		for(int n = 1; n < size6+1;n ++) {
		  c_mat3(i,j,k,l,m,n) = scale*c_mat2(i,j,k,l,m,n) + c_mat1(i,j,k,l,m,n);
		}
	       }
	      }
	     }
	    }
	   }
	  }
	  break;

     } //end switch cases

} //end CMatrix triad test

//~~~~~~begin FMatrix tests~~~~~~~

//copy
void cpu_test_class::test_copy_fmat(int ntimes, int dims)
{

  // auto start = std::chrono::high_resolution_clock::now();

   switch(dims) {
	
     //1D
     case(1):
      // auto start = std::chrono::high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 1; j < size1+1; j++) {
	  f_mat2(j) = f_mat1(j);
	  }
	}
	//std::cout<< "Elapsed time of an average of " <<ntimes<< "runs with FArray of size "<< length<< "is: "<< watch.ElapsedMilliseconds() << " ms!\n";
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	//std::cout<<"Average runtime for a copy of FArray of size "<< length " with" << ntimes << "runs is" << duration.count() <<"s\n";
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
		f_mat2(i,j) = f_mat1(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int k = 1; k < size3+1; k++) {
	   for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
	      f_mat2(i,j,k) = f_mat1(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int l = 1; l < size4+1; l++) {
	  for(int k = 1; k < size3+1; k++) {
	   for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
	 	f_mat2(i,j,k,l) = f_mat1(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int m = 1; m < size5+1; m++) {
	   for(int l = 1; l < size4+1; l++) {
	    for(int k = 1; k < size3+1; k++) {
	     for(int j = 1; j < size2+1; j++) {
	      for(int i = 1; i < size1+1; i++) {
		f_mat2(i,j,k,l,m) = f_mat1(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int n = 1; n < size6+1; n++) {
	  for(int m = 1; m < size5+1; m++) {
	   for(int l = 1; l < size4+1; l++) {
	    for(int k = 1; k < size3+1; k++) {
	     for(int j = 1; j < size2+1; j++) {
	      for(int i = 1; i < size1+1; i++) {
		f_mat2(i,j,k,l,m,n) = f_mat1(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end copy test

//scale
void cpu_test_class::test_scale_fmat(int ntimes, int dims, int scale)
{

  // auto start = std::chrono::high_resolution_clock::now();

   switch(dims) {
	
     //1D
     case(1):
      // auto start = std::chrono::high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 1; j < size1+1; j++) {
	  f_mat2(j) = scale*f_mat2(j);
	  }
	}
	//std::cout<< "Elapsed time of an average of " <<ntimes<< "runs with FArray of size "<< length<< "is: "<< watch.ElapsedMilliseconds() << " ms!\n";
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	//std::cout<<"Average runtime for a copy of FArray of size "<< length " with" << ntimes << "runs is" << duration.count() <<"s\n";
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
		f_mat2(i,j) = scale*f_mat2(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int k = 1; k < size3+1; k++) {
	   for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
	      f_mat2(i,j,k) = scale*f_mat2(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int l = 1; l < size4+1; l++) {
	  for(int k = 1; k < size3+1; k++) {
	   for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
	 	f_mat2(i,j,k,l) = scale*f_mat2(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int m = 1; m < size5+1; m++) {
	   for(int l = 1; l < size4+1; l++) {
	    for(int k = 1; k < size3+1; k++) {
	     for(int j = 1; j < size2+1; j++) {
	      for(int i = 1; i < size1+1; i++) {
		f_mat2(i,j,k,l,m) = scale*f_mat2(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int n = 1; n < size6+1; n++) {
	  for(int m = 1; m < size5+1; m++) {
	   for(int l = 1; l < size4+1; l++) {
	    for(int k = 1; k < size3+1; k++) {
	     for(int j = 1; j < size2+1; j++) {
	      for(int i = 1; i < size1+1; i++) {
		f_mat2(i,j,k,l,m,n) = scale*f_mat2(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end scale test

//triad
void cpu_test_class::test_triad_fmat(int ntimes, int dims, int scale)
{

  // auto start = std::chrono::high_resolution_clock::now();

   switch(dims) {
	
     //1D
     case(1):
      // auto start = std::chrono::high_resolution_clock::now();
       for(int i = 0; i < ntimes; i++) {
	for(int j = 1; j < size1+1; j++) {
	  f_mat3(j) = scale*f_mat2(j) + f_mat1(j);
	  }
	}
	//std::cout<< "Elapsed time of an average of " <<ntimes<< "runs with FArray of size "<< length<< "is: "<< watch.ElapsedMilliseconds() << " ms!\n";
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = finish - start;
	//std::cout<<"Average runtime for a copy of FArray of size "<< length " with" << ntimes << "runs is" << duration.count() <<"s\n";
	break;

     //2D
     case(2):
	//run the example ntimes
	for(int run = 0; run < ntimes; run++) {
	  for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
		f_mat3(i,j) = scale*f_mat2(i,j) + f_mat1(i,j);
		}
	     }
	  }
	break;

      
     //3D
     case(3):

	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin dimension loops
	 for(int k = 1; k < size3+1; k++) {
	   for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
	      f_mat3(i,j,k) = scale*f_mat2(i,j,k) + f_mat1(i,j,k);
	     }
	    }
	   }
	  }
	 break;

     //4D
     case(4):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int l = 1; l < size4+1; l++) {
	  for(int k = 1; k < size3+1; k++) {
	   for(int j = 1; j < size2+1; j++) {
	    for(int i = 1; i < size1+1; i++) {
	 	f_mat3(i,j,k,l) = scale*f_mat2(i,j,k,l) + f_mat1(i,j,k,l);
	       }
	      }
	     }
	    }
	   }
	  break;

     //5D
     case(5):
	
	//run ntimes
	for(int run = 0; run < ntimes; run++) {
	  //begin iterating over dims
	  for(int m = 1; m < size5+1; m++) {
	   for(int l = 1; l < size4+1; l++) {
	    for(int k = 1; k < size3+1; k++) {
	     for(int j = 1; j < size2+1; j++) {
	      for(int i = 1; i < size1+1; i++) {
		f_mat3(i,j,k,l,m) = scale*f_mat2(i,j,k,l,m) + f_mat1(i,j,k,l,m);
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

     //6D
     case(6):

	//run times
	for(int run = 0; run < ntimes; run++) {
	 //begin iterating over dims
	 for(int n = 1; n < size6+1; n++) {
	  for(int m = 1; m < size5+1; m++) {
	   for(int l = 1; l < size4+1; l++) {
	    for(int k = 1; k < size3+1; k++) {
	     for(int j = 1; j < size2+1; j++) {
	      for(int i = 1; i < size1+1; i++) {
		f_mat3(i,j,k,l,m,n) = scale*f_mat2(i,j,k,l,m,n) + f_mat1(i,j,k,l,m,n);
	       }
	      }
	     }
	    }
	   }
	  }
	}
	break;  

    } //end switch

} //end triad










#endif
