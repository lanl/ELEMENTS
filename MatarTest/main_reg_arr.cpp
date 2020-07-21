/* Stream Benchmard for regular arrays
 * Note, unlike the test for the matar data,
 * I don't have a nice class for the test. The 
 * initialization, copy, scale, add is all shown in the main. 
 */

#include <stdio.h>
#include <chrono>
#include <iostream>



//sizes 
#define size1 4000000
#define size2 2000
#define size3 150
#define size4 45
#define size5 21
#define size6 13
#define ntimes 20

using namespace std;
using namespace std::chrono;


//begin!
int main() 
{

  printf("~~~~~~~1D Array~~~~~\n");
  
  //1D
  //int arr1_1d [size1];
  //int arr2_1d [size1];
  //int arr3_1d [size1];

  int * arr1_1d = (int*)malloc(size1 * sizeof(int*));
  int * arr2_1d = (int*)malloc(size1 * sizeof(int*));
  int * arr3_1d = (int*)malloc(size1 * sizeof(int*));

  //initialize 1D
  for(int i = 0; i < size1; i++) {
   arr1_1d[i] = 1;
  }

  //copy
  auto start_1d_copy = high_resolution_clock::now();

   for(int runs = 0; runs < ntimes; runs++) {
    for(int i = 0; i < size1; i++) {
	arr2_1d[i] = arr1_1d[i];
     }
    }

  auto end_1d_copy = high_resolution_clock::now();
  auto tot_1d_copy = duration_cast<microseconds>(end_1d_copy - start_1d_copy);
  cout<< "Time for regular 1D array of size "<<size1<<" is " << tot_1d_copy.count()/ntimes <<" microseconds!"<<endl;

  //scale
  auto start_1d_scale = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size1; i++) {
 	arr2_1d[i] = 1*arr2_1d[i];
      }
     }
 
  auto end_1d_scale = high_resolution_clock::now();
  auto tot_1d_scale = duration_cast<microseconds>(end_1d_scale - start_1d_scale);
  cout<<"Avg scale time for regular 1D array of size "<<size1<<" is "<<tot_1d_scale.count()/ntimes<<" microseconds!"<<endl;

  //triad
  auto start_1d_triad = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size1; i++) {
	arr3_1d[i] = 1*arr2_1d[i] + arr1_1d[i];
 	}
      }

  auto end_1d_triad = high_resolution_clock::now();
  auto tot_1d_triad = duration_cast<microseconds>(end_1d_triad - start_1d_triad);
  cout<<"Avg triad time for regular 1D array of size "<<size1<<" is "<<tot_1d_triad.count()/ntimes<<" microseconds!"<<endl;

  //free 1d array
  free(arr1_1d);
  free(arr2_1d);
  free(arr3_1d);


  //2D
  printf("~~~~~2D Array~~~~~\n");
  
  //int arr1_2d[size2][size2];
  //int arr2_2d[size2][size2];
  //int arr3_2d[size2][size2];

 
  int **arr1_2d = (int**)malloc(size2 * sizeof(int*));
  int **arr2_2d = (int**)malloc(size2 * sizeof(int*));
  int **arr3_2d = (int**)malloc(size2 * sizeof(int*));
  
  for(int i = 0; i < size2; i++) {
    arr1_2d[i] = (int*)malloc(size2 * sizeof(int));
    arr2_2d[i] = (int*)malloc(size2 * sizeof(int));
    arr3_2d[i] = (int*)malloc(size2 * sizeof(int));
  }
 

  //this was giving me an error, come back to it
 /*
  //create usin gnew
  int** arr1_2d = new *int[size2];
  int** arr2_2d = new *int[size2];
  int** arr3_2d = new *int[size2];

  for(int i = 0; i < size2; i++) {
     arr1_2d[i] = new int[size2];
     arr2_2d[i] = new int[size2];
     arr3_2d[i] = new int[size2];
   }
 */


  ///initialize
  for(int i = 0; i < size2; i++){
   for(int j = 0; j < size2; j++) {
    arr1_2d[i][j] = 1;
    }
   }


  //copy
  auto start_2d_copy = high_resolution_clock::now();

   for(int runs = 0; runs < ntimes; runs++) {
    for(int i = 0; i < size2; i++) {
     for(int j = 0; j< size2; j++) {
	arr2_2d[i][j] = arr1_2d[i][j];
     }
    }
   }

  auto end_2d_copy = high_resolution_clock::now();
  auto tot_2d_copy = duration_cast<microseconds>(end_2d_copy - start_2d_copy);
  cout<< "Time for regular 2D array of size "<<size2<<" is " << tot_2d_copy.count()/ntimes <<" microseconds!"<<endl;

  //scale
  auto start_2d_scale = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size2; i++) {
      for(int j = 0; j < size2; j++) {
 	arr2_2d[i][j] = 1*arr2_2d[i][j];
      }
     }
    }
 
  auto end_2d_scale = high_resolution_clock::now();
  auto tot_2d_scale = duration_cast<microseconds>(end_2d_scale - start_2d_scale);
  cout<<"Avg scale time for regular 2D array of size "<<size2<<" is "<<tot_2d_scale.count()/ntimes<<" microseconds!"<<endl;

  //triad
  auto start_2d_triad = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size2; i++) {
      for(int j = 0; j < size2; j++) {
	arr3_2d[i][j] = 1*arr2_2d[i][j] + arr1_2d[i][j];
 	}
      }
     }
  auto end_2d_triad = high_resolution_clock::now();
  auto tot_2d_triad = duration_cast<microseconds>(end_2d_triad - start_2d_triad);
  cout<<"Avg triad time for regular 2D array of size "<<size2<<" is "<<tot_2d_triad.count()/ntimes<<" microseconds!"<<endl;

  //3D
  printf("~~~~~~~3D~~~~~~~\n");

  //declare
  //int arr1_3d[size3][size3][size3];
  //int arr2_3d[size3][size3][size3];
  //int arr3_3d[size3][size3][size3];

   int ***arr1_3d = new int**[size3];
   int ***arr2_3d = new int**[size3];
   int ***arr3_3d = new int**[size3];

   for(int i = 0; i < size3; i++) {
     arr1_3d[i] = new int*[size3];
     arr2_3d[i] = new int*[size3];
     arr3_3d[i] = new int*[size3];
 
     for(int j = 0; j < size3; j++) {
      arr1_3d[i][j] = new int[size3];
      arr2_3d[i][j] = new int[size3];
      arr3_3d[i][j] = new int[size3];
     }
    }


  ///initialize
  for(int i = 0; i < size3; i++){
   for(int j = 0; j < size3; j++) {
    for(int k = 0; k < size3; k++) {
      arr1_3d[i][j][k] = 1;
    }
   }
  }
  //copy
  auto start_3d_copy = high_resolution_clock::now();

   for(int runs = 0; runs < ntimes; runs++) {
    for(int i = 0; i < size3; i++) {
     for(int j = 0; j< size3; j++) {
      for(int k = 0; k < size3; k++) {
	arr2_3d[i][j][k] = arr1_3d[i][j][k];
     }
    }
   }
  }
  auto end_3d_copy = high_resolution_clock::now();
  auto tot_3d_copy = duration_cast<microseconds>(end_3d_copy - start_3d_copy);
  cout<< "Time for regular 3D array of size "<<size3<<" is " << tot_3d_copy.count()/ntimes <<" microseconds!"<<endl;

  //scale
  auto start_3d_scale = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size3; i++) {
      for(int j = 0; j < size3; j++) {
       for(int k = 0; k < size3; k++) {
 	arr2_3d[i][j][k] = 1*arr2_3d[i][j][k];
      }
     }
    }
   } 
  auto end_3d_scale = high_resolution_clock::now();
  auto tot_3d_scale = duration_cast<microseconds>(end_3d_scale - start_3d_scale);
  cout<<"Avg scale time for regular 3D array of size "<<size3<<" is "<<tot_3d_scale.count()/ntimes<<" microseconds!"<<endl;

  //triad
  auto start_3d_triad = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size3; i++) {
      for(int j = 0; j < size3; j++) {
       for(int k = 0; k< size3; k++) {
	arr3_3d[i][j][k] = 1*arr2_3d[i][j][k] + arr1_3d[i][j][k];
 	}
      }
     }
    }
  auto end_3d_triad = high_resolution_clock::now();
  auto tot_3d_triad = duration_cast<microseconds>(end_3d_triad - start_3d_triad);
  cout<<"Avg triad time for regular 3D array of size "<<size3<<" is "<<tot_3d_triad.count()/ntimes<<" microseconds!"<<endl;

    //free arrays
    //code Eappen sent me
   for(int i = 0; i < size3; i++) {
     for(int j = 0; j < size3; j++) {
       delete [] arr1_3d[i][j];
       delete [] arr2_3d[i][j];
       delete [] arr3_3d[i][j];
      }
     delete [] arr1_3d[i];
     delete [] arr2_3d[i];
     delete [] arr3_3d[i];
    }
  
    delete [] arr1_3d;
    delete [] arr2_3d;
    delete [] arr3_3d;



  //4D
  printf("~~~~~~~~~~~~~4D~~~~~~~~~~~~~~~\n");


  //segmentation fault below
 /*
  int ****arr1_4d = new int***[size4];
  int ****arr2_4d = new int***[size4];
  int ****arr3_4d = new int***[size4];

  for(int i = 0; i < size4; i++) {
  
   arr1_4d[i] = new int**[size4];
   arr2_4d[i] = new int**[size4];
   arr3_4d[i] = new int**[size4];

   for(int j = 0; j < size4; j++) {
    
     arr1_4d[i][j] = new int*[size4];
     arr2_4d[i][j] = new int*[size4];
     arr2_4d[i][j] = new int*[size4];

     for(int k = 0; k < size4; k++) {

      arr1_4d[i][j][k] = new int[size4];
      arr2_4d[i][j][k] = new int[size4];
      arr3_4d[i][j][k] = new int[size4];
    
     }
    }
   }
  */
 
  int arr1_4d[size4][size4][size4][size4];
  int arr2_4d[size4][size4][size4][size4];
  int arr3_4d[size4][size4][size4][size4];
 
  //initialize
  for(int i = 0; i < size4; i++){
   for(int j = 0; j < size4; j++) {
    for(int k = 0; k < size4; k++) {
     for(int l = 0; l < size4; l++) {
       arr1_4d[i][j][k][l] = 1;
    }
   }
  }
 }
  //copy
  auto start_4d_copy = high_resolution_clock::now();

   for(int runs = 0; runs < ntimes; runs++) {
    for(int i = 0; i < size4; i++) {
     for(int j = 0; j< size4; j++) {
      for(int k = 0; k < size4; k++) {
	for(int l = 0; l < size4; l++) {
         arr2_4d[i][j][k][l] = arr1_4d[i][j][k][l];
     }
    }
   }
  }
 }

  auto end_4d_copy = high_resolution_clock::now();
  auto tot_4d_copy = duration_cast<microseconds>(end_4d_copy - start_4d_copy);
  cout<< "Time for regular 4D array of size "<<size4<<" is " << tot_4d_copy.count()/ntimes <<" microseconds!"<<endl;

  //scale
  auto start_4d_scale = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size4; i++) {
      for(int j = 0; j < size4; j++) {
       for(int k = 0; k < size4; k++) {
	for(int l = 0; l < size4; l++) { 
	  arr2_4d[i][j][k][l] = 1*arr2_4d[i][j][k][l];
	}
      }
     }
    }
   } 
  auto end_4d_scale = high_resolution_clock::now();
  auto tot_4d_scale = duration_cast<microseconds>(end_4d_scale - start_4d_scale);
  cout<<"Avg scale time for regular 4D array of size "<<size4<<" is "<<tot_4d_scale.count()/ntimes<<" microseconds!"<<endl;

  //triad
  auto start_4d_triad = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size4; i++) {
      for(int j = 0; j < size4; j++) {
       for(int k = 0; k< size4; k++) {
	for(int l = 0; l < size4; l++) {
	  arr3_4d[i][j][k][l] = 1*arr2_4d[i][j][k][l] + arr1_4d[i][j][k][l];
	}
 	}
      }
     }
    }
  auto end_4d_triad = high_resolution_clock::now();
  auto tot_4d_triad = duration_cast<microseconds>(end_4d_triad - start_4d_triad);
  cout<<"Avg triad time for regular 4D array of size "<<size4<<" is "<<tot_4d_triad.count()/ntimes<<" microseconds!"<<endl;


  //5D
  printf("~~~~~~~~~~~5D~~~~~~~~~~~~\n");

  int arr1_5d[size5][size5][size5][size5][size5];
  int arr2_5d[size5][size5][size5][size5][size5];
  int arr3_5d[size5][size5][size5][size5][size5];

  //initialize
  for(int i = 0; i < size5; i++){
   for(int j = 0; j < size5; j++) {
    for(int k = 0; k < size5; k++) {
     for(int l = 0; l < size5; l++) {
      for(int m = 0; m < size5; m++) {
	arr1_5d[i][j][k][l][m] = 1;
	}	
       }
      }
     }
    }


  //copy
  auto start_5d_copy = high_resolution_clock::now();

   for(int runs = 0; runs < ntimes; runs++) {
    for(int i = 0; i < size5; i++) {
     for(int j = 0; j < size5; j++) {
      for(int k = 0; k < size5; k++) {
	for(int l = 0; l < size5; l++) {
	 for(int m = 0; m < size5; m++){
           arr2_5d[i][j][k][l][m] = arr1_5d[i][j][k][l][m];
	}
     }
    }
   }
  }
 }

  auto end_5d_copy = high_resolution_clock::now();
  auto tot_5d_copy = duration_cast<microseconds>(end_5d_copy - start_5d_copy);
  cout<< "Time for regular 5D array of size "<<size5<<" is " << tot_5d_copy.count()/ntimes <<" microseconds!"<<endl;

  //scale
  auto start_5d_scale = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size5; i++) {
      for(int j = 0; j < size5; j++) {
       for(int k = 0; k < size5; k++) {
	for(int l = 0; l < size5; l++) { 
	 for(int m = 0; m < size5; m++) {
	  arr2_5d[i][j][k][l][m] = 1*arr2_5d[i][j][k][l][m];
	 }
	}
      }
     }
    }
   } 
  auto end_5d_scale = high_resolution_clock::now();
  auto tot_5d_scale = duration_cast<microseconds>(end_5d_scale - start_5d_scale);
  cout<<"Avg scale time for regular 5D array of size "<<size5<<" is "<<tot_5d_scale.count()/ntimes<<" microseconds!"<<endl;

  //triad
  auto start_5d_triad = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size5; i++) {
      for(int j = 0; j < size5; j++) {
       for(int k = 0; k < size5; k++) {
	for(int l = 0; l < size5; l++) {
	 for(int m = 0; m < size5; m++) {
	   arr3_5d[i][j][k][l][m] = 1*arr2_5d[i][j][k][l][m] + arr1_5d[i][j][k][l][m];
	}
	}
 	}
      }
     }
    }
  auto end_5d_triad = high_resolution_clock::now();
  auto tot_5d_triad = duration_cast<microseconds>(end_5d_triad - start_5d_triad);
  cout<<"Avg triad time for regular 5D array of size "<<size5<<" is "<<tot_5d_triad.count()/ntimes<<" microseconds!"<<endl;

  //6D
  printf("~~~~~~~~~6D~~~~~~~~~\n");

  int arr1_6d[size6][size6][size6][size6][size6][size6];
  int arr2_6d[size6][size6][size6][size6][size6][size6];
  int arr3_6d[size6][size6][size6][size6][size6][size6];

  //initialize
  
  for(int i = 0; i < size6; i++){
   for(int j = 0; j < size6; j++) {
    for(int k = 0; k < size6; k++) {
     for(int l = 0; l < size6; l++) {
      for(int m = 0; m < size6; m++) {
	for(int n = 0; n < size6; n++){
	 arr1_6d[i][j][k][l][m][n] = 1;
	}
	}	
       }
      }
     }
    }


  //copy
  auto start_6d_copy = high_resolution_clock::now();

   for(int runs = 0; runs < ntimes; runs++) {
    for(int i = 0; i < size6; i++) {
     for(int j = 0; j < size6; j++) {
      for(int k = 0; k < size6; k++) {
	for(int l = 0; l < size6; l++) {
	 for(int m = 0; m < size6; m++){
	  for(int n = 0; n < size6; n++){
           arr2_6d[i][j][k][l][m][n] = arr1_6d[i][j][k][l][m][n];
	  }
	}
     }
    }
   }
  }
 }

  auto end_6d_copy = high_resolution_clock::now();
  auto tot_6d_copy = duration_cast<microseconds>(end_6d_copy - start_6d_copy);
  cout<< "Time for regular 6D array of size "<<size6<<" is " << tot_6d_copy.count()/ntimes <<" microseconds!"<<endl;

  //scale
  auto start_6d_scale = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size6; i++) {
      for(int j = 0; j < size6; j++) {
       for(int k = 0; k < size6; k++) {
	for(int l = 0; l < size6; l++) { 
	 for(int m = 0; m < size6; m++) {
	  for(int n = 0; n < size6; n++) {
       	    arr2_6d[i][j][k][l][m][n] = 1*arr2_6d[i][j][k][l][m][n];
	   }
	 }
	}
      }
     }
    }
   } 
  auto end_6d_scale = high_resolution_clock::now();
  auto tot_6d_scale = duration_cast<microseconds>(end_6d_scale - start_6d_scale);
  cout<<"Avg scale time for regular 6D array of size "<<size6<<" is "<<tot_6d_scale.count()/ntimes<<" microseconds!"<<endl;

  //triad
  auto start_6d_triad = high_resolution_clock::now();

    for(int runs = 0; runs < ntimes; runs++) {
     for(int i = 0; i < size6; i++) {
      for(int j = 0; j < size6; j++) {
       for(int k = 0; k < size6; k++) {
	for(int l = 0; l < size6; l++) {
	 for(int m = 0; m < size6; m++) {
	  for(int n = 0; n < size6; n++){
	   arr3_6d[i][j][k][l][m][n] = 1*arr2_6d[i][j][k][l][m][n] + arr1_6d[i][j][k][l][m][n];
         }
	}
       }
      }
     }
    }
   }
  auto end_6d_triad = high_resolution_clock::now();
  auto tot_6d_triad = duration_cast<microseconds>(end_6d_triad - start_6d_triad);
  cout<<"Avg triad time for regular 6D array of size "<<size6<<" is "<<tot_6d_triad.count()/ntimes<<" microseconds!"<<endl;


} //DONE! end main

/*
//~~~~test1: copy~~~~~
void std_copy(int arr1[], int arr2[], int dim)
{
  switch(dim) {

    //1D
    case(1):
     
      for(int runs = 0; runs < ntimes; runs++) {
	for(int i = 0; i < size1; i++) {
	  arr2[i] = arr1[i];
	}
       }
      break;

    //2D
    case(2):

      for(int runs = 0; runs < ntimes; runs++) {
	for(int i = 0; i < size2; i++) {
	 for(int j = 0; j < size2; j++) {
	  arr2[i][j] = arr1[i][j];
	  }
	 }
	}
	break;

     //3D
     case(3):

	for(int runs = 0; runs < ntimes; runs++){
	 for(int i = 0; i < size3; i++) {
	  for(int j = 0; j < size3; j++) {
	   for(int k = 0; k < size3; k++) {
		arr2[i][j][k] = arr1[i][j][k];
		}
	       }
	      }
	     }
	  break;

    //4D
    case(4):

     for(int runs = 0; runs < ntimes; runs++) {
      for(int i = 0; i < size4; i++) {
	for(int j = 0; j < size4; j++) {
	 for(int k = 0; k < size4; k++) {
	  for(int l = 0; l < size4; l++) {
		arr2[i][j][k][l] = arr1[i][j][k][l];
		}
	      }
	     }
	    }
	   }
	 break;

    //5D
    case(5):

      for(int runs = 0; runs < ntimes; runs++) {
	for(int i = 0; i < size5; i++) {
	 for(int j = 0; j < size5; j++) {
	  for(int k = 0; k < size5; k++) {
	   for(int l = 0; l < size5; l++) {
	    for(int m = 0; m < size5; m++){
		arr2[i][j][k][l][m] = arr1[i][j][k][l][m];
		}
	       }
	      }
	     }
	    }
	   }
	  break;

     //6D
     case(6):

      for(int runs = 0; runs < ntimes; runs++) {
	for(int i = 0; i < size6; i++) {
	 for(int j = 0; j < size6; j++) {
	  for(int k = 0; k < size6; k++) {
	   for(int l = 0; l < size6; l++) {
	    for(int m = 0; m < size6; m++){
	     for(int n = 0; n < size6; n++) {
		arr2[i][j][k][l][m][n] = arr1[i][j][k][l][m][n];
		}
	       }
	      }
	     }
	    }
	   }
	  }
	 break;

  } //end switch



} //end copy
*/















