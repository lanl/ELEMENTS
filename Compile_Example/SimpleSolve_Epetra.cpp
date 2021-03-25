// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
   \file   SimpleSolve.cpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Sat Jul 17 10:35:39 2010

   \brief  Simple example of Amesos2 usage.

   This example solves a simple sparse system of linear equations using the
   Amesos2 interface to the Superlu solver.
*/

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <EpetraExt_CrsMatrixIn.h>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  typedef double Scalar;

  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;
  int nrows = 6;
  int *entries_per_row = new int[nrows];
  double **row_entries = new double*[nrows];
  int **row_indices = new int*[nrows];
  entries_per_row[0] = 3;
  entries_per_row[1] = 2;
  entries_per_row[2] = 1;
  entries_per_row[3] = 2;
  entries_per_row[4] = 2;
  entries_per_row[5] = 2;

  for(int i = 0; i < nrows; i++){
    row_entries[i] = new double[entries_per_row[i]];
    row_indices[i] = new int[entries_per_row[i]];
  }

  //fill values to test problem
  row_indices[0][0] = 0; row_indices[0][1] = 2; row_indices[0][2] = 4;
  row_indices[1][0] = 0; row_indices[1][1] = 1;
  row_indices[2][0] = 2;
  row_indices[3][0] = 0; row_indices[3][1] = 3;
  row_indices[4][0] = 1; row_indices[4][1] = 4;
  row_indices[5][0] = 3; row_indices[5][1] = 5;

  row_entries[0][0] = 7; row_entries[0][1] = -3; row_entries[0][2] = -1;
  row_entries[1][0] = 2; row_entries[1][1] = 8;
  row_entries[2][0] = 1;
  row_entries[3][0] = -3; row_entries[3][1] = 5;
  row_entries[4][0] = -1; row_entries[4][1] = 4;
  row_entries[5][0] = -2; row_entries[5][1] = 6;

  using Teuchos::RCP;
  using Teuchos::rcp;

  // Before we do anything, check that SuperLU is enabled
  if( !Amesos2::query("KLU2") ){
    std::cerr << "KLU2 not enabled in this run.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  #ifdef HAVE_MPI
  const Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  const Epetra_SerialComm comm;
#endif
  size_t myRank = comm.MyPID();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  *fos << myRank << " : " << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  //Parallel map
  Epetra_Map map = Epetra_Map(nrows, 0, comm);

  RCP<MAT> A = rcp( new MAT(View, map, 3) ); // max of three entries in a row

  /*
   * We will solve a system with a known solution, for which we will be using
   * the following matrix:
   *
   * [ [ 7,  0,  -3, 0,  -1, 0 ]
   *   [ 2,  8,  0,  0,  0,  0 ]
   *   [ 0,  0,  1,  0,  0,  0 ]
   *   [ -3, 0,  0,  5,  0,  0 ]
   *   [ 0,  -1, 0,  0,  4,  0 ]
   *   [ 0,  0,  0,  -2, 0,  6 ] ]
   *
   */
  // Construct matrix
  if(myRank==0){
    A->InsertGlobalValues(0,entries_per_row[0], row_entries[0], row_indices[0]);
    A->InsertGlobalValues(1,entries_per_row[1], row_entries[1], row_indices[1]);
    A->InsertGlobalValues(2,entries_per_row[2], row_entries[2], row_indices[2]);
    A->InsertGlobalValues(3,entries_per_row[3], row_entries[3], row_indices[3]);
    A->InsertGlobalValues(4,entries_per_row[4], row_entries[4], row_indices[4]);
    A->InsertGlobalValues(5,entries_per_row[5], row_entries[5], row_indices[5]);
  }
  A->FillComplete();
  
  A->Print(*(fos->getOStream()));
  // Create random X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->Random();

  /* Create B
   *
   * Use RHS:
   *
   *  [[-7]
   *   [18]
   *   [ 3]
   *   [17]
   *   [18]
   *   [28]]
   */
  RCP<MV> B = rcp(new MV(map,numVectors));
  int data[6] = {-7,18,3,17,18,28};
  for( int i = 0; i < 6; ++i ){
      B->ReplaceGlobalValue(i,0,data[i]);
  }
  
  B->Print(*(fos->getOStream()));

  // Create solver interface to Superlu with Amesos2 factory method
  RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("KLU2", A, X, B);

  solver->symbolicFactorization().numericFactorization().solve();


  /* Print the solution
   *
   * Should be:
   *
   *  [[1]
   *   [2]
   *   [3]
   *   [4]
   *   [5]
   *   [6]]
   */

  *fos << "Solution :" << std::endl;
  // Print the solution
  X->Print(*(fos->getOStream()));
  *fos << std::endl;

  // We are done.
  return 0;
}
