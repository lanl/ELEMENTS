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

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Kokkos::View<Scalar*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  typedef Kokkos::View<SizeType*, array_layout, device_type, memory_traits> row_pointers;
 
  using vec_map_type = Tpetra::Map<LO, GO>;
  using vec_device_type = typename vec_map_type::device_type;
  typedef Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
  //typedef Tpetra::dual_view_type::t_dev vec_array;

  


  // Before we do anything, check that SuperLU is enabled
  if( !Amesos2::query("KLU2") ){
    std::cerr << "KLU2 not enabled in this run.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Tpetra::getDefaultComm();

  size_t myRank = comm->getRank();

  std::ostream &out = std::cout;

  out << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  // create a Map
  global_size_t nrows = 6;
  global_size_t *entries_per_row = new global_size_t[nrows];
  row_pointers row_offsets = row_pointers("row_offsets", nrows+1);
  
  row_offsets(0) = 0;
  entries_per_row[0] = 3;
  row_offsets(1) = 3;
  entries_per_row[1] = 2;
  row_offsets(2) = 5;
  entries_per_row[2] = 1;
  row_offsets(3) = 6;
  entries_per_row[3] = 2;
  row_offsets(4) = 8;
  entries_per_row[4] = 2;
  row_offsets(5) = 10;
  entries_per_row[5] = 2;
  row_offsets(6) = 12;

  global_size_t nnz = 0;

  //count number of non-zeros
  for(int i = 0; i < nrows; i++) nnz+=row_offsets(i);

  indices_array all_indices = indices_array("indices_array",nnz);
  values_array all_values = values_array("values_array",nnz);

  //set Kokkos view data

  all_indices(0) = 0; all_indices(1) = 2; all_indices(2) = 4;
  all_indices(3) = 0; all_indices(4) = 1;
  all_indices(5) = 2;
  all_indices(6) = 0; all_indices(7) = 3;
  all_indices(8) = 1; all_indices(9) = 4;
  all_indices(10) = 3; all_indices(11) = 5;

  all_values(0) = 7; all_values(1) = -3; all_values(2) = -1;
  all_values(3) = 2; all_values(4) = 8;
  all_values(5) = 1;
  all_values(6) = -3; all_values(7) = 5;
  all_values(8) = -1; all_values(9) = 4;
  all_values(10) = -2; all_values(11) = 6;

  RCP<Tpetra::Map<LO,GO> > map
    = rcp( new Tpetra::Map<LO,GO>(nrows,0,comm) );

  RCP<MAT> A = rcp( new MAT(map, map, row_offsets, all_indices, all_values) ); // max of three entries in a row

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
   RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  // Construct matrix
  if( myRank == 0 ){
    //A->insertGlobalValues(0,tuple<GO>(0,2,4),tuple<Scalar>(7,-3,-1));
    //A->insertGlobalValues(1,tuple<GO>(0,1),tuple<Scalar>(2,8));
    //A->insertGlobalValues(2,tuple<GO>(2),tuple<Scalar>(1));
    //A->insertGlobalValues(3,tuple<GO>(0,3),tuple<Scalar>(-3,5));
    //A->insertGlobalValues(4,tuple<GO>(1,4),tuple<Scalar>(-1,4));
    //A->insertGlobalValues(5,tuple<GO>(3,5),tuple<Scalar>(-2,6));
  }
  A->fillComplete();
  
  *fos << "Matrix :" << std::endl;
  A->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;

  // Create random X
  vec_array Xview = vec_array("Xview", nrows,1);
  RCP<MV> X = rcp(new MV(map, Xview));
  X->randomize();

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
  vec_array Bview = vec_array("Bview", nrows,1);
  RCP<MV> B = rcp(new MV(map, Bview));
  int data[6] = {-7,18,3,17,18,28};
  for( int i = 0; i < 6; ++i ){
    if( B->getMap()->isNodeGlobalElement(i) ){
      B->replaceGlobalValue(i,0,data[i]);
    }
  }

  *fos << "RHS :" << std::endl;
  B->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;

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
  X->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;

  // We are done.
  return 0;
}
