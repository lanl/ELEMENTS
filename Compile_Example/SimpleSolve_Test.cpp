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
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include <set>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#include "matar.h"


int main(int argc, char *argv[]) {
  //Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  //initialize MPI
  MPI_Init(&argc,&argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Tpetra::getDefaultComm();

  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  typedef Tpetra::Details::DefaultTypes::node_type node_type;
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
  typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  typedef Kokkos::View<SizeType*, array_layout, device_type, memory_traits> row_pointers;
 
  using vec_map_type = Tpetra::Map<LO, GO>;
  using vec_device_type = typename vec_map_type::device_type;
  typedef Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
  //typedef Tpetra::dual_view_type::t_dev vec_array;
  

  // Before we do anything, check that KLU2 is enabled
  if( !Amesos2::query("KLU2") ){
    std::cerr << "KLU2 not enabled in this run.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }


  size_t myRank = comm->getRank();

  std::ostream &out = std::cout;

  out << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;
  
  //construct global data for example; normally this would be from file input and distributed
  //according to the row map at that point
  global_size_t nrows = 6;
  int count_data[nrows];
  double rhs_data[nrows];
  int index_data[nrows][nrows];
  double value_data[nrows][nrows];
  count_data[0] = 3; count_data[1] = 2; count_data[2] = 1;
  count_data[3] = 2; count_data[4] = 2; count_data[5] = 2;

  index_data[0][0] = 0; index_data[0][1] = 2; index_data[0][2] = 4;
  index_data[1][0] = 0; index_data[1][1] = 1;
  index_data[2][0] = 2;
  index_data[3][0] = 0; index_data[3][1] = 3;
  index_data[4][0] = 1; index_data[4][1] = 4;
  index_data[5][0] = 3; index_data[5][1] = 5;

  value_data[0][0] = 7; value_data[0][1] = -3; value_data[0][2] = -1;
  value_data[1][0] = 2; value_data[1][1] = 8;
  value_data[2][0] = 1;
  value_data[3][0] = -3; value_data[3][1] = 5;
  value_data[4][0] = -1; value_data[4][1] = 4;
  value_data[5][0] = -2; value_data[5][1] = 6;

  rhs_data[0] = -7;
  rhs_data[1] = 18;
  rhs_data[2] = 3;
  rhs_data[3] = 17;
  rhs_data[4] = 18;
  rhs_data[5] = 28;
  
  // create a Map
  RCP<Tpetra::Map<LO,GO,node_type> > map
    = rcp( new Tpetra::Map<LO,GO,node_type>(nrows,0,comm) );

  global_size_t local_nrows = map->getNodeNumElements();
  global_size_t *entries_per_row = new global_size_t[local_nrows];
  row_pointers row_offsets = row_pointers("row_offsets", local_nrows+1);
  
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();

  //init row_offsets
  for(int i=0; i < local_nrows+1; i++)
    row_offsets(i) = 0;

  for(LO i=0; i < local_nrows; i++){
    entries_per_row[i] = count_data[map->getGlobalElement(i)];
    row_offsets(i+1) = entries_per_row[i] + row_offsets(i);
  }
  
  //row_offsets(1) = 3;
  //entries_per_row[1] = 2;
  //row_offsets(2) = 5;
  //entries_per_row[2] = 1;
  //row_offsets(3) = 6;
  //entries_per_row[3] = 2;
  //row_offsets(4) = 8;
  //entries_per_row[4] = 2;
  //row_offsets(5) = 10;
  //entries_per_row[5] = 2;
  //row_offsets(6) = 12;

  global_size_t nnz = row_offsets(local_nrows);

  //indices_array all_indices = indices_array("indices_array",nnz);
  //values_array all_values = values_array("values_array",nnz);
  CArrayKokkos<GO, array_layout, device_type, memory_traits> all_global_indices(nnz);
  CArrayKokkos<LO, array_layout, device_type, memory_traits> all_indices(nnz);
  CArrayKokkos<Scalar, Kokkos::LayoutRight, device_type, memory_traits> all_values(nnz);
  CArrayKokkos<Scalar, Kokkos::LayoutLeft, device_type, memory_traits> Bview(local_nrows);
  CArrayKokkos<Scalar, Kokkos::LayoutLeft, device_type, memory_traits> Xview(local_nrows);

  //set Kokkos view data
  LO entrycount = 0;
  for(LO i=0; i < local_nrows; i++){
    for(LO j=0; j < entries_per_row[i]; j++){
    all_global_indices(entrycount) = index_data[map->getGlobalElement(i)][j];
    all_values(entrycount) = value_data[map->getGlobalElement(i)][j];
    entrycount++;
    }
  }
  
  //build column map
  RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const RCP<const Tpetra::Map<LO,GO,node_type> > dommap = map;

  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,all_global_indices.get_kokkos_view(), nullptr);

  //convert global indices to local indices using column map
  entrycount = 0;
  for(LO i=0; i < local_nrows; i++){
    for(LO j=0; j < entries_per_row[i]; j++){
    all_indices(entrycount) = colmap->getLocalElement(all_global_indices(entrycount));
    entrycount++;
    }
  }

  RCP<MAT> A = rcp( new MAT(map, colmap, row_offsets, all_indices.get_kokkos_view(), all_values.get_kokkos_view()) ); // max of three entries in a row

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
  vec_array Xview_pass = vec_array("Xview_pass", local_nrows,1);
  Xview_pass.assign_data(Xview.pointer());
  RCP<MV> X = rcp(new MV(map, Xview_pass));
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
  vec_array Bview_pass = vec_array("Bview", local_nrows,1);
  for(LO i=0; i < local_nrows; i++)
    Bview(i) = rhs_data[map->getGlobalElement(i)];
  
  Bview_pass.assign_data(Bview.pointer());

  RCP<MV> B = rcp(new MV(map, Bview_pass));

  //for( int i = 0; i < 6; ++i ){
    //if( B->getMap()->isNodeGlobalElement(i) ){
      //B->replaceGlobalValue(i,0,data[i]);
    //}
  //}

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

  //*fos << "Xview :" << std::endl;
  //for(int print = 0; print < local_nrows; print++)
    //*fos << Xview(print) << "; ";
  //*fos << std::endl;

  // We are done.
  return 0;
}
