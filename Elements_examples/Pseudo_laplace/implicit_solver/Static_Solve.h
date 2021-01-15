#ifndef STATIC_SOLVE_H
#define STATIC_SOLVE_H  

#include "utilities.h"

using namespace utils;

namespace Static_Solve{

  void run(int argc, char *argv[]);

  void read_mesh(char *MESH);

  void init_global();

  void generate_bcs();

  void allocate_state();

  void initialize_state();

  void calculate_ref_elem();

  void apply_boundary();

  void vtk_writer();

  void ensight();

  void smooth_cells();

  void smooth_element();

  void get_nodal_jacobian();

  //sparse data structures for global mass matrix
  CArray <int> Global_Mass_Matrix_Sparse_Map;
  CArray <real_t> Global_Mass_Matrix;
  CArray <int> Global_Mass_Matrix_ids;
};

#endif // end HEADER_H
