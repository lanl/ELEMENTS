#ifndef MATH_UTILITY_H
#define MATH_UTILITY_H  

#include "utilities.h"

using namespace utils;

namespace StaticSolve{

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

  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh state declarations ---
  node_t          node;
  mat_pt_t        mat_pt;


  material_t  * material;
  mat_fill_t  * mat_fill;
  boundary_t  * boundary;


  // --- Mesh regions and material fills ---
  int NR = 0; // number of Regions
  int NC = 0; // number of contours
  int NF = 0; // number of fill
  int NB = 0; // number of boundary patch sets to tag


  // --- Graphics output variables ---
  int graphics_id = 0;
  int graphics_cyc_ival = 0;

  real_t graphics_times[250];
  real_t graphics_dt_ival = 1.0e8;
  real_t graphics_time = graphics_dt_ival;  // the times for writing graphics dump


  // --- Time and cycling variables ---
  real_t TIME = 0.0;
  real_t TFINAL = 1.e16;
  real_t dt = 1.e-8;
  real_t dt_max = 1.0e-2;
  real_t dt_min = 1.0e-8;
  real_t dt_cfl = 0.3;
  real_t dt_start = 1.0e-8;

  int rk_num_stages = 1;
  int rk_storage = 1;
  int rk_stage = 0;

  int cycle = 0;
  int cycle_stop = 1000000000;
  int stop_calc = 0;    // a flag to end the calculation when = 1

  real_t percent_comp = 0.0;

  // --- Precision variables ---
  real_t fuzz = 1.0e-16;  // machine precision
  real_t tiny = 1.0e-12;  // very very small (between real_t and single)
  real_t small= 1.0e-8;   // single precision


  // --- Dimensional and mesh constants ---
  int num_dim = 3;
  int p_order = 0;

  //sparse data structures for global mass matrix
  auto Global_Mass_Matrix_Sparse_Map;
  auto Global_Mass_Matrix;
  auto Global_Mass_Matrix_ids;
};

#endif // end HEADER_H
