#include "utilities.h"
#include "Simulation_Parameters.h"
#include "state.h"
#include "header.h"
#include "elements.h"
#include "swage.h"

using namespace utils;

Simulation_Parameters::Simulation_Parameters(){
  node = new node_t();
  mat_pt = new mat_pt_t();
  material = new material_t();
}

Simulation_Parameters::~Simulation_Parameters(){
  delete node;
  delete mat_pt;
  delete material;
}

void Simulation_Parameters::input(){
  // ---- time varaibles and cycle info ---- //

  //simulation spatial dimension
  num_dim = 3;

  // time step
  TFINAL = 0.5;
  
  dt = 0.001;
  dt_min = 1.e-8;
  dt_max = 1.e-2;
  dt_start = 1.e-5;
    
  cycle_stop = 100000000;

  rk_num_stages = 1;

  rk_storage = 2;

  dt_cfl = 0.5;

  //polynomial interpolation order
  p_order = 0;

  // ---- graphics information ---- //
  graphics_cyc_ival = 10;
  graphics_dt_ival  = 0.1;    // every 0.1 time units

  // ---- fill instructions and intial conditions ---- //

  NF = 2; // number of fills
    
  mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
  
  //Static isotropic parameters to move into a child class later
  Elastic_Modulus = 10000;
  Poisson_Ratio = 0;
  words_per_line = 1;
  elem_words_per_line = 8;

  //ensight file readin for Isotropical elasticity

  num_gauss_points = 2;
    
  // Global instructions
  mat_fill[0].volume = region::global;    // fill everywhere
  mat_fill[0].mat_id = 0;                 // material id
  mat_fill[0].field1 = 1.0;               // some field
  mat_fill[0].field2 = 0.0;               // some other field

    
  // Specific instructions
  mat_fill[1].volume = region::box;   // fill a sphere
  mat_fill[1].mat_id = 1;             // material id
  mat_fill[1].x1 = 0.0;
  mat_fill[1].x2 = 0.7;
  mat_fill[1].y1 = 0.0;
  mat_fill[1].y2 = 2.0;
  mat_fill[1].z1 = 0.0;
  mat_fill[1].z2 = 2.0;
  mat_fill[1].field1 = 10.0;  // some field
  mat_fill[1].field2 = 0.0;   // some other field

  // ---- boundary conditions ---- //
  NB = 2; // number of boundaries
  NBSF = 1; //number of surface density force conditions
  NBD = 1; //number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces
           //note this only implies a fixed displacement on the surface if no other basis functions have support on the surface
    
  // allocate boundary memory
  boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
    
  // Tag X=0 plane
  boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
  boundary[0].value = 0.0;
  boundary[0].thermal_bc = bdy::isothermal;
    
  // Tag X=2 plane
  boundary[1].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
  boundary[1].value = 2.0;
  boundary[1].thermal_bc = bdy::isothermal;
}