#ifndef STATIC_SOLVER_H
#define STATIC_SOLVER_H  

#include "utilities.h"
#include "../Solver.h"

namespace swage{
  class mesh_t;
}

namespace elements{
  class element_selector;
  class Element3D;
  class Element2D;
  class ref_element;
}

class Static_Solver: public Solver{

public:
  Static_Solver();
  ~Static_Solver();

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
  
  swage::mesh_t *init_mesh;
  swage::mesh_t *mesh;
  
  elements::element_selector *element_select;
  elements::Element3D *elem;
  elements::Element2D *elem2D;
  elements::ref_element  *ref_elem;
  

  class Simulation_Parameters *simparam;
  
};

#endif // end HEADER_H
