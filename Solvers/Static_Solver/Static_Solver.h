#ifndef STATIC_SOLVER_H
#define STATIC_SOLVER_H  

#include "utilities.h"
#include "../Solver.h"
#include "matar.h"

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

  void assemble();

  int solve();

  void local_matrix(int ielem, CArray <real_t> &Local_Matrix);

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

  void Element_Material_Properties(size_t, real_t &Element_Modulus, real_t &Poisson_Ratio);

  void Force_Vector_Construct();
  
  swage::mesh_t *init_mesh;
  swage::mesh_t *mesh;
  
  elements::element_selector *element_select;
  elements::Element3D *elem;
  elements::Element2D *elem2D;
  elements::ref_element  *ref_elem;
  

  class Simulation_Parameters *simparam;
  
  CArray <size_t> Element_Types;
  CArray <size_t> Global_Mass_Matrix_Assembly_Map;
  RaggedRightArray <size_t> Graph_Matrix;
  RaggedRightArray <size_t> DOF_Graph_Matrix;
  RaggedRightArray <real_t> Mass_Matrix;
  CArray <real_t> Nodal_Forces;
  CArray <size_t> Mass_Matrix_strides;
  CArray <size_t> Graph_Matrix_strides;
  
};

#endif // end HEADER_H
