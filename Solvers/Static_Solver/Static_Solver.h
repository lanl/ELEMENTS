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

  void Displacement_Boundary_Conditions();

  void Force_Vector_Construct();
  
  swage::mesh_t *init_mesh;
  swage::mesh_t *mesh;
  
  elements::element_selector *element_select;
  elements::Element3D *elem;
  elements::Element2D *elem2D;
  elements::ref_element  *ref_elem;
  

  class Simulation_Parameters *simparam;
  
  CArray <size_t> Element_Types;
  CArray <size_t> Global_Stiffness_Matrix_Assembly_Map;
  RaggedRightArray <size_t> Graph_Matrix;
  RaggedRightArray <size_t> DOF_Graph_Matrix;
  RaggedRightArray <real_t> Stiffness_Matrix;
  CArray <real_t> Nodal_Forces;
  CArray <size_t> Stiffness_Matrix_strides;
  CArray <size_t> Graph_Matrix_strides;

  //types of boundary conditions
  enum bc_type {NONE,DISPLACEMENT_CONDITION, X_DISPLACEMENT_CONDITION,
   Y_DISPLACEMENT_CONDITION, Z_DISPLACEMENT_CONDITION, LOADING_CONDITION};

  //lists what kind of boundary condition the nodal DOF is subjected to if any
  CArray <int> Node_DOF_Boundary_Condition_Type;
  //stores the displacement value for the boundary condition on this nodal DOF
  CArray <real_t> Node_DOF_Displacement_Boundary_Conditions;
  //stores applied point forces on nodal DOF
  CArray <real_t> Node_DOF_Force_Boundary_Conditions;
  //lists what kind of boundary condition each boundary set is assigned to
  CArray <int> Boundary_Condition_Type_List;
  //constant surface force densities corresponding to each boundary set (provide varying field later)
  CArray <real_t> Boundary_Surface_Force_Densities;
  //constant displacement condition applied to all nodes on a boundary surface (convenient option to avoid specifying nodes)
  CArray <real_t> Boundary_Surface_Displacements;
  
  //number of displacement boundary conditions acting on nodes; used to size the reduced global stiffness map
  size_t Number_DOF_BCS;
  
};

#endif // end HEADER_H
