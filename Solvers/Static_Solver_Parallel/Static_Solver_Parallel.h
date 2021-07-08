#ifndef STATIC_SOLVER_PARALLEL_H
#define STATIC_SOLVER_PARALLEL_H  

#include "utilities.h"
#include "../Solver.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
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

namespace swage{
  class mesh_t;
}

namespace elements{
  class element_selector;
  class Element3D;
  class Element2D;
  class ref_element;
}

class Static_Solver_Parallel: public Solver{

public:
  Static_Solver_Parallel();
  ~Static_Solver_Parallel();

  //Trilinos type definitions
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<real_t,LO,GO> MAT;
  typedef const Tpetra::CrsMatrix<real_t,LO,GO> const_MAT;
  typedef Tpetra::MultiVector<real_t,LO,GO> MV;

  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  typedef Tpetra::Details::DefaultTypes::node_type node_type;
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  using global_size_t = Tpetra::global_size_t;
  
  typedef Kokkos::View<real_t*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
  typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  typedef Kokkos::View<SizeType*, array_layout, device_type, memory_traits> row_pointers;
 
  using vec_map_type = Tpetra::Map<LO, GO>;
  using vec_device_type = typename vec_map_type::device_type;
  //typedef Kokkos::DualView<real_t**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
  typedef MV::dual_view_type::t_dev vec_array;
  typedef MV::dual_view_type::t_host host_vec_array;
  typedef MV::dual_view_type dual_vec_array;

  void run(int argc, char *argv[]);

  void read_mesh(char *MESH);

  void init_global();

  void assemble();

  int solve();

  void local_matrix(int ielem, CArray <real_t> &Local_Matrix);

  void local_matrix_multiply(int ielem, CArray <real_t> &Local_Matrix);

  void generate_bcs();

  void Get_Boundary_Patches();

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

  void init_boundary_sets(int num_boundary_sets);

  void tag_boundaries(int this_bc_tag, real_t val, int bdy_set);

  int check_boundary(Node_Combination &Patch_Nodes, int this_bc_tag, real_t val);\

  //debug and performance functions/variables
  double CPU_Time();
  void init_clock();
  double initial_CPU_time;
  
  swage::mesh_t *init_mesh;
  swage::mesh_t *mesh;
  
  elements::element_selector *element_select;
  elements::Element3D *elem;
  elements::Element2D *elem2D;
  elements::ref_element  *ref_elem;
  

  class Simulation_Parameters *simparam;
  
  //Local FEA data
  size_t nlocal_nodes;
  dual_vec_array dual_node_data; //first three indices of second dim should be positions
  dual_vec_array dual_nodal_forces;
  CArray<elements::elem_types::elem_type> Element_Types;
  CArray<size_t> Nodes_Per_Element_Type;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Stiffness_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Stiffness_Matrix;
  //CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Forces;
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Results; //result of linear solve; typically displacements and densities
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Stiffness_Matrix_Strides;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;

  //Ghost data on this MPI rank
  size_t nghost_nodes;
  CArrayKokkos<GO, Kokkos::LayoutLeft, node_type::device_type> ghost_nodes;
  CArrayKokkos<int, array_layout, device_type, memory_traits> ghost_node_ranks;
  Teuchos::ArrayView<const GO> ghost_nodes_pass;
  Teuchos::ArrayView<int> ghost_node_ranks_pass;

  //Local FEA data including ghosts
  size_t nall_nodes;
  size_t rnum_elem;
  dual_vec_array dual_all_node_data; //first three indices of second dim should be positions

  //Global FEA data
  size_t num_nodes, num_elem;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > map; //map of node indices
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_node_map; //map of node indices with ghosts on each rank
  Teuchos::RCP<MV> node_data_distributed;
  Teuchos::RCP<MV> all_node_data_distributed;
  Teuchos::RCP<MAT> Global_Stiffness_Matrix;
  Teuchos::RCP<MV> Global_Nodal_Forces;
  Teuchos::RCP<MAT> A;
  Teuchos::RCP<MV> B;
  Teuchos::RCP<MV> X;
  //Parallel map for the global set of nodes (before removing BCS)
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_dof_map; //map of local dofs (typically num_node_local*num_dim)
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_dof_map; //map of local and ghost dofs (typically num_node_all*num_dim)
  
  //Boundary Conditions Data
  //CArray <Nodal_Combination> Patch_Nodes;
  size_t nboundary_patches;
  size_t num_boundary_conditions;
  CArrayKokkos<Node_Combination, array_layout, device_type, memory_traits> Boundary_Patches;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches; //set of patches corresponding to each boundary condition
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> NBoundary_Condition_Patches;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches_strides;


  //element selection parameters and data
  size_t max_nodes_per_element;

  //types of boundary conditions
  enum bc_type {NONE,DISPLACEMENT_CONDITION, X_DISPLACEMENT_CONDITION,
   Y_DISPLACEMENT_CONDITION, Z_DISPLACEMENT_CONDITION, LOADING_CONDITION};

  //lists what kind of boundary condition the nodal DOF is subjected to if any
  CArrayKokkos<int, array_layout, device_type, memory_traits> Node_DOF_Boundary_Condition_Type;
  //stores the displacement value for the boundary condition on this nodal DOF
  CArray<real_t> Node_DOF_Displacement_Boundary_Conditions;
  //stores applied point forces on nodal DOF
  CArray<real_t> Node_DOF_Force_Boundary_Conditions;
  //lists what kind of boundary condition each boundary set is assigned to
  CArray<int> Boundary_Condition_Type_List;
  //constant surface force densities corresponding to each boundary set (provide varying field later)
  CArray<real_t> Boundary_Surface_Force_Densities;
  //constant displacement condition applied to all nodes on a boundary surface (convenient option to avoid specifying nodes)
  CArray<real_t> Boundary_Surface_Displacements;
  
  //number of displacement boundary conditions acting on nodes; used to size the reduced global stiffness map
  size_t Number_DOF_BCS;

  //MPI data
  int myrank; //index of this mpi rank in the world communicator
  int nranks; //number of mpi ranks in the world communicator
  MPI_Comm world; //stores the default communicator object (MPI_COMM_WORLD)

  //! mapping used to get local ghost index from the global ID.
  typedef ::Tpetra::Details::FixedHashTable<GO, LO, Kokkos::HostSpace::device_type>
    global_to_local_table_host_type;

  global_to_local_table_host_type global2local_map;
  CArrayKokkos<int, Kokkos::LayoutLeft, Kokkos::HostSpace::device_type> active_ranks;

  //Pertains to local mesh information being stored as prescribed by the row map
  global_size_t local_nrows;
  global_size_t min_gid;
  global_size_t max_gid;
  global_size_t index_base;

  //file readin variables
  std::ifstream *in;
  int words_per_line, elem_words_per_line;
  
};

#endif // end HEADER_H
