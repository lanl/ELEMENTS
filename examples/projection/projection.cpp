#include "elements.h"
#include "vtk_io.h"
#include "error.h"
#include "common.h"

void project_function_lagrange(
    std::function<Real(Real, Real, Real)> func, 
    SwageMesh &mesh, CArray<Real> &projection);

void project_function_legendre(
    std::function<Real(Real, Real, Real)> func, 
    SwageMesh &mesh, CArray<Real> &projection);

//ErrorCode initialize_solution(
//    std::function<Real(Real, Real, Real)> init_cond, 
//    SwageMesh &mesh, CArray<Real> &solution);

// Definition of global element variable to be able to safely link elements
// library
//elements::HexN elem;

int main() {
  const SizeType num_dim  = 3;
  SizeType elem_order     = 1; 
  std::string input_name("simple_2nd_order_mesh.vtk");
  std::string solution_name("u");
  std::string output_name("output.vtu");

  // Try to read input VTK grid and initialize a SWAGE mesh from it
  SwageMesh initial_mesh;
  try {
    VtkGrid input_grid = swage2vtk::read_grid(input_name);
    swage2vtk::init_swage_mesh(elem_order, input_grid, initial_mesh);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  // Interpolate quadrature points and evaluate Jacobian determinants
  SwageMesh mesh;
  refine_high_order_mesh(initial_mesh, mesh);
  evaluate_jacobian_determinants(mesh);

  // Project a function onto the mesh
  SizeType num_verts_1d = elem_order + 1;
  SizeType num_verts    = pow(num_verts_1d, num_dim);
  SizeType num_elems    = mesh.num_elems();
  CArray<Real> projection(num_elems, num_verts);
  auto func = [](Real x, Real y, Real z) -> Real {
    //return cos(M_PI*x);
    return 1.0;
  };

  try {
    project_function_lagrange(func, mesh, projection);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  // Try to initialize a VTK grid from the SWAGE mesh and the projection data
  try {
    VtkGrid output_grid = swage2vtk::init_vtk_grid(initial_mesh, solution_name, 
        projection);
    swage2vtk::write_grid(output_grid, output_name);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
