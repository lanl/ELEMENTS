#pragma once

#include "common/utilities.h"
#include "common/error.h"
#include "common/common.h"

#include "swage/swage.h"

typedef swage::mesh_t SwageMesh;

namespace swage2vtk {
  void init_swage_mesh_from_vtk_grid_file(const std::string &filename, 
      const int &elem_order, SwageMesh &mesh);

  void write_vtk_grid_file_from_swage_mesh(SwageMesh &mesh, 
      const std::string &solution_name, const CArray<double> &solution, 
      std::string &filename);
}
