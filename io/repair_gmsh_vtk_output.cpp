#include "vtk_io.h"

#include "pointsGenerators.h"

/*
 * Generate a map from the local vertex ordering in Gmsh to the local vertex
 * ordering in VTK for high-order elements
 */
void generate_gmsh_map(int elem_order, int *map, bool reverse) {
  int Nvd = Np + 1;

  fullMatrix<RealScalar> ref = gmshGenerateMonomialsHexahedron(Np);

  for(int msh_vid = 0; msh_vid < ref.size1(); ++msh_vid) {
    const RealScalar u = ref(msh_vid, 0);
    const RealScalar v = ref(msh_vid, 1);
    const RealScalar w = ref(msh_vid, 2);

    int I = int(std::round(u));
    int J = int(std::round(v));
    int K = int(std::round(w));

    int swg_vid;
    convert_ijk_coords_to_local_id(Nvd, I, J, K, swg_vid);

    if (!reverse) {
      map[swg_vid] = msh_vid;
    } else {
      map[msh_vid] = swg_vid;
    }
  }

  return 0;
}


      // Loop over Gmsh vertex IDs 
      for (int msh_vid = 0; msh_vid < num_verts_per_elem; msh_vid++) {
        int global_vert_id = global_vert_ids->GetId(msh_vid);
        int swg_vert_id = gmsh_reverse_map[msh_vid];
        mesh.nodes_in_elem(elem_id, swg_vert_id) = global_vert_id;
      }

  int *gmsh_reverse_map = new int[Nve];
  generate_gmsh_map(elem_order, gmsh_reverse_map, true);
