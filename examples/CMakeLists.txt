# "Average" example, an overview of ELEMENTS
add_subdirectory(average)

# "VTK" example, demonstration of visualization using VTK/Paraview
add_subdirectory(vtk)

# "hexref" example, demo of MATAR-Kokkoified hexahedral reference element
if (MATAR_WITH_KOKKOS)
  add_subdirectory(hexref)
endif()
