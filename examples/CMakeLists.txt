find_package(Elements)
if(Matar_ENABLE_KOKKOS)
    find_package(Kokkos)
endif()

# "Average" example, an overview of ELEMENTS
add_subdirectory(average)

# "VTK" example, demonstration of visualization using VTK/Paraview
add_subdirectory(vtk)
