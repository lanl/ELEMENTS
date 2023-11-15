project(Examples)

# CMake version > 3.8 required for C++17 standard
cmake_minimum_required(VERSION 3.8)

find_package(Elements)
if(Matar_ENABLE_KOKKOS)
    find_package(Kokkos)
endif()

# "Average" example, an overview of ELEMENTS
add_subdirectory(average)

# "VTK" example, demonstration of visualization using VTK/Paraview
add_subdirectory(vtk)
