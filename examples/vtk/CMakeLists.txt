# Demonstrates initialization of VTK grid from SWAGE mesh and output
if (VTK_FOUND)
  add_executable(
    plotting_jacobian_determinants 
      plotting_jacobian_determinants.cpp
  )

  target_link_libraries(
    plotting_jacobian_determinants 
      common
      elements 
      swage 
      swage2vtk 
  )
endif()
