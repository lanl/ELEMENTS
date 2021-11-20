# Demonstrates: 
# (1) reading input VTK grid from file; 
# (2) initializing SWAGE mesh from VTK grid; 
# (3) projection of a function onto the mesh; 
# (4) initialization of an output VTK grid file from the SWAGE mesh and the
# projection data; and 
# (5) writing output VTK grid to file.
if (ENABLE_VTK AND VTK_FOUND AND ENABLE_BLAS_LAPACK)
  add_executable(
    projection
    projection.cpp
    project_function.cpp
  )
  target_link_libraries(
    projection
    elements 
    swage 
    swage2vtk 
    common
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  )

  file(GLOB GRIDS *.vtk)
  file(COPY ${GRIDS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
