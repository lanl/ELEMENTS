if (VTK_FOUND AND ENABLE_VTK)
  add_library(swage2vtk vtk_io.cpp)
  target_link_libraries(swage2vtk ${VTK_LIBRARIES})
endif()
