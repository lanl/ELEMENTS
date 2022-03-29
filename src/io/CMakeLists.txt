if (WITH_VTK AND VTK_FOUND)
  add_library(swage2vtk vtk_io.cpp)
  target_link_libraries(swage2vtk ${VTK_LIBRARIES})

  target_include_directories(
    swage2vtk
      PUBLIC
        ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
        ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
  )
endif()
