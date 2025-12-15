# Examples directory
# Build examples when ELEMENTS_BUILD_EXAMPLES is enabled

if(ELEMENTS_BUILD_EXAMPLES)

  # Decomp example
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/decomp_example/mesh_decomp_example.cpp)
    add_executable(mesh_decomp_example
      decomp_example/mesh_decomp_example.cpp
    )

    # Include directories for example headers
    # Note: The example includes headers like "mesh.h", "state.h", "mesh_io.h"
    # which are in the decomp_example directory, and also needs access to
    # library headers from decomp, swage, elements, and matar
    target_include_directories(mesh_decomp_example
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/decomp_example
        ${CMAKE_SOURCE_DIR}/decomp
        ${CMAKE_SOURCE_DIR}/swage
        ${CMAKE_SOURCE_DIR}/elements
        ${CMAKE_SOURCE_DIR}/matar/src/include
    )

    # Link against ELEMENTS libraries
    # These INTERFACE libraries will propagate include directories and dependencies
    target_link_libraries(mesh_decomp_example
      PRIVATE
        Elements::decomp
        Elements::swage
        Elements::elements
        matar
    )

    # MPI is required for this example
    if(ELEMENTS_ENABLE_MPI)
      target_link_libraries(mesh_decomp_example
        PRIVATE
          MPI::MPI_CXX
      )
    else()
      message(WARNING "mesh_decomp_example requires MPI. Set ELEMENTS_ENABLE_MPI=ON")
      set_target_properties(mesh_decomp_example PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endif()

    # Kokkos is required (MATAR needs it)
    if(ELEMENTS_ENABLE_KOKKOS)
      target_link_libraries(mesh_decomp_example
        PRIVATE
          Kokkos::kokkos
      )
    else()
      message(WARNING "mesh_decomp_example may require Kokkos. Set ELEMENTS_ENABLE_KOKKOS=ON")
    endif()

    # Install example executable
    install(TARGETS mesh_decomp_example
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
  endif()

endif()

