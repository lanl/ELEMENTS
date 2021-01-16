add_library(elements elements.cpp elements.h)
add_library(swage swage.cpp swage.h)

include_directories(MATAR)
add_subdirectory(MATAR)

target_link_libraries(elements matar)
target_link_libraries(swage matar)

if (KOKKOS)
  add_definitions(-DHAVE_KOKKOS=1)
  #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)
  if (CUDA)
    add_definitions(-DHAVE_CUDA=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
  elseif (OPENMP)
    add_definitions(-DHAVE_OPENMP=1)
    #set_target_properties(matar PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
  endif()

  #target_link_libraries(matar INTERFACE Kokkos::kokkos)
  target_link_libraries(elements  Kokkos::kokkos)
  target_link_libraries(swage  Kokkos::kokkos)
endif()
