add_library(
  swage 
  swage.cpp 
  refine_high_order_mesh.cpp
  evaluate_jacobian_determinants.cpp
)
target_link_libraries(
  swage
  elements
)

target_include_directories(swage PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}
)
