# examples/average/CMakeLists.txt

# 1. Define the executable
#    Point directly to the source file inside 'src/'
add_executable(point_cloud_example src/point_cloud_example.cpp)
add_executable(bin_mesh_test src/bin_mesh_test.cpp)
add_executable(multi_point_test src/multi_point_test.cpp)
add_executable(point_cloud_test src/point_cloud_test.cpp)
add_executable(points_in_box_test src/points_in_box_test.cpp)

# 2. Add this example's specific include path
#    This allows main.cpp to find headers in examples/point_connectivity/include/
target_include_directories(point_cloud_example PRIVATE include)
target_include_directories(bin_mesh_test PRIVATE include)
target_include_directories(multi_point_test PRIVATE include)
target_include_directories(point_cloud_test PRIVATE include)
target_include_directories(points_in_box_test PRIVATE include)

# 3. Link against the main library (ELEMENTS)
#    This pulls in Kokkos, MATAR, and the main library headers automatically.
target_link_libraries(point_cloud_example PRIVATE ELEMENTS)
target_link_libraries(bin_mesh_test PRIVATE ELEMENTS)
target_link_libraries(multi_point_test PRIVATE ELEMENTS)
target_link_libraries(point_cloud_test PRIVATE ELEMENTS)
target_link_libraries(points_in_box_test PRIVATE ELEMENTS)