# examples/average/CMakeLists.txt

# 1. Define the executable
#    Point directly to the source file inside 'src/'
add_executable(gradient_test src/gradient_test.cpp)
add_executable(interpolation_test src/interpolation_test.cpp)
add_executable(kronecker_delta_test src/kronecker_delta_test.cpp)
add_executable(partition_unity_test src/partition_unity.cpp)

# 2. Add this example's specific include path
#    This allows main.cpp to find headers in examples/point_connectivity/include/
target_include_directories(gradient_test PRIVATE include)
target_include_directories(interpolation_test PRIVATE include)
target_include_directories(kronecker_delta_test PRIVATE include)
target_include_directories(partition_unity_test PRIVATE include)

# 3. Link against the main library (ELEMENTS)
#    This pulls in Kokkos, MATAR, and the main library headers automatically.
target_link_libraries(gradient_test PRIVATE ELEMENTS)
target_link_libraries(interpolation_test PRIVATE ELEMENTS)
target_link_libraries(kronecker_delta_test PRIVATE ELEMENTS)
target_link_libraries(partition_unity_test PRIVATE ELEMENTS)