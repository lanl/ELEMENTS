# examples/average/CMakeLists.txt

# 1. Define the executable
#    Point directly to the source file inside 'src/'
add_executable(point_connectivity src/point_connectivity.cpp)

# 2. Add this example's specific include path
#    This allows main.cpp to find headers in examples/point_connectivity/include/
target_include_directories(point_connectivity PRIVATE include)

# 3. Link against the main library (ELEMENTS)
#    This pulls in Kokkos, MATAR, and the main library headers automatically.
target_link_libraries(point_connectivity PRIVATE ELEMENTS)