# examples/average/CMakeLists.txt

# 1. Define the executable
#    Point directly to the source file inside 'src/'
add_executable(mesh_decomp src/mesh_decomp_example.cpp)

# 2. Add this example's specific include path
#    This allows main.cpp to find headers in examples/average/include/
target_include_directories(mesh_decomp PRIVATE include)

# 3. Link against the main library (ELEMENTS)
#    This pulls in Kokkos, MATAR, and the main library headers automatically.
target_link_libraries(mesh_decomp PRIVATE ELEMENTS)