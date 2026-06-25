# examples/mesh_form/CMakeLists.txt

# 1. Define the executable
#    Point directly to the source file inside 'src/'
add_executable(mesh_form src/mesh_form.cpp)

# 2. Add this example's specific include path
#    This allows main.cpp to find headers in examples/mesh_form/include/
target_include_directories(mesh_form PRIVATE include)

# 3. Link against the main library (ELEMENTS)
#    This pulls in Kokkos, MATAR, and the main library headers automatically.
target_link_libraries(mesh_form PRIVATE ELEMENTS)