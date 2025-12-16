# examples/average/CMakeLists.txt

# 1. Define the executable
#    Point directly to the source file inside 'src/'
add_executable(average_example src/average.cpp)

# 2. Add this example's specific include path
#    This allows main.cpp to find headers in examples/average/include/
target_include_directories(average_example PRIVATE include)

# 3. Link against the main library (ELEMENTS)
#    This pulls in Kokkos, MATAR, and the main library headers automatically.
target_link_libraries(average_example PRIVATE ELEMENTS)