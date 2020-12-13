cmake_minimum_required (VERSION 3.2)
project (Laplace)

set (CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(elements)


add_subdirectory(elements)
add_subdirectory(src)
add_subdirectory(driver)
add_subdirectory(implicit_solver)