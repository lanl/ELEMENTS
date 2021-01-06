cmake_minimum_required (VERSION 3.2)
project (Laplace)

set (CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(elements)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(elements)
add_subdirectory(src)
add_subdirectory(driver)
add_subdirectory(implicit_solver)