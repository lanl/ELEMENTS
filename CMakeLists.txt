cmake_minimum_required (VERSION 3.2)
project (Elements)

set (CMAKE_CXX_STANDARD 14)

include_directories(Include)
include_directories(Elements)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(Elements)
add_subdirectory(Solvers)