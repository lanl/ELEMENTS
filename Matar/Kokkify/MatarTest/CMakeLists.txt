cmake_minimum_required(VERSION 3.1)
project(mtest CXX C)
#project(unittest_esn CXX C)
#set(CMAKE_C_FLAGS "-g -G")
#set(CMAKE_CXX_FLAGS "-g -G")
# Enable verbose
set(CMAKE_C_FLAGS "-g -G") # VERBOSE=1")
set(CMAKE_CXX_FLAGS "-g -G") # VERBOSE=1")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

add_executable(mtest main.cpp pseudo_mesh.hpp pseudo_mesh.cpp matar.h kokkos_alias.h)
#add_executable(dummy dummy.cpp matar_test.h kokkos_alias.h)
#add_executable(mtest main.cpp pseudo_mesh.hpp pseudo_mesh.cpp pseudo_var.hpp pseudo_var.cpp matar.h)
#add_executable(unittest_esn main_unittest_esn.cpp pseudo_var.hpp pseudo_var.cpp matar.h)

include_directories(./)

add_subdirectory(Kokkos/kokkos)

target_link_libraries(mtest Kokkos::kokkos)
#target_link_libraries(dummy Kokkos::kokkos)
#target_link_libraries(unittest_esn Kokkos::kokkos)
