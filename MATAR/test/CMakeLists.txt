
add_executable(mtest main.cpp)

target_link_libraries(mtest matar)

include_directories(virtualFcnKokkos)
add_subdirectory(virtualFcnKokkos)
add_definitions(-DHAVE_KOKKOS=1 -DHAVE_OPENMP=1)

