# Includes
include_directories(include)

set(Average_SRC_CXX src/average.cpp)

add_executable(Average ${Average_SRC_CXX})

# Link libraries including MPI (required for mesh_io.h which includes mpi.h)
# and Kokkos (required for Kokkos macros and types in swage/mesh.h)
target_link_libraries(Average 
    PRIVATE
        Elements
        input
        matar
        MPI::MPI_CXX
        Kokkos::kokkos
)
