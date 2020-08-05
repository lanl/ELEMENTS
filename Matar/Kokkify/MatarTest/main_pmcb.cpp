#include <iostream>
#include <chrono>
#include "pseudo_mesh_cpu_benchmark.hpp"

int main(int argc, char** argv) {
    // The following code runs a (basic) version of John  Mcalpin's STREAM
    // benchmark on MATAR's Kokkos-specific Array and Matrix types
    size_t size_i_1D = 128;

    size_t size_i_2D = 128;
    size_t size_j_2D = 128;

    size_t size_j = 128;
    size_t size_k = 128;
    size_t size_l = 128;
    size_t size_m = 128;
    size_t size_n = 128;

    Kokkos::initialize();
    {
        using policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        using policy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
        using policy4D = Kokkos::MDRangePolicy<Kokkos::Rank<4>>;
        using policy5D = Kokkos::MDRangePolicy<Kokkos::Rank<5>>;
        using policy6D = Kokkos::MDRangePolicy<Kokkos::Rank<6>>;
      
        // Outline of what this code should do:
        // 1. Create 1D, 2D, etc. instances of MATAR's Kokkos-specific
        //    Array and Matrix types
        // 2. Make 1D, 2D, etc. Kokkos views
        // 3. Run the copy, scale, sum, and triad benchmarks on each of the
        //    above objects, and capture the results with chrono calipers
        //    and likwid markers (not sure how likwid will play with OpenMP)
        // 4. Perform 2D matrix-matrix multiplication with the 2D array and
        //    matrix objects, along with the appropriate Kokkos view, and
        //    capture the memory bandwidth information via likwid
        // 5. Run performance benchmark on a RaggedRightArrayKokkos object
        //    (follow the code in main.cpp)

        pmcb cb_1D;

        pmcb cb_2D;

        pmcb cb_3D;
        pmcb cb_4D;
        pmcb cb_5D;
        pmcb cb_6D;

    }
    Kokkos::finalize();

    return 0;
}
