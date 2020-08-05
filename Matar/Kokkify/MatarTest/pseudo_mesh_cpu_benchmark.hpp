#ifndef PSEUDO_MESH_CPU_BENCHMARK
#define PSEUDO_MESH_CPU_BENCHMARK

#include "matar.h"

class pmcb {
    public:
        size_t size1;
        size_t size2;
        size_t size3;
        size_t size4;
        size_t size5;
        size_t size6;

        FArrayKokkos<real_t>  fak;
        FMatrixKokkos<real_t> fmk;
        CArrayKokkos<real_t>  cak;
        CMatrixKokkos<real_t> cmk;

        // Default constructor
        pmcb();

        // 1D initialization function
        void init(size_t dim1);

        // 2D initialization function
        void init(size_t dim1, size_t dim2);

        // 3D initialization function
        void init(size_t dim1, size_t dim2, size_t dim3);

        // 4D initialization function
        void init(size_t dim1, size_t dim2, size_t dim3, size_t dim4);

        // 5D initialization function
        void init(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                  size_t dim5);

        // 6D initialization function
        void init(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                  size_t dim5, size_t dim6);

        // Destructor
        ~pmcb();
};

#endif
