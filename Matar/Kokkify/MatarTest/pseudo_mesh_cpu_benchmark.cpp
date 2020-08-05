#include "pseudo_mesh_cpu_benchmark.hpp"

// Default constructor
pmcb::pmcb() {}

// 1D initialization function
void pmcb::init(size_t dim1) {
    size1 = dim1;

    cak = CArrayKokkos<real_t>  (size1);
    cmk = CMatrixKokkos<real_t> (size1);
    fak = FArrayKokkos<real_t>  (size1);
    fmk = FMatrixKokkos<real_t> (size1);
}

// 2D initialization function
void pmcb::init(size_t dim1, size_t dim2) {
    size1 = dim1;
    size2 = dim2;

    cak = CArrayKokkos<real_t>  (size1, size2);
    cmk = CMatrixKokkos<real_t> (size1, size2);
    fak = FArrayKokkos<real_t>  (size1, size2);
    fmk = FMatrixKokkos<real_t> (size1, size2);
}

// 3D initialization function
void pmcb::init(size_t dim1, size_t dim2, size_t dim3) {
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;

    cak = CArrayKokkos<real_t>  (size1, size2, size3);
    cmk = CMatrixKokkos<real_t> (size1, size2, size3);
    fak = FArrayKokkos<real_t>  (size1, size2, size3);
    fmk = FMatrixKokkos<real_t> (size1, size2, size3);
}

// 4D initialization function
void pmcb::init(size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    size4 = dim4;

    cak = CArrayKokkos<real_t>  (size1, size2, size3, size4);
    cmk = CMatrixKokkos<real_t> (size1, size2, size3, size4);
    fak = FArrayKokkos<real_t>  (size1, size2, size3, size4);
    fmk = FMatrixKokkos<real_t> (size1, size2, size3, size4);
}

// 5D initialization function
void pmcb::init(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                size_t dim5) {
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    size4 = dim4;
    size5 = dim5;

    cak = CArrayKokkos<real_t>  (size1, size2, size3, size4, size5);
    cmk = CMatrixKokkos<real_t> (size1, size2, size3, size4, size5);
    fak = FArrayKokkos<real_t>  (size1, size2, size3, size4, size5);
    fmk = FMatrixKokkos<real_t> (size1, size2, size3, size4, size5);
}

// 6D initialization function
void pmcb::init(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                size_t dim5, size_t dim6) {
    size1 = dim1;
    size2 = dim2;
    size3 = dim3;
    size4 = dim4;
    size5 = dim5;
    size6 = dim6;

    cak = CArrayKokkos<real_t>  (size1, size2, size3, size4, size5, size6);
    cmk = CMatrixKokkos<real_t> (size1, size2, size3, size4, size5, size6);
    fak = FArrayKokkos<real_t>  (size1, size2, size3, size4, size5, size6);
    fmk = FMatrixKokkos<real_t> (size1, size2, size3, size4, size5, size6);
}

// Destructor
pmcb::~pmcb() {}
