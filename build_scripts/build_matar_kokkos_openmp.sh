#!/bin/sh

# Build configuration (modify this!)
ELEMENTS_DIR=/path/to/elements
MATAR_DIR=/path/to/matar/within/elements
PROJECT_BUILD_DIR=/path/to/project/build/directory
export OMP_NUM_THREADS=sensible_number_of_openmp_threads_for_your_machine
export OMP_PROC_BIND=spread


# Build and install Kokkos with OpenMP
KOKKOS_SOURCE_DIR=${MATAR_DIR}/src/Kokkos/kokkos
KOKKOS_BUILD_DIR=${PROJECT_BUILD_DIR}/kokkos
KOKKOS_INSTALL_DIR=${KOKKOS_BUILD_DIR}

OPTIONS=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_DIR}
  -DCMAKE_CXX_STANDARD=17
  -DKokkos_ENABLE_SERIAL=ON
  -DKokkos_ENABLE_OPENMP=ON
  -DKokkos_ENABLE_TESTS=OFF
  -DBUILD_TESTING=OFF
)

cmake "${OPTIONS[@]}" -S ${KOKKOS_SOURCE_DIR} -B ${KOKKOS_BUILD_DIR}
cmake --build ${KOKKOS_BUILD_DIR}
cmake --install ${KOKKOS_BUILD_DIR}


# Build and install ELEMENTS with MATAR/Kokkos
ELEMENTS_SOURCE_DIR=${ELEMENTS_DIR}
ELEMENTS_BUILD_DIR=${PROJECT_BUILD_DIR}/elements
ELEMENTS_INSTALL_DIR=${ELEMENTS_BUILD_DIR}

OPTIONS=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX=${ELEMENTS_INSTALL_DIR}
  -DMATAR_WITH_KOKKOS=ON
  -DKokkos_DIR=${KOKKOS_INSTALL_DIR}/lib/cmake/Kokkos
  -DKOKKOS=ON
  -DOPENMP=ON
)

cmake "${OPTIONS[@]}" -S ${ELEMENTS_SOURCE_DIR} -B ${ELEMENTS_BUILD_DIR}
cmake --build ${ELEMENTS_BUILD_DIR} --verbose
cmake --install ${ELEMENTS_BUILD_DIR}
