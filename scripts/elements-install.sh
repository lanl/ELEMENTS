#!/bin/bash -e

kokkos_build_type=${1}

rm -rf ${ELEMENTS_INSTALL_DIR}
mkdir -p ${ELEMENTS_BUILD_DIR} 

cmake_options=(
    -D CMAKE_INSTALL_PREFIX="${ELEMENTS_INSTALL_DIR}"
    -D CMAKE_PREFIX_PATH="${KOKKOS_INSTALL_DIR}"
)

if [ "$kokkos_build_type" = "none" ]; then
    cmake_options+=(
        -D Matar_ENABLE_KOKKOS=OFF
    )
else
    cmake_options+=(
        -D Matar_ENABLE_KOKKOS=ON
    )
fi

# Print CMake options for reference
echo "CMake Options: ${cmake_options[@]}"

# Configure Elements
cmake "${cmake_options[@]}" -B "${ELEMENTS_BUILD_DIR}" -S "${ELEMENTS_SOURCE_DIR}"

# Build Elements
echo "Building Elements..."
make -C ${ELEMENTS_BUILD_DIR} -j${ELEMENTS_BUILD_CORES}

# Install Elements
echo "Installing Elements..."
make -C ${ELEMENTS_BUILD_DIR} install

echo "Elements installation complete"
