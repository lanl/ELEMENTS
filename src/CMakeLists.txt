# --- swage/CMakeLists.txt ---

# Create the Interface Library
add_library(ELEMENTS INTERFACE)

# Link Dependencies
# These are inherited by anyone who links against ELEMENTS
target_link_libraries(ELEMENTS 
    INTERFACE
        Kokkos::kokkos 
        matar::matar
        MPI::MPI_CXX
        $<BUILD_INTERFACE:ptscotch>
        $<BUILD_INTERFACE:scotch>
        $<BUILD_INTERFACE:ptscotcherrexit>
)

# Define Include Directories
target_include_directories(ELEMENTS INTERFACE
    # 1. Include the current directory ('src')
    #    This allows users to resolve "ELEMENTS.h", "swage/mesh.h", and "elements/material.h"
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>

    # 2. MATAR's solvers headers (cramers_rule.hpp, lu_solver.hpp, qr_solver.hpp).
    #    MATAR's own target exports only src/include; the examples that need
    #    these solvers depend on ELEMENTS, not matar::matar directly, so this
    #    is where the dependency belongs.
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/matar/solvers>

    # 3. Standard install location
    $<INSTALL_INTERFACE:include>
)

# Install Rules: The Target
install(TARGETS ELEMENTS
    EXPORT ElementsTargets
)

# Install Rules: The Headers
# Copy headers from 'src' (current dir) to the install 'include' folder.
# This single command recursively copies 'ELEMENTS.h', the 'swage/' folder,
# and the 'elements/' folder into the install 'include' directory.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ 
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Install MATAR's solvers headers alongside, so the $<INSTALL_INTERFACE:include>
# path above resolves the same includes the examples use in the build tree.
install(DIRECTORY ${PROJECT_SOURCE_DIR}/matar/solvers/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)
