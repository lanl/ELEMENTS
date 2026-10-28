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
        ptscotch           
        scotch             
        ptscotcherrexit    
)

# Define Include Directories
target_include_directories(ELEMENTS INTERFACE
    # 1. Include the current directory ('src')
    #    This allows users to resolve "ELEMENTS.h", "swage/mesh.h", and "elements/material.h"
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    
    # 2. MATAR's header-only solvers (cramers_rule.hpp, lu_solver.hpp, ...),
    #    which the matar target no longer exposes
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/matar/solvers>

    # 3. Standard install location
    $<INSTALL_INTERFACE:include>
)

# Kokkos 5 requires C++20; propagate it to anyone linking ELEMENTS
target_compile_features(ELEMENTS INTERFACE cxx_std_20)

# Install Rules: The Target
install(TARGETS ELEMENTS
    EXPORT ElementsTargets
)

# Install Rules: The Headers
# Since we are in 'src', this single command recursively copies 
# 'ELEMENTS.h', the 'swage/' folder, and the 'elements/' folder 
# into the install 'include' directory.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ 
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)