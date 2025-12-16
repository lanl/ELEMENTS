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
    
    # 2. Standard install location
    $<INSTALL_INTERFACE:include>
)

# Install Rules: The Target
install(TARGETS ELEMENTS
    EXPORT ElementsTargets
)

# Install Rules: The Headers
# Copy headers from 'swage' (current dir) to the install 'include' folder
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ 
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Install Rules: The Headers
# Since we are in 'src', this single command recursively copies 
# 'ELEMENTS.h', the 'swage/' folder, and the 'elements/' folder 
# into the install 'include' directory.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ 
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)