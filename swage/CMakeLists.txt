# 1. List your compiled sources (REMOVE mesh.cpp)
# set(SOURCES
#     # basis.cpp
#     # utility.cpp (if you have other cpp files)
# )

# 2. (Optional but Recommended) List headers for IDE visibility
#    This makes them show up in Visual Studio/Xcode project trees,
#    even though they are not "compiled".
# set(HEADERS
#     mesh.h
#     ../elements/ref_elem.h
#     ../elements/ref_surf_elem.h
#     ../elements/ref_quadrature.h
# )

# Create the library
#    You can include ${HEADERS} here; CMake ignores them for compilation
#    but uses them for IDE display.
add_library(ELEMENTS INTERFACE)#${SOURCES} ${HEADERS})

# Link Dependencies
#    PUBLIC ensures that anyone using 'ELEMENTS' also gets the include paths
#    and compiler flags for Kokkos and MATAR automatically.
target_link_libraries(ELEMENTS 
    INTERFACE
        Kokkos::kokkos 
        matar::matar
)

# Define Include Directories
#    This is the most important part for header-only classes.
#    It tells CMake: "Look in ../include when building this lib,
#    and look in include/ when this lib is installed."
target_include_directories(ELEMENTS INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)

# Install Rules
#    (Interface libraries only need the EXPORT and header installation)
install(TARGETS ELEMENTS
    EXPORT ElementsTargets
    # No LIBRARY/ARCHIVE destination needed for header-only
)

install(DIRECTORY ../include/ DESTINATION include)