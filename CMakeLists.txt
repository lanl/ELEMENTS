project(ELEMENTS)
cmake_minimum_required(VERSION 3.12)  # same version as in MATAR CMakeLists

# User must provide information to locate MATAR dependency
# (Here I'm asking the user to define two separate variables because, as it is
# now, the MATAR library could be built anywhere in the filesystem but the
# header will remain behind with the source. It would be nice if this were
# fixed so that the header is installed/copied into an include/ subdirectory of
# the build directory and the library is installed/copied into a lib/
# subdirectory of the build directory. That would standardize the build
# somewhat and make external use of the library, like here, simpler and
# easier.)
# (Also might want to make a CMake find_package script for MATAR.)
set(MATAR_INC_DIR "MATAR_INC_DIR" CACHE PATH "Directory containing MATAR include")
set(MATAR_INC_DIR "MATAR_LIB_DIR" CACHE PATH "Directory containing MATAR include")

# Includes
include_directories(Utilities)
include_directories(${MATAR_INC_DIR})

# Libraries to be made
add_library(elements Elements/elements.cpp)
add_library(swage Elements/swage.cpp)
target_link_libraries(swage ${MATAR_LIB_DIR}/libmatar.a)

# Install headers in include/ subdirectory of install directory
file(GLOB ELEMENTS_INCLUDE ${PROJECT_SOURCE_DIR}/Elements/*.h)
install(FILES ${ELEMENTS_INCLUDE} DESTINATION include)
install(FILES Utilities/utilities.h DESTINATION include)

# Install libraries in lib/ subdirectory of install directory
install(TARGETS elements DESTINATION lib)
install(TARGETS swage DESTINATION lib)
