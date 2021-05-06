project(ELEMENTS)
cmake_minimum_required(VERSION 3.2)

set (CMAKE_CXX_STANDARD 14)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

# Packages & modules
find_package(Vector)
if (CMAKE_VECTOR_NOVEC)
  set(VECTOR_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}")
  set(VECTOR_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}")
endif (CMAKE_VECTOR_NOVEC)
if (CMAKE_VECTOR_VERBOSE)
  set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
  set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
  set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} ${VECTOR_Fortran_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)

INCLUDE(cmake/Modules/TargetDistclean.cmake OPTIONAL)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# Compiler flags
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")

# Includes
include_directories(elements)
include_directories(geometry)
include_directories(swage)
include_directories(utilities)
include_directories(matar)

# Libraries
add_subdirectory(elements)
add_subdirectory(geometry)
add_subdirectory(swage)

# Examples
add_subdirectory(examples)

# Install headers in include/ subdirectory of install directory
#file(GLOB ELEMENTS_INCLUDE ${PROJECT_SOURCE_DIR}/Elements/*.h)
#install(FILES ${ELEMENTS_INCLUDE} DESTINATION include)
#install(FILES Utilities/utilities.h DESTINATION include)
file(GLOB_RECURSE ALL_HEADERS *.h)
install(FILES ${ALL_HEADERS} DESTINATION include)

# Install libraries in lib/ subdirectory of install directory
install(TARGETS elements DESTINATION lib)
install(TARGETS geometry DESTINATION lib)
install(TARGETS swage DESTINATION lib)
