cmake_minimum_required (VERSION 3.2)
project (fierro)

set (CMAKE_CXX_STANDARD 14)

include_directories(../matar)
include_directories(../elements)
include_directories(../swage)
include_directories(../utilities)
include_directories(../geometry)
include_directories(src)

include_directories(include)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

message( ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
message( ${CMAKE_MODULE_PATH} )
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

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
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${VECTOR_Fortran_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${VECTOR_CXX_FLAGS}")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../elements ${CMAKE_CURRENT_BINARY_DIR}/elements)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../swage ${CMAKE_CURRENT_BINARY_DIR}/swage)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../geometry ${CMAKE_CURRENT_BINARY_DIR}/geometry)
add_subdirectory(src)
add_subdirectory(test)

#if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
#  set_target_properties(file.cpp PROPERTIES COMPILE_FLAGS "${VECTOR_CXX_FLAGS_RELEASE}")
#endif (${CMAKE_BUILD_TYPE} STREQUAL "Release”)
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECTOR_CXX_FLAGS_RELEASE}")