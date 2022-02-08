#set(Laplace_SRC_H ../include/header.h ../include/state.h ../../utili/utilities.h)

set(Average_SRC_CXX average.cpp)

add_executable(Average ${Average_SRC_CXX})

target_include_directories(
  Average
  PUBLIC
    ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
    ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)

target_link_libraries(
  Average 
    swage
    elements
    geometry
    input
)
