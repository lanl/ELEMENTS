set(Laplace_SRC_H ../include/header.h ../include/state.h ../include/utilities.h)
set(Laplace_SRC_CXX Pseudo_Laplacian.cpp)

add_executable(Laplace ${Laplace_SRC_CXX} ${Laplace_SRC_H})

target_link_libraries (Laplace elements)
target_link_libraries (Laplace input)