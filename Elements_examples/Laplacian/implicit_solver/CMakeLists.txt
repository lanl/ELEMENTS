set(Implicit_SRC_H ../include/header.h ../include/state.h ../include/utilities.h)
set(Implicit_SRC_CXX Implicit_Solver.cpp)

add_executable(Implicit_Solver ${Implicit_SRC_CXX} ${Implicit_SRC_H})

target_link_libraries (Implicit_Solver elements)
target_link_libraries (Implicit_Solver input)