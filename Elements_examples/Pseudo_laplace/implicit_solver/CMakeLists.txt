set(Implicit_Solver_SRC_H ../include/header.h ../include/state.h ../include/utilities.h)
set(Implicit_Solver_SRC_CXX Static_Solve.cpp)

add_executable(Implicit_Solver ${Implicit_Solver_SRC_CXX} ${Implicit_Solver_SRC_H})

target_link_libraries (Implicit_Solver swage)
target_link_libraries (Implicit_Solver elements)
target_link_libraries (Implicit_Solver input)