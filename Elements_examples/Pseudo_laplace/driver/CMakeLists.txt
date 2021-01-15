set(Laplace_SRC_H ../include/header.h ../include/state.h ../include/utilities.h)
set(Class_SRC_H ../include/Simulation_Parameters.h Solver.h Pseudo_Laplacian.h Pseudo_Laplacian.cpp)
set(Laplace_SRC_CXX main.cpp Solver.cpp ../include/Simulation_Parameters.cpp)

add_executable(main ${Laplace_SRC_CXX} ${Laplace_SRC_H} ${Class_SRC_H})

target_link_libraries (main swage)
target_link_libraries (main elements)
target_link_libraries (main input)