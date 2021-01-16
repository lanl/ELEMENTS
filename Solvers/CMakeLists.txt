set(Solver_SRC_H ../Include/header.h ../Include/state.h ../Include/utilities.h ../Include/math_util.h)
set(Class_SRC Solver.h ../Include/Simulation_Parameters.h ../include/Simulation_Parameters.cpp)
set(Laplace_SRC Pseudo_Laplace/Pseudo_Laplacian.h Pseudo_Laplace/Pseudo_Laplacian.cpp)
set(Static_SRC Static_Solver/Static_Solver.h Static_Solver/Static_Solver.cpp)
set(Solver_SRC_CXX main.cpp Solver.cpp ../Include/math_util.cpp)

include_directories(Static_Solver)
include_directories(Pseudo_Laplace)

add_executable(main ${Solver_SRC_CXX} ${Solver_SRC_H} ${Class_SRC} ${Laplace_SRC} ${Static_SRC})

target_link_libraries (main swage)
target_link_libraries (main elements)
target_link_libraries (main input)