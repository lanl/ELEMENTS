add_library(input input.cpp)

target_include_directories(
  input
  PUBLIC
    ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
    ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)
