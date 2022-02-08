add_library(geometry geometry.cpp)

target_include_directories(
  geometry
    PUBLIC
      ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
      ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)
