add_library(slam slam.cpp)

target_include_directories(
  slam
    PUBLIC
      ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
      ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)
