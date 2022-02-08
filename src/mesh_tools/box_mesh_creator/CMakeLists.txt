message(${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES})
add_executable(Meshbox meshbox.cpp)

target_include_directories(
  Meshbox
    PUBLIC
      ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
      ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)
