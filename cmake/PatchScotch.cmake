# This script performs a cross-platform sed-like replacement
file(READ "${FileToPatch}" content)
string(REPLACE "CMAKE_SOURCE_DIR" "CMAKE_CURRENT_SOURCE_DIR" content "${content}")
file(WRITE "${FileToPatch}" "${content}")