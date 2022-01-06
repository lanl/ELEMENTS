# List the header files to be included in the documentation
file(GLOB_RECURSE DOC_HEADERS ../include/*.h)

# Set up variables for use in automating Doxygen documentation generation
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/include)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)

# Replace the placeholders in the Doxygen configuration file template with the
# values of the variables defined above
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile) 
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Make the output directory (Doxygen won't do so automatically)
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# Add rule for building the documentation (running the Doxygen command)
add_custom_command(
  OUTPUT ${DOXYGEN_INDEX_FILE}
  DEPENDS ${DOC_HEADERS}
  COMMAND doxygen ${DOXYFILE_OUT}
  MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
  COMMENT "Generating docs"
)

# Add target that ensures the rule above is applied (command is executed)
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
