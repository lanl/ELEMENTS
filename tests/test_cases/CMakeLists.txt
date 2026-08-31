# tests/test_cases/CMakeLists.txt
#
# One executable, elements_tests, bundles every non-MPI test_*.cpp in this
# directory. A second executable, elements_mpi_tests, runs under mpirun
# separately since it needs multiple ranks.

include(GoogleTest)

file(GLOB ELEMENTS_TEST_SOURCES CONFIGURE_DEPENDS "test_*.cpp")
list(REMOVE_ITEM ELEMENTS_TEST_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/test_mpi_decomp.cpp")

add_executable(elements_tests test_main.cpp ${ELEMENTS_TEST_SOURCES})
target_link_libraries(elements_tests PRIVATE ELEMENTS gtest)
# PRE_TEST: enumerate tests when ctest runs, not by executing the binary
# right after linking (avoids Kokkos initialization during builds, which
# matters for the CUDA/HIP backends).
gtest_discover_tests(elements_tests DISCOVERY_MODE PRE_TEST)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_mpi_decomp.cpp" AND MPIEXEC_EXECUTABLE)
    add_executable(elements_mpi_tests mpi_test_main.cpp test_mpi_decomp.cpp)
    target_link_libraries(elements_mpi_tests PRIVATE ELEMENTS gtest)

    # --oversubscribe: GitHub Actions / small VMs often expose fewer slots
    # than ranks (Open MPI).
    add_test(
        NAME elements_mpi_tests
        COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:elements_mpi_tests>
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
