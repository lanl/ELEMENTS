# tests/CMakeLists.txt
#
# GoogleTest unit tests for ELEMENTS, modeled on MATAR's test suite
# (matar/test/CMakeLists.txt). ELEMENTS' own build forces
# MATAR_BUILD_TESTS=OFF, so googletest is not already fetched and this
# subdirectory fetches its own copy.

include(FetchContent)

if(NOT TARGET gtest)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
    )
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    # On Windows: keep gtest from overriding the parent project's runtime settings.
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

add_subdirectory(test_cases)
