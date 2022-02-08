# The simple MATAR to BLAS/LAPACK interface is only compiled if the optional
# BLAS/LAPACK dependency specified by the user
if (WITH_BLAS_LAPACK)
  add_library(
    common 
      common.cpp
      error.cpp
      matar_blas_lapack_interface.cpp
  )
else()
  add_library(
    common 
      common.cpp
      error.cpp
  )
endif()

target_include_directories(
  common
    PUBLIC
      ${ELEMENTS_PUBLIC_INCLUDE_DIRECTORIES}
      ${MATAR_PUBLIC_INCLUDE_DIRECTORIES}
)
