cmake_minimum_required(VERSION 3.5)

project(CATChem VERSION 0.1.0 LANGUAGES Fortran)

include(CMakePrintHelpers)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# If GNU compiler and Debug build type, add extra debug flags
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
cmake_print_variables(
  CMAKE_Fortran_COMPILER
  CMAKE_Fortran_COMPILER_ID
  CMAKE_BUILD_TYPE
)
set(_catchem_compiler_options "")
if(
  CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  list(
    APPEND
    _catchem_compiler_options
    -Og
    -Wall
    -Wextra
    -pedantic
    -fcheck=all
    -ffpe-trap=invalid,zero,overflow
    -fbacktrace
    -fmax-errors=0
  )
elseif(
  CMAKE_Fortran_COMPILER_ID MATCHES "Intel"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  list(
    APPEND
    _catchem_compiler_options
    -O0
    -warn
    all
    -check
    all
    -implicitnone
    -fpe0
    -traceback
    -error-limit
    0
    -diag-disable
    6379,10448
  )
endif()

add_subdirectory(src)
add_subdirectory(tests)
