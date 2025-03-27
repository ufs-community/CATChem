#
# GOCART-2G Process Library
#

set(
  _srcs
  GOCART/Process_Library/GOCART2G_Process.F90
  GOCART/Process_Library/GOCART2G_MieMod.F90
)
set(_lib GOCART2G)

add_library(${_lib} ${_srcs})
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
  target_compile_options(${_lib} PRIVATE -ffree-line-length-none)
endif()

find_package(NetCDF REQUIRED COMPONENTS Fortran)

target_link_libraries(${_lib} NetCDF::NetCDF_Fortran)

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

#
# MICM Fortran API
#

include(FetchContent)

FetchContent_Declare(
  musica_fortran
  GIT_REPOSITORY https://github.com/NCAR/musica.git
  GIT_TAG v0.9.0
)

set(MICM_ENABLE_TESTS OFF)
set(MICM_ENABLE_EXAMPLES OFF)

set(MUSICA_BUILD_C_CXX_INTERFACE ON)
set(MUSICA_BUILD_FORTRAN_INTERFACE ON)
set(MUSICA_ENABLE_INSTALL OFF)
set(MUSICA_ENABLE_TESTS OFF)
set(MUSICA_ENABLE_MICM ON)
set(MUSICA_ENABLE_TUVX OFF)

FetchContent_MakeAvailable(musica_fortran)
