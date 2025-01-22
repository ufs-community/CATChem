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
# target_link_libraries(${_lib} PUBLIC CATChem_core)

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
