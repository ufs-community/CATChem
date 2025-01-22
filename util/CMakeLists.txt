set(_srcs reademissions.F90 prepmetvars.F90)

set(_lib reademissions)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
target_link_libraries(${_lib} PUBLIC GOCART2G)

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
