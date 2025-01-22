set(
  _srcs
  CCPr_SeaSalt_Mod.F90
  ccpr_seasalt_common_mod.F90
  ccpr_scheme_gong03_mod.F90
  ccpr_scheme_gong97_mod.F90
  ccpr_scheme_geos12_mod.F90
)

set(_lib CATChem_process_seasalt)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
