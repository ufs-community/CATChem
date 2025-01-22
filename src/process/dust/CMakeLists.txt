set(
  _srcs
  CCPr_Dust_Mod.F90
  ccpr_dust_common_mod.F90
  ccpr_scheme_fengsha_mod.F90
  ccpr_scheme_ginoux_mod.F90
)

set(_lib CATChem_process_dust)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
