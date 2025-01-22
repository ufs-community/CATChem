set(
  _srcs
  charpak_mod.F90
  chemstate_mod.F90
  cmn_size_mod.F90
  config_mod.F90
  config_opt_mod.F90
  constants.F90
  diagstate_mod.F90
  emisstate_mod.F90
  error_mod.F90
  gridstate_mod.F90
  init_mod.F90
  metstate_mod.F90
  precision_mod.F90
  qfyaml_mod.F90
  species_mod.F90
  state_mod.F90
)

add_library(CATChem_core ${_srcs})
set_target_properties(
  CATChem_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_compile_options(CATChem_core PRIVATE ${_catchem_compiler_options})
