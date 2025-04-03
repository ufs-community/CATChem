set(_lib CATChem_process_chem)

set(_srcs CCPr_Chem_Mod.F90)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
target_link_libraries(${_lib} PRIVATE musica::musica-fortran)
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
