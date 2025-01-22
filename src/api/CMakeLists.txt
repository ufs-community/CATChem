set(_srcs catchem.F90)

set(_lib CATChem)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
target_link_libraries(${_lib} PUBLIC CATChem_process_dust)
target_link_libraries(${_lib} PUBLIC CATChem_process_seasalt)
target_link_libraries(${_lib} PUBLIC CATChem_process_plumerise)
target_link_libraries(${_lib} PUBLIC CATChem_process_drydep)
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_compile_options(${_lib} PRIVATE ${_catchem_compiler_options})
