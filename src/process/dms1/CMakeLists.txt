set(_srcs CCPr_DMS_Mod.F90 ccpr_scheme_GOCART_DMS_mod.F90)

set(_lib CATChem_process_DMS)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
target_link_libraries(${_lib} PUBLIC GOCART2G)

set_target_properties(${_lib} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
