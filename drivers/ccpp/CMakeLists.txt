set(_srcs ccpp_catchem_interface.F90)

set(_lib CCPP_CATChem_interface)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem)

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_compile_options(${_lib} PRIVATE ${_catchem_compiler_options})
