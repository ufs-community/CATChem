set(_srcs CATChem_API.F90)
add_library(CATChem ${_srcs})
target_link_libraries(CATChem PUBLIC CATChem_core)
set_target_properties(
  CATChem
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(CATChem PUBLIC ${CMAKE_BINARY_DIR}/include)
