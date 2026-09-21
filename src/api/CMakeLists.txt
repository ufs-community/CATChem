set(_srcs CATChem_API.F90)
add_library(CATChem ${_srcs})
target_link_libraries(
  CATChem
  PUBLIC
    CATChem_core
    catchem_validation_options
    CATChem_core_cpp
    CATChem_process_seasalt
    CATChem_process_dust
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_settling
    CATChem_process_so4chem
    CATChem_process_carbchem
)
set_target_properties(
  CATChem
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(CATChem PUBLIC ${CMAKE_BINARY_DIR}/include)
