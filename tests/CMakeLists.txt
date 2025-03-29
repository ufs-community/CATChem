include(CTest)

add_compile_options(${_catchem_compiler_options})

add_library(testing testing_mod.f90)
target_link_libraries(testing PRIVATE CATChem_core)
set_target_properties(
  testing
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

add_executable(test_main test_main.F90)
target_link_libraries(test_main PRIVATE CATChem_core)
target_link_libraries(test_main PRIVATE CATChem)
target_link_libraries(test_main PRIVATE testing)
set_target_properties(
  test_main
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
add_test(
  NAME test_main
  COMMAND test_main
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_test(
  NAME test_main_missing_config
  COMMAND test_main
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_property(TEST test_main_missing_config PROPERTY WILL_FAIL TRUE)

add_executable(test_dust test_dust.f90)
target_link_libraries(test_dust PRIVATE CATChem_core)
target_link_libraries(test_dust PRIVATE CATChem)
target_link_libraries(test_dust PRIVATE CATChem_process_dust)
target_link_libraries(test_dust PRIVATE testing)
set_target_properties(
  test_dust
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
add_test(
  NAME test_ccpr_dust
  COMMAND test_dust
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_seasalt test_seasalt.F90)
target_link_libraries(test_seasalt PRIVATE CATChem_core)
target_link_libraries(test_seasalt PRIVATE CATChem)
target_link_libraries(test_seasalt PRIVATE CATChem_process_seasalt)
target_link_libraries(test_seasalt PRIVATE testing)
set_target_properties(
  test_seasalt
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
add_test(
  NAME test_ccpr_seasalt
  COMMAND test_seasalt
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_prise test_plumerise.F90)
target_link_libraries(test_prise PRIVATE CATChem_core)
target_link_libraries(test_prise PRIVATE CATChem)
target_link_libraries(test_prise PRIVATE testing)
set_target_properties(
  test_prise
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
add_test(
  NAME test_plumerise
  COMMAND test_prise
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_drydep test_drydep.f90)
target_link_libraries(test_drydep PRIVATE CATChem_core)
target_link_libraries(test_drydep PRIVATE CATChem)
target_link_libraries(test_drydep PRIVATE CATChem_process_drydep)
target_link_libraries(test_drydep PRIVATE testing)
set_target_properties(
  test_drydep
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
add_test(
  NAME test_drydep
  COMMAND test_drydep
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_volcanic test_volcanic.F90)
target_link_libraries(test_volcanic PRIVATE CATChem_core)
target_link_libraries(test_volcanic PRIVATE CATChem)
target_link_libraries(test_volcanic PRIVATE CATChem_process_volcanic)
target_link_libraries(test_volcanic PRIVATE testing)
set_target_properties(
  test_volcanic
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

add_test(
  NAME test_volcanic
  COMMAND test_volcanic
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_dms test_dms.F90)
target_link_libraries(test_dms PRIVATE CATChem_core)
target_link_libraries(test_dms PRIVATE CATChem)
target_link_libraries(test_dms PRIVATE CATChem_process_DMS)
target_link_libraries(test_dms PRIVATE testing)
set_target_properties(
  test_dms
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

add_test(
  NAME test_dms
  COMMAND test_dms
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "Fortran module directory: ${CMAKE_BINARY_DIR}/include")

add_executable(test_chem test_chem.f90)
target_link_libraries(test_chem PRIVATE CATChem)
target_link_libraries(test_chem PRIVATE testing)
set_target_properties(
  test_chem
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
add_test(NAME test_chem COMMAND test_chem)
