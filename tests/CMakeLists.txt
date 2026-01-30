include(CTest)

add_compile_options(${_catchem_compiler_options})

add_library(testing testing_mod.f90)
target_link_libraries(testing PRIVATE CATChem_core)
set_target_properties(
  testing
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(testing PRIVATE ${CMAKE_BINARY_DIR}/include)

# Core framework tests
add_executable(test_CATChemCore test_CATChemCore.f90)
target_link_libraries(test_CATChemCore PRIVATE CATChem_core)
target_link_libraries(test_CATChemCore PRIVATE testing)
set_target_properties(
  test_CATChemCore
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(test_CATChemCore PRIVATE ${CMAKE_BINARY_DIR}/include)
add_dependencies(test_CATChemCore CATChem_core)
add_test(
  NAME test_CATChemCore
  COMMAND test_CATChemCore
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_StateManager test_StateManager.f90)
target_link_libraries(test_StateManager PRIVATE CATChem_core)
target_link_libraries(test_StateManager PRIVATE testing)
set_target_properties(
  test_StateManager
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(
  test_StateManager
  PRIVATE ${CMAKE_BINARY_DIR}/include
)
add_dependencies(test_StateManager CATChem_core)
add_test(
  NAME test_StateManager
  COMMAND test_StateManager
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_MetState test_MetState.f90)
target_link_libraries(test_MetState PRIVATE CATChem_core)
target_link_libraries(test_MetState PRIVATE testing)
set_target_properties(
  test_MetState
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(test_MetState PRIVATE ${CMAKE_BINARY_DIR}/include)
add_dependencies(test_MetState CATChem_core)
add_test(
  NAME test_MetState
  COMMAND test_MetState
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_GridManager test_GridManager.f90)
target_link_libraries(test_GridManager PRIVATE CATChem_core)
target_link_libraries(test_GridManager PRIVATE testing)
set_target_properties(
  test_GridManager
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(test_GridManager PRIVATE ${CMAKE_BINARY_DIR}/include)
add_dependencies(test_GridManager CATChem_core)
add_test(
  NAME test_GridManager
  COMMAND test_GridManager
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ConfigManager test_ConfigManager.f90)
target_link_libraries(test_ConfigManager PRIVATE CATChem_core)
target_link_libraries(test_ConfigManager PRIVATE testing)
set_target_properties(
  test_ConfigManager
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(
  test_ConfigManager
  PRIVATE ${CMAKE_BINARY_DIR}/include
)
add_dependencies(test_ConfigManager CATChem_core)
add_test(
  NAME test_ConfigManager
  COMMAND test_ConfigManager
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ProcessManager test_ProcessManager.f90)
target_link_libraries(test_ProcessManager PRIVATE CATChem_core)
target_link_libraries(test_ProcessManager PRIVATE testing)
set_target_properties(
  test_ProcessManager
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_ProcessManager
  COMMAND test_ProcessManager
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_DiagnosticManager test_DiagnosticManager.f90)
target_link_libraries(test_DiagnosticManager PRIVATE CATChem_core)
target_link_libraries(test_DiagnosticManager PRIVATE testing)
set_target_properties(
  test_DiagnosticManager
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_DiagnosticManager
  COMMAND test_DiagnosticManager
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_Error test_Error.f90)
target_link_libraries(test_Error PRIVATE CATChem_core)
target_link_libraries(test_Error PRIVATE testing)
set_target_properties(
  test_Error
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_Error
  COMMAND test_Error
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_Precision test_Precision.f90)
target_link_libraries(test_Precision PRIVATE CATChem_core)
target_link_libraries(test_Precision PRIVATE testing)
set_target_properties(
  test_Precision
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_Precision
  COMMAND test_Precision
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ChemSpeciesUtils test_ChemSpeciesUtils.f90)
target_link_libraries(test_ChemSpeciesUtils PRIVATE CATChem_core)
target_link_libraries(test_ChemSpeciesUtils PRIVATE testing)
set_target_properties(
  test_ChemSpeciesUtils
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_ChemSpeciesUtils
  COMMAND test_ChemSpeciesUtils
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ProcessInterface test_ProcessInterface.f90)
target_link_libraries(test_ProcessInterface PRIVATE CATChem_core)
target_link_libraries(test_ProcessInterface PRIVATE testing)
set_target_properties(
  test_ProcessInterface
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_ProcessInterface
  COMMAND test_ProcessInterface
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_GridGeometry test_GridGeometry.f90)
target_link_libraries(test_GridGeometry PRIVATE CATChem_core)
target_link_libraries(test_GridGeometry PRIVATE testing)
set_target_properties(
  test_GridGeometry
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_GridGeometry
  COMMAND test_GridGeometry
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_TimeState test_TimeState.f90)
target_link_libraries(test_TimeState PRIVATE CATChem_core)
target_link_libraries(test_TimeState PRIVATE testing)
set_target_properties(
  test_TimeState
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_TimeState
  COMMAND test_TimeState
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_Constants test_Constants.f90)
target_link_libraries(test_Constants PRIVATE CATChem_core)
target_link_libraries(test_Constants PRIVATE testing)
set_target_properties(
  test_Constants
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_Constants
  COMMAND test_Constants
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_VirtualColumn test_VirtualColumn.f90)
target_link_libraries(test_VirtualColumn PRIVATE CATChem_core)
target_link_libraries(test_VirtualColumn PRIVATE testing)
set_target_properties(
  test_VirtualColumn
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_VirtualColumn
  COMMAND test_VirtualColumn
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ColumnInterface test_ColumnInterface.f90)
target_link_libraries(test_ColumnInterface PRIVATE CATChem_core)
target_link_libraries(test_ColumnInterface PRIVATE testing)
set_target_properties(
  test_ColumnInterface
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_ColumnInterface
  COMMAND test_ColumnInterface
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# add_executable(test_StateInterface test_StateInterface.f90)
# target_link_libraries(test_StateInterface PRIVATE CATChem_core)
# target_link_libraries(test_StateInterface PRIVATE testing)
# set_target_properties(
#   test_StateInterface
#   PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
# )
# add_test(
#   NAME test_StateInterface
#   COMMAND test_StateInterface
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )

add_executable(test_DiagnosticInterface test_DiagnosticInterface.f90)
target_link_libraries(test_DiagnosticInterface PRIVATE CATChem_core)
target_link_libraries(test_DiagnosticInterface PRIVATE testing)
set_target_properties(
  test_DiagnosticInterface
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_DiagnosticInterface
  COMMAND test_DiagnosticInterface
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ProcessRegistry test_ProcessRegistry.f90)
target_link_libraries(test_ProcessRegistry PRIVATE CATChem_core)
target_link_libraries(test_ProcessRegistry PRIVATE testing)
set_target_properties(
  test_ProcessRegistry
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_ProcessRegistry
  COMMAND test_ProcessRegistry
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(test_ProcessFactory test_ProcessFactory.f90)
target_link_libraries(test_ProcessFactory PRIVATE CATChem_core)
target_link_libraries(test_ProcessFactory PRIVATE testing)
set_target_properties(
  test_ProcessFactory
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
add_test(
  NAME test_ProcessFactory
  COMMAND test_ProcessFactory
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# New CATChem API test
add_executable(test_catchem_api test_catchem_api.F90)
target_link_libraries(test_catchem_api PRIVATE CATChem)
target_link_libraries(test_catchem_api PRIVATE testing)
set_target_properties(
  test_catchem_api
  PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
    LINKER_LANGUAGE Fortran
)
target_include_directories(test_catchem_api PRIVATE ${CMAKE_BINARY_DIR}/include)
add_dependencies(test_catchem_api CATChem)
add_test(
  NAME test_catchem_api
  COMMAND test_catchem_api
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# # Legacy tests (to be updated)
# add_executable(test_main test_main.F90)
# target_link_libraries(test_main PRIVATE CATChem_core)
# target_link_libraries(test_main PRIVATE CATChem)
# target_link_libraries(test_main PRIVATE testing)
# set_target_properties(
#   test_main
#   PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
# )
# add_test(
#   NAME test_main
#   COMMAND test_main
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )
# add_test(
#   NAME test_main_missing_config
#   COMMAND test_main
#   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
# )
# set_property(TEST test_main_missing_config PROPERTY WILL_FAIL TRUE)

# add_executable(test_dust test_dust.f90)
# target_link_libraries(test_dust PRIVATE CATChem_core)
# target_link_libraries(test_dust PRIVATE CATChem)
# target_link_libraries(test_dust PRIVATE CATChem_process_dust)
# target_link_libraries(test_dust PRIVATE testing)
# set_target_properties(
#   test_dust
#   PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
# )
# add_test(
#   NAME test_ccpr_dust
#   COMMAND test_dust
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )

# add_executable(test_seasalt test_seasalt.F90)
# target_link_libraries(test_seasalt PRIVATE CATChem_core)
# target_link_libraries(test_seasalt PRIVATE CATChem)
# target_link_libraries(test_seasalt PRIVATE CATChem_process_seasalt)
# target_link_libraries(test_seasalt PRIVATE testing)
# set_target_properties(
#   test_seasalt
#   PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
# )
# add_test(
#   NAME test_ccpr_seasalt
#   COMMAND test_seasalt
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )

# add_executable(test_prise test_plumerise.F90)
# target_link_libraries(test_prise PRIVATE CATChem_core)
# target_link_libraries(test_prise PRIVATE CATChem)
# target_link_libraries(test_prise PRIVATE testing)
# set_target_properties(
#   test_prise
#   PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
# )
# add_test(
#   NAME test_plumerise
#   COMMAND test_prise
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )

# add_executable(test_drydep test_drydep.f90)
# target_link_libraries(test_drydep PRIVATE CATChem_core)
# target_link_libraries(test_drydep PRIVATE CATChem)
# target_link_libraries(test_drydep PRIVATE testing)
# set_target_properties(
#   test_drydep
#   PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
# )
# add_test(
#   NAME test_drydep
#   COMMAND test_drydep
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
# )

#add_executable(test_chem test_chem.f90)
#target_link_libraries(test_chem PRIVATE CATChem)
#target_link_libraries(test_chem PRIVATE testing)
#set_target_properties(
#  test_chem
#  PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
#     LINKER_LANGUAGE Fortran
#)
#add_test(NAME test_chem COMMAND test_chem)
