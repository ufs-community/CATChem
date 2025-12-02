set(
  _base_srcs
  constants.F90
  utilities_mod.F90
  met_utilities_mod.F90
  error_mod.F90
)

set(
  _state_srcs
  TimeState_Mod.F90
  species_mod.F90
  GridGeometry_Mod.F90
  metstate_mod.F90
  chemstate_mod.F90
  ExtEmisData_Mod.F90
  DiagnosticInterface_Mod.F90
)

set(
  _core_srcs
  ConfigManager_Mod.F90
  VirtualColumn_Mod.F90
  ColumnInterface_Mod.F90
  GridManager_Mod.F90
  DiagnosticManager_Mod.F90
  StateManager_Mod.F90
  ChemSpeciesUtils_Mod.F90
  EmissionConfigValidator_Mod.F90
  ProcessInterface_Mod.F90
  ProcessFactory_Mod.F90
  ProcessRegistry_Mod.F90
  ProcessManager_Mod.F90
  CATChemCore_Mod.F90
)

set(_srcs ${_base_srcs} ${_state_srcs} ${_core_srcs})
# add Python3 package for generating Fortran module macros
find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_library(CATChem_core ${_srcs})
set_target_properties(
  CATChem_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
# Force serial build for this target to avoid Fortran module race conditions
set_property(TARGET CATChem_core PROPERTY JOB_POOLS compile=serial)
set_property(TARGET CATChem_core PROPERTY JOB_POOL_COMPILE serial)

target_compile_options(CATChem_core PRIVATE ${_catchem_compiler_options})

target_link_libraries(CATChem_core PRIVATE yaml_interface catchem_utilities)

# Link GOCART2G library for external process dependencies
if(TARGET GOCART2G)
  target_link_libraries(CATChem_core PUBLIC GOCART2G)
endif()

set(METSTATE_MOD_F90 ${CMAKE_CURRENT_SOURCE_DIR}/metstate_mod.F90)
set(
  METSTATE_ACCESSOR_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_column_accessor.inc
)
set(
  METSTATE_ACCESSOR_INT_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_column_accessor_int.inc
)
set(
  METSTATE_ACCESSOR_LOGICAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_column_accessor_logical.inc
)
set(
  METSTATE_ALLOCATE_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_allocate_fields.inc
)
set(
  METSTATE_DEALLOCATE_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_deallocate_fields.inc
)
set(
  METSTATE_2D_SCALAR_ACCESSOR_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_2d_scalar_accessor.inc
)
set(
  METSTATE_2D_SCALAR_ACCESSOR_INT_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_2d_scalar_accessor_int.inc
)
set(
  METSTATE_2D_SCALAR_ACCESSOR_LOGICAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_2d_scalar_accessor_logical.inc
)
set(
  METSTATE_SCALAR_ACCESSOR_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_scalar_accessor.inc
)
set(
  METSTATE_SCALAR_ACCESSOR_INT_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_scalar_accessor_int.inc
)
set(
  METSTATE_SCALAR_ACCESSOR_LOGICAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_scalar_accessor_logical.inc
)
set(
  METSTATE_SET_FIELD_2D_REAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_2d_real.inc
)
set(
  METSTATE_SET_FIELD_2D_INT_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_2d_int.inc
)
set(
  METSTATE_SET_FIELD_2D_LOGICAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_2d_logical.inc
)
set(
  METSTATE_SET_FIELD_3D_REAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_3d_real.inc
)
set(
  METSTATE_SET_FIELD_3D_INT_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_3d_int.inc
)
set(
  METSTATE_SET_FIELD_3D_LOGICAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_3d_logical.inc
)
set(
  METSTATE_SET_FIELD_SCALAR_REAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_scalar_real.inc
)
set(
  METSTATE_SET_FIELD_SCALAR_INT_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_scalar_int.inc
)
set(
  METSTATE_SET_FIELD_SCALAR_LOGICAL_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_set_field_scalar_logical.inc
)
set(
  METSTATE_MULTIPLE_FIELDS_INTERFACE_INC
  ${CMAKE_CURRENT_BINARY_DIR}/metstate_multiple_fields_interface.inc
)
set(VIRTUALMET_TYPE_INC ${CMAKE_CURRENT_BINARY_DIR}/virtualmet_type.inc)
set(VIRTUALMET_POPULATE_INC ${CMAKE_CURRENT_BINARY_DIR}/virtualmet_populate.inc)
set(VIRTUALMET_CLEANUP_INC ${CMAKE_CURRENT_BINARY_DIR}/virtualmet_cleanup.inc)

add_custom_command(
  OUTPUT ${METSTATE_ACCESSOR_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_column_accessor.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    column_accessor ${METSTATE_MOD_F90} ${METSTATE_ACCESSOR_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_column_accessor.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_ALLOCATE_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_allocate_fields.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py allocate
    ${METSTATE_MOD_F90} ${METSTATE_ALLOCATE_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_allocate_fields.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_DEALLOCATE_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_deallocate_fields.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    deallocate ${METSTATE_MOD_F90} ${METSTATE_DEALLOCATE_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_deallocate_fields.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_2D_SCALAR_ACCESSOR_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_2d_scalar_accessor.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    2d_scalar_accessor ${METSTATE_MOD_F90} ${METSTATE_2D_SCALAR_ACCESSOR_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_2d_scalar_accessor.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SCALAR_ACCESSOR_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_scalar_accessor.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    scalar_accessor ${METSTATE_MOD_F90} ${METSTATE_SCALAR_ACCESSOR_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_scalar_accessor.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_ACCESSOR_INT_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_column_accessor_int.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    column_accessor_int ${METSTATE_MOD_F90} ${METSTATE_ACCESSOR_INT_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_column_accessor_int.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_ACCESSOR_LOGICAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_column_accessor_logical.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    column_accessor_logical ${METSTATE_MOD_F90} ${METSTATE_ACCESSOR_LOGICAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_column_accessor_logical.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_2D_SCALAR_ACCESSOR_INT_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_2d_scalar_accessor_int.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    2d_scalar_accessor_int ${METSTATE_MOD_F90}
    ${METSTATE_2D_SCALAR_ACCESSOR_INT_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_2d_scalar_accessor_int.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_2D_SCALAR_ACCESSOR_LOGICAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_2d_scalar_accessor_logical.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    2d_scalar_accessor_logical ${METSTATE_MOD_F90}
    ${METSTATE_2D_SCALAR_ACCESSOR_LOGICAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_2d_scalar_accessor_logical.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SCALAR_ACCESSOR_INT_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_scalar_accessor_int.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    scalar_accessor_int ${METSTATE_MOD_F90} ${METSTATE_SCALAR_ACCESSOR_INT_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_scalar_accessor_int.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SCALAR_ACCESSOR_LOGICAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_scalar_accessor_logical.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    scalar_accessor_logical ${METSTATE_MOD_F90}
    ${METSTATE_SCALAR_ACCESSOR_LOGICAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_scalar_accessor_logical.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_2D_REAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_2d_real.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_2d_real ${METSTATE_MOD_F90} ${METSTATE_SET_FIELD_2D_REAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_2d_real.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_2D_INT_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_2d_int.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_2d_int ${METSTATE_MOD_F90} ${METSTATE_SET_FIELD_2D_INT_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_2d_int.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_2D_LOGICAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_2d_logical.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_2d_logical ${METSTATE_MOD_F90}
    ${METSTATE_SET_FIELD_2D_LOGICAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_2d_logical.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_3D_REAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_3d_real.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_3d_real ${METSTATE_MOD_F90} ${METSTATE_SET_FIELD_3D_REAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_3d_real.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_3D_INT_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_3d_int.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_3d_int ${METSTATE_MOD_F90} ${METSTATE_SET_FIELD_3D_INT_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_3d_int.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_3D_LOGICAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_3d_logical.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_3d_logical ${METSTATE_MOD_F90}
    ${METSTATE_SET_FIELD_3D_LOGICAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_3d_logical.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_SCALAR_REAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_scalar_real.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_scalar_real ${METSTATE_MOD_F90}
    ${METSTATE_SET_FIELD_SCALAR_REAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_scalar_real.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_SCALAR_INT_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_scalar_int.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_scalar_int ${METSTATE_MOD_F90}
    ${METSTATE_SET_FIELD_SCALAR_INT_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_scalar_int.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_SET_FIELD_SCALAR_LOGICAL_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_set_field_scalar_logical.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    set_field_scalar_logical ${METSTATE_MOD_F90}
    ${METSTATE_SET_FIELD_SCALAR_LOGICAL_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_set_field_scalar_logical.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${METSTATE_MULTIPLE_FIELDS_INTERFACE_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating metstate_multiple_fields_interface.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    multiple_fields_interface ${METSTATE_MOD_F90}
    ${METSTATE_MULTIPLE_FIELDS_INTERFACE_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating metstate_multiple_fields_interface.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${VIRTUALMET_TYPE_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating virtualmet_type.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    virtualmet_type ${METSTATE_MOD_F90} ${VIRTUALMET_TYPE_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating virtualmet_type.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${VIRTUALMET_POPULATE_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating virtualmet_populate.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    virtualmet_populate ${METSTATE_MOD_F90} ${VIRTUALMET_POPULATE_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating virtualmet_populate.inc"
  VERBATIM
)

add_custom_command(
  OUTPUT ${VIRTUALMET_CLEANUP_INC}
  COMMAND
    ${CMAKE_COMMAND} -E echo
    "Generating virtualmet_cleanup.inc from MetState_Mod.F90"
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
    ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    virtualmet_cleanup ${METSTATE_MOD_F90} ${VIRTUALMET_CLEANUP_INC}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/generate_metstate_macros.py
    ${METSTATE_MOD_F90}
  COMMENT "Auto-generating virtualmet_cleanup.inc"
  VERBATIM
)

add_custom_target(
  metstate_macros
  ALL
  DEPENDS
    ${METSTATE_ACCESSOR_INC}
    ${METSTATE_ACCESSOR_INT_INC}
    ${METSTATE_ACCESSOR_LOGICAL_INC}
    ${METSTATE_ALLOCATE_INC}
    ${METSTATE_DEALLOCATE_INC}
    ${METSTATE_2D_SCALAR_ACCESSOR_INC}
    ${METSTATE_2D_SCALAR_ACCESSOR_INT_INC}
    ${METSTATE_2D_SCALAR_ACCESSOR_LOGICAL_INC}
    ${METSTATE_SCALAR_ACCESSOR_INC}
    ${METSTATE_SCALAR_ACCESSOR_INT_INC}
    ${METSTATE_SCALAR_ACCESSOR_LOGICAL_INC}
    ${METSTATE_SET_FIELD_INC}
    ${METSTATE_SET_FIELD_2D_REAL_INC}
    ${METSTATE_SET_FIELD_2D_INT_INC}
    ${METSTATE_SET_FIELD_2D_LOGICAL_INC}
    ${METSTATE_SET_FIELD_3D_REAL_INC}
    ${METSTATE_SET_FIELD_3D_INT_INC}
    ${METSTATE_SET_FIELD_3D_LOGICAL_INC}
    ${METSTATE_SET_FIELD_SCALAR_REAL_INC}
    ${METSTATE_SET_FIELD_SCALAR_INT_INC}
    ${METSTATE_SET_FIELD_SCALAR_LOGICAL_INC}
    ${METSTATE_MULTIPLE_FIELDS_INTERFACE_INC}
    ${VIRTUALMET_TYPE_INC}
    ${VIRTUALMET_POPULATE_INC}
    ${VIRTUALMET_CLEANUP_INC}
)

add_dependencies(CATChem_core metstate_macros)

set_source_files_properties(
  metstate_mod.F90
  PROPERTIES
    COMPILE_FLAGS "-cpp"
    OBJECT_DEPENDS
      "${METSTATE_ACCESSOR_INC};${METSTATE_ACCESSOR_INT_INC};${METSTATE_ACCESSOR_LOGICAL_INC};${METSTATE_ALLOCATE_INC};${METSTATE_DEALLOCATE_INC};${METSTATE_2D_SCALAR_ACCESSOR_INC};${METSTATE_2D_SCALAR_ACCESSOR_INT_INC};${METSTATE_2D_SCALAR_ACCESSOR_LOGICAL_INC};${METSTATE_SCALAR_ACCESSOR_INC};${METSTATE_SCALAR_ACCESSOR_INT_INC};${METSTATE_SCALAR_ACCESSOR_LOGICAL_INC};${METSTATE_SET_FIELD_INC};${METSTATE_SET_FIELD_2D_REAL_INC};${METSTATE_SET_FIELD_2D_INT_INC};${METSTATE_SET_FIELD_2D_LOGICAL_INC};${METSTATE_SET_FIELD_3D_REAL_INC};${METSTATE_SET_FIELD_3D_INT_INC};${METSTATE_SET_FIELD_3D_LOGICAL_INC};${METSTATE_SET_FIELD_SCALAR_REAL_INC};${METSTATE_SET_FIELD_SCALAR_INT_INC};${METSTATE_SET_FIELD_SCALAR_LOGICAL_INC};${METSTATE_MULTIPLE_FIELDS_INTERFACE_INC}"
)

set_source_files_properties(
  VirtualColumn_Mod.F90
  PROPERTIES
    COMPILE_FLAGS "-cpp"
    OBJECT_DEPENDS
      "${VIRTUALMET_TYPE_INC};${VIRTUALMET_POPULATE_INC};${VIRTUALMET_CLEANUP_INC}"
)

target_include_directories(CATChem_core PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
