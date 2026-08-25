set(
  _base_srcs
  constants.F90
  error_mod.F90
  UnitConversion_Mod.F90
  met_utilities_mod.F90
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
  ConfigManager_Mod.F90
  StateManager_Mod.F90
  DiagnosticManager_Mod.F90
  ChemSpeciesUtils_Mod.F90
  Interop_Mod.F90
)

set(_srcs ${_base_srcs} ${_state_srcs})

add_library(CATChem_core ${_srcs})
set_target_properties(
  CATChem_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
# Force serial build for this target to avoid Fortran module race conditions
set_property(TARGET CATChem_core PROPERTY JOB_POOLS compile=serial)
set_property(TARGET CATChem_core PROPERTY JOB_POOL_COMPILE serial)

foreach(opt IN LISTS _catchem_compiler_options)
  target_compile_options(
    CATChem_core
    PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${opt}>"
  )
endforeach()

target_link_libraries(CATChem_core PRIVATE yaml_interface catchem_utilities)

# Link GOCART2G library for external process dependencies
if(TARGET GOCART2G)
  target_link_libraries(CATChem_core PUBLIC GOCART2G)
endif()

# Conditional Kokkos linking
if(CATCHEM_ENABLE_KOKKOS)
  set(
    _cpp_core_srcs
    catchem_core.cpp
    catchem_api.cpp
    catchem_diagnostic.cpp
    catchem_config_manager.cpp
    catchem_grid_manager.cpp
    catchem_state_manager.cpp
    catchem_logger.cpp
  )
  add_library(CATChem_core_cpp STATIC ${_cpp_core_srcs})
  target_link_libraries(
    CATChem_core_cpp
    PUBLIC Kokkos::kokkos yaml-cpp::yaml-cpp
  )
  target_include_directories(
    CATChem_core_cpp
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_compile_definitions(CATChem_core_cpp PUBLIC CATCHEM_ENABLE_KOKKOS)

  target_link_libraries(CATChem_core PUBLIC CATChem_core_cpp)
  target_compile_definitions(CATChem_core PRIVATE CATCHEM_ENABLE_KOKKOS)
endif()

if(TARGET musica)
  target_link_libraries(CATChem_core PUBLIC musica::musica-fortran)
endif()

set_source_files_properties(metstate_mod.F90 PROPERTIES COMPILE_FLAGS "-cpp")

target_include_directories(CATChem_core PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
