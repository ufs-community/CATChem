foreach(opt IN LISTS _catchem_cxx_warning_options)
  add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:${opt}>")
endforeach()

set(
  _cpp_core_srcs
  catchem_core.cpp
  catchem_api.cpp
  catchem_api_lifecycle.cpp
  catchem_api_state.cpp
  catchem_api_diagnostics.cpp
  catchem_api_config.cpp
  catchem_api_utilities.cpp
  catchem_diagnostic.cpp
  catchem_config_manager.cpp
  catchem_grid_manager.cpp
  catchem_state_manager.cpp
  catchem_emission_manager.cpp
  catchem_logger.cpp
  catchem_runtime_lease.cpp
  catchem_handle_registry.hpp
  catchem_boundary_guard.hpp
  catchem_physical_validation.cpp
  catchem_physical_validation.hpp
  catchem_runtime_lease.hpp
)

# Focused architecture sources are added as they are introduced. Keeping this
# list explicit makes subsystem ownership visible without breaking incremental
# feature work when a source has not landed yet.
set(
  _catchem_architecture_sources
  catchem_validation.hpp
  catchem_execution_plan.cpp
  catchem_execution_plan.hpp
  catchem_field_contract.hpp
  catchem_handle_registry.hpp
)
foreach(_source IN LISTS _catchem_architecture_sources)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_source}")
    list(APPEND _cpp_core_srcs "${_source}")
  endif()
endforeach()
add_library(CATChem_core_cpp STATIC ${_cpp_core_srcs})
target_link_libraries(
  CATChem_core_cpp
  PUBLIC yaml-cpp::yaml-cpp catchem_validation_options
)
target_include_directories(CATChem_core_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Build provenance for diagnostic global attributes (feature 013, FR-011).
target_compile_definitions(
  CATChem_core_cpp
  PUBLIC CATCHEM_BUILD_VERSION="${CATCHEM_VERSION_FULL}" CATCHEM_BUILD_COMMIT="${CATCHEM_GIT_COMMIT}"
)

if(CATCHEM_ENABLE_KOKKOS)
  target_link_libraries(CATChem_core_cpp PUBLIC Kokkos::kokkos)
  target_compile_definitions(CATChem_core_cpp PUBLIC CATCHEM_ENABLE_KOKKOS)
else()
  target_link_libraries(CATChem_core_cpp PUBLIC mdspan::mdspan)
endif()

# Instead of removing CATChem_core, we make it an INTERFACE library that depends on CATChem_core_cpp
# and external libraries, so downstream targets still build without modification.
add_library(CATChem_core INTERFACE)
target_link_libraries(CATChem_core INTERFACE CATChem_core_cpp catchem_utilities)

if(TARGET GOCART2G)
  target_link_libraries(CATChem_core INTERFACE GOCART2G)
endif()
if(TARGET musica)
  target_link_libraries(CATChem_core INTERFACE musica::musica-fortran)
endif()
