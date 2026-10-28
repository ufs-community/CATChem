set(
  generic_core_files
  "${CATCHEM_SOURCE_DIR}/src/core/catchem_core.cpp"
  "${CATCHEM_SOURCE_DIR}/src/core/catchem_execution_plan.cpp"
  "${CATCHEM_SOURCE_DIR}/src/core/catchem_execution_plan.hpp"
  "${CATCHEM_SOURCE_DIR}/src/core/catchem_state_manager.hpp"
)

# Generic orchestration may consume descriptor roles/capabilities, but must not
# recognize particular species or select from a built-in process allowlist.
set(
  forbidden_pattern
  "(O3|NO2|SO2|SO4|DMS|gaschem|carbchem|seasalt|drydep|wetdep|settling)"
)
foreach(source IN LISTS generic_core_files)
  file(READ "${source}" contents)
  if(contents MATCHES "${forbidden_pattern}")
    message(
      FATAL_ERROR
      "Mechanism/process literal '${CMAKE_MATCH_1}' found in generic orchestration file ${source}"
    )
  endif()
endforeach()
