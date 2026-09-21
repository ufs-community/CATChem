foreach(opt IN LISTS _catchem_compiler_options)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${opt}>")
endforeach()
foreach(opt IN LISTS _catchem_cxx_warning_options)
  add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:${opt}>")
endforeach()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/catchem_test_config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/catchem_test_config.hpp"
)

# Attach one or more semicolon-separated integrity labels to tests while
# rejecting misspelled test names at configure time.
function(catchem_label_tests labels)
  foreach(test_name IN LISTS ARGN)
    if(NOT TEST ${test_name})
      message(FATAL_ERROR "Cannot label unknown CATChem test '${test_name}'")
    endif()
    set_property(TEST ${test_name} APPEND PROPERTY LABELS "${labels}")
  endforeach()
endfunction()

# Architecture tests use one helper so labels and core linkage cannot drift.
function(catchem_add_architecture_test test_name source labels)
  add_executable(${test_name} ${source})
  target_link_libraries(${test_name} PRIVATE CATChem_core_cpp)
  add_test(NAME ${test_name} COMMAND ${test_name})
  catchem_label_tests("${labels}" ${test_name})
endfunction()

# Built only when explicitly requested; timing never affects normal CTest pass/fail.
add_executable(
  benchmark_core_orchestration
  EXCLUDE_FROM_ALL
  benchmark_core_orchestration.cpp
)
target_link_libraries(benchmark_core_orchestration PRIVATE CATChem_core_cpp)

catchem_add_architecture_test(test_mechanism_independence test_mechanism_independence.cpp "mechanism;configuration")
catchem_add_architecture_test(test_multi_mechanism_instances test_multi_mechanism_instances.cpp "mechanism;lifecycle")
catchem_add_architecture_test(test_process_contracts test_process_contracts.cpp "process-contract;dataflow")
target_link_libraries(
  test_process_contracts
  PRIVATE
    CATChem_process_dust
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
)
target_include_directories(
  test_process_contracts
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/process/dust
    ${PROJECT_SOURCE_DIR}/src/process/seasalt
    ${PROJECT_SOURCE_DIR}/src/process/drydep
    ${PROJECT_SOURCE_DIR}/src/process/wetdep
    ${PROJECT_SOURCE_DIR}/src/process/so4chem
)
catchem_add_architecture_test(test_execution_plan test_execution_plan.cpp "process-contract;dataflow")
catchem_add_architecture_test(test_selective_sync test_selective_sync.cpp "process-contract;dataflow;performance")
catchem_add_architecture_test(test_process_lifecycle test_process_lifecycle.cpp "lifecycle;concurrency")
catchem_add_architecture_test(test_process_registry_concurrency test_process_registry_concurrency.cpp "lifecycle;concurrency")
catchem_add_architecture_test(test_handle_concurrency test_handle_concurrency.cpp "lifecycle;concurrency;boundary")
catchem_add_architecture_test(test_physical_validation test_physical_validation.cpp "physical-validation;dataflow")
catchem_add_architecture_test(test_boundary_errors test_boundary_errors.cpp "boundary;lifecycle;concurrency")
target_link_libraries(
  test_boundary_errors
  PRIVATE
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
catchem_add_architecture_test(test_c_api_compatibility test_c_api_compatibility.cpp "boundary;compatibility")
catchem_add_architecture_test(test_extension_process test_extension_process.cpp "mechanism;process-contract;lifecycle;dataflow")
catchem_add_architecture_test(test_core_construction test_core_construction.cpp "mechanism;configuration;lifecycle")
# Source paths come from the configured catchem_test_config.hpp rather than a
# quoted -D define: the Intel MPI compiler wrappers strip the quotes when the
# define is passed on the command line (seen in the ci-intel host profile),
# and test_extension_process never used the define.
target_include_directories(
  test_core_construction
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(
  NAME check_core_mechanism_independence
  COMMAND
    ${CMAKE_COMMAND} -DCATCHEM_SOURCE_DIR=${PROJECT_SOURCE_DIR} -P
    ${CMAKE_CURRENT_SOURCE_DIR}/check_core_mechanism_independence.cmake
)
catchem_label_tests("mechanism;configuration" check_core_mechanism_independence)
target_link_libraries(
  test_c_api_compatibility
  PRIVATE
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
target_link_libraries(
  test_handle_concurrency
  PRIVATE
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# Automatically link the musica configs folder so that CTest can find relative NetCDF files on GHA CI
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/configs")
  file(
    CREATE_LINK
      "${PROJECT_SOURCE_DIR}/src/external/musica/configs"
      "${CMAKE_CURRENT_BINARY_DIR}/configs"
    SYMBOLIC
  )
endif()

if(CATCHEM_BUILD_CCPP)
  add_executable(test_ccpp_integration test_ccpp_integration.f90)
  target_link_libraries(test_ccpp_integration PRIVATE CCPP_CATChem_interface)
  target_include_directories(
    test_ccpp_integration
    PRIVATE ${CMAKE_BINARY_DIR}/include
  )
  set_target_properties(
    test_ccpp_integration
    PROPERTIES LINKER_LANGUAGE Fortran
  )
  add_test(NAME test_ccpp_integration COMMAND test_ccpp_integration)
  catchem_label_tests("integration;dataflow;boundary;performance" test_ccpp_integration)
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/host_conformance
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

# Foundational Fortran testing helper library.  Restored for the spec-011 US4
# constants-parity gate (test_Constants).
add_library(testing testing_mod.f90)
set_target_properties(
  testing
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_link_libraries(testing PUBLIC catchem_utilities)

# Restored Low-Level Foundational Unit Tests
macro(add_fortran_test name_test)
  add_executable(${name_test} ${name_test}.f90)
  target_link_libraries(${name_test} PRIVATE CATChem_core testing)
  set_target_properties(
    ${name_test}
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(${name_test} PRIVATE ${CMAKE_BINARY_DIR}/include)
  add_dependencies(${name_test} CATChem_core)
  add_test(
    NAME ${name_test}
    COMMAND ${name_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endmacro()

# Constants parity is an explicit spec-011 US4 gate.  The other low-level
# Fortran core unit tests were retired: their behavior is covered by
# test_catchem_properties and the C++ suite.
add_fortran_test(test_Constants)

# High-level API init-sequence test — mirrors the NUOPC cap's usage of
# CATChem_Model (facade construction, met binding, timestep)
if(TARGET CATChem)
  add_executable(test_CATChem_API test_CATChem_API.f90)
  target_link_libraries(test_CATChem_API PRIVATE CATChem CATChem_core)
  set_target_properties(
    test_CATChem_API
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(
    test_CATChem_API
    PRIVATE ${CMAKE_BINARY_DIR}/include
  )
  add_dependencies(test_CATChem_API CATChem)
  add_test(
    NAME test_CATChem_API
    COMMAND test_CATChem_API
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()

# C++ Interop and Integration Suite — runs in both build flavors; the
# CATChem_core_cpp link carries Kokkos (or standalone mdspan) transitively.
add_executable(test_catchem_logger test_catchem_logger.cpp)
target_link_libraries(test_catchem_logger PRIVATE CATChem_core_cpp)
add_test(NAME test_catchem_logger COMMAND test_catchem_logger)

add_executable(test_catchem_error test_catchem_error.cpp)
target_link_libraries(test_catchem_error PRIVATE CATChem_core_cpp)
add_test(NAME test_catchem_error COMMAND test_catchem_error)

add_executable(test_catchem_interop test_catchem_interop.cpp)
target_link_libraries(
  test_catchem_interop
  PRIVATE
    CATChem_core_cpp
    CATChem_core
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
set_target_properties(test_catchem_interop PROPERTIES LINKER_LANGUAGE CXX)
add_test(NAME test_catchem_interop COMMAND test_catchem_interop)

add_executable(test_catchem_dataflow test_catchem_dataflow.cpp)
target_link_libraries(
  test_catchem_dataflow
  PRIVATE
    CATChem_core_cpp
    CATChem_core
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
set_target_properties(
  test_catchem_dataflow
  PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
)
add_test(NAME test_catchem_dataflow COMMAND test_catchem_dataflow)
set_tests_properties(test_catchem_dataflow PROPERTIES LABELS "dataflow;parity")
if(CATCHEM_ENABLE_KOKKOS)
  set_tests_properties(
    test_catchem_dataflow
    PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=4;KOKKOS_NUM_THREADS=4;OMP_PROC_BIND=false"
  )
endif()

add_executable(test_c_api_negative test_c_api_negative.cpp)
target_link_libraries(
  test_c_api_negative
  PRIVATE
    CATChem_core_cpp
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
add_test(NAME test_c_api_negative COMMAND test_c_api_negative)
catchem_label_tests("boundary" test_c_api_negative)
catchem_label_tests("dataflow;parity" test_catchem_dataflow)

add_executable(test_runtime_lifecycle test_runtime_lifecycle.cpp)
target_link_libraries(
  test_runtime_lifecycle
  PRIVATE
    CATChem_core_cpp
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
add_test(NAME test_runtime_lifecycle COMMAND test_runtime_lifecycle)
catchem_label_tests("lifecycle;dataflow" test_runtime_lifecycle)

add_executable(test_config_validation test_config_validation.cpp)
target_link_libraries(
  test_config_validation
  PRIVATE
    CATChem_core_cpp
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
# Fixture paths come from the configured catchem_test_config.hpp (see
# test_core_construction above for why quoted -D defines are avoided).
target_include_directories(
  test_config_validation
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(NAME test_config_validation COMMAND test_config_validation)
catchem_label_tests("configuration" test_config_validation)

add_executable(test_diagnostic_lifecycle test_diagnostic_lifecycle.cpp)
target_link_libraries(
  test_diagnostic_lifecycle
  PRIVATE
    CATChem_core_cpp
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
add_test(NAME test_diagnostic_lifecycle COMMAND test_diagnostic_lifecycle)
catchem_label_tests("diagnostics;dataflow" test_diagnostic_lifecycle)

add_executable(test_timestep_failure test_timestep_failure.cpp)
target_link_libraries(
  test_timestep_failure
  PRIVATE
    CATChem_core_cpp
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
add_test(NAME test_timestep_failure COMMAND test_timestep_failure)
catchem_label_tests("lifecycle;boundary;dataflow" test_timestep_failure)

add_test(
  NAME test_build_options
  COMMAND
    ${CMAKE_COMMAND} -DCATCHEM_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DCATCHEM_BINARY_ROOT=${CMAKE_CURRENT_BINARY_DIR}/build-option-tests -P
    ${CMAKE_CURRENT_SOURCE_DIR}/test_build_options.cmake
)
catchem_label_tests("configuration;build" test_build_options)

add_executable(test_separate_memory test_separate_memory.cpp)
target_link_libraries(test_separate_memory PRIVATE CATChem_core_cpp)
add_test(NAME test_separate_memory COMMAND test_separate_memory)
catchem_label_tests("dataflow" test_separate_memory)
if(
  CATCHEM_ENABLE_KOKKOS
  AND (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL)
)
  catchem_label_tests("separate-memory" test_separate_memory)
endif()

add_executable(
  test_host_integration_conformance
  test_host_integration_conformance.cpp
)
target_link_libraries(
  test_host_integration_conformance
  PRIVATE
    CATChem_core_cpp
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
# Fixture paths come from the configured catchem_test_config.hpp (see
# test_core_construction above for why quoted -D defines are avoided).
target_include_directories(
  test_host_integration_conformance
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)
add_test(
  NAME test_host_integration_conformance
  COMMAND test_host_integration_conformance
)
catchem_label_tests("boundary;integration;dataflow;conformance" test_host_integration_conformance)

# Existing process/property tests are the unchanged science-regression gate.
set(
  _catchem_science_regression_tests
  test_catchem_seasalt
  test_catchem_dust
  test_catchem_drydep
  test_catchem_wetdep
  test_catchem_settling
  test_catchem_so4chem
  test_catchem_carbchem
  test_catchem_properties
  test_catchem_scenarios
)

# Individual C++ Process Unit Tests
macro(add_process_unit_test process_name target_lib)
  set(test_name "test_catchem_${process_name}")
  add_executable(${test_name} "${test_name}.cpp")
  target_link_libraries(
    ${test_name}
    PRIVATE CATChem_core_cpp CATChem_core ${target_lib}
  )
  set_target_properties(
    ${test_name}
    PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
  )
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endmacro()

add_process_unit_test(seasalt CATChem_process_seasalt)
add_process_unit_test(dust CATChem_process_dust)
add_process_unit_test(drydep CATChem_process_drydep)
add_process_unit_test(wetdep CATChem_process_wetdep)
add_process_unit_test(settling CATChem_process_settling)
add_process_unit_test(so4chem CATChem_process_so4chem)
add_process_unit_test(carbchem CATChem_process_carbchem)

# Runtime YAML scheme options must reach the science layer and unknown keys
# must fail validation (settling scale_factor round-trip).
add_process_unit_test(option_propagation CATChem_process_settling)

# legacy optics-table (simple_scheme: true) settling path and the
# top-level "mie:" configuration parsing.  Both read the generated
# catchem_test_config.hpp (optics fixture directory), so the binary dir is
# added to the include path explicitly.
add_executable(test_catchem_settling_optics test_catchem_settling_optics.cpp)
target_link_libraries(
  test_catchem_settling_optics
  PRIVATE CATChem_core_cpp CATChem_core CATChem_process_settling
)
target_include_directories(
  test_catchem_settling_optics
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)
set_target_properties(
  test_catchem_settling_optics
  PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
)
add_test(
  NAME test_catchem_settling_optics
  COMMAND test_catchem_settling_optics
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
catchem_label_tests("process;settling;optics" test_catchem_settling_optics)

add_executable(test_catchem_config_mie test_catchem_config_mie.cpp)
target_link_libraries(
  test_catchem_config_mie
  PRIVATE CATChem_core_cpp CATChem_core
)
target_include_directories(
  test_catchem_config_mie
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)
set_target_properties(
  test_catchem_config_mie
  PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
)
add_test(NAME test_catchem_config_mie COMMAND test_catchem_config_mie)
catchem_label_tests("configuration" test_catchem_config_mie)

add_executable(test_catchem_derived_met test_catchem_derived_met.cpp)
target_link_libraries(
  test_catchem_derived_met
  PRIVATE CATChem_core_cpp CATChem_core
)
set_target_properties(
  test_catchem_derived_met
  PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
)
add_test(
  NAME test_catchem_derived_met
  COMMAND test_catchem_derived_met
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
catchem_label_tests("unit;cpp;derived_met" test_catchem_derived_met)

# Indexing-contract regression lock: CAP->core->bridge flat column-major order
# and interface extents (spec 011 US5).
add_executable(test_catchem_index_contract test_catchem_index_contract.cpp)
target_link_libraries(
  test_catchem_index_contract
  PRIVATE CATChem_core_cpp CATChem_core
)
set_target_properties(
  test_catchem_index_contract
  PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
)
add_test(
  NAME test_catchem_index_contract
  COMMAND test_catchem_index_contract
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
catchem_label_tests("unit;cpp;dataflow;index_contract" test_catchem_index_contract)

# Modern Randomized Property-Based Testing
add_executable(test_catchem_properties test_catchem_properties.cpp)
target_link_libraries(
  test_catchem_properties
  PRIVATE
    CATChem_core_cpp
    CATChem_core
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
set_target_properties(test_catchem_properties PROPERTIES LINKER_LANGUAGE CXX)
file(
  COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_config.yml
    ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_new_config.yml
    ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_species.yml
    ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_emission.yml
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/Configs/Default
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Configs
)
add_test(
  NAME test_catchem_properties
  COMMAND test_catchem_properties
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# The Default configuration is the production-like full process schedule used
# by the legacy-vs-C++ numerical parity harness.  Keep a focused contract test
# separate from the smaller root-level configuration used by API unit tests.
add_executable(test_default_config_contract test_default_config_contract.cpp)
target_link_libraries(test_default_config_contract PRIVATE CATChem_core_cpp)
add_test(
  NAME test_default_config_contract
  COMMAND test_default_config_contract
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
catchem_label_tests("configuration;parity" test_default_config_contract)

# Standalone candidate executable used by run_default_numerical_parity.py.
# It is intentionally not a normal CTest: the matching legacy executable is
# built in a separate upstream/develop worktree.
add_executable(
  default_parity_candidate_runner
  default_parity_candidate_runner.cpp
)
target_link_libraries(
  default_parity_candidate_runner
  PRIVATE
    CATChem_core_cpp
    CATChem_process_seasalt
    CATChem_process_dust
    CATChem_process_carbchem
    CATChem_process_settling
    CATChem_process_drydep
    CATChem_process_so4chem
    CATChem_process_wetdep
)

# The comparator is intentionally independent of either core: legacy and C++
# runners write normalized JSON snapshots, so they cannot collide on Fortran
# module names or symbols.  This fixture validates the comparator in normal
# CI; use it with outputs from separately-built worktrees for science parity.
find_program(PYTHON_EXECUTABLE NAMES python3 python)
if(PYTHON_EXECUTABLE)
  add_test(
    NAME test_default_parity_met_profile
    COMMAND
      ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/build_parity_met_profile.py --profile
      ${CMAKE_CURRENT_SOURCE_DIR}/MetProfiles/Profile_NCWCP.csv --levels 64
      --columns 1 --output
      ${CMAKE_CURRENT_BINARY_DIR}/default_parity_met_profile.json
  )
  catchem_label_tests("parity" test_default_parity_met_profile)

  add_test(
    NAME test_parity_snapshot_comparator
    COMMAND
      ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/compare_parity_snapshots.py
      ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/parity/legacy.json
      ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/parity/candidate.json
  )
  catchem_label_tests("parity" test_parity_snapshot_comparator)

  # C++ equivalents of the legacy per-process integration tests.  Each case
  # enters through the C++ core, binds the shared MET/species fixture, and
  # runs exactly one process from the Default configuration.
  foreach(
    process_name
    seasalt
    dust
    carbchem
    settling
    drydep
    so4chem
    wetdep
  )
    add_test(
      NAME test_cpp_${process_name}_integration
      COMMAND
        ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/run_cpp_process_integration.py
        --candidate-runner $<TARGET_FILE:default_parity_candidate_runner>
        --process ${process_name} --source-root ${CMAKE_SOURCE_DIR}
    )
    catchem_label_tests("integration;cpp;${process_name}" test_cpp_${process_name}_integration)
  endforeach()

  # Match the three scheme variants exercised by the upstream Fortran
  # sea-salt integration test.  They are separate YAML-configured C++ cases.
  foreach(seasalt_scheme gong97 gong03 geos12)
    add_test(
      NAME test_cpp_seasalt_${seasalt_scheme}_integration
      COMMAND
        ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/run_cpp_process_integration.py
        --candidate-runner $<TARGET_FILE:default_parity_candidate_runner>
        --process seasalt --seasalt-scheme ${seasalt_scheme} --source-root
        ${CMAKE_SOURCE_DIR}
    )
    catchem_label_tests("integration;cpp;seasalt;${seasalt_scheme}" test_cpp_seasalt_${seasalt_scheme}_integration)
  endforeach()

  # Missing-host-input guard: dropping PMID from the staged MET profile must
  # abort the settling process with a named error (spec 011 US1).
  add_test(
    NAME test_cpp_settling_missing_pmid
    COMMAND
      ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/run_cpp_process_integration.py
      --candidate-runner $<TARGET_FILE:default_parity_candidate_runner>
      --process settling --drop-field PMID --source-root ${CMAKE_SOURCE_DIR}
  )
  catchem_label_tests("integration;cpp;settling" test_cpp_settling_missing_pmid)
endif()

# Multi-Atmospheric Scenario Property Test
add_executable(test_catchem_scenarios test_catchem_scenarios.cpp)
target_link_libraries(
  test_catchem_scenarios
  PRIVATE
    CATChem_core_cpp
    CATChem_core
    CATChem_process_settling
    CATChem_process_seasalt
    CATChem_process_drydep
    CATChem_process_wetdep
    CATChem_process_so4chem
    CATChem_process_dust
    CATChem_process_carbchem
)
set_target_properties(
  test_catchem_scenarios
  PROPERTIES CXX_STANDARD 20 LINKER_LANGUAGE CXX
)
add_test(
  NAME test_catchem_scenarios
  COMMAND test_catchem_scenarios
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

foreach(_science_test IN LISTS _catchem_science_regression_tests)
  if(TEST ${_science_test})
    set_property(
      TEST ${_science_test}
      APPEND
      PROPERTY LABELS "science-regression"
    )
  endif()
endforeach()

# Static Fortran Core Import Check
find_program(PYTHON_EXECUTABLE NAMES python3 python)
if(PYTHON_EXECUTABLE)
  add_test(
    NAME test_compare_nuopc_tracer_extrema
    COMMAND
      ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/test_compare_nuopc_tracer_extrema.py
  )
  catchem_label_tests("regression;boundary" test_compare_nuopc_tracer_extrema)
endif()
if(PYTHON_EXECUTABLE)
  add_test(
    NAME test_check_fortran_core_imports
    COMMAND
      ${PYTHON_EXECUTABLE}
      "${PROJECT_SOURCE_DIR}/tools/check_fortran_core_imports.py"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )
endif()

if(CATCHEM_BUILD_MUSICA)
  # Photolysis Process Unit/Integration Test
  add_executable(test_catchem_photolysis test_catchem_photolysis.cpp)
  target_link_libraries(
    test_catchem_photolysis
    PRIVATE
      CATChem_core_cpp
      CATChem_core
      catchem_process_photolysis
      yaml-cpp::yaml-cpp
  )
  set_target_properties(test_catchem_photolysis PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME test_catchem_photolysis COMMAND test_catchem_photolysis)

  # Coupled Process Integration Test
  add_executable(test_catchem_gaschem test_catchem_gaschem.cpp)
  target_link_libraries(
    test_catchem_gaschem
    PUBLIC
      CATChem_core_cpp
      CATChem_core
      catchem_process_photolysis
      catchem_process_gaschem
  )
  set_target_properties(test_catchem_gaschem PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME test_catchem_gaschem COMMAND test_catchem_gaschem)
endif()

add_executable(test_catchem_gaschem_units test_catchem_gaschem_units.cpp)
add_test(NAME test_catchem_gaschem_units COMMAND test_catchem_gaschem_units)

# ESMF-backed NUOPC transform test — builds and runs ONLY when the cap is
# built (CATCHEM_BUILD_NUOPC=ON). cc_nuopc links ESMF/NetCDF/MPI PUBLIC,
# so its usage requirements arrive transitively; its Fortran modules live
# in the cap's module directory.
if(CATCHEM_BUILD_NUOPC)
  # libesmf.a references the MPI C++ bindings; an executable link needs
  # them explicitly (the cap static library never resolves them itself).
  find_package(MPI COMPONENTS Fortran CXX REQUIRED)
  add_executable(test_nuopc_transform test_nuopc_transform.f90)
  target_link_libraries(
    test_nuopc_transform
    PRIVATE cc_nuopc CATChem CATChem_core MPI::MPI_Fortran MPI::MPI_CXX
  )
  set_target_properties(
    test_nuopc_transform
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(
    test_nuopc_transform
    PRIVATE
      ${CMAKE_BINARY_DIR}/include
      $<TARGET_PROPERTY:NUOPC_CATChem_cap,Fortran_MODULE_DIRECTORY>
  )
  # Stage the production field mapping and the test model configs into the
  # test's working directory.
  file(
    COPY
      ${PROJECT_SOURCE_DIR}/drivers/nuopc/CATChem_field_mapping.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_new_config.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_species.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_emission.yml
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/host_conformance
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_test(
    NAME test_nuopc_transform
    COMMAND test_nuopc_transform
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  catchem_label_tests("integration;dataflow;boundary" test_nuopc_transform)

  add_executable(nuopc_contract_harness nuopc_contract_harness.F90)
  target_link_libraries(
    nuopc_contract_harness
    PRIVATE cc_nuopc CATChem CATChem_core MPI::MPI_Fortran MPI::MPI_CXX
  )
  set_target_properties(
    nuopc_contract_harness
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(
    nuopc_contract_harness
    PRIVATE
      ${CMAKE_BINARY_DIR}/include
      $<TARGET_PROPERTY:NUOPC_CATChem_cap,Fortran_MODULE_DIRECTORY>
  )
  # Give the lifecycle harness its own run directory.  The cap resolves its
  # YAML paths relative to the embedding application's working directory, as
  # it does under UFS; staging the process-free test YAML under the cap's
  # conventional filename keeps this test hermetic without changing the
  # production configuration used by test_nuopc_transform.
  set(NUOPC_CONTRACT_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/nuopc_contract_run)
  file(MAKE_DIRECTORY ${NUOPC_CONTRACT_RUN_DIR})
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_nuopc_lifecycle_config.yml
    ${NUOPC_CONTRACT_RUN_DIR}/CATChem_new_config.yml
    COPYONLY
  )
  file(
    COPY
      ${PROJECT_SOURCE_DIR}/drivers/nuopc/CATChem_field_mapping.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_species.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_emission.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/nuopc_contract_run.yml
    DESTINATION ${NUOPC_CONTRACT_RUN_DIR}
  )
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/host_conformance
    DESTINATION ${NUOPC_CONTRACT_RUN_DIR}
  )
  add_test(
    NAME nuopc_contract_valid
    COMMAND nuopc_contract_harness
    WORKING_DIRECTORY ${NUOPC_CONTRACT_RUN_DIR}
  )
  # A lifecycle regression must fail rather than leave CI or a developer's
  # test run waiting indefinitely in a NUOPC phase.
  set_tests_properties(nuopc_contract_valid PROPERTIES TIMEOUT 30)
  catchem_label_tests("integration;contract;boundary" nuopc_contract_valid)

  # Feature 013: generic process-diagnostic writer round-trip.  This links the
  # NUOPC cap directly (not the self-contained cc_nuopc archive) so the build
  # graph does not duplicate the drivers/nuopc Fortran modules; the cap already
  # propagates ESMF/NetCDF/MPI usage requirements, and its module directory is
  # the same include/ the cap writes into.
  add_executable(test_nuopc_diag_output test_nuopc_diag_output.f90)
  target_link_libraries(
    test_nuopc_diag_output
    PRIVATE NUOPC_CATChem_cap CATChem CATChem_core MPI::MPI_Fortran MPI::MPI_CXX
  )
  set_target_properties(
    test_nuopc_diag_output
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
      LINKER_LANGUAGE Fortran
  )
  target_include_directories(
    test_nuopc_diag_output
    PRIVATE
      ${CMAKE_BINARY_DIR}/include
      $<TARGET_PROPERTY:NUOPC_CATChem_cap,Fortran_MODULE_DIRECTORY>
  )
  # Stage the all-process diagnostic config and the species/emission inventories
  # it references into the test's working directory.
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_diag_output_config.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_diag_narrow_config.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_species.yml
      ${CMAKE_CURRENT_SOURCE_DIR}/CATChem_emission.yml
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_test(
    NAME test_nuopc_diag_output
    COMMAND test_nuopc_diag_output
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  catchem_label_tests("integration;diagnostics" test_nuopc_diag_output)
endif()
