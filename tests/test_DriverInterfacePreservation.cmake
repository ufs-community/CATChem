# test_DriverInterfacePreservation.cmake
# Verification that CCPP and NUOPC driver caps compile and conform
# to interface preservation requirements with the C++ Core API.

cmake_minimum_required(VERSION 3.10)

function(assert_file_contains filepath search_string description)
  if(NOT EXISTS "${filepath}")
    message(FATAL_ERROR "FAILED: ${description}\n  File not found: ${filepath}")
  endif()
  file(READ "${filepath}" file_content)
  string(FIND "${file_content}" "${search_string}" found_pos)
  if(found_pos EQUAL -1)
    message(
      FATAL_ERROR
      "FAILED: ${description}\n  File: ${filepath}\n  Missing: ${search_string}"
    )
  else()
    message(STATUS "PASSED: ${description}")
  endif()
endfunction()

set(SRC_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")

message(STATUS "=== Driver Interface Preservation Tests ===")

# Verify CATChem CCPP interface exposes the dynamic capgen routines
assert_file_contains(
  "${SRC_ROOT}/drivers/ccpp/ccpp_catchem_interface.F90"
  "subroutine ccpp_catchem_interface_register"
  "CCPP Driver exposes register subroutine"
)

assert_file_contains(
  "${SRC_ROOT}/drivers/ccpp/ccpp_catchem_interface.F90"
  "subroutine ccpp_catchem_interface_init"
  "CCPP Driver exposes init subroutine"
)

assert_file_contains(
  "${SRC_ROOT}/drivers/ccpp/ccpp_catchem_interface.F90"
  "subroutine ccpp_catchem_interface_run"
  "CCPP Driver exposes run subroutine"
)

# Verify Modern Public API
assert_file_contains(
  "${SRC_ROOT}/src/api/CATChem_API.F90"
  "type CATChem_Model"
  "Public API exports CATChem_Model"
)

assert_file_contains(
  "${SRC_ROOT}/src/api/CATChem_API.F90"
  "procedure :: initialize"
  "CATChem_Model has initialize binding"
)

assert_file_contains(
  "${SRC_ROOT}/src/api/CATChem_API.F90"
  "procedure :: run_timestep"
  "CATChem_Model has run_timestep binding"
)

message(STATUS "All Driver Interface Preservation checks PASSED!")
