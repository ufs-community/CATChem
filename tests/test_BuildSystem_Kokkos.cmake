# test_BuildSystem_Kokkos.cmake
# Build system verification tests for the CATCHEM_ENABLE_KOKKOS option.
#
# Architecture under test: the C++ core and all process C++ sources compile
# in BOTH build flavors. Kokkos is an acceleration option — with
# CATCHEM_ENABLE_KOKKOS=OFF the standalone kokkos/mdspan library provides
# the Kokkos::mdspan data layer and kernels run serially on the host.
#
# Test 1: CATCHEM_ENABLE_KOKKOS option exists, defaults to OFF, and the
#         deprecated ENABLE_KOKKOS spelling is honored via a shim
# Test 2: ON configures C++20 and locates/fetches Kokkos
# Test 3: OFF fetches the standalone mdspan library
# Test 4: The C++ core is built unconditionally; only the Kokkos link and
#         compile definition are conditional
# Test 5: Process C++ sources are unconditional and carry no Kokkos guards
# Test 6: The Kokkos compat header provides the host-only fallbacks
# Test 7: The C++ test suite is not gated on Kokkos
#
# Usage:
#   cmake -P tests/test_BuildSystem_Kokkos.cmake

cmake_minimum_required(VERSION 3.10)

# Helper: check that a file contains a specific string
function(assert_file_contains filepath search_string description)
  file(READ "${filepath}" file_content)
  string(FIND "${file_content}" "${search_string}" found_pos)
  if(found_pos EQUAL -1)
    message(
      FATAL_ERROR
      "FAILED: ${description}\n  File: ${filepath}\n  Expected to find: ${search_string}"
    )
  else()
    message(STATUS "PASSED: ${description}")
  endif()
endfunction()

# Helper: check that a file does NOT contain a specific string
function(assert_file_not_contains filepath search_string description)
  file(READ "${filepath}" file_content)
  string(FIND "${file_content}" "${search_string}" found_pos)
  if(NOT found_pos EQUAL -1)
    message(
      FATAL_ERROR
      "FAILED: ${description}\n  File: ${filepath}\n  Should NOT contain: ${search_string}"
    )
  else()
    message(STATUS "PASSED: ${description}")
  endif()
endfunction()

message(STATUS "")
message(STATUS "=== Build System Tests: CATCHEM_ENABLE_KOKKOS ===")
message(STATUS "")

# Determine source root (this script is in tests/)
get_filename_component(SRC_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

# --- Test 1: CATCHEM_ENABLE_KOKKOS option configuration ---
message(STATUS "--- Test 1: CATCHEM_ENABLE_KOKKOS option configuration ---")

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "option(CATCHEM_ENABLE_KOKKOS"
  "Top-level CMakeLists.txt defines CATCHEM_ENABLE_KOKKOS option"
)

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "CATCHEM_ENABLE_KOKKOS \"Enable Kokkos GPU/parallel support\" OFF"
  "CATCHEM_ENABLE_KOKKOS defaults to OFF"
)

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "ENABLE_KOKKOS is deprecated; use CATCHEM_ENABLE_KOKKOS"
  "Deprecated ENABLE_KOKKOS spelling is mapped to CATCHEM_ENABLE_KOKKOS"
)

# --- Test 2: CATCHEM_ENABLE_KOKKOS=ON configuration ---
message(STATUS "--- Test 2: CATCHEM_ENABLE_KOKKOS=ON configuration ---")

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "CMAKE_CXX_STANDARD 20"
  "MUSICA builds use C++20"
)

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "CMAKE_CXX_STANDARD 17"
  "Non-MUSICA builds use C++17"
)

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "find_package(Kokkos"
  "ON path looks for a system Kokkos"
)

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "FetchContent_MakeAvailable(kokkos)"
  "ON path falls back to fetching a pinned Kokkos"
)

# --- Test 3: OFF uses the vendored mdspan (no network at configure) ---
message(STATUS "--- Test 3: CATCHEM_ENABLE_KOKKOS=OFF configuration ---")

if(NOT EXISTS "${SRC_ROOT}/src/external/mdspan/include/mdspan/mdspan.hpp")
  message(
    FATAL_ERROR
    "FAILED: vendored mdspan single header is missing\n  Expected: src/external/mdspan/include/mdspan/mdspan.hpp"
  )
else()
  message(STATUS "PASSED: vendored mdspan single header is present")
endif()

assert_file_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "src/external/mdspan/include"
  "OFF path points the mdspan target at the vendored header"
)

assert_file_not_contains(
  "${SRC_ROOT}/CMakeLists.txt"
  "FetchContent_MakeAvailable(mdspan)"
  "OFF path performs no mdspan fetch (network-restricted UFS nodes)"
)

# --- Test 4: C++ core is unconditional; Kokkos link is conditional ---
message(STATUS "--- Test 4: C++ core built in both flavors ---")

assert_file_contains(
  "${SRC_ROOT}/src/core/CMakeLists.txt"
  "add_library(CATChem_core_cpp STATIC"
  "CATChem_core_cpp is declared unconditionally"
)

assert_file_contains(
  "${SRC_ROOT}/src/core/CMakeLists.txt"
  "target_compile_definitions(CATChem_core_cpp PUBLIC CATCHEM_ENABLE_KOKKOS)"
  "Kokkos compile definition is exported when ON"
)

assert_file_contains(
  "${SRC_ROOT}/src/core/CMakeLists.txt"
  "mdspan::mdspan"
  "OFF path links the standalone mdspan target"
)

# --- Test 5: process C++ sources unconditional, no Kokkos guards ---
message(STATUS "--- Test 5: Process C++ sources unconditional ---")

foreach(
  process
  seasalt
  dust
  drydep
  wetdep
  settling
  so4chem
  carbchem
)
  set(_process_cml "${SRC_ROOT}/src/process/${process}/CMakeLists.txt")
  assert_file_contains(
    "${_process_cml}"
    "catchem_process_${process}.cpp"
    "${process} compiles its C++ source in both flavors"
  )
  assert_file_contains(
    "${_process_cml}"
    "CATChem_core_cpp"
    "${process} links CATChem_core_cpp for transitive Kokkos/mdspan usage"
  )
  assert_file_not_contains(
    "${_process_cml}"
    "if(CATCHEM_ENABLE_KOKKOS)"
    "${process} has no CMake-level Kokkos conditionality"
  )
  assert_file_not_contains(
    "${_process_cml}"
    "if(ENABLE_KOKKOS)"
    "${process} does not use the deprecated ENABLE_KOKKOS guard"
  )
endforeach()

# --- Test 6: compat header fallbacks ---
message(STATUS "--- Test 6: Kokkos compat header ---")

assert_file_contains(
  "${SRC_ROOT}/src/core/catchem_kokkos_compat.hpp"
  "#define KOKKOS_INLINE_FUNCTION inline"
  "Compat header masks KOKKOS_INLINE_FUNCTION for host-only builds"
)

assert_file_contains(
  "${SRC_ROOT}/src/core/catchem_kokkos_compat.hpp"
  "#include <mdspan/mdspan.hpp>"
  "Compat header provides mdspan for host-only builds"
)

assert_file_contains(
  "${SRC_ROOT}/src/core/catchem_kokkos_compat.hpp"
  "#define MDSPAN_IMPL_STANDARD_NAMESPACE Kokkos"
  "Compat header pins the vendored mdspan namespace to Kokkos"
)

# --- Test 7: C++ tests not gated on Kokkos ---
message(STATUS "--- Test 7: C++ test suite runs in both flavors ---")

assert_file_contains(
  "${SRC_ROOT}/tests/CMakeLists.txt"
  "add_executable(test_catchem_interop"
  "Interop test suite is registered"
)

assert_file_not_contains(
  "${SRC_ROOT}/tests/CMakeLists.txt"
  "if(CATCHEM_ENABLE_KOKKOS)"
  "tests/CMakeLists.txt does not gate the C++ suite on Kokkos"
)

message(STATUS "")
message(STATUS "=== All build system tests PASSED ===")
message(STATUS "")
