if(NOT DEFINED CATCHEM_SOURCE_DIR OR NOT DEFINED CATCHEM_BINARY_ROOT)
  message(FATAL_ERROR "CATCHEM_SOURCE_DIR and CATCHEM_BINARY_ROOT are required")
endif()

execute_process(
  COMMAND
    ${CMAKE_COMMAND} -S "${CATCHEM_SOURCE_DIR}" -B
    "${CATCHEM_BINARY_ROOT}/canonical" -DCATCHEM_BUILD_TESTING=OFF
    -DCATCHEM_BUILD_NUOPC=OFF -DCATCHEM_BUILD_CCPP=OFF
    -DCATCHEM_ENABLE_KOKKOS=OFF
  RESULT_VARIABLE canonical_result
  OUTPUT_QUIET
  ERROR_VARIABLE canonical_error
)
if(NOT canonical_result EQUAL 0)
  message(FATAL_ERROR "Canonical options failed: ${canonical_error}")
endif()

execute_process(
  COMMAND
    ${CMAKE_COMMAND} -S "${CATCHEM_SOURCE_DIR}" -B
    "${CATCHEM_BINARY_ROOT}/conflict" -DCATCHEM_BUILD_TESTING=OFF
    -DCATCHEM_BUILD_NUOPC=ON -DBUILD_NUOPC=OFF
  RESULT_VARIABLE conflict_result
  OUTPUT_QUIET
  ERROR_VARIABLE conflict_error
)
if(conflict_result EQUAL 0 OR NOT conflict_error MATCHES "conflicts with")
  message(
    FATAL_ERROR
    "Contradictory NUOPC options were not rejected with migration guidance"
  )
endif()

execute_process(
  COMMAND
    ${CMAKE_COMMAND} -S "${CATCHEM_SOURCE_DIR}" -B
    "${CATCHEM_BINARY_ROOT}/trace" -DCATCHEM_BUILD_TESTING=OFF
    -DCATCHEM_BUILD_NUOPC=OFF -DCATCHEM_TRACE_NUOPC=ON
  RESULT_VARIABLE trace_result
  OUTPUT_QUIET
  ERROR_VARIABLE trace_error
)
if(
  trace_result EQUAL 0
  OR NOT trace_error MATCHES "requires CATCHEM_BUILD_NUOPC"
)
  message(FATAL_ERROR "Unsupported NUOPC tracing combination was not rejected")
endif()

execute_process(
  COMMAND
    ${CMAKE_COMMAND} -S "${CATCHEM_SOURCE_DIR}" -B
    "${CATCHEM_BINARY_ROOT}/ccpp-missing-host" -DCATCHEM_BUILD_TESTING=OFF
    -DCATCHEM_BUILD_NUOPC=OFF -DCATCHEM_BUILD_CCPP=ON
  RESULT_VARIABLE ccpp_result
  OUTPUT_QUIET
  ERROR_VARIABLE ccpp_error
)
if(ccpp_result EQUAL 0 OR NOT ccpp_error MATCHES "CATCHEM_CCPP_HOST_MODULE_DIR")
  message(
    FATAL_ERROR
    "CCPP without host modules was not rejected with actionable guidance"
  )
endif()
