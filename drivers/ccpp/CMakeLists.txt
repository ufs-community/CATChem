set(_srcs ccpp_catchem_interface.F90)

if(CATCHEM_CCPP_TEST_STUBS)
  add_library(
    catchem_ccpp_test_stubs
    ${PROJECT_SOURCE_DIR}/tests/ccpp_test_stubs.F90
  )
  set_target_properties(
    catchem_ccpp_test_stubs
    PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/ccpp-test-modules
  )
elseif(NOT CATCHEM_CCPP_HOST_MODULE_DIR)
  message(
    FATAL_ERROR
    "CATCHEM_BUILD_CCPP requires CATCHEM_CCPP_HOST_MODULE_DIR containing the host CCPP modules, "
    "or CATCHEM_CCPP_TEST_STUBS=ON for the contract-only test harness"
  )
endif()

set(_lib CCPP_CATChem_interface)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem)
if(CATCHEM_CCPP_TEST_STUBS)
  target_link_libraries(${_lib} PUBLIC catchem_ccpp_test_stubs)
  target_include_directories(
    ${_lib}
    PRIVATE ${CMAKE_BINARY_DIR}/ccpp-test-modules
  )
else()
  target_include_directories(${_lib} PRIVATE ${CATCHEM_CCPP_HOST_MODULE_DIR})
endif()

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_compile_options(${_lib} PRIVATE ${_catchem_compiler_options})
