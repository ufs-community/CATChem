# CMakeLists.txt for seasalt process tests
# Generated on: 2025-11-14T23:01:21.936351

cmake_minimum_required(VERSION 3.12)

# Test executables
set(
  TEST_SOURCES
  unit/test_seasalt_unit.F90
  integration/test_seasalt_integration.F90
)

# Create test executables
foreach(test_source ${TEST_SOURCES})
  get_filename_component(test_name ${test_source} NAME_WE)

  add_executable(${test_name} ${test_source})

  # Link against process modules
  target_link_libraries(
    ${test_name}
    CATChem_process_seasalt
    CATChem_core
    testing
  )

  # Set Fortran module directory
  set_target_properties(
    ${test_name}
    PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )

  # Add include directories for module files
  target_include_directories(
    ${test_name}
    PRIVATE "${CMAKE_BINARY_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )

  # Add to test suite
  add_test(NAME ${test_name} COMMAND ${test_name})

  # Set test properties
  set_tests_properties(
    ${test_name}
    PROPERTIES TIMEOUT 60 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endforeach()

# Test data files
configure_file(
  "${CMAKE_SOURCE_DIR}/tests/Configs/Default/CATChem_new_config.yml"
  "${CMAKE_CURRENT_BINARY_DIR}/seasalt_config.yaml"
  COPYONLY
)

# Additional test data files - copy other config files that might be needed
file(GLOB TEST_DATA_FILES "${CMAKE_SOURCE_DIR}/tests/Configs/Default/*.yml")
foreach(test_data_file ${TEST_DATA_FILES})
  get_filename_component(filename ${test_data_file} NAME)
  configure_file(
    "${test_data_file}"
    "${CMAKE_CURRENT_BINARY_DIR}/${filename}"
    COPYONLY
  )
endforeach()

# Custom test targets with unique names
add_custom_target(
  run_seasalt_unit
  COMMAND test_seasalt_unit
  DEPENDS test_seasalt_unit
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Running unit tests for seasalt process"
)

add_custom_target(
  run_seasalt_integration
  COMMAND test_seasalt_integration
  DEPENDS test_seasalt_integration
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Running integration tests for seasalt process"
)

add_custom_target(
  run_seasalt_all
  DEPENDS run_seasalt_unit run_seasalt_integration
  COMMENT "Running all tests for seasalt process"
)

# Test coverage (if available)
if(CODE_COVERAGE)
  find_program(GCOV_PATH gcov)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)

  if(GCOV_PATH AND LCOV_PATH AND GENHTML_PATH)
    # Add coverage flags
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

    # Coverage target
    add_custom_target(
      coverage_seasalt
      COMMAND
        ${LCOV_PATH} --directory . --capture --output-file coverage_seasalt.info
      COMMAND
        ${LCOV_PATH} --remove coverage_seasalt.info '/usr/*' --output-file
        coverage_seasalt.info
      COMMAND ${GENHTML_PATH} -o coverage_seasalt_html coverage_seasalt.info
      DEPENDS run_seasalt_all
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      COMMENT "Generating test coverage report for seasalt process"
    )
  endif()
endif()

# Performance benchmarks
add_custom_target(
  benchmark_seasalt
  DEPENDS run_seasalt_integration
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Running performance benchmarks for seasalt process"
)

# Memory checking with valgrind (if available)
find_program(VALGRIND_PATH valgrind)
if(VALGRIND_PATH)
  add_custom_target(
    memcheck_seasalt
    COMMAND
      ${VALGRIND_PATH} --tool=memcheck --leak-check=full --show-leak-kinds=all
      --track-origins=yes ./test_seasalt_integration
    DEPENDS test_seasalt_integration
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Running memory check for seasalt process"
  )
endif()

# Install test executables (optional)
if(INSTALL_TESTS)
  install(TARGETS ${TEST_SOURCES} RUNTIME DESTINATION bin/tests)

  install(
    FILES "${CMAKE_SOURCE_DIR}/tests/Configs/Default/CATChem_new_config.yml"
    DESTINATION share/tests/seasalt
    RENAME seasalt_config.yaml
  )
endif()
