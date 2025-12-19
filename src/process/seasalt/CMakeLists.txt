# SeaSalt Process CMakeLists.txt
# Generated on: 2025-11-14T23:01:21.792441
# Author: Barry Baker & Wei Li
# Description: Process for computing sea salt aerosol emissions over ocean surfaces

# Define seasalt process sources
set(
  SEASALT_PROCESS_SOURCES
  SeaSaltCommon_Mod.F90
  ProcessSeaSaltInterface_Mod.F90
  SeaSaltProcessCreator_Mod.F90
)

# Define seasalt scheme sources
set(
  SEASALT_SCHEME_SOURCES
  schemes/SeaSaltScheme_GONG97_Mod.F90
  schemes/SeaSaltScheme_GONG03_Mod.F90
  schemes/SeaSaltScheme_GEOS12_Mod.F90
)

# Combine all sources
set(SEASALT_ALL_SOURCES ${SEASALT_PROCESS_SOURCES} ${SEASALT_SCHEME_SOURCES})

# Create the seasalt process library
set(_lib CATChem_process_seasalt)
add_library(${_lib} ${SEASALT_ALL_SOURCES})

# Link with required libraries
target_link_libraries(${_lib} PUBLIC CATChem_core)

# Additional dependencies

# Set module directory
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# Include directories
target_include_directories(${_lib} PRIVATE ${CMAKE_BINARY_DIR}/include)

# Compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(${_lib} PRIVATE -fdefault-real-8 -fdefault-double-8)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(${_lib} PRIVATE -r8)
endif()

# Add preprocessor definitions
target_compile_definitions(
  ${_lib}
  PRIVATE PROCESS_NAME="seasalt" PROCESS_VERSION="1.0.0"
)

# Install targets
install(
  TARGETS ${_lib}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

# Install module files
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.mod"
)

# Add subdirectories
if(BUILD_TESTING)
  # Add tests subdirectory - tests are located in tests/process/seasalt/
  set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests/process/seasalt")
  if(EXISTS "${TESTS_DIR}/CMakeLists.txt")
    add_subdirectory("${TESTS_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/tests")
  endif()
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Documentation
if(BUILD_DOCUMENTATION)
  find_package(Doxygen QUIET)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_NAME "SeaSalt Process")
    set(
      DOXYGEN_PROJECT_BRIEF
      "Process for computing sea salt aerosol emissions over ocean surfaces"
    )
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    set(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_FILE_PATTERNS "*.F90 *.md")
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build/* */tests/*")

    doxygen_add_docs(
      seasalt_docs
      ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating seasalt process documentation"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/
      DESTINATION share/doc/seasalt
      OPTIONAL
    )
  endif()
endif()

# Performance testing
if(BUILD_BENCHMARKS)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    add_executable(seasalt_benchmark benchmarks/benchmark_seasalt.F90)
    target_link_libraries(
      seasalt_benchmark
      PRIVATE ${_lib} benchmark::benchmark
    )
  endif()
endif()

# Process-specific targets
# Emission inventory preprocessing (if needed)
# add_custom_target(seasalt_preprocess_inventory
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/preprocess_inventory.py
#     COMMENT "Preprocessing emission inventory for seasalt"
# )

# Print configuration summary
message(STATUS "Configured SeaSalt process:")
message(STATUS "  Process type: emission")
message(STATUS "  Number of schemes: 3")
message(STATUS "  Default scheme: ")
message(STATUS "  Has size bins: No")
message(STATUS "  Multiphase: No")
message(STATUS "  Generate tests: Yes")
message(STATUS "  Generate docs: Yes")
message(STATUS "  Generate examples: Yes")
