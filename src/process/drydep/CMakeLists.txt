# DryDep Process CMakeLists.txt
# Generated on: 2025-11-14T22:58:26.596456
# Author: Wei Li
# Description: Process for computing dry deposition of gas and aerosol species

# Define drydep process sources
set(
  DRYDEP_PROCESS_SOURCES
  DryDepCommon_Mod.F90
  ProcessDryDepInterface_Mod.F90
  DryDepProcessCreator_Mod.F90
)

# Define drydep scheme sources
set(
  DRYDEP_SCHEME_SOURCES
  schemes/DryDepScheme_WESELY_Mod.F90
  schemes/DryDepScheme_GOCART_Mod.F90
  schemes/DryDepScheme_ZHANG_Mod.F90
)

# Combine all sources
set(DRYDEP_ALL_SOURCES ${DRYDEP_PROCESS_SOURCES} ${DRYDEP_SCHEME_SOURCES})

# Create the drydep process library
set(_lib CATChem_process_drydep)
add_library(${_lib} ${DRYDEP_ALL_SOURCES})

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
  PRIVATE PROCESS_NAME="drydep" PROCESS_VERSION="1.0.0"
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
  # Add tests subdirectory - tests are located in tests/process/drydep/
  set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests/process/drydep")
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
    set(DOXYGEN_PROJECT_NAME "DryDep Process")
    set(
      DOXYGEN_PROJECT_BRIEF
      "Process for computing dry deposition of gas and aerosol species"
    )
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    set(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_FILE_PATTERNS "*.F90 *.md")
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build/* */tests/*")

    doxygen_add_docs(
      drydep_docs
      ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating drydep process documentation"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/
      DESTINATION share/doc/drydep
      OPTIONAL
    )
  endif()
endif()

# Performance testing

# Process-specific targets

# Print configuration summary
message(STATUS "Configured DryDep process:")
message(STATUS "  Process type: deposition")
message(STATUS "  Number of schemes: 3")
message(STATUS "  Default scheme: ")
message(STATUS "  Has size bins: No")
message(STATUS "  Multiphase: No")
message(STATUS "  Generate tests: Yes")
message(STATUS "  Generate docs: Yes")
message(STATUS "  Generate examples: Yes")
