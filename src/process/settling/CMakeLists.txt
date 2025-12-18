# Settling Process CMakeLists.txt
# Generated on: 2025-12-18T14:12:33.107935
# Author: Wei Li
# Description: Process for computing gravitational settling of aerosol species

# Define settling process sources
set(
  SETTLING_PROCESS_SOURCES
  SettlingCommon_Mod.F90
  ProcessSettlingInterface_Mod.F90
  SettlingProcessCreator_Mod.F90
)

# Define settling scheme sources
set(SETTLING_SCHEME_SOURCES schemes/SettlingScheme_GOCART_Mod.F90)

# Combine all sources
set(SETTLING_ALL_SOURCES ${SETTLING_PROCESS_SOURCES} ${SETTLING_SCHEME_SOURCES})

# Create the settling process library
set(_lib CATChem_process_settling)
add_library(${_lib} ${SETTLING_ALL_SOURCES})

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
  PRIVATE PROCESS_NAME="settling" PROCESS_VERSION="1.0.0"
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
  # Add tests subdirectory - tests are located in tests/process/settling/
  set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests/process/settling")
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
    set(DOXYGEN_PROJECT_NAME "Settling Process")
    set(
      DOXYGEN_PROJECT_BRIEF
      "Process for computing gravitational settling of aerosol species"
    )
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    set(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_FILE_PATTERNS "*.F90 *.md")
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build/* */tests/*")

    doxygen_add_docs(
      settling_docs
      ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating settling process documentation"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/
      DESTINATION share/doc/settling
      OPTIONAL
    )
  endif()
endif()

# Performance testing

# Process-specific targets

# Print configuration summary
message(STATUS "Configured Settling process:")
message(STATUS "  Process type: deposition")
message(STATUS "  Number of schemes: 1")
message(STATUS "  Default scheme: ")
message(STATUS "  Has size bins: No")
message(STATUS "  Multiphase: No")
message(STATUS "  Generate tests: Yes")
message(STATUS "  Generate docs: Yes")
message(STATUS "  Generate examples: Yes")
