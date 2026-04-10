# CarbChem Process CMakeLists.txt
# Generated on: 2026-04-10T16:46:37.829974
# Author: Wei Li
# Description: Process for computing chemical production and loss of carbon species

# Define carbchem process sources
set(
  CARBCHEM_PROCESS_SOURCES
  CarbChemCommon_Mod.F90
  ProcessCarbChemInterface_Mod.F90
  CarbChemProcessCreator_Mod.F90
)

# Define carbchem scheme sources
set(CARBCHEM_SCHEME_SOURCES schemes/CarbChemScheme_GOCART_Mod.F90)

# Combine all sources
set(CARBCHEM_ALL_SOURCES ${CARBCHEM_PROCESS_SOURCES} ${CARBCHEM_SCHEME_SOURCES})

# Create the carbchem process library
set(_lib CATChem_process_carbchem)
add_library(${_lib} ${CARBCHEM_ALL_SOURCES})

# Link with required libraries
target_link_libraries(${_lib} PUBLIC CATChem_core)

# Additional dependencies
# Chemistry processes may need additional math libraries
target_link_libraries(${_lib} PUBLIC ${LAPACK_LIBRARIES})

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
  # Chemistry processes may benefit from aggressive optimization
  target_compile_options(${_lib} PRIVATE -O3 -funroll-loops)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(${_lib} PRIVATE -r8)
  target_compile_options(${_lib} PRIVATE -O3 -unroll)
endif()

# Add preprocessor definitions
target_compile_definitions(
  ${_lib}
  PRIVATE PROCESS_NAME="carbchem" PROCESS_VERSION="1.0.0"
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
# Note: tests are registered via tests/CMakeLists.txt (tests/process/carbchem/)

# Documentation
if(BUILD_DOCUMENTATION)
  find_package(Doxygen QUIET)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_NAME "CarbChem Process")
    set(
      DOXYGEN_PROJECT_BRIEF
      "Process for computing chemical production and loss of carbon species"
    )
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    set(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_FILE_PATTERNS "*.F90 *.md")
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build/* */tests/*")

    doxygen_add_docs(
      carbchem_docs
      ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating carbchem process documentation"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/
      DESTINATION share/doc/carbchem
      OPTIONAL
    )
  endif()
endif()

# Performance testing
if(BUILD_BENCHMARKS)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    add_executable(carbchem_benchmark benchmarks/benchmark_carbchem.F90)
    target_link_libraries(
      carbchem_benchmark
      PRIVATE ${_lib} benchmark::benchmark
    )
  endif()
endif()

# Process-specific targets

# Print configuration summary
message(STATUS "Configured CarbChem process:")
message(STATUS "  Process type: chemistry")
message(STATUS "  Number of schemes: 1")
message(STATUS "  Default scheme: ")
message(STATUS "  Number of species: 4")
message(STATUS "  Has size bins: No")
message(STATUS "  Multiphase: No")
message(STATUS "  Generate tests: Yes")
message(STATUS "  Generate docs: Yes")
message(STATUS "  Generate examples: No")
