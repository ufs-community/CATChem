# GasChem Process CMakeLists.txt
# Generated on: 2026-06-09T15:53:02.114257
# Author: Maggie Bruckner
# Description: Process for MICM gas phase chemical solver

# Define GasChem process sources
set(
  GASCHEM_PROCESS_SOURCES
  GasChemCommon_Mod.F90
  ProcessGasChemInterface_Mod.F90
  GasChemProcessCreator_Mod.F90
)

# Define GasChem scheme sources
set(GASCHEM_SCHEME_SOURCES schemes/GasChemScheme_MUSICA_Mod.F90)

# Combine all sources
set(GASCHEM_ALL_SOURCES ${GASCHEM_PROCESS_SOURCES} ${GASCHEM_SCHEME_SOURCES})

# Create the GasChem process library
set(_lib CATChem_process_GasChem)
add_library(${_lib} ${GASCHEM_ALL_SOURCES})

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
  PRIVATE PROCESS_NAME="GasChem" PROCESS_VERSION="1.0.0"
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
# Note: tests are registered via tests/CMakeLists.txt (tests/process/GasChem/)

# Documentation
if(BUILD_DOCUMENTATION)
  find_package(Doxygen QUIET)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_NAME "GasChem Process")
    set(DOXYGEN_PROJECT_BRIEF "Process for MICM gas phase chemical solver")
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    set(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_FILE_PATTERNS "*.F90 *.md")
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build/* */tests/*")

    doxygen_add_docs(
      GasChem_docs
      ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating GasChem process documentation"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/
      DESTINATION share/doc/GasChem
      OPTIONAL
    )
  endif()
endif()

# Performance testing
if(BUILD_BENCHMARKS)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    add_executable(GasChem_benchmark benchmarks/benchmark_GasChem.F90)
    target_link_libraries(
      GasChem_benchmark
      PRIVATE ${_lib} benchmark::benchmark
    )
  endif()
endif()

# Process-specific targets

# Print configuration summary
message(STATUS "Configured GasChem process:")
message(STATUS "  Process type: chemistry")
message(STATUS "  Number of schemes: 1")
message(STATUS "  Default scheme: ")
message(STATUS "  Has size bins: No")
message(STATUS "  Multiphase: No")
message(STATUS "  Generate tests: Yes")
message(STATUS "  Generate docs: Yes")
message(STATUS "  Generate examples: No")
