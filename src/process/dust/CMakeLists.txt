# Dust Process CMakeLists.txt
# Generated on: 2026-04-17T13:57:10.272387
# Author: Wei Li & Barry Baker
# Description: Process for computing windblown dust emissions

# Define dust process sources
set(
  DUST_PROCESS_SOURCES
  DustCommon_Mod.F90
  catchem_process_dust.cpp
  DustScienceBridge.F90
)

# Define dust scheme sources
set(
  DUST_SCHEME_SOURCES
  schemes/DustScheme_FENGSHA_Mod.F90
  schemes/DustScheme_GINOUX_Mod.F90
)

# Combine all sources
set(DUST_ALL_SOURCES ${DUST_PROCESS_SOURCES} ${DUST_SCHEME_SOURCES})

# Create the dust process library
set(_lib CATChem_process_dust)
add_library(${_lib} ${DUST_ALL_SOURCES})

# Link with required libraries (CATChem_core_cpp carries the Kokkos or
# standalone-mdspan usage requirements for the C++ sources)
target_link_libraries(${_lib} PUBLIC CATChem_core CATChem_core_cpp)

# Additional dependencies

# Set module directory
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# Include directories
target_include_directories(
  ${_lib}
  PRIVATE ${CMAKE_BINARY_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../../core
)

# Compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(
    ${_lib}
    PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8;-fdefault-double-8>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(${_lib} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-r8>)
endif()

# Add preprocessor definitions
target_compile_definitions(
  ${_lib}
  PRIVATE PROCESS_NAME="dust" PROCESS_VERSION="1.0.0"
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
# Note: tests are registered via tests/CMakeLists.txt (tests/process/dust/)

# Documentation
if(BUILD_DOCUMENTATION)
  find_package(Doxygen QUIET)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_NAME "Dust Process")
    set(DOXYGEN_PROJECT_BRIEF "Process for computing windblown dust emissions")
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    set(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_RECURSIVE YES)
    set(DOXYGEN_FILE_PATTERNS "*.F90 *.md")
    set(DOXYGEN_EXCLUDE_PATTERNS "*/build/* */tests/*")

    doxygen_add_docs(
      dust_docs
      ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating dust process documentation"
    )

    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/
      DESTINATION share/doc/dust
      OPTIONAL
    )
  endif()
endif()

# Performance testing
if(BUILD_BENCHMARKS)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    add_executable(dust_benchmark benchmarks/benchmark_dust.F90)
    target_link_libraries(dust_benchmark PRIVATE ${_lib} benchmark::benchmark)
  endif()
endif()

# Process-specific targets
# Emission inventory preprocessing (if needed)
# add_custom_target(dust_preprocess_inventory
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/preprocess_inventory.py
#     COMMENT "Preprocessing emission inventory for dust"
# )

# Print configuration summary
message(STATUS "Configured Dust process:")
message(STATUS "  Process type: emission")
message(STATUS "  Number of schemes: 2")
message(STATUS "  Default scheme: ")
message(STATUS "  Has size bins: No")
message(STATUS "  Multiphase: No")
message(STATUS "  Generate tests: Yes")
message(STATUS "  Generate docs: Yes")
message(STATUS "  Generate examples: No")
