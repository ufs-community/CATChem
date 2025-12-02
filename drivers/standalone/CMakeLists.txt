# CMakeLists.txt for CATChem Standalone Column Driver
cmake_minimum_required(VERSION 3.12)

project(
  CATChem_Standalone_Driver
  VERSION 2.0.0
  DESCRIPTION "Standalone column driver for CATChem atmospheric chemistry"
  LANGUAGES Fortran
)

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)

# Check for NetCDF (optional for output)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(NETCDF QUIET netcdf-fortran)
endif()

# Set compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(
    CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Wall -Wextra -fdefault-real-8 -fdefault-double-8"
  )
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -r8")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -unroll")
endif()

# Include directories
include_directories(
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/api
  ${CMAKE_SOURCE_DIR}/src/process
)

# Define source files for the standalone driver
set(DRIVER_SOURCES standalone_column_driver.F90)

# Define CATChem library sources (would link to actual CATChem library in real build)
set(
  CATCHEM_CORE_SOURCES
  ${CMAKE_SOURCE_DIR}/src/core/precision_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/error_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/constants.F90
  ${CMAKE_SOURCE_DIR}/src/core/config_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/ConfigManager_Mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/state_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/metstate_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/chemstate_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/emisstate_mod.F90
  ${CMAKE_SOURCE_DIR}/src/core/diagstate_mod.F90
  ${CMAKE_SOURCE_DIR}/src/api/CATChemAPI_Mod.F90
  ${CMAKE_SOURCE_DIR}/src/process/ProcessManager_Mod.F90
  ${CMAKE_SOURCE_DIR}/src/process/ProcessInterface_Mod.F90
  ${CMAKE_SOURCE_DIR}/src/process/ProcessFactory_Mod.F90
)

# Create CATChem library
add_library(catchem_core ${CATCHEM_CORE_SOURCES})

# Set module output directory
set_target_properties(
  catchem_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
)

# Create standalone driver executable
add_executable(standalone_column_driver ${DRIVER_SOURCES})

# Link libraries
target_link_libraries(standalone_column_driver catchem_core ${NETCDF_LIBRARIES})

# Include module directory
target_include_directories(
  standalone_column_driver
  PRIVATE ${CMAKE_BINARY_DIR}/modules ${NETCDF_INCLUDE_DIRS}
)

# Set output directory
set_target_properties(
  standalone_column_driver
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Installation
install(
  TARGETS standalone_column_driver
  RUNTIME DESTINATION bin COMPONENT drivers
)

# Create configuration files for testing
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/column_driver_config.yml.in
  ${CMAKE_BINARY_DIR}/column_driver_config.yml
  @ONLY
)

# Add custom target for running the driver
add_custom_target(
  run_column_driver
  COMMAND ${CMAKE_BINARY_DIR}/bin/standalone_column_driver
  DEPENDS standalone_column_driver
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Running standalone column driver"
)

# Documentation target
find_program(DOXYGEN_EXECUTABLE doxygen)
if(DOXYGEN_EXECUTABLE)
  add_custom_target(
    driver_docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating driver documentation"
  )
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "CATChem Standalone Driver Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  Fortran flags: ${CMAKE_Fortran_FLAGS}")
if(NETCDF_FOUND)
  message(STATUS "  NetCDF support: YES")
else()
  message(STATUS "  NetCDF support: NO")
endif()
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
