cmake_minimum_required(VERSION 3.10)

project(CATChem VERSION 0.1.0 LANGUAGES Fortran C CXX)

# Set C++ standard for yaml-cpp compatibility
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to build NUOPC cap (defaults to OFF)
option(BUILD_NUOPC "Build NUOPC cap and drivers" OFF)

include(CMakePrintHelpers)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMakeModules/Modules")

# If GNU compiler and Debug build type, add extra debug flags
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
cmake_print_variables(
  CMAKE_Fortran_COMPILER
  CMAKE_Fortran_COMPILER_ID
  CMAKE_BUILD_TYPE
)
set(_catchem_compiler_options "")

# Handle Intel compiler deprecation warnings and incompatible flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # Always disable the deprecation warning for ifort
  list(APPEND _catchem_compiler_options -diag-disable=6379,10448)
  # Also disable problematic flags that might come from external dependencies
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable=6379,10448")
  # Disable memory sanitizer to avoid linking issues
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-sanitize=all")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-sanitize=all")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-sanitize=all")
endif()

if(
  CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  list(
    APPEND
    _catchem_compiler_options
    -Og
    -Wall
    -Wextra
    -pedantic
    -fcheck=all
    -ffpe-trap=invalid,zero,overflow
    -fbacktrace
    -fmax-errors=0
    -ffree-line-length-none
  )
elseif(
  CMAKE_Fortran_COMPILER_ID MATCHES "Intel"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  list(
    APPEND
    _catchem_compiler_options
    -O0
    -warn
    all
    -check
    all
    -implicitnone
    -fpe0
    -traceback
    -error-limit
    0
  )
elseif(
  CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  list(
    APPEND
    _catchem_compiler_options
    -O0
    -warn
    all
    -check
    all
    -implicitnone
    -fpe0
    -traceback
    -error-limit
    0
  )
endif()

# Find external dependencies if not already available (for standalone builds)
# Check if targets exist before trying to find packages
if(NOT TARGET NetCDF::NetCDF_Fortran)
  find_package(NetCDF REQUIRED COMPONENTS Fortran)
endif()

if(BUILD_NUOPC)
  if(NOT TARGET esmf)
    find_package(ESMF REQUIRED)
  endif()
endif()

# Enable testing for standalone builds
enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

if(BUILD_NUOPC)
  add_subdirectory(drivers/nuopc)

  # Copy module files from include/ to mod/ for UFS integration during build
  # add_custom_target(copy_modules_to_mod ALL
  #   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/../mod
  #   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/../mod
  #   DEPENDS NUOPC_CATChem_cap CATChem_core
  #   COMMENT "Copying module files from include/ to mod/ for UFS integration"
  # )

  # Create a self-contained catchem library by combining all necessary object files
  # This avoids export dependency issues by not linking to other CMake targets

  # Create catchem library directly from source files to avoid dependency exports
  add_library(
    cc_nuopc
    STATIC
    drivers/nuopc/aqmio.F90
    drivers/nuopc/catchem_nuopc_cap.F90
    drivers/nuopc/catchem_nuopc_interface.F90
  )

  # Link only to external dependencies (not internal CATChem targets)
  target_link_libraries(cc_nuopc PUBLIC esmf NetCDF::NetCDF_Fortran)

  # Manually handle the CATChem dependencies by linking to built libraries
  # Only need to depend on CATChem - it already depends on CATChem_core, yaml_interface, etc.
  add_dependencies(cc_nuopc CATChem)

  # Link to library files in dependency order (most dependent first)
  target_link_libraries(
    cc_nuopc
    PRIVATE
      $<TARGET_FILE:CATChem>
      $<TARGET_FILE:CATChem_process_seasalt>
      $<TARGET_FILE:CATChem_core>
      $<TARGET_FILE:yaml_interface>
      $<TARGET_FILE:catchem_utilities>
  )

  # Include directories for compilation
  target_include_directories(
    cc_nuopc
    PRIVATE ${CMAKE_BINARY_DIR}/include
    #$<TARGET_PROPERTY:CATChem,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:esmf,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:NetCDF::NetCDF_Fortran,INTERFACE_INCLUDE_DIRECTORIES>
  )

  set_target_properties(
    cc_nuopc
    PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../mod
      EXPORT_NAME cc_nuopc
  )

  target_compile_options(cc_nuopc PRIVATE ${_catchem_compiler_options})

  ###############################################################################
  ### Installation and Export for UFS integration
  ###############################################################################

  # Install only the self-contained catchem library - all dependencies are linked directly
  install(
    TARGETS cc_nuopc
    EXPORT ufs-config
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endif()

# Install module files to both locations for compatibility
# Keep include/ for standalone CATChem, add mod/ for UFS integration
# install(DIRECTORY ${CMAKE_BINARY_DIR}/include/
#         DESTINATION include/
#         COMPONENT Development
#         FILES_MATCHING PATTERN "*.mod")

# install(DIRECTORY ${CMAKE_BINARY_DIR}/include/
#         DESTINATION mod/
#         COMPONENT Development
#         FILES_MATCHING PATTERN "*.mod")

# Note: cc_nuopc target is exported via ufs-config export set in UFS CMakeLists.txt
