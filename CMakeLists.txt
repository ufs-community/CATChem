cmake_minimum_required(VERSION 3.21)

# Extract version from pyproject.toml
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" PYPROJECT_TOML)
string(REGEX MATCH "version = \"([^\"]+)\"" _ ${PYPROJECT_TOML})
set(CATCHEM_VERSION_FULL ${CMAKE_MATCH_1})

# CMake project() VERSION requires numeric format: major[.minor[.patch[.tweak]]]
if(CATCHEM_VERSION_FULL MATCHES "^([0-9]+(\\.[0-9]+)*)")
  set(CATCHEM_VERSION "${CMAKE_MATCH_1}")
else()
  message(
    WARNING
    "${CATCHEM_VERSION_FULL} includes strings or other characters (likely a release candidate), defaulting to 0.0.0"
  )
  set(CATCHEM_VERSION "0.0.0")
endif()

project(CATChem VERSION ${CATCHEM_VERSION} LANGUAGES Fortran C CXX)

# Global variables for path management
set(CATCHEM_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}")
set(CATCHEM_SRC_DIR "${CATCHEM_ROOT_DIR}/src")

# Option to build NUOPC cap (defaults to OFF)
option(CATCHEM_BUILD_NUOPC "Build NUOPC cap and drivers" OFF)
option(CATCHEM_BUILD_TESTING "Build CATChem testing framework" ON)
option(CATCHEM_TRACE_NUOPC "Enable ESMF trace regions in NUOPC cap" OFF)

if(CATCHEM_TRACE_NUOPC)
  add_compile_definitions(CATCHEM_TRACE_NUOPC)
endif()

# Option to build with MUSICA chemistry library (defaults to OFF)
option(CATCHEM_BUILD_MUSICA "Build with MUSICA chemistry library" OFF)

# Set C++ standard — MUSICA requires C++20, yaml-cpp needs at least C++11
if(CATCHEM_BUILD_MUSICA)
  set(CMAKE_CXX_STANDARD 20)
else()
  set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# Handle Intel compiler deprecation warnings
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  # Always disable the deprecation warning for ifort
  list(APPEND _catchem_compiler_options -diag-disable=6379,10448)
  # Also disable problematic flags that might come from external dependencies
  set(
    CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -diag-disable=6379,10448,10440"
  )
endif()

# No Fortran line length restrictions on all GNU build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
  list(APPEND _catchem_compiler_options -ffree-line-length-none)
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
    -fno-sanitize=all
    -ffpe-trap=invalid,zero,overflow
    -fbacktrace
    -fmax-errors=0
  )
elseif(
  CMAKE_Fortran_COMPILER_ID MATCHES "^Intel"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  list(
    APPEND
    _catchem_compiler_options
    -O0
    -warn
    all
    -check
    all,nouninit
    -implicitnone
    -fpe0
    -traceback
    -error-limit
    0
  )
endif()

# Disable memory sanitizer to avoid linking issues
if(
  CMAKE_C_COMPILER_ID STREQUAL "GNU"
  OR CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM"
)
  # gcc and icx both accept standard sanitizer flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-sanitize=all")
  # Note: Intel Classic C (icc) doesn't use -fsanitize, so we safely do nothing
endif()
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
  OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM"
)
  # g++ and icpx both accept standard sanitizer flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-sanitize=all")
  # Note: Intel Classic C++ (icpc) doesn't use -fsanitize, so we safely do nothing
endif()

# Disable Intel oneAPI C/C++ warning about -g without -O0
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rno-debug-disables-optimization")
endif()
if(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Rno-debug-disables-optimization")
endif()

# Find external dependencies if not already available (for standalone builds)
# Check if targets exist before trying to find packages
if(NOT TARGET NetCDF::NetCDF_Fortran)
  find_package(NetCDF REQUIRED COMPONENTS Fortran)
endif()

if(CATCHEM_BUILD_NUOPC)
  # Resolve ESMF target: prefer uppercase ESMF (FindESMF), fall back to
  # lowercase esmf (MPAS internal), or find_package as last resort.
  if(TARGET ESMF)
    set(CATCHEM_ESMF_TARGET ESMF)
  elseif(TARGET esmf)
    set(CATCHEM_ESMF_TARGET esmf)
  else()
    find_package(ESMF REQUIRED)
    set(CATCHEM_ESMF_TARGET esmf)
  endif()
  message(STATUS "CATChem: using ESMF target '${CATCHEM_ESMF_TARGET}'")
endif()

# Auto-detect if default real is promoted to 8-byte by the parent build
# or by per-target flags, and define USE_REAL8 so precision_mod sets fp=f8.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
  if(CMAKE_Fortran_FLAGS MATCHES "-fdefault-real-8")
    add_compile_definitions(USE_REAL8)
    message(STATUS "CATChem: detected -fdefault-real-8, enabling USE_REAL8")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
  if(CMAKE_Fortran_FLAGS MATCHES "-r8|-real-size 64|-fdefault-real-8")
    add_compile_definitions(USE_REAL8)
    message(
      STATUS
      "CATChem: detected real-8 promotion flag, enabling USE_REAL8"
    )
  endif()
endif()

# Find yaml-cpp for configuration parsing
if(NOT TARGET yaml-cpp::yaml-cpp AND NOT TARGET yaml-cpp)
  find_package(yaml-cpp QUIET)
endif()

add_subdirectory(src)

# Enable testing for standalone builds
if(CATCHEM_BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if(CATCHEM_BUILD_NUOPC)
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
    drivers/nuopc/catchem_regrid_mod.F90
    drivers/nuopc/catchem_latlon_output_mod.F90
    drivers/nuopc/aqmio.F90
    drivers/nuopc/catchem_emis_mod.F90
    drivers/nuopc/catchem_nuopc_cap.F90
    drivers/nuopc/catchem_nuopc_interface.F90
  )

  # Link only to external dependencies (not internal CATChem targets)
  target_link_libraries(
    cc_nuopc
    PUBLIC ${CATCHEM_ESMF_TARGET} NetCDF::NetCDF_Fortran
  )

  # Manually handle the CATChem dependencies by linking to built libraries
  # Only need to depend on CATChem - it already depends on CATChem_core, yaml_interface, etc.
  add_dependencies(cc_nuopc CATChem)

  # Link to library files in dependency order (most dependent first)
  target_link_libraries(
    cc_nuopc
    PRIVATE
      $<TARGET_FILE:CATChem>
      $<TARGET_FILE:CATChem_process_seasalt>
      $<TARGET_FILE:CATChem_process_dust>
      $<TARGET_FILE:CATChem_process_drydep>
      $<TARGET_FILE:CATChem_process_wetdep>
      $<TARGET_FILE:CATChem_process_settling>
      $<TARGET_FILE:CATChem_process_so4chem>
      $<TARGET_FILE:CATChem_process_carbchem>
      $<TARGET_FILE:CATChem_core>
      $<TARGET_FILE:yaml_interface>
      $<TARGET_FILE:catchem_utilities>
  )

  # Link GOCART2G library if available (needed for drydep process)
  if(TARGET GOCART2G)
    target_link_libraries(cc_nuopc PRIVATE $<TARGET_FILE:GOCART2G>)
  endif()

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
