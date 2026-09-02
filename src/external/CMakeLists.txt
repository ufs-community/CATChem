# External dependencies

# Create a basic utilities library that other modules can depend on
add_library(catchem_utilities ../core/Precision_Mod.F90)

set_target_properties(
  catchem_utilities
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

#
# GOCART-2G Process Library
#

set(
  _srcs
  GOCART/Process_Library/GOCART2G_Process.F90
  GOCART/Process_Library/GOCART2G_MieMod.F90
)
set(_lib GOCART2G)

add_library(${_lib} ${_srcs})
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
  target_compile_options(${_lib} PRIVATE -ffree-line-length-none)
  # Detect if parent build promotes default real to 8-byte and propagate to GOCART2G
  if(CMAKE_Fortran_FLAGS MATCHES "-fdefault-real-8")
    target_compile_options(${_lib} PRIVATE -fdefault-real-8 -fdefault-double-8)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # Detect if parent build promotes default real to 8-byte and propagate to GOCART2G
  # Intel supports: -r8, -real-size 64, and ifx also accepts -fdefault-real-8
  if(CMAKE_Fortran_FLAGS MATCHES "-r8|-real-size 64|-fdefault-real-8")
    target_compile_options(${_lib} PRIVATE -r8)
  endif()
endif()

# The fpp preprocessor shipped with IntelLLVM 2025.2.1 miscounts macro
# arguments in GOCART2G_MieMod.F90 ("#error: number of arguments doesn't
# match"); 2025.3.x is fixed. Where an unaffected fpp from an earlier oneAPI
# install exists (e.g. Gaea C6), preprocess this target with it instead —
# same workaround the ufs-weather-model carries for its GOCART targets:
# https://github.com/ufs-community/ufs-weather-model/issues/3213#issuecomment-4711120326
# Remove once affected platforms move to a fixed compiler.
set(_gocart_fpp "/opt/intel/oneapi/compiler/2025.0/bin/fpp")
if(
  CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM"
  AND CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 2025.2.1
  AND EXISTS "${_gocart_fpp}"
)
  message(STATUS "GOCART2G: using ${_gocart_fpp} (IntelLLVM 2025.2.1 fpp bug)")
  target_compile_options(
    ${_lib}
    PRIVATE -fpp-name=${_gocart_fpp} -Qoption,fpp,-P
  )
endif()

#find_package(NetCDF REQUIRED COMPONENTS Fortran)

# Use NetCDF from parent scope (already found by main UFS CMakeLists.txt or standalone CMakeLists.txt)
if(TARGET NetCDF::NetCDF_Fortran)
  target_link_libraries(${_lib} NetCDF::NetCDF_Fortran)
endif()

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

#
# yaml-cpp
#

if(NOT yaml-cpp_FOUND AND NOT TARGET yaml-cpp::yaml-cpp AND NOT TARGET yaml-cpp)
  #  message(
  #    FATAL_ERROR
  #    "yaml-cpp NOT FOUND BY find_package (Temporary error for testing)"
  #  )
  message(STATUS "yaml-cpp not found on system, using internal version")

  include(FetchContent)

  # Configure yaml-cpp build options
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Build yaml-cpp tests")
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Build yaml-cpp tools")
  set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Build yaml-cpp contrib")

  # Declare yaml-cpp as already available via FetchContent
  # This prevents MUSICA from fetching their own copy
  FetchContent_Declare(yaml-cpp SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp)
  FetchContent_MakeAvailable(yaml-cpp)
endif()

# Ensure both yaml-cpp and yaml-cpp::yaml-cpp targets exist regardless of how it was found
if(TARGET yaml-cpp AND NOT TARGET yaml-cpp::yaml-cpp)
  add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
endif()
if(TARGET yaml-cpp::yaml-cpp AND NOT TARGET yaml-cpp)
  add_library(yaml-cpp ALIAS yaml-cpp::yaml-cpp)
endif()

#
# MUSICA Fortran API (optional, controlled by CATCHEM_BUILD_MUSICA)
#

if(CATCHEM_BUILD_MUSICA)
  # Set MUSICA build options before adding the subdirectory
  set(MICM_ENABLE_TESTS OFF)
  set(MUSICA_BUILD_C_CXX_INTERFACE ON)
  set(MUSICA_BUILD_FORTRAN_INTERFACE ON)
  set(MUSICA_ENABLE_INSTALL OFF)
  set(MUSICA_ENABLE_TESTS OFF)
  set(MUSICA_ENABLE_MICM ON)
  set(MUSICA_ENABLE_TUVX ON)
  set(MUSICA_ENABLE_CARMA OFF)

  # Add MUSICA as a subdirectory
  add_subdirectory(musica)
endif()

# Add our yaml interface (yaml-cpp is already available from above)
add_subdirectory(yaml_interface)
