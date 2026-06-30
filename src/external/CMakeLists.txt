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
# yaml-cpp — declare via FetchContent so MUSICA and its dependencies reuse it
#

include(FetchContent)

# Configure yaml-cpp build options
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Build yaml-cpp tests")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Build yaml-cpp tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Build yaml-cpp contrib")

# Declare yaml-cpp as already available via FetchContent
# This prevents MUSICA from fetching their own copy
FetchContent_Declare(yaml-cpp SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp)
FetchContent_MakeAvailable(yaml-cpp)

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
