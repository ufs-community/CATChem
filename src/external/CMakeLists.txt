# External dependencies

# Create a basic utilities library that other modules can depend on
add_library(catchem_utilities ../core/Precision_Mod.F90)

set_target_properties(
  catchem_utilities
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# TODO: GOCART-2G Process Library - currently disabled as source files are missing
# Uncomment when GOCART source files are available
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
# MUSICA Fortran API
#

# Set MUSICA build options before adding the subdirectory
# set(MICM_ENABLE_TESTS OFF)
# set(MICM_ENABLE_EXAMPLES OFF)
# set(MUSICA_BUILD_C_CXX_INTERFACE ON)
# set(MUSICA_BUILD_FORTRAN_INTERFACE ON)
# set(MUSICA_ENABLE_INSTALL OFF)
# set(MUSICA_ENABLE_TESTS OFF)
# set(MUSICA_ENABLE_MICM ON)
# set(MUSICA_ENABLE_TUVX ON)

# # Add MUSICA as a subdirectory
# add_subdirectory(musica)

#
# YAML-CPP standalone build (since we're not using MUSICA)
#

# Configure yaml-cpp build options
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Build yaml-cpp tests")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Build yaml-cpp tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Build yaml-cpp contrib")

# Use the local submodule
add_subdirectory(yaml-cpp)

# Add our yaml interface (now that yaml-cpp is available independently)
add_subdirectory(yaml_interface)
