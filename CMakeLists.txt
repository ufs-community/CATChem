cmake_minimum_required(VERSION 3.5)

project(
  CATChem
  VERSION 0.1.0
  LANGUAGES Fortran)

# If GNU compiler and Debug build type, add extra debug flags
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95" AND CMAKE_BUILD_TYPE MATCHES "Debug")
  add_compile_options(
    -g
    -Wall
    -Wextra
    -Og
    -pedantic
    -fcheck=all
    -ffpe-trap=invalid,zero,overflow
    -fbacktrace
    -fmax-errors=5)
endif()

set(_param_src
    parameters/catchem_config.F90
    parameters/catchem_constants.F90
    parameters/ddep_data_mod.F90
    parameters/dust_data_mod.F90
    parameters/opt_data_mod.F90
    parameters/plume_data_mod.F90
    parameters/seas_data_mod.F90)
add_library(CATChem_params STATIC ${_param_src})

set(_scheme_src
    src/aero_soa_vbs_data_mod.F90
    src/dep_vertmx_mod.F90
    src/drydep_gocart_mod.F90
    src/drydep_wesely_mod.F90
    src/dust_afwa_mod.F90
    src/dust_fengsha_mod.F90
    src/gocart_aerosols_mod.F90
    src/gocart_chem_mod.F90
    src/gocart_diag_mod.F90
    src/gocart_dmsemis_mod.F90
    src/gocart_dust_default_mod.F90
    src/gocart_dust_simple_mod.F90
    src/gocart_seas_default_mod.F90
    src/gocart_seas_ngac_mod.F90
    src/gocart_seas_simple_mod.F90
    src/gocart_settling_mod.F90
    src/opt_aer_mod.F90
    src/opt_aer_out_mod.F90
    src/opt_aer_ra_mod.F90
    src/opt_averaging_mod.F90
    src/opt_driver_mod.F90
    src/opt_gocart_mod.F90
    src/opt_mod.F90
    src/plume_rise_mod.F90
    src/plume_scalar_mod.F90
    src/plume_zero_mod.F90
    src/vash_settling_mod.F90)
add_library(CATChem_schemes STATIC ${_scheme_src})
target_link_libraries(CATChem_schemes PUBLIC CATChem_params)
