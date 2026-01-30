# NUOPC Cap for CATChem
# Author: Barry Baker
# Date: November 2024

set(
  _srcs
  aqmio.F90
  catchem_nuopc_cap.F90
  catchem_nuopc_interface.F90
  catchem_emis_mod.F90
  # catchem_nuopc_netcdf_out.F90
  # catchem_nuopc_cf_input.F90
  # catchem_types.F90
)

set(_lib NUOPC_CATChem_cap)

# Note: ESMF and NetCDF are already found by parent CMakeLists.txt
# find_package(ESMF MODULE REQUIRED)  # Removed to avoid conflicts
# find_package(NetCDF REQUIRED)      # Removed to avoid conflicts

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem esmf NetCDF::NetCDF_Fortran)

# target_include_directories(${_lib} PUBLIC
#   $<TARGET_PROPERTY:esmf,INTERFACE_INCLUDE_DIRECTORIES>
#   $<TARGET_PROPERTY:NetCDF::NetCDF_Fortran,INTERFACE_INCLUDE_DIRECTORIES>
# )

set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

target_compile_options(${_lib} PRIVATE ${_catchem_compiler_options})

# # Build standalone driver
# add_executable(catchem_nuopc_driver catchem_nuopc_driver.F90)
# target_link_libraries(catchem_nuopc_driver PRIVATE ${_lib})
# target_include_directories(catchem_nuopc_driver PRIVATE ${CMAKE_BINARY_DIR}/include)

# # Installation
# install(TARGETS ${_lib}
#   EXPORT CATChemTargets
#   LIBRARY DESTINATION lib
#   ARCHIVE DESTINATION lib
# )

# install(TARGETS catchem_nuopc_driver
#   RUNTIME DESTINATION bin
# )

# # Install configuration files
# install(FILES
#   ../../parm/config/catchem_nuopc_metadata.yml
#   ../../parm/catchem_output_config.yml
#   ../../parm/catchem_input_config.yml
#   ../../parm/catchem_field_config.yml
#   ../../parm/catchem_nuopc_config.yml
#   DESTINATION share/catchem/configs
# )
