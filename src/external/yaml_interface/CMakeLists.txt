# YAML Interface Library
add_library(yaml_interface yaml_interface.cpp yaml_interface_mod.F90)

# Link with yaml-cpp (from MICM/MUSICA) and basic utilities
target_link_libraries(
  yaml_interface
  PRIVATE yaml-cpp::yaml-cpp catchem_utilities
)

# Set include directories
target_include_directories(
  yaml_interface
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Set Fortran module output directory
set_target_properties(
  yaml_interface
  PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include"
)

# Export the target (removed from CATChemTargets to avoid yaml-cpp export issues)
install(
  TARGETS yaml_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

# Install headers and modules
install(FILES yaml_interface.h DESTINATION include)

install(
  DIRECTORY "${CMAKE_BINARY_DIR}/modules/"
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.mod"
)
