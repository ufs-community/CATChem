# =====================================================================
# CATChem standalone ESMF/NUOPC driver and application
# =====================================================================
# Builds:
#   - catchem_standalone_driver : library with the grid utilities and the
#     NUOPC_Driver component (catchem_driver)
#   - catchem_app               : standalone executable (main program)
#
# Requires the CATChem NUOPC cap target `cc_nuopc` to have been defined by
# the top-level CMakeLists (i.e. CATCHEM_BUILD_NUOPC=ON).
# =====================================================================

if(NOT TARGET cc_nuopc)
  message(
    FATAL_ERROR
    "drivers/standalone requires the CATChem NUOPC cap (cc_nuopc). "
    "Enable it with -DCATCHEM_BUILD_NUOPC=ON."
  )
endif()

set(
  _srcs
  catchem_standalone_grid_mod.F90
  catchem_driver.F90
)

add_library(catchem_standalone_driver STATIC ${_srcs})
add_dependencies(catchem_standalone_driver cc_nuopc)

# cc_nuopc pulls in the CATChem core/process archives by raw file path
# ($<TARGET_FILE:...>), which deliberately hides them from CMake's dependency
# graph so the exported cc_nuopc target doesn't drag in the non-exported
# yaml-cpp. That raw-path trick also drops yaml_interface's transitive yaml-cpp
# dependency AND the fact that it is C++. By linking yaml_interface here as a
# real CMake target, CMake restores both usage requirements automatically: it
# pulls in yaml-cpp and, seeing a C++ library in the link closure, adds the C++
# standard library to the final (Fortran) link. No library names are hardcoded.
target_link_libraries(
  catchem_standalone_driver
  PUBLIC cc_nuopc yaml_interface ${CATCHEM_ESMF_TARGET} NetCDF::NetCDF_Fortran
)

# cc_nuopc writes its Fortran modules to <repo>/mod; CATChem core modules go
# to <build>/include. Make both visible.
target_include_directories(
  catchem_standalone_driver
  PRIVATE
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_BINARY_DIR}/../mod
)

set_target_properties(
  catchem_standalone_driver
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

target_compile_options(catchem_standalone_driver PRIVATE ${_catchem_compiler_options})

# ---------------------------------------------------------------------
# Standalone executable
# ---------------------------------------------------------------------
add_executable(catchem_app catchem_app.F90)
target_link_libraries(catchem_app PRIVATE catchem_standalone_driver)
target_include_directories(catchem_app PRIVATE ${CMAKE_BINARY_DIR}/include)
target_compile_options(catchem_app PRIVATE ${_catchem_compiler_options})

# Force the Fortran linker. The main program is Fortran, but the C++ yaml-cpp
# dependency in the link closure would otherwise make CMake pick the C++ linker
# (C++ has a higher linker preference than Fortran), which then rejects
# gfortran-only flags such as -fopenmp. With the Fortran linker, CMake still
# adds the C++ runtime automatically because yaml_interface is linked as a
# C++ target. This mirrors the existing CATChem process tests.
set_target_properties(
  catchem_app
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    LINKER_LANGUAGE Fortran
)

# Copy the example driver configure file next to the build for convenience
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/catchem_standalone.configure
  ${CMAKE_BINARY_DIR}/catchem_standalone.configure
  COPYONLY
)

# Copy the CATChem science configuration YAML files (config + species +
# emissions) next to the binary's working directory so a standalone run can
# find them. catchem_app resolves catchem_config_file: (and the species /
# emission paths it references) relative to the current working directory,
# which is the build root where catchem_standalone.configure lives. The
# canonical set lives under tests/Configs/Default, the same files the process
# tests use.
file(GLOB _catchem_standalone_yml "${PROJECT_SOURCE_DIR}/tests/Configs/Default/*.yml")
foreach(_yml ${_catchem_standalone_yml})
  get_filename_component(_yml_name ${_yml} NAME)
  configure_file(${_yml} ${CMAKE_BINARY_DIR}/${_yml_name} COPYONLY)
endforeach()

# ---------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------
install(
  TARGETS catchem_app
  RUNTIME DESTINATION bin COMPONENT drivers
)

install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/catchem_standalone.configure
  DESTINATION share/catchem COMPONENT drivers
)

message(STATUS "CATChem: standalone driver/app (catchem_app) enabled")
