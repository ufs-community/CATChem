# External dependencies

# Create a basic utilities library that other modules can depend on
add_library(catchem_utilities ../process/catchem_bridge_utilities.F90)

set_target_properties(
  catchem_utilities
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# The precision module MUST be compiled with the SAME real-8 promotion the
# process libraries use unconditionally, so its kind(0.0) is 8-byte and the
# compile-time precision guard validates against the real promotion.  Do NOT
# gate this on the global Fortran flags: the promotion is per-target.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
  target_compile_options(
    catchem_utilities
    PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8;-fdefault-double-8>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
  target_compile_options(
    catchem_utilities
    PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-r8>
  )
endif()

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
# GOCART2G holds the shared science kernels (Chem_Settling, wet/dry removal,
# emission helpers) called by the process schemes.  It MUST be promoted to
# real-8 unconditionally to match those schemes and the `fp` kind; gating on
# the global Fortran flags (which the per-target promotion does not set) would
# leave GOCART2G in single precision while its callers run double -- a silent
# precision/ABI mismatch across the call boundary.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
  target_compile_options(${_lib} PRIVATE -ffree-line-length-none)
  target_compile_options(
    ${_lib}
    PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8;-fdefault-double-8>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
  target_compile_options(${_lib} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-r8>)
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

  # MUSICA and the dependencies it fetches (micm, mechanism_configuration) are
  # ordinary in-tree targets, so their include directories are not treated as
  # system headers and CATChem's Debug warning flags would report findings
  # inside them. Re-export their include directories as SYSTEM.
  foreach(
    _dep
    musica
    musica::micm
    micm::micm
    micm
    musica::mechanism_configuration
    mechanism_configuration
  )
    if(TARGET ${_dep})
      get_target_property(_dep_real ${_dep} ALIASED_TARGET)
      if(NOT _dep_real)
        set(_dep_real ${_dep})
      endif()
      set_property(
        TARGET ${_dep_real}
        APPEND
        PROPERTY
          INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
            $<TARGET_PROPERTY:${_dep_real},INTERFACE_INCLUDE_DIRECTORIES>
      )
    endif()
  endforeach()

  # musica-fortran exports -ffree-line-length-none as a PUBLIC compile option
  # with no language guard (src/external/musica/fortran/CMakeLists.txt), so
  # every C++ translation unit that links it gets
  # "cc1plus: warning: command-line option '-ffree-line-length-none' is valid
  # for Fortran but not for C++". Re-scope its interface compile options to
  # Fortran here rather than patching the submodule.
  if(TARGET musica-fortran)
    get_target_property(_mf_opts musica-fortran INTERFACE_COMPILE_OPTIONS)
    if(_mf_opts)
      set(_mf_scoped "")
      foreach(_opt IN LISTS _mf_opts)
        if(_opt MATCHES "COMPILE_LANGUAGE")
          list(APPEND _mf_scoped "${_opt}")
        else()
          list(APPEND _mf_scoped "$<$<COMPILE_LANGUAGE:Fortran>:${_opt}>")
        endif()
      endforeach()
      set_target_properties(
        musica-fortran
        PROPERTIES INTERFACE_COMPILE_OPTIONS "${_mf_scoped}"
      )
    endif()
  endif()
endif()
