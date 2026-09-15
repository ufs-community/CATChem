cmake_minimum_required(VERSION 3.20.2)

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

# Deprecated: ENABLE_KOKKOS was renamed to CATCHEM_ENABLE_KOKKOS to match the
# other CATCHEM_* options. Honor the old name for out-of-tree callers.
if(DEFINED ENABLE_KOKKOS AND NOT DEFINED CATCHEM_ENABLE_KOKKOS)
  message(DEPRECATION "ENABLE_KOKKOS is deprecated; use CATCHEM_ENABLE_KOKKOS")
  set(CATCHEM_ENABLE_KOKKOS "${ENABLE_KOKKOS}")
endif()

# Option to enable Kokkos GPU/parallel support (defaults to OFF)
option(CATCHEM_ENABLE_KOKKOS "Enable Kokkos GPU/parallel support" OFF)

# Option to build with MUSICA chemistry library (defaults to OFF)
option(CATCHEM_BUILD_MUSICA "Build with MUSICA chemistry library" OFF)

# MUSICA's gaschem/photolysis processes build on the Kokkos-backed C++ core.
if(CATCHEM_BUILD_MUSICA AND NOT CATCHEM_ENABLE_KOKKOS)
  message(
    FATAL_ERROR
    "CATCHEM_BUILD_MUSICA=ON requires CATCHEM_ENABLE_KOKKOS=ON (gaschem and "
    "photolysis are built on the Kokkos-backed C++ core)"
  )
endif()

# C++ standard: C++17 everywhere except MUSICA builds — the mdspan data
# layer and the C++ core need only 17, and 17 has broader support on the
# UFS HPC toolchains; MUSICA requires C++20.
if(CATCHEM_BUILD_MUSICA)
  set(CMAKE_CXX_STANDARD 20)
else()
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CATCHEM_ENABLE_KOKKOS)
  if(CMAKE_VERSION VERSION_LESS "3.16")
    message(
      FATAL_ERROR
      "CATCHEM_ENABLE_KOKKOS requires CMake 3.16 or later (Kokkos minimum)"
    )
  endif()

  # CATChem's C++ core depends on mdspan support (Kokkos::mdspan and the
  # <mdspan/mdspan.hpp> headers), which is only available in modern Kokkos
  # (5.x). Older distro packages (e.g. Ubuntu's Kokkos 3.7) are not usable.
  # Prefer a sufficiently new system install and otherwise fetch and build a
  # pinned release.
  set(
    CATCHEM_KOKKOS_MIN_VERSION
    "5.0"
    CACHE STRING
    "Minimum acceptable Kokkos version"
  )
  set(
    CATCHEM_KOKKOS_VERSION
    "5.1.1"
    CACHE STRING
    "Kokkos release to fetch when a suitable system install is not found"
  )

  find_package(Kokkos ${CATCHEM_KOKKOS_MIN_VERSION} QUIET)

  if(Kokkos_FOUND)
    message(STATUS "Using system Kokkos ${Kokkos_VERSION}")

    # Work around a packaging bug in some Kokkos installs where the exported
    # Kokkos targets advertise a non-existent "/include" path in
    # INTERFACE_INCLUDE_DIRECTORIES (its _IMPORT_PREFIX resolves to "/"). CMake
    # treats this as a hard error at generate time, so strip any include
    # directories that do not actually exist on disk. This must be applied to
    # every imported Kokkos target because Kokkos::kokkos inherits include
    # directories transitively from its component targets.
    foreach(
      _kokkos_target
      Kokkos::kokkos
      Kokkos::kokkoscore
      Kokkos::kokkoscontainers
      Kokkos::kokkosalgorithms
      Kokkos::kokkossimd
    )
      if(TARGET ${_kokkos_target})
        get_target_property(
          _kokkos_include_dirs
          ${_kokkos_target}
          INTERFACE_INCLUDE_DIRECTORIES
        )
        if(_kokkos_include_dirs)
          set(_kokkos_include_dirs_clean "")
          foreach(_kokkos_inc IN LISTS _kokkos_include_dirs)
            if(_kokkos_inc MATCHES "\\$<" OR EXISTS "${_kokkos_inc}")
              list(APPEND _kokkos_include_dirs_clean "${_kokkos_inc}")
            else()
              message(
                STATUS
                "Dropping non-existent include directory '${_kokkos_inc}' from ${_kokkos_target}"
              )
            endif()
          endforeach()
          set_target_properties(
            ${_kokkos_target}
            PROPERTIES
              INTERFACE_INCLUDE_DIRECTORIES "${_kokkos_include_dirs_clean}"
          )
        endif()
      endif()
    endforeach()
  else()
    message(
      STATUS
      "Suitable Kokkos (>= ${CATCHEM_KOKKOS_MIN_VERSION}) not found; fetching and building Kokkos ${CATCHEM_KOKKOS_VERSION}"
    )
    include(FetchContent)
    # Enable the OpenMP host backend by default for the fetched build.
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
    FetchContent_Declare(
      kokkos
      URL
        "https://github.com/kokkos/kokkos/releases/download/${CATCHEM_KOKKOS_VERSION}/kokkos-${CATCHEM_KOKKOS_VERSION}.tar.gz"
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
      SYSTEM # treat the fetched headers as system headers (no warnings from them)
    )
    FetchContent_MakeAvailable(kokkos)
  endif()
else()
  # Host-only build: the C++ core's data layer is Kokkos::mdspan, provided
  # without Kokkos by a vendored single-header build of the kokkos/mdspan
  # reference implementation (src/external/mdspan — see its README for
  # provenance and regeneration). Vendored rather than fetched because UFS
  # compiles on network-restricted compute nodes.
  add_library(mdspan INTERFACE)
  add_library(mdspan::mdspan ALIAS mdspan)
  target_include_directories(
    mdspan
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src/external/mdspan/include"
  )
endif()
# Integration build options. CATCHEM_BUILD_NUOPC is authoritative; retain the
# historical BUILD_NUOPC cache variable only as a migration path.
set(_catchem_nuopc_option_was_defined FALSE)
if(DEFINED CATCHEM_BUILD_NUOPC)
  set(_catchem_nuopc_option_was_defined TRUE)
endif()
option(CATCHEM_BUILD_NUOPC "Build NUOPC cap and drivers" OFF)
if(DEFINED BUILD_NUOPC)
  if(
    _catchem_nuopc_option_was_defined
    AND NOT "${BUILD_NUOPC}" STREQUAL "${CATCHEM_BUILD_NUOPC}"
  )
    message(
      FATAL_ERROR
      "BUILD_NUOPC=${BUILD_NUOPC} conflicts with "
      "CATCHEM_BUILD_NUOPC=${CATCHEM_BUILD_NUOPC}; remove BUILD_NUOPC and "
      "use the CATCHEM_* option"
    )
  endif()
  message(DEPRECATION "BUILD_NUOPC is deprecated; use CATCHEM_BUILD_NUOPC")
  set(
    CATCHEM_BUILD_NUOPC
    "${BUILD_NUOPC}"
    CACHE BOOL
    "Build NUOPC cap and drivers"
    FORCE
  )
endif()
option(CATCHEM_BUILD_CCPP "Build the CCPP host interface" OFF)

# Debug builds treat compiler warnings as errors so the tree cannot drift back
# from zero warnings. OFF is the escape hatch for a toolchain that reports
# something CATChem cannot fix. Release builds are unaffected (the warning
# flag lists below are Debug-only).
option(
  CATCHEM_WARNINGS_AS_ERRORS
  "Treat Debug-build compiler warnings as errors"
  ON
)
option(
  CATCHEM_CCPP_TEST_STUBS
  "Build minimal CCPP framework stubs for contract testing"
  OFF
)
set(
  CATCHEM_CCPP_HOST_MODULE_DIR
  ""
  CACHE PATH
  "Directory containing host CCPP Fortran modules"
)
option(CATCHEM_BUILD_TESTING "Build CATChem testing framework" ON)
option(
  CATCHEM_ENABLE_SANITIZERS
  "Enable target-scoped address and undefined-behavior sanitizers"
  OFF
)
option(
  CATCHEM_ENABLE_BOUNDS_CHECK
  "Enable target-scoped Fortran bounds checking"
  OFF
)
option(CATCHEM_TRACE_NUOPC "Enable ESMF trace regions in NUOPC cap" OFF)
# Option to build the standalone ESMF/NUOPC driver + application (requires NUOPC)
option(
  CATCHEM_BUILD_STANDALONE
  "Build standalone ESMF/NUOPC driver and app"
  OFF
)

if(CATCHEM_TRACE_NUOPC AND NOT CATCHEM_BUILD_NUOPC)
  message(FATAL_ERROR "CATCHEM_TRACE_NUOPC=ON requires CATCHEM_BUILD_NUOPC=ON")
endif()

if(CATCHEM_TRACE_NUOPC)
  add_compile_definitions(CATCHEM_TRACE_NUOPC)
endif()

# CATCHEM_BUILD_MUSICA and the C++ standard are configured next to
# CATCHEM_ENABLE_KOKKOS near the top of this file.

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

add_library(catchem_validation_options INTERFACE)
if(CATCHEM_ENABLE_SANITIZERS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(
      catchem_validation_options
      INTERFACE
        "$<$<COMPILE_LANGUAGE:C,CXX>:-fsanitize=address,undefined;-fno-omit-frame-pointer>"
    )
    target_link_options(
      catchem_validation_options
      INTERFACE -fsanitize=address,undefined
    )
  else()
    message(
      FATAL_ERROR
      "CATCHEM_ENABLE_SANITIZERS is unsupported by ${CMAKE_CXX_COMPILER_ID}"
    )
  endif()
endif()
if(CATCHEM_ENABLE_BOUNDS_CHECK)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
    target_compile_options(
      catchem_validation_options
      INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-fcheck=bounds>"
    )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    # SHELL: keeps "-check bounds" together. Without it CMake de-duplicates
    # the repeated "-check" token against the Debug list's "-check all,nouninit"
    # and ifx receives a bare "bounds" ("error #10236: File not found: 'bounds'").
    target_compile_options(
      catchem_validation_options
      INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-check bounds>"
    )
  else()
    message(
      FATAL_ERROR
      "CATCHEM_ENABLE_BOUNDS_CHECK is unsupported by ${CMAKE_Fortran_COMPILER_ID}"
    )
  endif()
endif()

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

# C/C++ warning flags for CATChem-owned translation units (core, process
# wrappers, tests). Debug only. Applied per directory so that vendored code
# under src/external and fetched dependencies (Kokkos) are not affected.
set(_catchem_cxx_warning_options "")
if(
  CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang|IntelLLVM"
  AND CMAKE_BUILD_TYPE MATCHES "Debug"
)
  # -pedantic: reject GNU extensions and features from a newer standard than
  # the one selected (the host-only build is C++17 for UFS compatibility).
  list(APPEND _catchem_cxx_warning_options -Wall -Wextra -pedantic)
  if(CATCHEM_WARNINGS_AS_ERRORS)
    list(APPEND _catchem_cxx_warning_options -Werror)
  endif()
endif()

if(CATCHEM_WARNINGS_AS_ERRORS AND CMAKE_BUILD_TYPE MATCHES "Debug")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|G95")
    list(APPEND _catchem_compiler_options -Werror)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    # SHELL: keeps "-warn errors" as one unit; the list already holds
    # "-warn all" and CMake would otherwise de-duplicate the second "-warn".
    list(APPEND _catchem_compiler_options "SHELL:-warn errors")
  endif()
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

#find MPI for standalone builds
if(CATCHEM_BUILD_STANDALONE)
  find_package(MPI REQUIRED)
endif()

if(CATCHEM_BUILD_NUOPC)
  # Resolve ESMF target: prefer uppercase ESMF (FindESMF), fall back to
  # lowercase esmf (MPAS internal), or find_package as last resort.
  if(TARGET ESMF)
    set(CATCHEM_ESMF_TARGET ESMF)
  elseif(TARGET esmf)
    set(CATCHEM_ESMF_TARGET esmf)
  elseif(TARGET ESMF::ESMF)
    set(CATCHEM_ESMF_TARGET ESMF::ESMF)
  else()
    find_package(ESMF REQUIRED)
    # find_package may create the target under different names depending on the
    # provider: the repo's FindESMF.cmake creates `esmf`, while an installed
    # ESMF CMake config package (e.g. Homebrew/upstream) creates `ESMF::ESMF`
    # (or `ESMF`). Re-resolve the actual target instead of assuming `esmf`,
    # otherwise CMake passes a bogus `-lesmf` flag and the link fails.
    if(TARGET esmf)
      set(CATCHEM_ESMF_TARGET esmf)
    elseif(TARGET ESMF::ESMF)
      set(CATCHEM_ESMF_TARGET ESMF::ESMF)
    elseif(TARGET ESMF)
      set(CATCHEM_ESMF_TARGET ESMF)
    else()
      message(
        FATAL_ERROR
        "find_package(ESMF) succeeded but created no usable target "
        "(expected esmf, ESMF::ESMF, or ESMF)."
      )
    endif()
  endif()
  message(STATUS "CATChem: using ESMF target '${CATCHEM_ESMF_TARGET}'")
endif()

# Precision selection for the Fortran science bridges/schemes.
#
# DOUBLE precision is the default.  The bridge precision module
# (catchem_bridge_precision) sets `fp = f8` unless USE_REAL4 is explicitly
# defined, matching how GOCART/UFS build the science (-fdefault-real-8 / -r8)
# and matching the real-8 promotion every process library + the
# catchem_utilities target apply.  A compile-time guard in that module fails
# the build if `fp` ever diverges from the actual default-real promotion, so
# the precision cannot silently mismatch across the C++/Fortran boundary
# regardless of how a parent build (e.g. UFS `add_subdirectory(CATChem)`)
# configures flags.
#
# To build single precision, define USE_REAL4 AND drop the -fdefault-real-8 /
# -r8 promotion from the scheme + utilities targets so the two stay consistent.

# Find yaml-cpp for configuration parsing
if(NOT TARGET yaml-cpp::yaml-cpp AND NOT TARGET yaml-cpp)
  find_package(yaml-cpp QUIET)
endif()

add_subdirectory(src)

if(CATCHEM_BUILD_CCPP)
  add_subdirectory(drivers/ccpp)
endif()

# Enable testing for standalone builds
if(CATCHEM_BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if(CATCHEM_BUILD_NUOPC)
  find_package(MPI COMPONENTS Fortran REQUIRED)
  add_subdirectory(drivers/nuopc)

  # Fortran consumers must see generated .mod files while compiling, not only
  # when linking cc_nuopc.  UFSDriver uses cc_nuopc directly, so keep the
  # module directory on the target's public build interface.
  set(CATCHEM_NUOPC_MODULE_DIR "${CMAKE_BINARY_DIR}/include")

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
    drivers/nuopc/catchem_emis_data_mod.F90
    drivers/nuopc/catchem_emis_mod.F90
    drivers/nuopc/catchem_nuopc_cap.F90
    drivers/nuopc/catchem_nuopc_interface.F90
  )

  # UFS integration consumes this duplicate, self-contained NUOPC archive;
  # propagate the trace option directly to it rather than relying on the
  # directory-wide definition above.
  if(CATCHEM_TRACE_NUOPC)
    target_compile_definitions(cc_nuopc PRIVATE CATCHEM_TRACE_NUOPC)
  endif()

  # Link only to external dependencies (not internal CATChem targets)
  target_link_libraries(
    cc_nuopc
    PUBLIC ${CATCHEM_ESMF_TARGET} NetCDF::NetCDF_Fortran MPI::MPI_Fortran
  )

  # Manually handle the CATChem dependencies by linking to built libraries
  # Only need to depend on CATChem - it already depends on CATChem_core, etc.
  add_dependencies(cc_nuopc CATChem)

  # cc_nuopc compiles the same drivers/nuopc sources as NUOPC_CATChem_cap and
  # both emit their .mod files into ${CATCHEM_NUOPC_MODULE_DIR}. Without an
  # ordering between them, a parallel build lets two gfortran instances
  # rename the same <module>.mod0 -> <module>.mod file and one of them fails
  # with "Cannot rename module file". Serialise the two compilations.
  add_dependencies(cc_nuopc NUOPC_CATChem_cap)

  # Link to library files in dependency order (most dependent first).
  # NOTE: $<TARGET_FILE:> carries no transitive usage requirements, so this
  # list must mirror the full dependency closure of the CATChem targets by
  # hand. CATChem_core_cpp implements the C API the Fortran core and API
  # bind to (host-only or Kokkos flavor); yaml-cpp (below) is its config
  # parser. A Kokkos-enabled UFS build will additionally need
  # Kokkos::kokkos here.
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
      $<TARGET_FILE:CATChem_core_cpp>
      $<TARGET_FILE:catchem_utilities>
  )

  # yaml-cpp for the C++ core's config parsing. Prefer target-based linking
  # here (not $<TARGET_FILE:...>) so both real and ALIAS/imported targets work
  # across vendored and system/spack provisioning modes.
  if(TARGET yaml-cpp::yaml-cpp)
    target_link_libraries(cc_nuopc PRIVATE yaml-cpp::yaml-cpp)
  elseif(TARGET yaml-cpp)
    target_link_libraries(cc_nuopc PRIVATE yaml-cpp)
  endif()

  # Link GOCART2G library if available (needed for drydep process)
  if(TARGET GOCART2G)
    target_link_libraries(cc_nuopc PRIVATE $<TARGET_FILE:GOCART2G>)
  endif()

  # Include directories for compilation
  target_include_directories(
    cc_nuopc
    PUBLIC $<BUILD_INTERFACE:${CATCHEM_NUOPC_MODULE_DIR}>
    #$<TARGET_PROPERTY:CATChem,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:esmf,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:NetCDF::NetCDF_Fortran,INTERFACE_INCLUDE_DIRECTORIES>
  )

  set_target_properties(
    cc_nuopc
    PROPERTIES
      # Keep generated module interfaces inside the build tree.  A path
      # relative to the CMake working directory resolves outside a container
      # build and makes downstream NUOPC applications unable to compile
      # against cc_nuopc in a clean build.
      Fortran_MODULE_DIRECTORY ${CATCHEM_NUOPC_MODULE_DIR}
      EXPORT_NAME cc_nuopc
  )

  # UFS configurations historically add <build>/mod to their Fortran module
  # search path.  Keep that layout usable while the public target interface
  # above remains the canonical dependency-safe path.
  add_custom_command(
    TARGET cc_nuopc
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/mod
    COMMAND
      ${CMAKE_COMMAND} -E copy_directory ${CATCHEM_NUOPC_MODULE_DIR}
      ${CMAKE_BINARY_DIR}/mod
    COMMENT
      "Publishing CATChem NUOPC Fortran modules to mod/ for UFS compatibility"
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

  # If CATChem builds an in-tree yaml-cpp target, export it with cc_nuopc so
  # parent install(EXPORT ufs-config ...) does not fail on missing export-set
  # dependencies. Imported/system yaml-cpp targets are excluded.
  if(TARGET yaml-cpp)
    get_target_property(_catchem_yamlcpp_is_imported yaml-cpp IMPORTED)
    if(NOT _catchem_yamlcpp_is_imported)
      install(
        TARGETS yaml-cpp
        EXPORT ufs-config
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
      )
    endif()
  endif()

  # Standalone ESMF/NUOPC driver + application (built on top of cc_nuopc)
  if(CATCHEM_BUILD_STANDALONE)
    add_subdirectory(drivers/standalone)
  endif()
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
