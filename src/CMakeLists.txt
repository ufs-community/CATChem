# Always build external dependencies first (yaml-cpp is required for core)
add_subdirectory(external)

add_subdirectory(core)

# Make process modules optional (some depend on external libraries)
option(
  BUILD_PROCESSES
  "Build process modules (may require external dependencies)"
  ON
)
if(BUILD_PROCESSES)
  add_subdirectory(process)
endif()

# API depends on process modules and C++ core
option(BUILD_API "Build API (requires process modules)" ON)
if(BUILD_API)
  add_subdirectory(api)
endif()
