# WetDep Schemes CMakeLists.txt
# Generated on: 2025-12-15T16:30:33.761046

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - WetDepScheme_JACOB_Mod.F90

# Scheme-specific configurations can be added here if needed

# JACOB scheme configuration
# Required meteorological fields for JACOB:
#   - T
#   - TSTEP
#   - AIRDEN_DRY
#   - MAIRDEN
#   - PFLLSAN
#   - PFILSAN
#   - PEDGE
#   - REEVAPLS

# Scheme validation targets
add_custom_target(
  validate_jacob_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating jacob scheme implementation"
  COMMENT "Validating jacob scheme"
)

# Combined validation target
add_custom_target(
  validate_all_schemes
  DEPENDS validate_jacob_scheme
  COMMENT "Validating all wetdep schemes"
)
