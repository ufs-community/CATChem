# WetDep Schemes CMakeLists.txt
# Generated on: 2026-09-22T13:06:55.613528

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - WetDepScheme_JACOB_Mod.F90
# - WetDepScheme_GOCART_Mod.F90

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

# GOCART scheme configuration
# Required meteorological fields for GOCART:
#   - T
#   - TSTEP
#   - MAIRDEN
#   - PEDGE
#   - PFLLSAN
#   - PFILSAN
#   - PRECCON
#   - PRECLSC

# Scheme validation targets
add_custom_target(
  validate_jacob_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating jacob scheme implementation"
  COMMENT "Validating jacob scheme"
)
add_custom_target(
  validate_gocart_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating gocart scheme implementation"
  COMMENT "Validating gocart scheme"
)

# Combined validation target
add_custom_target(
  validate_all_schemes
  DEPENDS validate_jacob_scheme validate_gocart_scheme
  COMMENT "Validating all wetdep schemes"
)
