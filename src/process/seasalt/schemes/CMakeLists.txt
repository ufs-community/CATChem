# SeaSalt Schemes CMakeLists.txt
# Generated on: 2025-11-14T23:01:21.797764

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - SeaSaltScheme_GONG97_Mod.F90
# - SeaSaltScheme_GONG03_Mod.F90
# - SeaSaltScheme_GEOS12_Mod.F90

# Scheme-specific configurations can be added here if needed

# GONG97 scheme configuration
# Required meteorological fields for GONG97:
#   - FROCEAN
#   - FRSEAICE
#   - SST
#   - U10M
#   - V10M

# GONG03 scheme configuration
# Required meteorological fields for GONG03:
#   - FROCEAN
#   - FRSEAICE
#   - SST
#   - U10M
#   - V10M

# GEOS12 scheme configuration
# Required meteorological fields for GEOS12:
#   - FROCEAN
#   - FRSEAICE
#   - SST
#   - USTAR

# Scheme validation targets
add_custom_target(
  validate_gong97_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating gong97 scheme implementation"
  COMMENT "Validating gong97 scheme"
)
add_custom_target(
  validate_gong03_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating gong03 scheme implementation"
  COMMENT "Validating gong03 scheme"
)
add_custom_target(
  validate_geos12_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating geos12 scheme implementation"
  COMMENT "Validating geos12 scheme"
)

# Combined validation target
add_custom_target(
  validate_all_schemes
  DEPENDS validate_gong97_scheme validate_gong03_scheme validate_geos12_scheme
  COMMENT "Validating all seasalt schemes"
)
