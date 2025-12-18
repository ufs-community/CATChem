# Settling Schemes CMakeLists.txt
# Generated on: 2025-12-18T14:12:33.125087

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - SettlingScheme_GOCART_Mod.F90

# Scheme-specific configurations can be added here if needed

# GOCART scheme configuration
# Required meteorological fields for GOCART:
#   - T
#   - TSTEP
#   - AIRDEN
#   - RH
#   - Z
#   - PMID
#   - DELP

# Scheme validation targets
add_custom_target(
  validate_gocart_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating gocart scheme implementation"
  COMMENT "Validating gocart scheme"
)

# Combined validation target
add_custom_target(
  validate_all_schemes
  DEPENDS validate_gocart_scheme
  COMMENT "Validating all settling schemes"
)
