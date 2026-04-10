# CarbChem Schemes CMakeLists.txt
# Generated on: 2026-04-10T16:46:37.844032

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - CarbChemScheme_GOCART_Mod.F90

# Scheme-specific configurations can be added here if needed

# GOCART scheme configuration
# Required meteorological fields for GOCART:
#   - DELP
#   - TSTEP
#   - AIRDEN
#   - PMID

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
  COMMENT "Validating all carbchem schemes"
)
