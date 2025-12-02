# DryDep Schemes CMakeLists.txt
# Generated on: 2025-11-25T22:20:02.442271

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - DryDepScheme_WESELY_Mod.F90
# - DryDepScheme_GOCART_Mod.F90
# - DryDepScheme_ZHANG_Mod.F90

# Scheme-specific configurations can be added here if needed

# WESELY scheme configuration
# Required meteorological fields for WESELY:
#   - USTAR
#   - TSTEP
#   - TS
#   - SWGDN
#   - SUNCOSmid
#   - OBK
#   - CLDFRC
#   - BXHEIGHT
#   - Z0
#   - PS
#   - FRLAI
#   - ILAND
#   - SALINITY
#   - FRLANDUSE
#   - TSKIN
#   - LON
#   - LAT
#   - LUCNAME
#   - IsSnow
#   - IsIce
#   - IsLand

# GOCART scheme configuration
# Required meteorological fields for GOCART:
#   - USTAR
#   - TSTEP
#   - T
#   - AIRDEN
#   - Z
#   - LWI
#   - PBLH
#   - HFLUX
#   - Z0H
#   - U10M
#   - V10M
#   - FRLAKE
#   - GWETTOP

# ZHANG scheme configuration
# Required meteorological fields for ZHANG:
#   - USTAR
#   - TSTEP
#   - TS
#   - OBK
#   - BXHEIGHT
#   - Z0
#   - RH
#   - PS
#   - U10M
#   - V10M
#   - FRLANDUSE
#   - ILAND
#   - LUCNAME
#   - IsSnow
#   - IsIce

# Scheme validation targets
add_custom_target(
  validate_wesely_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating wesely scheme implementation"
  COMMENT "Validating wesely scheme"
)
add_custom_target(
  validate_gocart_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating gocart scheme implementation"
  COMMENT "Validating gocart scheme"
)
add_custom_target(
  validate_zhang_scheme
  COMMAND ${CMAKE_COMMAND} -E echo "Validating zhang scheme implementation"
  COMMENT "Validating zhang scheme"
)

# Combined validation target
add_custom_target(
  validate_all_schemes
  DEPENDS validate_wesely_scheme validate_gocart_scheme validate_zhang_scheme
  COMMENT "Validating all drydep schemes"
)
