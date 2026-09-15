# GasChem Schemes CMakeLists.txt
# Generated on: 2026-08-24T17:51:00.709076

# This file is included by the parent CMakeLists.txt
# Schemes are built as part of the main process library

# Scheme source files are defined in the parent CMakeLists.txt:
# - GasChemScheme_MUSICA_Mod.F90

# Scheme-specific configurations can be added here if needed

# MUSICA scheme configuration
# Required meteorological fields for MUSICA:
#   - TSTEP
#   - T
#   - AIRDEN
#   - PMID
#   - SPHU


# Scheme validation targets
add_custom_target(validate_no_phot_scheme
    COMMAND ${CMAKE_COMMAND} -E echo "Validating no_phot scheme implementation"
    COMMENT "Validating no_phot scheme"
)

# Combined validation target
add_custom_target(validate_all_schemes
    DEPENDS
        validate_no_phot_scheme
    COMMENT "Validating all GasChem schemes"
)