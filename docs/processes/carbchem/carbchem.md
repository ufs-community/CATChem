# Carbon Chemistry Process

**Process Type:** Chemistry
**Description:** Process for computing chemical production and loss of carbon species
**Author:** Wei Li
**Generated:** 2026-04-10T16:46:38.030407

## Overview

The CarbChem process implements Process for computing chemical production and loss of carbon species. This process provides a modular, extensible framework for chemistry calculations within the CATChem chemical transport model.

## Available Schemes

### GOCART Scheme

**Name:** `gocart`
**Description:** GOCART carbon species chemical production and loss scheme
**Author:** Wei Li
**Reference:** GOCART2G process library carbonChemLoss function; GOCART-2G described in Collow et al. [2024], <https://doi.org/10.5194/gmd-17-1443-2024>
#### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|--------|-------------|
| `time_days_hydrophobic_to_hydrophilic` | 2.5 |  -  | Rate of conversion of hydrophobic to hydrophilic [days] |

#### Required Meteorological Fields

- `DELP` - Delta-P (wet) across box [Pa]
- `TSTEP` - Model time step [s]
- `AIRDEN` - Wet air density [kg/m3]
- `PMID` - Average wet air pressure [Pa] defined as arithmetic average of edge pressures



## Process Interface

### Species

The carbchem process operates on the following chemical species:

- `oc1`
- `oc2`
- `bc1`
- `bc2`

### Required Inputs



### Process Diagnostics

| Diagnostic | Units | Description |
|------------|-------|-------------|
| `Production_mass_per_species_per_level` | kg/kg | Production mass (negative for loss) per species per level |

## Usage

### Basic Integration

The process is a C++ class, `catchem::CarbChemProcess` (`src/process/carbchem/catchem_process_carbchem.cpp`), registered with the `ProcessRegistry` at start-up; a host model does not instantiate it directly. It is enabled and configured in the runtime YAML and runs when the core advances the configured `run_phases`:

```yaml
processes:
  carbchem:
    activate: true
    scheme: 'gocart'
```

Hosts drive the core through the Fortran proxy `CATChem_Model` (`src/api/CATChem_API.F90`) or `catchem::Core` in C++; see the [Process Interface API](../../api/process-interface.md). The science kernel itself stays in Fortran and is called through `CarbchemScienceBridge.F90`.

### Scheme Selection

Available schemes: `gocart`. Select one in the YAML configuration; scheme-specific parameters live under a key with the scheme's name (see the parameter tables above):

```yaml
processes:
  carbchem:
    scheme: 'gocart'
```

## Implementation Details

### Pure Science Kernels

Each scheme is implemented as a pure science kernel with no infrastructure dependencies:

```fortran
! GOCART scheme
pure subroutine compute_gocart( &
   num_layers, num_species, params, &
   DELP, &   TSTEP, &   AIRDEN, &   PMID, &
   species_conc, emission_flux)
```

### Host Model Responsibilities

The host model (CATChem infrastructure) handles:

- Parameter initialization and validation
- Input array validation and error handling
- Memory management and array allocation
- Integration with model time stepping
- Diagnostic output management

## Configuration

### YAML Configuration Example

```yaml
processes:
  carbchem:
    enabled: true
    scheme: "gocart"
    parameters:
      time_days_hydrophobic_to_hydrophilic: 2.5
    diagnostics:
      enabled: true
      output_frequency: "daily"
```

## Technical Specifications

- **Parallelization:** Column
- **Memory Requirements:** Low
- **Timestep Dependency:** Independent
- **Multiphase Support:** No
- **Size Bin Support:** No
- **Vectorization:** Supported

## Files

### Source Code
- `src/process/carbchem/CarbChemCommon_Mod.F90` - Scheme options, constants, and shared types used by the schemes
- `src/process/carbchem/CarbChemScienceBridge.F90` - `bind(C)` bridge that unpacks the C++ state and calls the Fortran scheme
- `src/process/carbchem/catchem_process_carbchem.cpp` - C++ process class `catchem::CarbChemProcess`: configuration, species resolution, diagnostics registration, and the per-step call into the bridge
- `src/process/carbchem/catchem_process_carbchem.hpp` - C++ process class `catchem::CarbChemProcess` (declaration)
- `src/process/carbchem/schemes/CarbChemScheme_GOCART_Mod.F90` - GOCART carbon species chemical production and loss scheme

### Tests
- `tests/test_catchem_carbchem.cpp`

### Documentation
- `docs/processes/carbchem/carbchem.md` - This documentation

## Contributing

When modifying or extending this process:

1. **Science Changes:** Modify the scheme modules in `schemes/`
2. **Interface Changes:** Update the main interface module
3. **New Schemes:** Add new scheme modules and update the creator
4. **Tests:** Add corresponding unit and integration tests
5. **Documentation:** Update this documentation file

## References

- GOCART: GOCART2G process library carbonChemLoss function; GOCART-2G described in Collow, A. B., et al. [2024], Geosci. Model Dev., 17, 1443–1468, <https://doi.org/10.5194/gmd-17-1443-2024>

---
*This documentation was automatically generated by the CATChem Process Generator on 2026-04-10T16:46:38.030407*
