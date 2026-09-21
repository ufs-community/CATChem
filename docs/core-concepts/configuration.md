---
type: explanation
category: explanation
tags: [configuration, yaml, cpp-core, config-manager]
---

# Configuration System

This section describes the C++ configuration system in CATChem (`catchem::ConfigManager`), which provides type-safe, hierarchical YAML configuration parsing.

## Overview

The configuration system in CATChem is based on YAML, parsed natively in C++ via `yaml-cpp`. The system is centered around `catchem::ConfigManager`, which loads, validates, and exposes typed structures for simulation parameters, grid specifications, process activation settings, diagnostic settings, species properties, and emission category mappings.

## Core Concepts

### Hierarchical Structure & Type Safety

Configuration files (`CATChem_config.yml`, `CATChem_species.yml`, `CATChem_emission.yml`) organize model settings hierarchically:

- **Simulation**: Name, species file path, emission mapping file path, verbosity.
- **Grid**: Vertical levels (`number_of_levels`), soil layers.
- **Timesteps**: Transport timestep, chemistry timestep.
- **Diagnostics**: Output frequency, output directory, output variable lists.
- **Processes**: Per-process activation (`activate`), scheme choices (`scheme`), and scheme parameters.
- **Species**: Molecular weights, gas/aerosol flags, dry/wet deposition parameters, particle radii, densities, and Mie parameters.
- **Emission Mappings**: Source categories, variable mapping names, regridding metadata, and scale factors.

## C++ Configuration Ownership

The C++ `ConfigManager` loads configuration files via:

```cpp
catchem::ConfigManager config_mgr;
config_mgr.load_from_file("CATChem_config.yml");
config_mgr.load_species_file("CATChem_species.yml");
config_mgr.load_emission_mapping_file("CATChem_emission.yml");
```

Host models and NUOPC query configuration data directly through C API bindings (`catchem_config_get_yaml_bool`, `catchem_config_get_output_frequency`, etc.) without Fortran configuration objects.
