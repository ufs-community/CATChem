# Configuration API

This section covers the modernized configuration management APIs that handle YAML-based setup and runtime configuration in CATChem.

## Overview

The configuration system has been fully elevated to the C++ core to unify runtime configuration and species metadata loading:

- **catchem::ConfigManager**: C++ class managing YAML-based settings using the `yaml-cpp` library.
- **Unified Metadata Parsing**: Species configurations (`CATChem_species.yml`) parsed natively in C++ to populate species attributes (e.g. molecular weights).
- **C-API Interface**: Exposes BIND(C) parameters so host models can configure dimensions and timestep thresholds seamlessly.

---

## Core Components

### catchem::ConfigManager

The C++ class responsible for parsing the main configuration file:

```cpp
#pragma once
#include <string>
#include <string_view>
#include <yaml-cpp/yaml.h>

namespace catchem {

    class ConfigManager {
    public:
        ConfigData data;
        YAML::Node root_node;
        bool is_loaded = false;
        std::string config_file_path;

        ConfigManager() = default;
        void load_from_file(const std::string& filename);

        YAML::Node get_process_config(std::string_view process_name) const;
    };
```
} // namespace catchem
```

---

## Configuration Files Format

### 1. Main Simulation Config (`CATChem_config.yml`)
Configures the run timeline, referenced species/emission files, and scheduled
processes (see `configs/Default/CATChem_config.yml` for a production example):

```yaml
simulation:
  name: test
  start_date: 20240501 0000
  end_date: 20240501 0100
  species_filename: ./Configs/Default/CATChem_species.yml
  emission_filename: ./Configs/Default/CATChem_emission.yml

grid:
  nx: 100
  ny: 100
  nz: 50
```

### 2. Species Metadata Config (`CATChem_species.yml`)
Configures active chemical species, molecular weights, and physics flags
(see `configs/Default/CATChem_species.yml`):

```yaml
- name: so2
  __description: Sulfur dioxide
  __is_gas: true
  __is_drydep: true
  __is_wetdep: true

- name: dust1
  __description: Dust bin 1
  __is_dust: true
  __radius: 5.0e-07
  __density: 2500.0
```

---

## Parameter Access Patterns

### 1. Native C++ Parameter Retrieval
Developers fetch parameters directly inside process initialization blocks using template-based accessors:

```cpp
void DryDepProcess::init(std::shared_ptr<StateManager> state) {
    auto config_mgr = state->get_config_manager();

    // Type-safe parameter retrieval with default value
    double scale_factor = config_mgr->get_parameter<double>("processes.dry_dep.scale_factor", 1.0);
    bool check_mass = config_mgr->get_parameter<bool>("processes.dry_dep.check_mass_conservation", false);
}
```

### 2. Species Metadata Lookup
Active species metadata is pre-loaded during StateManager initialization, allowing C++ physics kernels to query properties efficiently:

```cpp
// Query molecular weight directly from StateManager
double o3_mw = state->species_metadata["O3"].molecular_weight;
```

---

## Best Practices

### Performance
1.  **Cache in `init`**: Never fetch config parameters inside the `run()` timestep loops. Read and cache settings to local variables during `init()`.
2.  **Avoid Redundant Parsing**: Let the C++ `ConfigManager` act as the single parsing point. Upstream models can retrieve dimensions through the flat BIND(C) API `catchem_get_grid_dimensions` rather than re-parsing files.

### Safety
1.  **Define Defaults**: Always provide sensible default values to `get_parameter` calls to prevent execution crashes on missing parameters.
2.  **Shield Exceptions**: Ensure YAML parsing errors are shielded in BIND(C) boundaries. Catch `YAML::Exception` during load phases.

---

## See Also

- [State Management API](state-management.md) - Chemical Species Views setup
- [Process Interface API](process-interface.md) - YAML settings inside processes
- [Column Interface API](column-interface.md) - Grid Manager layout

---
