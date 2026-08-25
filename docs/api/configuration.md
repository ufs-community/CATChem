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
Configures physical grid dimensions, model timelines, and scheduled processes:

```yaml
# CATChem Simulation Configuration
model:
  name: "CATChem Simulation"
  timestep: 300.0  # seconds

grid:
  nx: 100
  ny: 100
  nz: 50

processes:
  - name: "photolysis"
    enabled: true
    parameters:
      config_file: "src/external/musica/configs/tuvx/tuv_5_4.yml"
  - name: "gaschem"
    enabled: true
    parameters:
      config_dir: "src/external/musica/configs/tuvx/from_host/"
```

### 2. Species Metadata Config (`CATChem_species.yml`)
Configures active chemical species, molecular weights, and properties:

```yaml
# Chemical Species Configuration
species:
  - name: "O3"
    molecular_weight: 0.0479982  # kg/mol
    description: "Ozone"
    category: "gas"

  - name: "NO2"
    molecular_weight: 0.0460055  # kg/mol
    description: "Nitrogen Dioxide"
    category: "gas"
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
