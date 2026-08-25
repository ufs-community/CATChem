# Dust Examples

This document provides usage examples for the Dust process.

## Basic Configuration

```yaml
processes:
  - name: dust
    enabled: true
    scheme: Fengsha
    species: [O3, NO2, SO2]
    diagnostics: [process_rate, tendency]
```

## C++ Usage

```cpp
#include <catchem_core.hpp>

extern "C" void catchem_register_dust_cpp();

catchem_register_dust_cpp();
auto core = std::make_shared<catchem::Core>("tests/Configs/Default/CATChem_config.yml");
auto process = catchem::ProcessRegistry::get_instance().create("dust");
process->init(core->get_state_manager());
core->add_process(process);

// Execute timestep
core->run_timestep(1800.0);
```
