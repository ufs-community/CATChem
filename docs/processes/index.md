# CATChem Processes

This section contains documentation for all available CATChem processes.

## Overview

CATChem processes are modular components that implement specific atmospheric physics or chemistry mechanisms. Each process follows a standardized interface and can be configured independently.

## Process Categories

### Emission Processes
Source processes that add species to the atmosphere:
- Anthropogenic emissions
- Biogenic emissions
- Dust emissions
- Sea salt emissions

### Transformation Processes
Processes that convert species through chemical or physical mechanisms:
- Gas-phase chemistry
- Aerosol chemistry
- Particle settling
- Coagulation

### Loss Processes
Removal processes that remove species from the atmosphere:
- Dry deposition
- Radioactive decay
- Photolysis

### Transport Processes
Processes that move species spatially:
- Advection
- Turbulent mixing
- Convective transport

### Multi-Phase Processes
Processes involving multiple atmospheric phases:
- Aqueous chemistry (gas-liquid)
- Heterogeneous chemistry (gas-solid)
- Phase transitions

## Available Processes

- **[SeaSalt](seasalt/index.md)** - SeaSalt atmospheric process
- **[Dust](dust/index.md)** - Dust atmospheric process
- **[Dust](dust/index.md)** - Emission process (sources)
- **[CleanTest](cleantest/index.md)** - Emission process (sources)
- **[SeaSalt](seasalt/index.md)** - Emission process (sources)
- **[TestProcess](testprocess/index.md)** - Emission process (sources)
*This section will be automatically updated as new processes are generated.*

## Using Processes

### Basic Configuration

All processes are configured through YAML files:

```yaml
processes:
  - name: process_name
    enabled: true
    scheme: scheme_name
    species: [species_list]
    parameters:
      param1: value1
      param2: value2
    diagnostics: [diagnostic_list]
```

### Process Interface

All processes implement the standard interface:

```fortran
type, extends(ProcessInterface) :: MyProcessType
contains
  procedure :: init => my_process_init
  procedure :: run => my_process_run
  procedure :: finalize => my_process_finalize
end type
```

## Creating New Processes

See the [Process Generator Tutorial](../developer-guide/processes/generator-tutorial.md) for detailed instructions on creating new processes.

## Process Development

For information on developing processes manually or extending generated processes:

- [Process Architecture](../developer-guide/processes/architecture.md)
- [Creating Processes](../developer-guide/processes/creating.md)
- [Testing Guide](../developer-guide/processes/testing.md)
