# CATChem Processes

This section contains documentation for all available CATChem physical and chemical processes.

## Overview

CATChem processes are modular components that implement specific atmospheric physics or chemistry mechanisms. Each process follows a standardized interface and can be configured independently.

Under the modernized C++ core, processes are implemented as native C++ classes extending `catchem::ProcessInterface`, allowing them to execute in parallel on CPU host and GPU devices via Kokkos.

---

## Process Categories

### Chemistry & Transformation Processes
Processes that convert species through chemical or physical mechanisms:
- **[GasChem](gaschem/index.md)** - Gas-phase chemistry utilizing NCAR's C++ MICM solver via the musica library.
- **[Photolysis](photolysis/index.md)** - Photolysis rates determination utilizing NCAR's C++ TUV-x engine.
- **Particle settling** - Gravitational settling with Cunningham slip correction and Sutherland viscosity.
- Aerosol chemistry, coagulation, and phase transitions.

### Emission Processes
Source processes that add species to the atmosphere:
- **[Dust](dust/dust.md)** - Mineral dust emission and transport.
- **[SeaSalt](seasalt/seasalt.md)** - Marine aerosol processes.
- Anthropogenic and biogenic emissions.

### Loss Processes
Removal processes that remove species from the atmosphere:
- Dry deposition
- Wet deposition
- Radioactive decay

---

## Available Processes

**Emission**

- **[Dust](dust/dust.md)** - Windblown mineral dust emission (FENGSHA, GINOUX)
- **[Sea Salt](seasalt/seasalt.md)** - Marine aerosol emission (GEOS12, GONG03, GONG97)

**Chemistry**

- **[Carbon Chemistry](carbchem/carbchem.md)** - Carbonaceous aerosol chemistry (GOCART)
- **[Sulfate Chemistry](so4chem/so4chem.md)** - Sulfate aerosol chemistry (GOCART)
- **[Gas-phase Chemistry](gaschem/index.md)** - Gas-phase chemistry (MICM via MUSICA)
- **[Photolysis](photolysis/index.md)** - Photolysis rates (TUV-x via MUSICA)

**Deposition**

- **[Dry Deposition](drydep/drydep.md)** - Surface dry deposition (GOCART, WESELY, ZHANG)
- **[Wet Deposition](wetdep/wetdep.md)** - Precipitation scavenging (JACOB)

**Transport**

- **[Settling](settling/settling.md)** - Gravitational settling (GOCART)

---

## Using Processes

### Configuration

Set `activate: true` under the process name in `processes` and list the process
in `run_phases` where it should execute. See the
[configuration reference](../user-guide/configuration.md) for the YAML layout and
the [photolysis guide](photolysis/index.md#configure-photolysis-and-gas-phase-chemistry)
for a coupled TUV-x and MICM schedule.

### Process Interface (C++)

All native processes extend the standard C++ interface:

```cpp
namespace catchem {

    class ProcessInterface {
    public:
        virtual ~ProcessInterface() = default;
        virtual std::string get_name() const = 0;
        virtual void init(std::shared_ptr<StateManager> state) = 0;
        virtual void run(std::shared_ptr<StateManager> state) = 0;
        virtual void finalize() = 0;
    };

} // namespace catchem
```

---

## Process Development & Architecture Guides

For information on developing and integrating processes:

- **[Developer Architecture Guide](../developer-guide/architecture.md)** - C++ Core and memory layout.
- **[Process Interface API Reference](../api/process-interface.md)** - Details on ProcessRegistry and linker callbacks.
- **[Process Generator Tutorial](../developer-guide/processes/process-generator.md)** - Creating new processes with the generator.

---
