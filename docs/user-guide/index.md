# User Guide

Welcome to the CATChem User Guide. This comprehensive guide covers everything you need to know to effectively use the CATChem library and modeling component for atmospheric chemistry modeling.

## Getting Started

- **[HPC Installation](hpc-installation.md)** - Learn how to install CATChem on HPC systems
- **[Build System](build-system.md)** - Learn how to build CATChem
- **[Standalone Driver](standalone-driver.md)** - Run CATChem as a self-contained ESMF/NUOPC application
- **[Overview](overview.md)** - Understanding CATChem's architecture and capabilities
- **[Configuration System](configuration.md)** - How to configure CATChem for your needs
- **[Input Files](input-files.md)** - Required input data and formats

## Process Documentation

CATChem uses a modular process-based architecture. Each atmospheric process is implemented as a separate module. Process descriptions are automatically generated based on the notes in the code files to ensure efficient and accurate documentation. For information on how to add a new process to CATChem see the **[Process Development](../developer-guide/processes/index.md)** section of the CATChem Developer Guide.

- **[Process Overview](../processes/index.md)** - Introduction to CATChem processes

### Transport Processes
- **[Settling](../processes/settling/settling.md)** - Gravitational settling processes

### Chemical Processes
- **[Carbon Chemistry](../processes/carbchem/carbchem.md)** - Carbonaceous aerosol chemistry
- **[Sulfate Chemistry](../processes/so4chem/so4chem.md)** - Sulfate aerosol chemistry

### Emission Processes
- **[Dust Emissions](../processes/dust/dust.md)** - Windblown dust emissions
- **[Sea Salt Emissions](../processes/seasalt/seasalt.md)** - Marine aerosol processes

### Loss Processes
- **[Dry Deposition](../processes/drydep/drydep.md)** - Surface dry deposition
- **[Wet Deposition](../processes/wetdep/wetdep.md)** - Precipitation scavenging

## Core Concepts

For a deep dive into the foundational infrastructure of CATChem, including state management, configuration, diagnostics, and error handling, please see the **[Core Concepts](../core-concepts/index.md)** section.

## Next Steps

- **Developers**: See the **[Developer Guide](../developer-guide/index.md)**
- **API Users**: Check the **[API Reference](../api/index.md)**
