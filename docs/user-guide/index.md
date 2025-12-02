# User Guide

Welcome to the CATChem User Guide. This comprehensive guide covers everything you need to know to effectively use the CATChem library and modeling component for atmospheric chemistry modeling.

## Getting Started

- **[HPC Installation](hpc-installation.md)** - Learn how to install CATChem on HPC systems
- **[Build System](build-system.md)** - Learn how to build CATChem
- **[Overview](overview.md)** - Understanding CATChem's architecture and capabilities
- **[Configuration System](configuration.md)** - How to configure CATChem for your needs
- **[Input Files](input-files.md)** - Required input data and formats

## Process Documentation

CATChem uses a modular process-based architecture. Each atmospheric process is implemented as a separate module. Process descriptions are automatically generated based on the notes in the code files to ensure efficient and accurate documentation. For information on how to add a new process to CATChem see the **[Process Development](../developer-guide/processes/index.md)** section of the CATChem Developer Guide.

- **[Process Overview](../processes/index.md)** - Introduction to CATChem processes

### Transport Processes
- **[Settling](../processes/settling_process.md)** - Gravitational settling processes
- **[Plumerise](../processes/plume_rise/plume_rise.md)** - Plume rise process for fires, industrial stacks, and other elevated point sources

### Chemical Processes

More information coming soon!

### Emission Processes
- **[Biogenic Emissions](../processes/biogenic_emission/biogenic_emission.md)** - Biogenic emissions from vegetation
- **[Dust Emissions](../processes/dust/dust.md)** - Windblown dust emissions
- **[Sea Salt Emissions](../processes/seasalt/index.md)** - Marine aerosol processes
- **[External Emission Plume Rise](../processes/EXTERNAL_DATA_PLUME_RISE_SOLUTION.md)** - Complete solution for external emission data and plume rise
- **[External Emissions](../processes/external_emission_data/external_emission_data.md)** - Data management for external emissions in CATChem
     -  **[Integration Guide](../processes/external_emission_data/INTEGRATION_GUIDE.md)** - Integration Guide for external emissions

### Loss Processes

More information coming soon!

## Core Concepts

For a deep dive into the foundational infrastructure of CATChem, including state management, configuration, diagnostics, and error handling, please see the **[Core Concepts](../core-concepts/index.md)** section.

## Next Steps

- **Developers**: See the **[Developer Guide](../developer-guide/index.md)**
- **API Users**: Check the **[API Reference](../api/index.md)**
