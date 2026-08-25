# State Management

This section covers the core state management concepts in CATChem, providing unified data handling across all atmospheric processes.

## Overview

The state management system in CATChem is designed to provide a centralized and efficient way to manage all model data. It is built around a central data repository called the `StateManager` (`catchem::StateManager`), which holds all the different states of the model, such as the chemical state, meteorological state, time state, and diagnostic state.

The key design goals of the state management system are:

- **Unified Data Access**: Provide a single, consistent interface for all processes to access and modify model data.
- **Performance**: Enable efficient data access patterns, such as column-based processing via Kokkos subviews, to optimize performance.
- **Safety**: Ensure data integrity and provide robust error handling.
- **Flexibility**: Support a wide range of use cases, from simple standalone models to complex, fully coupled Earth system models.

## Core Components

### StateManager

The `StateManager` is the heart of the state management system. It is a C++ container holding Kokkos memory views and structured sub-states:

- **`ChemState`**: Manages chemical species concentrations (`conc`), metadata properties, and species indices.
- **`MetState`**: Manages meteorological fields, such as temperature, pressure, density, and surface properties.
- **`TimeState`**: Manages simulation timelines, Julian dates, solar zenith angles, and calendar utilities.
- **`DiagnosticManager`**: Manages diagnostic variables and registries for process outputs.

The `StateManager` is responsible for initializing sub-states, managing zero-copy pointer bindings, and synchronizing host (CPU) and device (GPU) memory.

### State Types

Each state type is a specialized data container designed to manage a specific type of data.

- **`ChemState`**: Manages chemical concentrations and species metadata attributes.
- **`MetState`**: Manages 2D and 3D meteorological arrays (`T`, `PMID`, `PEDGE`, `AIRDEN`, `USTAR`, etc.).
- **`DiagnosticManager`**: Manages diagnostic fields, registries, and output buffers across processes.

## Data Access Patterns

The state management system supports several data access patterns, each with its own trade-offs in terms of performance and ease of use.

- **Direct Access**: This is the simplest data access pattern, where data is accessed directly from the state containers. This pattern is suitable for simple operations, but it can be inefficient for complex calculations.
- **Column Processing**: This is a more efficient data access pattern, where data is processed one column at a time. This pattern is well-suited for atmospheric calculations, which are often performed on vertical columns of air.
- **Batch Operations**: This is the most efficient data access pattern, where data is processed in batches. This pattern is suitable for bulk data handling, such as initializing the model with data from a file.

## Error Handling

The state management system includes a comprehensive error handling mechanism that is designed to ensure data integrity and provide robust error handling. All state operations return a return code that indicates whether the operation was successful. If an error occurs, the return code will be non-zero, and an error message will be logged.

## Thread Safety

The state management system is designed to be thread-safe, which means that it can be used in parallel processing environments. The state containers provide methods for safely accessing and modifying data from multiple threads.

## Memory Management

The state management system includes an automatic memory management mechanism that is designed to minimize memory usage and prevent memory leaks. The state containers are responsible for allocating and deallocating memory as needed.
