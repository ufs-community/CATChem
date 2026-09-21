# API Reference

Welcome to the CATChem API documentation. This section provides comprehensive reference documentation for all CATChem modules, types, and procedures.

!!! note "Documentation Structure"
    Our API documentation combines hand-written guides with auto-generated reference material from the source code using [MkDoxy](https://mkdoxy.kubaandrysek.cz/). This ensures both comprehensive coverage and up-to-date accuracy.

## � API Organization

### Core Modules
High-level APIs for the main CATChem systems:

- **[State Management](state-management.md)** - StateContainer, ChemState, MetState, & DiagState data handling
- **[Process Interface](process-interface.md)** - Process development and integration APIs
- **[Column Interface](column-interface.md)** - Column virtualization and 1D processing
- **[Configuration Manager](configuration.md)** - YAML configuration system

## 🔗 Auto-Generated Documentation

**[→ Complete Auto-Generated API Reference](../CATChem/group__catchem.md)**

The complete API documentation includes:

- **[Modules Index](../CATChem/modules.md)** - List of all modules
- **[Namespaces Index](../CATChem/namespaces.md)** - List of all namespaces
- **[Files Index](../CATChem/files.md)** - List of all files
- **[Functions Index](../CATChem/functions.md)** - List of all functions

### Key Auto-Generated Sections

**Core System Modules and Processes:**

- **[CATChem Directory](../CATChem/group__catchem.md)** - Main CATChem Directory
- **[Core CATChem API](../CATChem/group__catchem__api.md)** - Core CATChem API functions and data types
- **[Core Modules](../CATChem/group__core__modules.md)** - Core modules and data types for CATChem
- **[Processes](../CATChem/group__process__modules.md)** - All atmospheric chemistry processes

**Utilities:**

- **[Constants](../CATChem/constants_8_f90.md)** - Physical and mathematical constants
- **[Utilities](../CATChem/utilities__mod_8_f90.md)** - Common utility functions and tools


## Quick Reference

### Key Types

| Type | Module / Namespace | Description |
|------|---------------------|-------------|
| `catchem::Core` | C++ Namespace `catchem` | Central orchestration engine |
| `catchem::StateManager` | C++ Namespace `catchem` | Central memory state (using Kokkos Views) |
| `catchem::ProcessInterface` | C++ Namespace `catchem` | Base virtual interface class for physics/chemistry processes |
| `catchem::ProcessRegistry` | C++ Namespace `catchem` | Creator-lambda process registry |
| `catchem::CoreCreateOptions` | C++ Namespace `catchem` | Unified direct, configured, or host-grid Core construction inputs |
| `CATChemType` | Fortran `CATChemAPI_Mod` | BIND(C) wrapper API delegate |
| `StateContainerType` | Fortran `state_mod` | Fortran delegate container wrapping C++ StateManager |

### Common Patterns

=== "Process Initialization (C++)"

    ```cpp
    #include <catchem_process_registry.hpp>

    // Register implementations at link/load time; YAML activation determines
    // which registered processes are instantiated for this Core.
    auto core = catchem::Core("CATChem_config.yml");
    ```

=== "Process Initialization (Fortran)"

    ```fortran
    ! Modern Fortran delegates to C++ ProcessRegistry
    use ProcessName_Mod
    type(ProcessNameType) :: process
    call process%init(container, rc)
    ```

=== "Diagnostic Access (C++)"

    ```cpp
    #include <catchem_diagnostic_manager.hpp>

    // Query and write diagnostic midpoint photolysis rate
    auto diag_mgr = core->get_diagnostic_manager();
    double* jrate_ptr = diag_mgr->get_field_pointer("photolysis_rate_jfoo");
    jrate_ptr[cell_idx] = calculated_jrate;
    ```

=== "Diagnostic Access (Fortran)"

    ```fortran
    ! Direct bind pointers to C++ DiagnosticManager buffers
    use DiagnosticInterface_Mod
    type(DiagnosticFieldType), pointer :: field
    field => diag_mgr%get_field('field_name', rc)
    call field%get_data(data_array, rc)
    ```

## Search Tips

- Use the search box above to find specific C++ namespaces, classes, procedures, or Fortran wrapper types
- Browse by module or namespace for related functionality
- Check the inheritance hierarchy for C++ `ProcessInterface` subclasses
- Look at usage examples in the `tests/` directory (e.g. `tests/test_catchem_gaschem.cpp`)

## Conventions

### Naming Conventions
- **C++ Namespaces**: `catchem`
- **C++ Classes**: `camelCase` (with leading upper letter, e.g. `StateManager`)
- **C++ Methods**: `snake_case` (e.g. `run_timestep`, `sync_to_device`)
- **Fortran Modules**: `ModuleName_Mod`
- **Fortran Types**: `TypeNameType`
- **Fortran Procedures**: `snake_case`
- **Constants**: `UPPER_CASE`

### Return Codes
All C-bound procedures use integer return codes following the convention:
- `CC_SUCCESS = 0` - Successful operation
- `CC_FAILURE = -1` - Generic failure
- C++ methods utilize standard exceptions (e.g., `std::runtime_error`) shielded inside BIND(C) boundaries.

### Memory Management
- Central memory managed entirely in C++ `StateManager` using Kokkos Views
- Dual-space capabilities: synchronized dynamically between Host (CPU) and Device (GPU) memory layout
- Fortran pointer variables are dynamically bound at runtime to raw C++ host pointers without any duplicate allocations.

## Contributing

Found an issue with the documentation? The API docs are generated automatically, so:

1. **Source Code Issues**: Update the source code comments and docstrings
2. **Organization Issues**: Modify the MkDoxy configuration in `mkdocs.yml`
3. **Missing Documentation**: Add Doxygen-style comments to the source code

For details on documentation standards, see the [developer guide section on documentation](../developer-guide/documentation.md).
