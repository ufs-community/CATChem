# Diagnostic System

This section describes the diagnostic system in CATChem, which allows processes to register and manage their own diagnostic outputs at runtime.

## Overview

The diagnostic system in CATChem is designed to be flexible and extensible. It allows processes to define their own diagnostic variables and to control how and when they are output. The system is built around two key components: the `DiagnosticInterface` and the `DiagnosticManager`.

## Core Components

### The Diagnostic Interface

The `DiagnosticInterface` is a module that defines the interface for the diagnostic system. It provides the following key data types:

- **`DiagnosticFieldType`**: This type is used to define a single diagnostic field. It contains metadata about the field, such as its name, description, and units, as well as the data itself.
- **`DiagnosticRegistryType`**: This type is used to manage a collection of diagnostic fields for a single process. It provides methods for registering, unregistering, and querying diagnostic fields.

### The Diagnostic Manager

The `DiagnosticManager` is a module that provides a central manager for the diagnostic system. It is responsible for:

- Managing the diagnostic registries for all processes.
- Collecting diagnostic data from all processes.
- Writing diagnostic data to output files.
- Providing a centralized way to configure the diagnostic system.

## Usage

To use the diagnostic system, a process must first create a `DiagnosticRegistryType` object and register it with the `DiagnosticManager`. Then, the process can create `DiagnosticFieldType` objects and register them with its diagnostic registry.

Once a diagnostic field is registered, the process can update its value at any time. The `DiagnosticManager` will then automatically collect the data and write it to an output file at the specified frequency.

## Features

The diagnostic system in CATChem provides a number of advanced features, including:

- **Multiple data types**: The system supports a wide range of data types, including scalars, 1D, 2D, and 3D arrays of real, integer, and logical values.
- **Flexible metadata**: The system allows processes to specify a wide range of metadata for each diagnostic field, including its name, description, units, and output frequency.
- **Process-specific diagnostics**: Each process has its own diagnostic registry, which allows it to manage its own diagnostic outputs independently of other processes.
- **Runtime query and collection**: The system allows users to query and collect diagnostic data at runtime, without having to restart the model.
- **Optional diagnostic output**: The system allows users to enable or disable diagnostic output for each process and for each diagnostic field.
# Diagnostic shape and writer rules

Diagnostic registration is idempotent only for an identical name, type, units, and shape. Incompatible re-registration fails without replacing live storage. Pointer retrieval validates rank and every extent; a destination that cannot represent all semantic axes must request an explicit selection instead of receiving an implicit first-species or first-level slice.

Host and device diagnostic writers declare the current side. Final timestep synchronization copies only from the latest writer, preserving diagnostics produced by host-side Fortran bridges as well as execution-space kernels.

## Registering with the axes + labels contract

Processes register diagnostic fields through `DiagnosticManager::register_field_contract`,
which pairs every dimension with a **semantic axis** so the NUOPC writer can decide *what
the variable means* instead of guessing from names or storage rank:

| `SemanticAxis` | Meaning | Written as |
|---|---|---|
| `Column` | Flattened local columns (`col = i + (j-1)*nx`), leading axis, always | grid `(x, y)` |
| `Singleton` | Extent-1 trailing axis (a per-column total) | 2D variable |
| `Level` | Vertical level axis | 3D variable, `lev` dimension kept intact |
| `Species` / `Category` | **Packed** axis: one slot per species/bin | one variable per slot, named `<field>_<label>` |
| `Interface`, `SoilLayer` | Reserved surface axes | not yet emitted |

A packed axis must carry one **unpack label** per slot — typically the species
`short_name` or bin label — supplied at registration and taken from the same
resolved species/bin list the scheme iterates, so output names always match the
data layout:

```cpp
state->diagnostic_manager()->register_field_contract(
    "dust_emission_bin", "Dust Emission Per Bin", "kg/m2/s",
    {ncol, nbins}, {SemanticAxis::Column, SemanticAxis::Category},
    dust_bin_labels);   // e.g. {"dust1", "dust3", ...}
```

The contract is validated strictly at registration (fail fast, fail loudly):

- at most **one** packed (`Species`/`Category`) axis per field;
- label count **equals** the packed extent — and labels are required exactly when
  the field has a packed axis;
- every label is NetCDF-name-safe (`[A-Za-z_][A-Za-z0-9_]*`) and unique;
- the leading axis is `Column` with extent equal to the run's column count.

A packed field whose slot has no label is a hard writer error naming the field —
never an anonymous `column_N` and never a silent drop. See
[NUOPC Integration](../developer-guide/integration/nuopc.md#process-diagnostic-output)
for how the contract maps to the output file.
