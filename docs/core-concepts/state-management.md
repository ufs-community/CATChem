---
type: explanation
category: explanation
tags: [state-management, cpp-core, interop-field, zero-copy]
---

# State Management

This section covers state management concepts in CATChem (`catchem::StateManager`), providing zero-copy memory interop across C++ processes, host models (NUOPC), and Fortran science scheme bridges.

## Overview

The state management system in CATChem is centralized inside C++ `catchem::StateManager`. It manages 1D, 2D, and 3D memory fields using non-owning `InteropField` views, eliminating memory churn and pointer re-allocation during host model transfers.

Key features of `catchem::StateManager`:

- **Zero-Copy Host Binding**: Binds host meteorological and chemical arrays directly to non-owning `double*` views.
- **Dynamic Rebinding**: Updates host memory locations in $O(1)$ time without reallocating buffers.
- **On-the-Fly Meteorological Derivations**: Computes dry air density (`derive_airden_dry`) and layer thickness (`derive_bxheight`) when not directly supplied by the host driver.
- **Non-Negative Bounds Enforcement**: Clips small numerical solver underflows ($C < 0.0$) at the end of each timestep.

## Core Components

### `catchem::StateManager`
The central manager storing pointers and views:
- **Meteorology (`met`)**: Holds non-owning views for 3D fields (`T`, `PMID`, `PEDGE`, `AIRDEN`, `DELP`, `BXHEIGHT`, `CLDF`, `QV`, `RH`) and 2D fields (`PS`, `TS`, `LAT`, `LON`, `USTAR`, `HFLUX`, `PBLH`, `FROCEAN`, `FRSEAICE`).
- **Unified Chemistry (`chem`)**: Stores the unified 3D concentration view (`conc`) of shape `n_cols x n_levels x n_species` along with species metadata lists.

---

## Data Access Patterns

1. **Host Binding (`NUOPC` / Host Driver)**:
   Host drivers pass raw pointers to `bind_met_field_2d`, `bind_met_field_3d`, and `bind_unified_chemistry`.

2. **C++ Process Execution**:
   C++ processes extract raw pointers from `StateManager` and metadata from species definitions.

3. **Fortran Science Bridge Dispatch**:
   C++ processes pass raw pointers to flat BIND(C) Fortran science bridges, which execute vectorized numerical kernels.
