# CATChem Modernized Processes Overview

This document provides an overview of the modernized physical and chemical processes in CATChem, highlighting the improvements in physics, performance, and native C++ software architecture.

---

## Modernized Processes

### 1. Gas-Phase Chemistry (GasChem)

**Location**: `src/process/gaschem/`  
**Type**: Atmospheric Transformation (chemical kinetics)  
**Status**: ✅ Completed

#### Key Improvements

*   **Native C++ MICM Solver**: Direct wrapping of NCAR's C++ Model Independent Chemistry Module (**MICM**) using the `musica` library, completely bypassing legay Fortran proxies.
*   **Automatic Zero-Overhead Photolysis Coupling**: Automatically scans rate parameters for `"PHOTO.<label>"` strings and binds them dynamically in-place to midpoint J-rates stored in the `DiagnosticManager` under `"photolysis_rate_<label>"`.
*   **Bidirectional Density Scaling**: Enforces exact Volume Mixing Ratio (VMR, ppmv) to molar density ($\text{mol/m}^3$) scaling using Sutherland-derived dry air molar density.
*   **Physical Safeguards**: Strict clamping floor of `1.0e-20` on concentration conversion loops to eliminate numerical singularities and NaN propagation.

#### Diagnostics

*   Dynamic photolysis diagnostic rates (`photolysis_rate_*`), species concentrations.

---

### 2. Photolysis

**Location**: `src/process/photolysis/`  
**Type**: Atmospheric Transformation / Radiation  
**Status**: ✅ Completed

#### Key Improvements

*   **Native C++ TUV-x Engine**: Direct execution of NCAR's Tropospheric Ultraviolet and Visible (**TUV-x**) solver on the CPU host.
*   **Column-Wise Input Extraction**: Slices independent 1D atmospheric columns (temperature, pressure, SZA, $O_3$ profile) using unmanaged Kokkos host Views.
*   **Edge-to-Midpoint Interpolation**: Automatically translates calculated edge-level J-rates to layer midpoints (cell centers) matching the chemical solver grid layout.

#### Diagnostics

*   `photolysis_rate_<rx_name>` (e.g. `photolysis_rate_jfoo`), solar_zenith_angle.

---

### 3. Gravitational Settling

**Location**: `src/process/settling/`  
**Type**: Transport (gravitational settling)  
**Status**: ✅ Completed

#### Key Improvements

*   **Advanced Stokes Scheme**:
    *   Temperature-dependent dynamic viscosity using Sutherland's law.
    *   Cunningham slip correction for small particles.
    *   Support for non-spherical particles via shape factors.
    *   CFL-stable subcycling for numerical stability.
*   **C++ Kokkos Parallelization**: Ported calculation loops to native C++ Kokkos parallel functors, facilitating execution on both multi-core CPUs and GPU devices.

#### Available Schemes

*   `StokesScheme`: Advanced Stokes settling with slip correction.
*   `IntermediateReynoldsScheme`: For larger particles (intermediate Reynolds numbers).

#### Diagnostics

*   `settling_velocity`, `settling_flux`, `cfl_number`.

---

### 4. YSU Vertical Dispersion

**Location**: `src/process/ysuverticaldispersion/`  
**Type**: Transport (vertical mixing)  
**Status**: ✅ Completed

#### Key Improvements

*   **Scale-Aware YSU Scheme**:
    *   Enhanced entrainment calculations.
    *   Grid-resolution dependent mixing coefficients.
    *   Improved turbulence parameterization.
*   **Architecture Enhancements**:
    *   Unified C++ Kokkos memory view mapping.
    *   Optimized column processing loop.

#### Diagnostics

*   `mixing_coefficients`, `entrainment_rate`, `boundary_layer_height`.

---

## Shared C++ Code Quality Standards

All modernized C++ processes adhere to strict software engineering standards:

1.  **Kokkos Device Portability**: Execution loops are written as parallel lambda kernels, ensuring compilation safety on CUDA, HIP, and OpenMP backends.
2.  **No Duplicate Allocations**: Raw pointers are fetched dynamically via C-API boundaries and wrapped as unmanaged Views, achieving zero-copy performance.
3.  **Comprehensive TDD Verification**: Supported by robust CTest unit and coupled integration tests asserting mathematical invariants and convergence safety.

---

## See Also

- **[GasChem Process Detailed Documentation](gaschem/index.md)**
- **[Photolysis Process Detailed Documentation](photolysis/index.md)**
- **[Developer Architecture Guide](../developer-guide/architecture.md)**
- **[Process Interface API Reference](../api/process-interface.md)**

---
