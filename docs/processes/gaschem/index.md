# Modernized C++ Gas-Phase Chemistry (GasChem) Process Documentation

## Overview

The modernized `GasChemProcess` is a native C++ physics/chemistry process that acts as a core atmospheric transformation node in CATChem. It wraps the C++ Model Independent Chemistry Module (**MICM**) solver using NCAR's **MUSICA** library, completely replacing the legacy Fortran proxy layers and double-buffering.

*   **Namespace / Class**: `catchem::GasChemProcess : public catchem::ProcessInterface`
*   **Header Location**: `src/process/gaschem/catchem_process_gaschem.hpp`
*   **Source Location**: `src/process/gaschem/catchem_process_gaschem.cpp`

---

## Architecture & How It Works

During a chemistry simulation timestep, `GasChemProcess::run` executes a highly-vectorized, multi-step pipeline on the CPU host:

1.  **Device-to-Host Synchronization**: Retrieves active meteorological and concentration fields from Kokkos device memories:
    ```cpp
    state->sync_to_host();
    ```
2.  **Column Batching (`columns_per_batch`)**: Rather than allocating a single `musica::State` across the entire MPI rank decomposition (`column_count * level_count`), `GasChemProcess::init` allocates one reusable `musica::State` of size `columns_per_batch * level_count` (default `columns_per_batch = 1`), and `GasChemProcess::run` loops over local columns in batches:
    $$i_{\text{cell}} = i_{\text{lev}} \times N_{\text{batch\_cols}} + i_{\text{batch\_col}}$$
    Any trailing unused columns in a partial final batch ($N_{\text{cols}} \bmod N_{\text{batch\_cols}} \neq 0$) are padded with a copy of the last valid column in the batch during the solve and ignored on copy-back.
3.  **Meteorology Guard**: Throws `std::domain_error` if temperature, pressure, or dry air density is non-finite or not positive in any cell. This check is not affected by `negative_concentration_policy`.
4.  **Pre-Resolved Species & Environmental State Mapping**:
    *   During `GasChemProcess::init`, CATChem species indices are mapped once (case-insensitively) to MICM variable indices so the cell loops perform zero string allocations or hash-map lookups.
    *   Converts dry air density (`AIRDEN_DRY` in $\text{kg/m}^3$) to dry air molar density ($\text{mol/m}^3$) using the molecular weight of dry air ($0.0289644\text{ kg/mol}$):
        $$\rho_{\text{molar}} = \frac{\rho_{\text{dry}}}{0.0289644}$$
    *   Converts gas concentrations from Volume Mixing Ratios (VMR, $\text{ppmv}$) to molar densities ($\text{mol/m}^3$) using dry air molar density:
        $$C_{\text{molar}} = C_{\text{ppmv}} \times 10^{-6} \times \rho_{\text{molar}}$$
5.  **Concentration Guards (`negative_concentration_policy`)**: Checks every mapped species on input and again on solver output. A non-finite value always throws. A negative value throws by default, or is floored to zero when the policy is `clamp`. See [Negative concentrations](#negative-concentrations-negative_concentration_policy).
6.  **Solver Execution**: Runs the C++ MICM Rosenbrock standard-ordered solver for each batch and reports the batch index and column range if convergence is not reached:
    ```cpp
    micm_instance_->Solve(micm_state_.get(), tstep);
    ```
7.  **Back-Conversion & Synchronization**: Converts output concentrations back to VMR ($\text{ppmv}$), updates the host state views, and marks host concentrations modified.

---

## Automatic Zero-Overhead Photolysis Coupling

A standout feature of the C++ core is the dynamic, zero-hardcoded coupling channel with the C++ Photolysis process:

1.  **Pre-Scan Rate Parameters in `init()`**: During initialization, `GasChemProcess::init` iterates through the solver's rate parameter map once.
2.  **Prefix Extraction**: If a rate parameter matches the pattern `"PHOTO.<label>"` (e.g., `"PHOTO.jfoo"`):
    *   It strips the prefix to extract the raw photolysis reaction identifier: `jfoo`.
    *   It caches the diagnostic search key: `"photolysis_rate_jfoo"` alongside its MICM parameter index.
3.  **Per-Timestep Pointer Binding**: At the start of `GasChemProcess::run` (outside the batch and cell loops), it queries `DiagnosticManager` once per photolysis parameter to resolve raw host pointers, then copies the calculated photolysis J-rates directly into MICM's rate parameter array inside the batch loop.
4.  **Loss Coefficient Defaults**: Any mechanism loss coefficients matching `"LOSS.<label>"` are initialized to `1.0` to preserve the chemistry mechanism's built-in scaling factors.

---

## Configuration

The process is dynamically added and configured via YAML layout files. The MICM
mechanism is supplied either as a **single-file MUSICA v1 configuration** (`config_file`)
or as a **v0 CAMP multi-file directory** (`config_dir`). An explicit `config_file`
always takes precedence; otherwise `config_dir` is probed for a `config.yaml` /
`config.yml` v1 entry point before falling back to the v0 directory parse.
`columns_per_batch` (optional, default `1`) controls how many vertical columns are
packed into the reusable `musica::State` per `Solve()` call.
`negative_concentration_policy` (optional, default `reject`) selects what happens to a
negative concentration; see the next section.

```yaml
processes:
  gaschem:
    activate: true
    columns_per_batch: 1
    negative_concentration_policy: reject   # reject | clamp
    # Option A (preferred): MUSICA v1 single-file mechanism.
    config_file: "src/external/musica/configs/v1/chapman/config.yaml"
    # Option B: v0 CAMP mechanism directory (config.json + species.json + reactions.json).
    # A directory that contains config.yaml/config.yml is auto-routed to the v1 parser.
    # config_dir: "src/external/musica/configs/v0/chapman/"
```

> **Note:** when migrating a v0 mechanism to v1, the photolysis reaction names must
> match the TUV-x radiator labels: `gaschem` binds each `"PHOTO.<label>"` rate
> parameter to the `photolysis_rate_<label>` diagnostic produced by the photolysis
> process (e.g. v0 `jO3->O1D` vs v1 `jo3_b` are *different* labels).

### Negative concentrations (`negative_concentration_policy`)

A host model can hand back mixing ratios that are slightly below zero. Transport and
mixing undershoot next to sharp gradients, for example just above a surface emission
on a near-zero background. `negative_concentration_policy` decides what `gaschem` does
with such a value, for the species that map into the MICM mechanism.

| Value | Negative input from the host | Negative output from the solver |
|---|---|---|
| `reject` (default) | throws `std::domain_error` | throws `std::domain_error` |
| `clamp` | set to zero before the solve | set to zero before it is written back |

Under both values:

*   A non-finite concentration (NaN or infinity) throws.
*   Non-physical temperature, pressure, or dry air density throws.
*   Any other value for the key fails `GasChemProcess::init` with `std::invalid_argument`.

The key is independent of [`physical_validation.policy`](../../user-guide/configuration.md#physical_validation),
which governs meteorology. `gaschem` does not read that policy, so a run can floor
tracer noise from the host and still stop on a bad temperature.

A rejected value is reported in scientific notation with its species and zero-based,
PET-local location:

```text
GasChem received a negative or non-finite concentration (-1.000000e-12) for O at column=4, level=2, batch=1 (zero-based indices)
```

Under `clamp`, each PET logs one warning the first time it floors a value, with the
number of cells and the most negative value and where it was. Later timesteps log the
same line at debug level, and `finalize` logs the totals since `init`:

```text
[WARN ] GasChemProcess: floored negative concentrations to zero | inputs_floored=1 min_input=-1.000000e-12 min_input_species=O min_input_column=4 min_input_level=2 outputs_floored=0 min_output=0.000000e+00 min_output_species=none min_output_column=-1 min_output_level=-1
[INFO ] GasChemProcess: negative concentrations floored to zero since init | inputs_floored=1 min_input=-1.000000e-12 outputs_floored=0 min_output=0.000000e+00
```

Flooring adds mass equal to the sum of the floored values. Check `min_input` after
enabling `clamp`: values near machine noise are expected, while a large negative
points to a problem upstream that the floor would otherwise hide. The policy in effect
is printed on the `GasChemProcess: initialized MICM state` log line.

---

## Verification & Testing Suite

We employ a strict Test-Driven Development (TDD) harness containing two robust verification targets:

### 1. Mathematical Invariant Unit Tests (`test_catchem_gaschem_units`)
Verifies mathematical conversion reversibility and clamping safety across a wide physical envelope:
*   Asserts that $(\text{ppmv} \rightarrow \text{mol/m}^3 \rightarrow \text{ppmv})$ is strictly identity-preserving ($< 1.0\times 10^{-15}$ relative tolerance).
*   Validates conversion safety across a broad range of dry air densities ($0.5$ to $1.2\text{ kg/m}^3$) and mixing ratios ($10^{-6}$ to $100.0\text{ ppmv}$).
*   Enforces concentration flooring of `1.0e-20` on negative out-of-bound inputs.

### 2. End-to-End Simulation Tests (`test_catchem_gaschem`)
Runs a coupled multi-level atmospheric column simulation containing both the C++ `photolysis` and `gaschem` processes:
1.  Registers both processes dynamically inside a C++ `Core` instance.
2.  Synthesizes physical conditions modeling a noon summer atmospheric column.
3.  Executes `core->run_timestep(3600.0)`.
4.  Asserts that photolysis rate diagnostics are successfully populated and mapped to MICM, and that solver concentrations successfully converge with stable mass balance.

The same target runs a six-column case across `columns_per_batch` values of 1, 2, 4,
and 6, and a `negative_concentration_policy` suite: `reject` throws with the value and
location in the message, `clamp` solves bit-for-bit like a zero input and counts each
floored cell once, and neither accepts NaN or a non-physical temperature.
