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
    *   During `GasChemProcess::init`, every MICM variable is resolved once to a CATChem species through the active mechanism definition (`state->chemistry().mechanism`), so the cell loops perform zero string allocations or hash-map lookups. Resolution is case-insensitive and honours species `aliases`, so a mechanism that names a species differently from the host can be bridged from the species file:
        ```yaml
        - name: o3p
          aliases: [O]      # the MICM mechanism calls this species O
        ```
    *   Two MICM variables that resolve to the same CATChem species fail initialization.
    *   **Every species in the MICM mechanism must be in the CATChem species file.** A MICM variable with no CATChem counterpart fails initialization with a message naming it. MICM integrates every species in its mechanism, and the solver state is reused from one column batch to the next, so an unsupplied species would otherwise start from another column's result.
    *   Initialization logs one line listing CATChem gas species that the mechanism does not contain. Those are left untouched.
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

These four are the only options. Any other key under `gaschem`, including a misspelt
one such as `column_per_batch` or a nested block, stops the run at start-up with the
accepted names listed, as does an out-of-range value. The framework keys `activate`,
`diagnostics` and `diag_species` are not options and are always allowed. The check
runs before MICM is constructed: `GasChemProcess::parse_settings` is registered as the
process's settings validator, and `init()` calls the same function, so a process that
is initialized directly is checked as well. A legacy top-level `process:` block is
read the same way as `processes:`.

```text
Unknown option(s) for process 'gaschem': column_per_batch. Accepted options: columns_per_batch, config_dir, config_file, negative_concentration_policy.
```

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

### Species coverage

The species file must list every species of the MICM mechanism, by name or by alias.
Third-body species such as `M` are not MICM variables and need no entry. A species
that the host model should not transport still needs an entry, marked as not advected:

```yaml
- name: o1d
  __is_gas: true
  __is_advected: false     # solved by MICM, carried by CATChem, not transported by the host
  molecular weight [kg mol-1]: 0.01600
```

If any mechanism species is missing, initialization stops:

```text
GasChemProcess: MICM mechanism species missing from the CATChem species list: N2. Add them to the species file; use is_advected: false for species the host should not transport.
```

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

Under `clamp`, each PET logs up to two warnings: one the first time it floors an
input from the host, and one the first time it floors an output from the solver. The
two are tracked separately because a negative from the solver is a different problem
from noise in the host's tracers, and an earlier input warning must not hide it. The
`first_occurrence` field says which the line is for (`inputs`, `solver_outputs`, or
`inputs,solver_outputs` when both first happen on the same timestep). Each line gives
the number of cells and the most negative value and where it was. Later timesteps log
the same line at debug level with `first_occurrence=none`, and `finalize` logs the
totals since `init`:

```text
[WARN ] GasChemProcess: floored negative concentrations to zero | first_occurrence=inputs inputs_floored=1 min_input=-1.000000e-12 min_input_species=O min_input_column=4 min_input_level=2 outputs_floored=0 min_output=0.000000e+00 min_output_species=none min_output_column=-1 min_output_level=-1
[INFO ] GasChemProcess: negative concentrations floored to zero since init | inputs_floored=1 min_input=-1.000000e-12 outputs_floored=0 min_output=0.000000e+00
```

The counters and minima are updated on every timestep whatever the log level. The
per-timestep line is only formatted when it will be written.

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
and 6 and requires every batch size to match the single-column result within ten times
the solver's own tolerances (its relative tolerance, plus the mechanism's declared
absolute tolerance converted to ppmv). It also runs a `negative_concentration_policy` suite: `reject` throws with the value and
location in the message, `clamp` solves bit-for-bit like a zero input and counts each
floored cell once, and neither accepts NaN or a non-physical temperature. A further
suite runs three timesteps and checks that floored inputs and floored solver outputs
are each warned about once, that the counters do not depend on the log level, and that
`reject` stops on a negative from the solver. The Chapman mechanism never returns a
negative, so these cases use a `GasChemProcess` subclass that overrides the protected
`adjust_solver_output_for_testing` hook, which does nothing in the process itself.
A settings suite checks the option parser through the registered validator, through a
directly initialized process, through a legacy `process:` block, and through a `Core`
built from a runtime YAML.
