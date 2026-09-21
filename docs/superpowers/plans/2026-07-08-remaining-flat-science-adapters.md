# Modernize Remaining Processes Flat-Science Adapters Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply the exact same "Direct C++ to Flat-Fortran Science Adapter" pattern used in `DryDep` to the remaining three legacy Fortran processes: `SeaSalt` (emissions), `WetDep` (wet deposition), and `SO4chem` (chemistry).

**Architecture:** For each process, we bypass the legacy Fortran `StateManager` and `VirtualColumn` containers by writing a flat C-linkable Fortran bridge (`<Process>ScienceBridge.F90`). This bridge accepts raw host pointers from C++ and passes standard Fortran array slices to the unmodified science routines inside `src/process/<process>/schemes/`. The C++ caller (`catchem_process_<process>.cpp`) extracts host view pointers from `catchem::StateManager` (and `catchem::DiagnosticManager`) and executes the bridge. The legacy Fortran wrappers (`Process<Name>Interface_Mod.F90` and `<Name>ProcessCreator_Mod.F90`) are deleted.

**Tech Stack:** C++20, Fortran 2008, Kokkos, ISO_C_BINDING, CMake.

## Global Constraints
- Target C++20 utilizing standard-conforming Kokkos namespaces and mdspan.
- The unported flat Fortran science files under `src/process/<process>/schemes/` MUST remain completely untouched.
- Memory layouts across language boundaries must remain aligned (LayoutLeft, column-major) for zero-copy CPU executions.
- Compilation and verification tests must be executed inside the `cece-dev:latest` Docker environment.

---

### Task 1: Modernize SeaSalt Process

**Files:**
- Create: `src/process/seasalt/SeaSaltScienceBridge.F90`
- Modify: `src/process/seasalt/catchem_process_seasalt.hpp` and `.cpp`
- Modify: `src/process/seasalt/CMakeLists.txt`
- Delete: `src/process/seasalt/ProcessSeaSaltInterface_Mod.F90`, `src/process/seasalt/SeaSaltProcessCreator_Mod.F90`

**Interfaces:**
- Consumes: C++ pointers to Met views (e.g. `FROCEAN`, `FRSEAICE`, `SST`, `U10M`, `V10M`), Chem views, and Diagnostics.
- Produces: `run_seasalt_science_bridge` C-linkable symbol.

- [ ] **Step 1: Write `SeaSaltScienceBridge.F90`**
Look at the signatures of `compute_gong97`, `compute_gong03`, and `compute_geos12` in the `schemes/` folder. Write a `BIND(C)` routine `run_seasalt_science_bridge` that maps C pointers (`c_f_pointer`) to standard Fortran arrays, loops over columns, and dynamically calls the specified `active_scheme`.

- [ ] **Step 2: Rewrite `catchem_process_seasalt.cpp/hpp`**
Update `SeaSaltProcess::init` to register dynamic diagnostics in `state->diag_mgr`. Update `SeaSaltProcess::run` to extract raw view pointers (including from `mdspan` fields if necessary via `.data_handle()`), call `run_seasalt_science_bridge`, and fence with `sync_to_host()` and `sync_to_device()`.

- [ ] **Step 3: Update `CMakeLists.txt` and Delete Legacy Code**
Remove `ProcessSeaSaltInterface_Mod.F90` and `SeaSaltProcessCreator_Mod.F90`. Add `SeaSaltScienceBridge.F90` to the build.

- [ ] **Step 4: Verify in Docker**
Compile `CATChem_process_seasalt` to ensure syntax correctness.
`docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make CATChem_process_seasalt"`

---

### Task 2: Modernize WetDep Process

**Files:**
- Create: `src/process/wetdep/WetDepScienceBridge.F90`
- Modify: `src/process/wetdep/catchem_process_wetdep.hpp` and `.cpp`
- Modify: `src/process/wetdep/CMakeLists.txt`
- Delete: `src/process/wetdep/ProcessWetDepInterface_Mod.F90`, `src/process/wetdep/WetDepProcessCreator_Mod.F90`

**Interfaces:**
- Consumes: C++ pointers to Met views (e.g. `AIRDEN_DRY`, `MAIRDEN`, `PEDGE`, `PFILSAN`, `PFLLSAN`, `REEVAPLS`, `T`), Chem views, and Diagnostics.
- Produces: `run_wetdep_science_bridge` C-linkable symbol.

- [ ] **Step 1: Write `WetDepScienceBridge.F90`**
Look at the signature of `compute_jacob` in the `schemes/` folder. Write a `BIND(C)` routine `run_wetdep_science_bridge` that maps C pointers to standard Fortran arrays, loops over columns, and calls `compute_jacob`.

- [ ] **Step 2: Rewrite `catchem_process_wetdep.cpp/hpp`**
Update `WetDepProcess::init` and `run` matching the direct adapter pattern.

- [ ] **Step 3: Update `CMakeLists.txt` and Delete Legacy Code**
Remove legacy wrappers. Add `WetDepScienceBridge.F90`.

- [ ] **Step 4: Verify in Docker**
Compile `CATChem_process_wetdep` to ensure syntax correctness.
`docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make CATChem_process_wetdep"`

---

### Task 3: Modernize SO4chem Process

**Files:**
- Create: `src/process/so4chem/SO4chemScienceBridge.F90`
- Modify: `src/process/so4chem/catchem_process_so4chem.hpp` and `.cpp`
- Modify: `src/process/so4chem/CMakeLists.txt`
- Delete: `src/process/so4chem/ProcessSO4chemInterface_Mod.F90`, `src/process/so4chem/SO4chemProcessCreator_Mod.F90`

**Interfaces:**
- Consumes: C++ pointers to met fields, persistent column states (from `VirtualColumn` logic previously), Chem views, and Diagnostics.
- Produces: `run_so4chem_science_bridge` C-linkable symbol.

- [ ] **Step 1: Write `SO4chemScienceBridge.F90`**
Look at the signature of `compute_gocart` inside `SO4chemScheme_GOCART_Mod.F90`. Note the persistent state variables per column (e.g. `firsttime`, `nymd_last`, `nhms_last_recycle`, `xh2o2_init`, `PSO4_from_SO2_per_level`, etc.). These need to be allocated in C++ as 2D/3D views or vectors, passed via raw pointer, and mapped via `c_f_pointer` so the Fortran bridge can slice them `(icol, ...)` across timesteps.

- [ ] **Step 2: Rewrite `catchem_process_so4chem.cpp/hpp`**
In `SO4chemProcess`, explicitly manage the persistent column state buffers (like `firsttime`, `xh2o2_init`) as standard C++ buffers or `Kokkos::View`s so they are maintained between `run()` calls, and pass them to the bridge.

- [ ] **Step 3: Update `CMakeLists.txt` and Delete Legacy Code**
Remove legacy wrappers. Add `SO4chemScienceBridge.F90`.

- [ ] **Step 4: Verify in Docker**
Compile `CATChem_process_so4chem` to ensure syntax correctness.
`docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make CATChem_process_so4chem"`

---

### Task 4: Complete Integration Verification

**Files:**
- Modify: `tests/test_catchem_interop.cpp`

- [ ] **Step 1: Append Integration tests for SeaSalt, WetDep, and SO4chem**
Inside `test_catchem_interop.cpp`, append TEST 10, TEST 11, and TEST 12 to verify execution of each adapter inside the C++ Core pipeline, just like `TEST 9` for `DryDep`.

- [ ] **Step 2: Build and run test_catchem_interop in Docker**
```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_interop && ./tests/test_catchem_interop"
```

- [ ] **Step 3: Build and run test_catchem_api in Docker**
```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build-test cece-dev:latest bash -c "cmake .. -DENABLE_KOKKOS=ON -DENABLE_TESTING=ON && make test_catchem_api && cp ../tests/CATChem_species.yml ./ && cp ../tests/CATChem_new_config.yml ./ && ./tests/test_catchem_api"
```
Expected: PASS
