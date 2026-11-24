# Platform integrity migration

CATChem host integrations should use the checked C boundary declared in `catchem_api.hpp`. Checked calls return a stable `catchem_dataflow_status`; `catchem_get_last_error` adds operation-specific detail. Legacy pointer-returning calls remain compatibility wrappers, but cannot communicate the full failure category.

Opaque handles are typed, versioned, and live only until their owner is destroyed. A state handle obtained from a core becomes invalid when that core is destroyed. Hosts must not cache it across reinitialization. Each import transaction starts with `catchem_state_begin_import_generation`; every required meteorology, import, concentration, and diagnostic binding belongs to that generation.

Field bindings declare rank, extents, and semantic axes. Atmospheric layers, atmospheric interfaces, and soil layers are distinct. Read access synchronizes data but does not claim writer ownership; mutable access must use a write operation. Kokkos ownership is lease based: CATChem finalizes only a runtime it initialized, after the final CATChem lease is released. A host-owned Kokkos runtime remains host-owned.

Chemical indices must never be constants in an adapter. Resolve case-normalized species names against the active mechanism whenever the mechanism mapping is created. Mechanism changes therefore change cached indices without changing adapter source.

Configuration is validated before process construction and reports all independently discoverable issues with YAML paths. The old implicit 50-species fallback is removed. After a partial timestep failure, query the timestep outcome. `Reusable` permits retry, `RequiresReimport` requires a new complete import generation, and `RequiresReinitialize` requires rebuilding the model.

Build with `CATCHEM_BUILD_NUOPC` and `CATCHEM_BUILD_CCPP`. `BUILD_NUOPC` is deprecated and contradictory values are rejected. Target-scoped validation is enabled with `CATCHEM_ENABLE_SANITIZERS` and `CATCHEM_ENABLE_BOUNDS_CHECK`.
