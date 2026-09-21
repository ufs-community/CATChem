# Specification: Pure Proxy Core and Zero-Copy Coupled Exports

* **Status:** Draft (Awaiting Review)
* **Authors:** Gemini CLI Architect
* **Created:** July 9, 2026
* **Target Version:** 2.2.0
* **Pillars:** Perfect Stateless Wrappers, Zero-Leak Memory Architecture, Direct Export Buffer Alignment

---

## 1. Executive Summary & Core Bottlenecks

Through rigorous architectural trace analysis, we have uncovered a massive opportunity to simplify the CATChem core while addressing a colossal, silent memory leak in standard coupled executions.

### The Leaky Duplicate Core
Currently, when `metstate_mod.F90` or `chemstate_mod.F90` is initialized, standard Fortran invokes custom heap allocation routines (`allocate_metstate_arrays` and `allocate(unified_conc)`), allocating millions of elements on the heap.
However, during startup, `StateManager_Mod.F90` immediately overrides these pointer mappings with standard C++ pointers fetched from C++ StateManager views:
```fortran
call get_cpp_field(this%cpp_ptr, "T", ptr%T, [nx, ny, nz], rc)
```
In standard Fortran 2003/2008, re-associating an active pointer to another address via `c_f_pointer` or pointer assignment silently orphans the previously allocated heap segment, causing **massive memory leaks on every execution block**.

### The Solution: The Pure Proxy Core
We will transition both `MetStateType` and `ChemStateType` to **Perfect Stateless Proxies**.
1. **0% Fortran Allocations:** Fortran will allocate exactly zero bytes of local memory for physical variables. All pointer-member arrays remain null until bound dynamically to C++ host views.
2. **Zero-Copy Coupled Exports:** Because the ESMF/NUOPC Cap now binds the framework's output fields directly to the C++ core StateManager and chemical concentration views, we completely bypass the need to copy tracer outputs or final tendencies back and forth at the end of coupling steps. C++ writes outputs **directly to ESMF/CCPP export buffers in-place**.

---

## 2. Stateless Core Proxy Architecture

### MetStateType Simplified Proxy (`src/core/metstate_mod.F90`)
We will remove `allocate_metstate_arrays` and all case-by-case array allocations. The initialization step will simply associate the Opacity pointer handle:

```fortran
subroutine metstate_init(this, nx, ny, nlevs, nsoil, nsoiltype, nsurftype, error_mgr, rc)
   class(MetStateType), intent(inout) :: this
   integer, intent(in) :: nx, ny, nlevs
   ! Decoupled: Zero heap allocations occur here!
   call this%geometry%set(nx, ny, nlevs)
   this%NLEVS = nlevs
   this%State = 'MET'
   rc = CC_SUCCESS
end subroutine metstate_init
```

---

### ChemStateType Simplified Proxy (`src/core/chemstate_mod.F90`)
Similarly, we will purge local allocations for species descriptors inside the compatibility shell, mapping pointers once on synchronization:

```fortran
subroutine chemstate_init(this, max_species, error_mgr, rc, grid)
   class(ChemStateType), intent(inout) :: this
   type(GridGeometryType), pointer, optional, intent(in) :: grid

   this%State = 'Chem'
   this%nSpecies = 0
   if (present(grid)) this%Grid => grid
   rc = CC_SUCCESS
end subroutine chemstate_init
```

---

## 3. Zero-Copy Coupled Exports Specification

Since `chem_state%ChemSpecies(s)%conc` slices point directly to C++ unified chemistry layout views, which the ESMF Cap binds directly to ESMF standard array addresses at startup, **no copying of final concentration outputs is required**.

Inside `transform_catchem_to_field` inside `catchem_nuopc_interface.F90`:
* Standard trace-gas concentrations (`fptr4d`) are **automatically updated in-place** inside ESMF buffers.
* Diagnostic PM2.5 and PM10 categories are populated cleanly by fetching double precision references from C++ State / Diagnostic Managers once without allocating intermediate local arrays in Fortran.

---

## 4. Verification and Downstream Stability

All standard compiled targets and low-level unit tests (`test_MetState`, `test_TimeState`, `test_catchem_properties`) will be verified within the active GCC build container to confirm exact numerical consistency, leak-free memory footprints, and 100% green status.
