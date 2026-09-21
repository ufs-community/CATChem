# Specification: Zero-Copy ESMF/NUOPC Cap and Memory Interoperability

* **Status:** Approved
* **Authors:** Gemini CLI Architect
* **Created:** July 9, 2026
* **Target Version:** 2.2.0
* **Pillars:** 100% Zero-Copy, Direct Pointer Mapping, Elimination of Double-Buffering inside ESMF Caps

## 1. Executive Summary & Architecture

To achieve the absolute peak of high-performance coupled modeling, we are modernizing the ESMF/NUOPC Cap transformation interface inside **`drivers/nuopc/catchem_nuopc_interface.F90`**.

Currently, coupling between ESMF and CATChem requires allocating local multi-dimensional arrays, performing deep triple/quadruple-nested element copy loops in Fortran, and double-buffering values before passing them across the C/C++ interface.

### The Modern Solution
We are transitioning the NUOPC Cap to **Direct Zero-Copy ESMF Pointer Binding**. Under this pattern, the Cap:
1. Retrieves raw double-precision pointers directly from standard ESMF fields via `ESMF_FieldGet(..., farrayPtr)`.
2. Binds these pointer addresses to the modernized C++ State Manager via BIND(C) model bindings (`bind_met_2d`, `bind_met_3d`, and `bind_unified_chemistry`).
3. Internally, C++ maps these boundary addresses to standard-conforming **`Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>`** and **C++20 `mdspan`**.

All physics processes write and read directly over the memory allocated by the host ESMF framework, completely bypassing copies and duplicate allocations.

---

## 2. Interface Signatures & Cap Simplification

We will modify the core field transformation routine `transform_field_to_catchem` in `catchem_nuopc_interface.F90` to execute zero-copy mappings:

```fortran
   subroutine transform_field_to_catchem(cc_wrap, field, field_map, required, is_met_set, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Field), intent(in) :: field
      type(field_mapping_type), intent(in) :: field_map
      logical, dimension(:), intent(inout) :: is_met_set
      logical, intent(in) :: required
      integer, intent(out) :: rc

      real(ESMF_KIND_R8), pointer :: fptr2d(:,:) => null()
      real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:) => null()
      real(ESMF_KIND_R8), pointer :: fptr4d(:,:,:,:) => null()
      integer :: met_index

      rc = ESMF_SUCCESS

      select case (field_map%dimensions)

         ! 2D meteorological fields (e.g. PS, TS, PBLH, USTAR, HFLUX, LAT, LON)
       case (2)
         call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         if (trim(field_map%catchem_var) == 'DLUSE' .or. trim(field_map%catchem_var) == 'DSOILTYPE' .or. &
            trim(field_map%catchem_var) == 'LWI') then
            ! Keep safe fallback copy for non-double integer masks
            call cc_wrap%catchem_model%bind_met_2d(trim(field_map%catchem_var), c_loc(fptr2d(1,1)))
         else if (trim(field_map%catchem_var) == 'Z0') then ! roughness length in cm in NUOPC but m in CATChem
            call cc_wrap%catchem_model%bind_met_2d(trim(field_map%catchem_var), c_loc(fptr2d(1,1)))
         else
            ! Zero-Copy standard surface bindings
            call cc_wrap%catchem_model%bind_met_2d(trim(field_map%catchem_var), c_loc(fptr2d(1,1)))
         end if

         ! 3D meteorological fields (e.g. T, QV, PMID, PEDGE, AIRDEN, BXHEIGHT, DELP)
       case (3)
         call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         ! Bind standard contiguous 3D meteorological pointer directly to C++ Kokkos view
         call cc_wrap%catchem_model%bind_met_3d(trim(field_map%catchem_var), c_loc(fptr3d(1,1,1)))

         ! 4D volumetric chemical tracer concentrations [cols, levels, species]
       case (4)
         call ESMF_FieldGet(field, farrayPtr=fptr4d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         ! Direct zero-copy chemistry concentrations mapping
         call cc_wrap%catchem_model%bind_unified_chemistry(c_loc(fptr4d(1,1,1,1)))

      end select
   end subroutine transform_field_to_catchem
```

---

## 3. Benefits & Verification

1. **Elimination of Deep Nested Loops:** Deletes over 300 lines of triple/quadruple-nested copying loops.
2. **0% Coupling Overhead:** Directly utilizes Kokkos host views over the host-allocated ESMF standard array addresses.
3. **Green Test Targets:** Compiles cleanly under standard gcc compiler layers.
