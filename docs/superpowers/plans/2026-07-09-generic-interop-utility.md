# Generic Fortran-C-C++ Interop Pointer Utility Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eradicate all generated `.inc` macros and thousands of lines of select-case Fortran boilerplate by implementing a dynamic, generic runtime pointer-association utility (`Interop_Mod.F90`).

**Architecture:** We use Fortran 2018 standard ISO_C_BINDING to map raw C/C++ heap addresses directly to Fortran array pointer slices at runtime via string-name lookup. All metadata and physics memory allocations remain strictly inside the C++ Core orchestrator, maintaining a single source of truth while keeping the unported ESMF/NUOPC Caps compile-safe and binary-compatible.

**Tech Stack:** C++20, Kokkos, Fortran 2018 (iso_c_binding).

## Global Constraints
- Language Target: Target C++20 utilizing Kokkos host views and C++20 standard-conforming mdspan backports, avoiding C++23 direct `<mdspan>` dependencies.
- Layout Alignment: Retain Fortran column-major storage layout (Kokkos::LayoutLeft) for zero-copy CPU executions.
- Single Source of Truth: Core C++ handles 100% of orchestration, memory management, configuration loading, and diagnostics, completely bypassing legacy Fortran orchestration.
- Flat-Science Adapter Pattern: Cleanly map C++ unmanaged double views to Fortran array slices inside flat BIND(C) bridges (*ScienceBridge.F90) via standard c_f_pointer, avoiding duplicate physics code.
- Language Boundary Exception Checks: Wrap all BIND(C) export endpoints in robust C++ try-catch blocks to prevent escaping C++ exceptions.
- Docker Build Environment: All compilations and tests MUST be executed inside the `cece-dev:latest` Docker container, explicitly installing `python3` via standard `apt-get` updates.

---

### Task 1: Create the Generic Interoperability Module (`Interop_Mod.F90`)

**Files:**
- Create: `src/core/Interop_Mod.F90`

**Interfaces:**
- Consumes: None (raw standard ISO_C_BINDING)
- Produces: `get_cpp_field` generic module procedure with signatures:
  - `subroutine get_field_1d(cpp_ptr, name, f_ptr, dims, rc)`
  - `subroutine get_field_2d(cpp_ptr, name, f_ptr, dims, rc)`
  - `subroutine get_field_3d(cpp_ptr, name, f_ptr, dims, rc)`

- [ ] **Step 1: Write `src/core/Interop_Mod.F90`**

Write the complete generic pointer interop code:
```fortran
module Interop_Mod
   use iso_c_binding
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   use Precision_Mod, only: fp

   implicit none
   private

   public :: get_cpp_field

   interface
      type(c_ptr) function catchem_state_get_pointer_1d(state_ptr, name) bind(C, name="catchem_state_get_pointer_1d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function

      type(c_ptr) function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function

      type(c_ptr) function catchem_state_get_pointer_3d(state_ptr, name) bind(C, name="catchem_state_get_pointer_3d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function
   end interface

   interface get_cpp_field
      module procedure get_field_1d
      module procedure get_field_2d
      module procedure get_field_3d
   end interface

contains

   subroutine get_field_1d(cpp_ptr, name, f_ptr, dims, rc)
      type(c_ptr), intent(in) :: cpp_ptr
      character(len=*), intent(in) :: name
      real(fp), pointer, intent(out) :: f_ptr(:)
      integer, intent(in) :: dims(1)
      integer, intent(out) :: rc

      type(c_ptr) :: raw_ptr

      rc = CC_FAILURE
      nullify(f_ptr)

      if (.not. c_associated(cpp_ptr)) return

      raw_ptr = catchem_state_get_pointer_1d(cpp_ptr, trim(name) // c_null_char)
      if (c_associated(raw_ptr)) then
         call c_f_pointer(raw_ptr, f_ptr, dims)
         rc = CC_SUCCESS
      end if
   end subroutine get_field_1d

   subroutine get_field_2d(cpp_ptr, name, f_ptr, dims, rc)
      type(c_ptr), intent(in) :: cpp_ptr
      character(len=*), intent(in) :: name
      real(fp), pointer, intent(out) :: f_ptr(:,:)
      integer, intent(in) :: dims(2)
      integer, intent(out) :: rc

      type(c_ptr) :: raw_ptr

      rc = CC_FAILURE
      nullify(f_ptr)

      if (.not. c_associated(cpp_ptr)) return

      raw_ptr = catchem_state_get_pointer_2d(cpp_ptr, trim(name) // c_null_char)
      if (c_associated(raw_ptr)) then
         call c_f_pointer(raw_ptr, f_ptr, dims)
         rc = CC_SUCCESS
      end if
   end subroutine get_field_2d

   subroutine get_field_3d(cpp_ptr, name, f_ptr, dims, rc)
      type(c_ptr), intent(in) :: cpp_ptr
      character(len=*), intent(in) :: name
      real(fp), pointer, intent(out) :: f_ptr(:,:,:)
      integer, intent(in) :: dims(3)
      integer, intent(out) :: rc

      type(c_ptr) :: raw_ptr

      rc = CC_FAILURE
      nullify(f_ptr)

      if (.not. c_associated(cpp_ptr)) return

      raw_ptr = catchem_state_get_pointer_3d(cpp_ptr, trim(name) // c_null_char)
      if (c_associated(raw_ptr)) then
         call c_f_pointer(raw_ptr, f_ptr, dims)
         rc = CC_SUCCESS
      end if
   end subroutine get_field_3d

end module Interop_Mod
```

- [ ] **Step 2: Temporarily add to `src/core/CMakeLists.txt` for compilation verification**

We add `Interop_Mod.F90` to the compiled state sources `_state_srcs` to verify it compiles.

```cmake
set(
  _state_srcs
  TimeState_Mod.F90
  species_mod.F90
  GridGeometry_Mod.F90
  metstate_mod.F90
  chemstate_mod.F90
  ExtEmisData_Mod.F90
  DiagnosticInterface_Mod.F90
  ConfigManager_Mod.F90
  StateManager_Mod.F90
  DiagnosticManager_Mod.F90
  ChemSpeciesUtils_Mod.F90
  Interop_Mod.F90
)
```

- [ ] **Step 3: Run Docker compilation to ensure the new module compiles cleanly**

Command:
```bash
docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace cece-dev:latest bash -c "apt-get update && apt-get install -y python3 && cmake -B build-cece -DCATCHEM_BUILD_NUOPC=OFF -DENABLE_KOKKOS=ON -DCATCHEM_BUILD_TESTING=ON && cmake --build build-cece --parallel 4"
```
Expected output: Compiles cleanly with zero errors.

---

### Task 2: Store `cpp_ptr` within `MetStateType` and Bind Pointers

**Files:**
- Modify: `src/core/metstate_mod.F90`
- Modify: `src/core/StateManager_Mod.F90`

**Interfaces:**
- Consumes: `StateManagerType` and `MetStateType`
- Produces: `MetStateType%cpp_ptr` populated at initialization, and all array fields declared as `POINTER` instead of `ALLOCATABLE` to support zero-copy dynamic assignment.

- [ ] **Step 1: Convert `ALLOCATABLE` fields to standard `POINTER` in `MetStateType`**

Change the variable declarations inside `TYPE, PUBLIC :: MetStateType` to use `POINTER` with default nullification. For example:
```fortran
      ! Grid flags (2D: nx, ny)
      LOGICAL, POINTER         :: IsLand(:,:) => null()
      LOGICAL, POINTER         :: IsWater(:,:) => null()
      LOGICAL, POINTER         :: IsIce(:,:) => null()
      LOGICAL, POINTER         :: IsSnow(:,:) => null()
```
And similarly for all REAL and INTEGER arrays in the type definition (including `T`, `QV`, `PS`, `DELP` etc.). Add `type(c_ptr) :: cpp_ptr = c_null_ptr` to `MetStateType`.

- [ ] **Step 2: Bind pointers dynamically in StateManager_Mod when `get_met_state_ptr` is resolved**

In `StateManager_Mod.F90`, update `state_mgr_get_met_state_ptr` to automatically bind the pointers of `MetStateType` to the raw C++ memory if `cpp_ptr` is associated:
```fortran
   function state_mgr_get_met_state_ptr(this) result(ptr)
      use Interop_Mod, only: get_cpp_field
      class(StateManagerType), intent(in) :: this
      type(MetStateType), pointer :: ptr
      integer :: rc, nx, ny, nz, nl

      ptr => this%met_state
      if (associated(ptr) .and. c_associated(this%cpp_ptr)) then
         ptr%cpp_ptr = this%cpp_ptr
         call ptr%geometry%get_dimensions(nx, ny, nz)

         ! Associate key meteorological pointers dynamically!
         call get_cpp_field(this%cpp_ptr, "T", ptr%T, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "QV", ptr%QV, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "RH", ptr%RH, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "PMID", ptr%PMID, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "PEDGE", ptr%PEDGE, [nx, ny, nz+1], rc)
         call get_cpp_field(this%cpp_ptr, "AIRDEN", ptr%AIRDEN, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "AIRDEN_DRY", ptr%AIRDEN_DRY, [nx, ny, nz], rc)
         call get_cpp_field(this%cpp_ptr, "BXHEIGHT", ptr%BXHEIGHT, [nx, ny, nz], rc)

         ! Associate 2D pointers
         call get_cpp_field(this%cpp_ptr, "PS", ptr%PS, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "TS", ptr%TS, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "PBLH", ptr%PBLH, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "USTAR", ptr%USTAR, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "HFLUX", ptr%HFLUX, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "OBK", ptr%OBK, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "LAT", ptr%LAT, [nx, ny], rc)
         call get_cpp_field(this%cpp_ptr, "LON", ptr%LON, [nx, ny], rc)
      end if
   end function state_mgr_get_met_state_ptr
```

---

### Task 3: Restructure Accessor Procedures in `metstate_mod.F90`

**Files:**
- Modify: `src/core/metstate_mod.F90`

**Interfaces:**
- Consumes: `Interop_Mod` pointer retrievals
- Produces: Backward-compatible public interfaces for `get_column_ptr`, `get_2Dto0D_value`, and `get_scalar_value`.

- [ ] **Step 1: Refactor type-bound query routines to directly return pointer slices**

Refactor `metstate_get_column_ptr_subroutine` to dynamically retrieve pointers or slice our pointers directly:
```fortran
   subroutine metstate_get_column_ptr_subroutine(this, field_name, i, j, col_ptr, rc)
      use Interop_Mod, only: get_cpp_field
      use Error_Mod, only: CC_SUCCESS, CC_FAILURE

      class(MetStateType), intent(inout), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      real(fp), pointer :: col_ptr(:)
      integer, intent(out) :: rc

      real(fp), pointer :: full_3d(:,:,:)
      integer :: nx, ny, nz, nl

      rc = CC_FAILURE
      nullify(col_ptr)

      call this%geometry%get_dimensions(nx, ny, nz)

      nl = nz
      if (trim(field_name) == "PEDGE" .or. trim(field_name) == "PFILSAN" .or. trim(field_name) == "PFLLSAN") then
         nl = nz + 1
      end if

      ! If C++ is associated, query dynamically. Else, slice local pointers.
      if (c_associated(this%cpp_ptr)) then
         call get_cpp_field(this%cpp_ptr, field_name, full_3d, [nx, ny, nl], rc)
         if (rc == CC_SUCCESS) then
            col_ptr => full_3d(i, j, :)
         end if
      else
         ! Standalone/Fallback Local pointer slicing (replaces generated cases)
         select case (trim(field_name))
          case ('T')
            if (associated(this%T)) col_ptr => this%T(i, j, :)
          case ('QV')
            if (associated(this%QV)) col_ptr => this%QV(i, j, :)
          case ('RH')
            if (associated(this%RH)) col_ptr => this%RH(i, j, :)
          case ('PMID')
            if (associated(this%PMID)) col_ptr => this%PMID(i, j, :)
          case ('PEDGE')
            if (associated(this%PEDGE)) col_ptr => this%PEDGE(i, j, :)
          case ('AIRDEN')
            if (associated(this%AIRDEN)) col_ptr => this%AIRDEN(i, j, :)
          case ('AIRDEN_DRY')
            if (associated(this%AIRDEN_DRY)) col_ptr => this%AIRDEN_DRY(i, j, :)
          case ('BXHEIGHT')
            if (associated(this%BXHEIGHT)) col_ptr => this%BXHEIGHT(i, j, :)
         end select
         if (associated(col_ptr)) rc = CC_SUCCESS
      end if
   end subroutine metstate_get_column_ptr_subroutine
```

---

### Task 4: Eradicate Code Generation Tools and Include Files in Build System

**Files:**
- Delete: `cmake/generate_metstate_macros.py`
- Modify: `src/core/CMakeLists.txt`

**Interfaces:** None

- [ ] **Step 1: Delete `cmake/generate_metstate_macros.py`**

Command:
```bash
git rm cmake/generate_metstate_macros.py
```

- [ ] **Step 2: Clean up `src/core/CMakeLists.txt`**

Remove all generated INCLUDE macro custom commands (`add_custom_command` rules for `METSTATE_ACCESSOR_INC`, `METSTATE_ALLOCATE_INC`, `METSTATE_DEALLOCATE_INC`, `METSTATE_2D_SCALAR_ACCESSOR_INC` etc.) and the `metstate_macros` custom target.

Ensure clean linkage against `Interop_Mod.F90`.

---

### Task 3: Build & Interop Verification

**Files:**
- Test: All 9 CATChem unit tests

- [ ] **Step 1: Build from scratch inside Docker**

Verify there are absolutely no file-not-found compile errors during macro generation bypass.

Command:
```bash
docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace cece-dev:latest bash -c "apt-get update && apt-get install -y python3 && cmake -B build-cece -DCATCHEM_BUILD_NUOPC=OFF -DENABLE_KOKKOS=ON -DCATCHEM_BUILD_TESTING=ON && cmake --build build-cece --parallel 4"
```

- [ ] **Step 2: Execute CTest Suite**

Verify all test cases pass cleanly.

Command:
```bash
docker run --rm -v /Users/barry/Documents/CATChem:/workspace -w /workspace cece-dev:latest bash -c "ctest --test-dir build-cece --output-on-failure"
```
Expected output: `100% tests passed, 0 tests failed out of 9`.
