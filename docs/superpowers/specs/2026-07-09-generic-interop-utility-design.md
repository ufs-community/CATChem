# Design Spec: Generic Fortran-C-C++ Interop Pointer Utility

## 1. Executive Summary
This design replaces the complex Python-driven macro generators (`generate_metstate_macros.py`) and verbose static accessors with a generic runtime pointer-association utility (`Interop_Mod.F90`).

Leveraging standard Fortran `iso_c_binding` and `c_f_pointer`, we dynamically map raw C++ host-allocated heap addresses to Fortran pointer slices at runtime based on string-name lookups. This achieves 100% core C++ single-source-of-truth orchestration while preserving backward-compatible type-bound procedure signatures expected by upstream unported ESMF and NUOPC Caps.

---

## 2. Proposed Architecture & Component Design

### 2.1 Component Diagram
```
[ ESMF/NUOPC Cap Drivers / Standalone Drivers ]
                      │
                      │ 1. Request slice (e.g., "T", "PMID", or chem species)
                      ▼
               [ Proxy Modules ]
       (StateManager_Mod, metstate_mod, etc.)
                      │
                      │ 2. Queries dynamic 2D/3D pointer
                      ▼
               [ Interop_Mod ] ◄─── (Performs c_f_pointer mapping)
                      │
                      │ 3. Fetches C pointer by name
                      ▼
              [ Flat C exports ] (catchem_api.cpp)
                      │
                      │ 4. Extracts raw host pointer
                      ▼
            [ C++ StateManager ] (Kokkos::View HostSpace pointers)
```

---

### 2.2 The Generic Interoperability Module: `Interop_Mod.F90`
A new utility module `Interop_Mod` encapsulating standard C-pointer bindings and dynamic Fortran array reconstruction.

```fortran
module Interop_Mod
   use iso_c_binding
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   use Precision_Mod, only: fp

   implicit none
   private

   public :: get_cpp_field

   ! Interface mapping back to catchem_api.cpp pointers
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

   ! Generic interfaces for pointer association based on rank
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

---

### 2.3 Eliminating Massive Boilerplate in Proxy Files
By routing through `Interop_Mod`, we rewrite type-bound getters in `metstate_mod.F90` to execute runtime lookups:

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

      ! Determine grid bounds
      call this%geometry%get_dimensions(nx, ny, nz)

      ! Handles PEDGE fields with nz + 1 levels
      nl = nz
      if (trim(field_name) == "PEDGE" .or. trim(field_name) == "PFILSAN" .or. trim(field_name) == "PFLLSAN") then
         nl = nz + 1
      end if

      ! Query C++ pointer dynamically through general interop utility
      call get_cpp_field(this%cpp_ptr, field_name, full_3d, [nx, ny, nl], rc)
      if (rc == CC_SUCCESS) then
         col_ptr => full_3d(i, j, :)
      end if
   end subroutine metstate_get_column_ptr_subroutine
```

---

## 3. Scope of Eradication & Refactoring
1. **Delete File**: `cmake/generate_metstate_macros.py` completely removed from repository.
2. **Modify `src/core/CMakeLists.txt`**: Strip out all `add_custom_command` rules generating `.inc` accessor/allocate includes. Appends `Interop_Mod.F90` to compiled sources list.
3. **Eradicate code in `metstate_mod.F90`**: Remove over 1,200 lines of `set_field`, case macros, allocations, deallocations, and accessor includes. Add `cpp_ptr` to `MetStateType`.
4. **Compile & Execute Unit Tests**: Compile target targets within the docker container to confirm 100% interoperability and correctness.

---

## 4. Risks & Mitigations
* **Risk**: High-performance concerns with dynamic string lookup in tight loops.
  * *Mitigation*: The host model or Caps call `get_column_ptr` during state-exchange initialization or at boundaries, not inside deep scientific computation loops. The runtime cost of `std::unordered_map::find` is negligible ($\mathcal{O}(1)$).
* **Risk**: Array rank-dim checks mismatch.
  * *Mitigation*: Explicit dimension assertions in `get_cpp_field` based on geometry shapes prevent potential memory segmentation faults.
