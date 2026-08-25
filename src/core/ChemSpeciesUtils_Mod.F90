!> \file ChemSpeciesUtils_Mod.F90
!! \brief Lightweight backward-compatible Fortran wrapper for ChemSpeciesUtils procedures
!!
module ChemSpeciesUtils_Mod
   use iso_c_binding
   use StateManager_Mod, only: StateManagerType
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE

   implicit none
   private

   public :: create_species_mapping

   ! Interface mapping back to catchem_api.cpp
   interface
      integer(c_int) function catchem_state_get_species_index(state_ptr, name) bind(C, name="catchem_state_get_species_index")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function
   end interface

contains

   ! Maps an incoming set of tracer names to 1-based CATChem core species indices
   subroutine create_species_mapping(state_mgr, names, mapping, rc)
      type(StateManagerType), intent(in) :: state_mgr
      character(len=*), intent(in) :: names(:)
      integer, intent(out) :: mapping(:)
      integer, intent(out) :: rc

      integer :: i, j, f_len, cc_idx
      character(kind=c_char) :: c_name(64)

      rc = CC_SUCCESS
      mapping = -1

      if (.not. c_associated(state_mgr%cpp_ptr)) then
         rc = CC_FAILURE
         return
      end if

      do i = 1, size(names)
         f_len = len_trim(names(i))
         if (f_len == 0 .or. f_len > 63) then
            mapping(i) = -1
            cycle
         end if

         ! Convert name to null-terminated standard C-string
         do j = 1, f_len
            c_name(j) = names(i)(j:j)
         end do
         c_name(f_len+1) = c_null_char

         ! Call standard C-API to find the 1-based index (matching Fortran arrays)
         cc_idx = int(catchem_state_get_species_index(state_mgr%cpp_ptr, c_name))
         mapping(i) = cc_idx
      end do

   end subroutine create_species_mapping

end module ChemSpeciesUtils_Mod
