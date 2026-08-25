!> \file Interop_Mod.F90
!! \brief Standard conforming ISO_C_BINDING dynamic pointer association utility
!!
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
