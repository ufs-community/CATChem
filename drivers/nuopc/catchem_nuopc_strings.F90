!> \brief String conversion helpers for log messages in the NUOPC drivers.
!!
!! Trace and diagnostic messages are built by concatenation and handed to
!! ESMF_LogWrite directly, so no routine needs a formatting buffer for them.
!! The formats match the edit descriptors the former internal writes used
!! (I0, L1, F0.2, ES12.4) with surrounding blanks removed.
!!
!! \ingroup catchem_nuopc_group
module catchem_nuopc_strings
   use ESMF, only: ESMF_KIND_R8
   use catchem_nuopc_lengths, only: NUOPC_KEYWORD_LEN
   implicit none
   private
   public :: int_to_string, logical_to_string, real_to_string, real_to_string_sci

contains

   !> Integer as decimal text with no padding (format I0).
   pure function int_to_string(val) result(str)
      integer, intent(in) :: val
      character(len=:), allocatable :: str
      character(len=NUOPC_KEYWORD_LEN) :: buf
      write(buf, '(I0)') val
      str = trim(buf)
   end function int_to_string

   !> Logical as 'T' or 'F' (format L1).
   pure function logical_to_string(val) result(str)
      logical, intent(in) :: val
      character(len=1) :: str
      str = merge('T', 'F', val)
   end function logical_to_string

   !> Real with two decimals (format F0.2), left-adjusted in a fixed-length
   !! result; callers trim it.  Used by the cap's clock messages.
   function real_to_string(val) result(str)
      real(ESMF_KIND_R8), intent(in) :: val
      character(len=NUOPC_KEYWORD_LEN) :: str
      write(str, '(f0.2)') val
      str = adjustl(str)
   end function real_to_string

   !> Real in scientific notation (format ES12.4) with surrounding blanks removed.
   pure function real_to_string_sci(val) result(str)
      real(ESMF_KIND_R8), intent(in) :: val
      character(len=:), allocatable :: str
      character(len=NUOPC_KEYWORD_LEN) :: buf
      write(buf, '(ES12.4)') val
      str = trim(adjustl(buf))
   end function real_to_string_sci

end module catchem_nuopc_strings
