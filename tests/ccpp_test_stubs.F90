module machine
   use iso_fortran_env, only: real64
   implicit none
   integer, parameter :: kind_phys = real64
end module machine

module ccpp_constituent_prop_mod
   use machine, only: kind_phys
   implicit none
   type :: ccpp_constituent_properties_t
      character(len=64) :: std_name = ''
   contains
      procedure :: instantiate
   end type
   type :: ccpp_constituent_prop_ptr_t
      type(ccpp_constituent_properties_t), pointer :: p => null()
   end type
contains
   subroutine instantiate(self, std_name, long_name, diag_name, units, vertical_dim, default_value, min_value, &
      molar_mass, advected, errcode, errmsg)
      class(ccpp_constituent_properties_t), intent(inout) :: self
      character(len=*), intent(in) :: std_name, long_name, diag_name, units, vertical_dim
      real(kind_phys), intent(in) :: default_value, min_value
      real(kind_phys), intent(in), optional :: molar_mass
      logical, intent(in) :: advected
      integer, intent(out) :: errcode
      character(len=*), intent(out) :: errmsg
      self%std_name = std_name
      errcode = 0
      errmsg = ''
   end subroutine instantiate
end module ccpp_constituent_prop_mod

module ccpp_const_utils
   use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
   implicit none
contains
   subroutine ccpp_const_get_idx(props, name, idx, errmsg, errflg)
      type(ccpp_constituent_prop_ptr_t), intent(in) :: props(:)
      character(len=*), intent(in) :: name
      integer, intent(out) :: idx
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg
      integer :: i
      idx = -1
      do i = 1, size(props)
         if (associated(props(i)%p)) then
            if (trim(props(i)%p%std_name) == trim(name)) idx = i
         end if
      end do
      if (idx > 0) then
         errflg = 0
         errmsg = ''
      else
         errflg = 1
         errmsg = 'species not found'
      end if
   end subroutine ccpp_const_get_idx
end module ccpp_const_utils
