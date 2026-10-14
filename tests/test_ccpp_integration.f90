program test_ccpp_integration
   use iso_c_binding, only: c_int, c_char, c_null_char, c_ptr
   use CATChem_API, only: CATChem_Model
   use ccpp_catchem_interface, only: ccpp_catchem_interface_get_staging_metrics
   implicit none
   interface
      integer(c_int) function catchem_state_get_species_count_checked(state, count) &
         bind(C, name='catchem_state_get_species_count_checked')
         import :: c_int, c_ptr
         type(c_ptr), value :: state
         integer(c_int), intent(out) :: count
      end function
      integer(c_int) function catchem_state_get_species_name_at_checked(state, index, name, length) &
         bind(C, name='catchem_state_get_species_name_at_checked')
         import :: c_int, c_char, c_ptr
         type(c_ptr), value :: state
         integer(c_int), value :: index, length
         character(c_char), intent(out) :: name(*)
      end function
   end interface
   type(CATChem_Model) :: model
   integer :: allocations, gathers, scatters, rc, issues, i
   integer(c_int) :: count, status
   character(kind=c_char) :: c_name(64)
   character(len=64) :: name
   character(len=512) :: report
   character(len=64), parameter :: expected(3) = [character(len=64) :: &
      'unfamiliar_alpha', 'unfamiliar_beta', 'unfamiliar_gamma']

   call ccpp_catchem_interface_get_staging_metrics(allocations, gathers, scatters)
   if (allocations /= 0 .or. gathers /= 0 .or. scatters /= 0) error stop 'invalid initial staging metrics'
   call model%initialize('host_conformance/CATChem_config.yml', 2, 1, 3, rc=rc)
   if (rc /= 0) error stop 'CCPP parity fixture initialization failed'
   status = catchem_state_get_species_count_checked(model%state_mgr_ptr, count)
   if (status /= 0_c_int .or. count /= 3_c_int) error stop 'CCPP mechanism count parity failed'
   do i = 1, 3
      status = catchem_state_get_species_name_at_checked(model%state_mgr_ptr, int(i, c_int), c_name, 64_c_int)
      if (status /= 0_c_int) error stop 'CCPP mechanism name query failed'
      call from_c(c_name, name)
      if (trim(name) /= trim(expected(i))) error stop 'CCPP mechanism ordering parity failed'
   end do
   call model%set_physical_validation_policy(99, rc)
   if (rc /= 8) error stop 'CCPP invalid-policy category parity failed'
   call model%set_physical_validation_policy(0, rc)
   call model%get_physical_validation_report(issues, report, rc)
   if (rc /= 0 .or. issues /= 0) error stop 'CCPP physical-report parity failed'
   call model%finalize(rc)
   if (rc /= 0) error stop 'CCPP parity fixture finalization failed'
contains
   subroutine from_c(c_value, value)
      character(kind=c_char), intent(in) :: c_value(*)
      character(len=*), intent(out) :: value
      integer :: j
      value = ''
      do j = 1, len(value)
         if (c_value(j) == c_null_char) exit
         value(j:j) = c_value(j)
      end do
   end subroutine from_c
end program test_ccpp_integration
