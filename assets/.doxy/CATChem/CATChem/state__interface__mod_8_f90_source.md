

# File state\_interface\_mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**state\_interface\_mod.F90**](state__interface__mod_8_f90.md)

[Go to the documentation of this file](state__interface__mod_8_f90.md)


```Fortran

module state_interface_mod
   use precision_mod
   use error_mod, only: cc_success, cc_failure
   implicit none
   private

   ! Public interfaces - only basic validation utilities
   public :: statevalidatorutilstype

   ! Public enumerations
   public :: state_type_met, state_type_chem, state_type_emis, state_type_diag
   public :: state_status_uninitialized, state_status_initialized, state_status_valid, state_status_error

   ! Public utility procedures
   public :: get_state_type_name, allocate_met_field

   !=========================================================================
   ! Constants and Enumerations
   !=========================================================================

   integer, parameter :: STATE_TYPE_MET = 1
   integer, parameter :: STATE_TYPE_CHEM = 2
   integer, parameter :: STATE_TYPE_EMIS = 3
   integer, parameter :: STATE_TYPE_DIAG = 4
   integer, parameter :: STATE_TYPE_CONFIG = 5
   integer, parameter :: STATE_TYPE_GRID = 6

   integer, parameter :: STATE_STATUS_UNINITIALIZED = 0
   integer, parameter :: STATE_STATUS_INITIALIZED = 1
   integer, parameter :: STATE_STATUS_VALID = 2
   integer, parameter :: STATE_STATUS_ERROR = -1

   !=========================================================================
   ! Basic Utility Types
   !=========================================================================

   type :: statevalidatorutilstype

   contains
      procedure :: validate_dimensions => validator_validate_dimensions
      procedure :: validate_bounds => validator_validate_bounds
      procedure :: validate_consistency => validator_validate_consistency
      procedure :: check_nan_values => validator_check_nan_values
      procedure :: check_negative_values => validator_check_negative_values

   end type statevalidatorutilstype

contains

   !=========================================================================
   ! StateValidatorUtilsType Implementation
   !=========================================================================

   subroutine validator_validate_dimensions(this, array_shape, expected_shape, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      integer, intent(in) :: array_shape(:)
      integer, intent(in) :: expected_shape(:)
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      if (size(array_shape) /= size(expected_shape)) then
         rc = cc_failure
         return
      endif

      do i = 1, size(array_shape)
         if (array_shape(i) /= expected_shape(i)) then
            rc = cc_failure
            return
         endif
      enddo

   end subroutine validator_validate_dimensions

   subroutine validator_validate_bounds(this, values, min_val, max_val, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      real(fp), intent(in) :: min_val, max_val
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      do i = 1, size(values)
         if (values(i) < min_val .or. values(i) > max_val) then
            rc = cc_failure
            return
         endif
      enddo

   end subroutine validator_validate_bounds

   subroutine validator_validate_consistency(this, array1, array2, tolerance, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: array1(:), array2(:)
      real(fp), intent(in) :: tolerance
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      if (size(array1) /= size(array2)) then
         rc = cc_failure
         return
      endif

      do i = 1, size(array1)
         if (abs(array1(i) - array2(i)) > tolerance) then
            rc = cc_failure
            return
         endif
      enddo

   end subroutine validator_validate_consistency

   subroutine validator_check_nan_values(this, values, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      do i = 1, size(values)
         if (values(i) /= values(i)) then  ! NaN check
            rc = cc_failure
            return
         endif
      enddo

   end subroutine validator_check_nan_values

   subroutine validator_check_negative_values(this, values, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      do i = 1, size(values)
         if (values(i) < 0.0_fp) then
            rc = cc_failure
            return
         endif
      enddo

   end subroutine validator_check_negative_values

   !=========================================================================
   ! Utility Functions
   !=========================================================================

   function get_state_type_name(state_type) result(name)
      integer, intent(in) :: state_type
      character(len=32) :: name

      select case (state_type)
       case (state_type_met)
         name = 'Meteorology'
       case (state_type_chem)
         name = 'Chemistry'
       case (state_type_emis)
         name = 'Emissions'
       case (state_type_diag)
         name = 'Diagnostics'
       case (state_type_config)
         name = 'Configuration'
       case (state_type_grid)
         name = 'Grid'
       case default
         name = 'Unknown'
      end select

   end function get_state_type_name

   subroutine allocate_met_field(met_state, field_name, rc)
      use metstate_mod, only: metstatetype
      class(MetStateType), intent(inout) :: met_state
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc

      call met_state%allocate_field(field_name, rc)
   end subroutine allocate_met_field

end module state_interface_mod
```


