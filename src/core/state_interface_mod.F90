!> \file state_interface_mod.F90
!! \brief State interface module for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides utilities for working with CATChem state objects.
!! Most state management functionality is delegated to CATChemCore_Mod.
!!
module state_interface_mod
   use precision_mod
   use error_mod, only: CC_SUCCESS, CC_FAILURE
   implicit none
   private

   ! Public interfaces - only basic validation utilities
   public :: StateValidatorUtilsType

   ! Public enumerations
   public :: STATE_TYPE_MET, STATE_TYPE_CHEM, STATE_TYPE_EMIS, STATE_TYPE_DIAG
   public :: STATE_STATUS_UNINITIALIZED, STATE_STATUS_INITIALIZED, STATE_STATUS_VALID, STATE_STATUS_ERROR

   ! Public utility procedures
   public :: get_state_type_name, allocate_met_field

   !=========================================================================
   ! Constants and Enumerations
   !=========================================================================

   !> State type enumeration
   integer, parameter :: STATE_TYPE_MET = 1
   integer, parameter :: STATE_TYPE_CHEM = 2
   integer, parameter :: STATE_TYPE_EMIS = 3
   integer, parameter :: STATE_TYPE_DIAG = 4
   integer, parameter :: STATE_TYPE_CONFIG = 5
   integer, parameter :: STATE_TYPE_GRID = 6

   !> State status enumeration
   integer, parameter :: STATE_STATUS_UNINITIALIZED = 0
   integer, parameter :: STATE_STATUS_INITIALIZED = 1
   integer, parameter :: STATE_STATUS_VALID = 2
   integer, parameter :: STATE_STATUS_ERROR = -1

   !=========================================================================
   ! Basic Utility Types
   !=========================================================================

   !> \brief State validation utilities
   !!
   !! This type provides common validation routines that can be used
   !! across different state types.
   !!
   type :: StateValidatorUtilsType

   contains
      procedure :: validate_dimensions => validator_validate_dimensions
      procedure :: validate_bounds => validator_validate_bounds
      procedure :: validate_consistency => validator_validate_consistency
      procedure :: check_nan_values => validator_check_nan_values
      procedure :: check_negative_values => validator_check_negative_values

   end type StateValidatorUtilsType

contains

   !=========================================================================
   ! StateValidatorUtilsType Implementation
   !=========================================================================

   !> \brief Validate array dimensions
   subroutine validator_validate_dimensions(this, array_shape, expected_shape, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      integer, intent(in) :: array_shape(:)
      integer, intent(in) :: expected_shape(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (size(array_shape) /= size(expected_shape)) then
         rc = CC_FAILURE
         return
      endif

      do i = 1, size(array_shape)
         if (array_shape(i) /= expected_shape(i)) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_validate_dimensions

   !> \brief Validate value bounds
   subroutine validator_validate_bounds(this, values, min_val, max_val, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      real(fp), intent(in) :: min_val, max_val
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      do i = 1, size(values)
         if (values(i) < min_val .or. values(i) > max_val) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_validate_bounds

   !> \brief Validate consistency between related arrays
   subroutine validator_validate_consistency(this, array1, array2, tolerance, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: array1(:), array2(:)
      real(fp), intent(in) :: tolerance
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (size(array1) /= size(array2)) then
         rc = CC_FAILURE
         return
      endif

      do i = 1, size(array1)
         if (abs(array1(i) - array2(i)) > tolerance) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_validate_consistency

   !> \brief Check for NaN values
   subroutine validator_check_nan_values(this, values, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      do i = 1, size(values)
         if (values(i) /= values(i)) then  ! NaN check
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_check_nan_values

   !> \brief Check for negative values where they shouldn't exist
   subroutine validator_check_negative_values(this, values, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      do i = 1, size(values)
         if (values(i) < 0.0_fp) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_check_negative_values

   !=========================================================================
   ! Utility Functions
   !=========================================================================

   !> \brief Get state type name from enumeration
   function get_state_type_name(state_type) result(name)
      integer, intent(in) :: state_type
      character(len=32) :: name

      select case (state_type)
       case (STATE_TYPE_MET)
         name = 'Meteorology'
       case (STATE_TYPE_CHEM)
         name = 'Chemistry'
       case (STATE_TYPE_EMIS)
         name = 'Emissions'
       case (STATE_TYPE_DIAG)
         name = 'Diagnostics'
       case (STATE_TYPE_CONFIG)
         name = 'Configuration'
       case (STATE_TYPE_GRID)
         name = 'Grid'
       case default
         name = 'Unknown'
      end select

   end function get_state_type_name

   !> \brief Allocate a specific field in a MetStateType object
   !!
   !! Direct utility function for field allocation in MetStateType.
   !! For more complex state management, use CATChemCore_Mod.
   !!
   !! \param[inout] met_state   MetStateType object
   !! \param[in]    field_name  Name of the field to allocate (or 'ALL')
   !! \param[out]   rc          Return code (CC_SUCCESS or error code)
   subroutine allocate_met_field(met_state, field_name, rc)
      use MetState_Mod, only: MetStateType
      class(MetStateType), intent(inout) :: met_state
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc

      call met_state%allocate_field(field_name, rc)
   end subroutine allocate_met_field

end module state_interface_mod
