!> \file StateManager_Mod.F90
!! \brief Unified state management module for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides unified state management for CATChem, including:
!! - State container and basic management
!! - State validation utilities
!! - State type constants and utilities
!! - Direct field access functions
!!
!! Most complex lifecycle management is delegated to CATChemCore_Mod.
!!
module StateManager_Mod
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use ConfigManager_Mod, only: ConfigManagerType
   use MetState_Mod, only: MetStateType
   use TimeState_Mod, only: TimeStateType
   use ChemState_Mod, only: ChemStateType
   use GridManager_Mod, only: GridManagerType
   use DiagnosticManager_Mod, only: DiagnosticManagerType
   use VirtualColumn_Mod, only: VirtualColumnType

   implicit none
   private

   ! Public types
   public :: StateManagerType
   public :: StateValidatorUtilsType

   ! Public enumerations and constants
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
   ! Types
   !=========================================================================

   !> \brief Simplified state container managed by CATChemCore
   !!
   !! This type provides basic state object management. Complex lifecycle
   !! operations are handled by CATChemCore_Mod.
   !!
   type :: StateManagerType
      private

      ! Core state objects
      type(MetStateType),   allocatable :: met_state   !< Meteorological fields
      type(TimeStateType),  allocatable :: time_state  !< Time state and solar calculations
      type(ChemStateType),  allocatable :: chem_state  !< Chemical species concentrations
      type(ErrorManagerType)            :: error_mgr   !< Error manager

      ! Manager pointers (owned by CATChemCore)
      type(ConfigManagerType), pointer :: config => null()  !< Configuration manager
      type(GridManagerType), pointer :: grid_mgr => null()  !< Grid manager
      type(DiagnosticManagerType), pointer :: diag_mgr => null()  !< Diagnostic manager

      ! Simple metadata
      logical :: is_initialized = .false.              !< Initialization status
      logical :: is_configured = .false.               !< Configuration status
      character(len=256) :: name = ''                  !< Container name
      real(fp), public :: tstep = 0.0_fp               !< time step for all processes

   contains
      ! Basic lifecycle (called by CATChemCore)
      procedure :: init => manager_init
      procedure :: cleanup => manager_cleanup
      procedure :: finalize => manager_finalize
      procedure :: is_ready => manager_is_ready
      procedure :: set_configured => manager_set_configured

      ! State object accessors
      procedure :: get_config_ptr => manager_get_config_ptr
      procedure :: set_config => manager_set_config
      procedure :: get_met_state_ptr => manager_get_met_state_ptr
      procedure :: get_time_state_ptr => manager_get_time_state_ptr
      procedure :: get_chem_state_ptr => manager_get_chem_state_ptr
      procedure :: get_error_manager => manager_get_error_manager
      procedure :: get_grid_manager => manager_get_grid_manager
      procedure :: set_grid_manager => manager_set_grid_manager
      procedure :: get_diagnostic_manager => manager_get_diagnostic_manager
      procedure :: set_diagnostic_manager => manager_set_diagnostic_manager
      procedure :: create_virtual_column => manager_create_virtual_column
      procedure :: apply_virtual_column => manager_apply_virtual_column
      procedure :: populate_virtual_column => populate_virtual_column

      ! Utilities
      procedure :: set_name => manager_set_name
      procedure :: print_info => manager_print_info
      procedure :: get_memory_usage => manager_get_memory_usage

   end type StateManagerType

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
   ! StateManagerType Implementation
   !=========================================================================

   !> \brief Initialize the state manager (called by CATChemCore)
   subroutine manager_init(this, name, rc)
      class(StateManagerType), intent(inout) :: this
      character(len=*), optional, intent(in) :: name
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Set manager name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'StateManager'
      endif

      ! Allocate and initialize state objects
      ! Note: config will be set by external call if needed

      if (.not. allocated(this%met_state)) allocate(this%met_state)

      if (.not. allocated(this%time_state)) allocate(this%time_state)

      if (.not. allocated(this%chem_state)) allocate(this%chem_state)

      this%is_initialized = .true.
      this%is_configured = .false.

   end subroutine manager_init

   !> \brief Clean up the state manager
   subroutine manager_cleanup(this, rc)
      class(StateManagerType), intent(inout), target :: this
      integer, intent(out) :: rc

      integer :: config_rc, met_rc, time_rc, chem_rc
      type(ErrorManagerType), pointer :: error_mgr_ptr

      rc = CC_SUCCESS

      ! Get pointer to error manager for cleanup calls
      error_mgr_ptr => this%get_error_manager()

      ! Clean up and deallocate state objects - call their cleanup procedures first!
      if (allocated(this%met_state)) then
         call this%met_state%cleanup('ALL', met_rc)
         if (met_rc /= CC_SUCCESS) rc = met_rc  ! Don't stop cleanup on error
         deallocate(this%met_state)
      end if

      if (allocated(this%time_state)) then
         call this%time_state%cleanup(error_mgr_ptr, time_rc)
         if (time_rc /= CC_SUCCESS) rc = time_rc  ! Don't stop cleanup on error
         deallocate(this%time_state)
      end if

      if (allocated(this%chem_state)) then
         call this%chem_state%cleanup(chem_rc)
         if (chem_rc /= CC_SUCCESS) rc = chem_rc  ! Don't stop cleanup on error
         deallocate(this%chem_state)
      end if

      ! Finalize and deallocate state objects
      if (associated(this%config)) then
         call this%config%finalize(config_rc)
         if (config_rc /= CC_SUCCESS) rc = config_rc  ! Don't stop cleanup on error
         nullify(this%config)  ! Just nullify pointer, don't deallocate (owned by CATChemCore)
      end if

      this%is_initialized = .false.
      this%is_configured = .false.
      this%name = ''
      this%tstep = 0.0_fp

   end subroutine manager_cleanup

   !> \brief Check if manager is ready
   function manager_is_ready(this) result(ready)
      class(StateManagerType), intent(in) :: this
      logical :: ready

      ready = this%is_initialized .and. this%is_configured .and. &
         associated(this%config) .and. &
         allocated(this%met_state) .and. &
         allocated(this%time_state) .and. &
         allocated(this%chem_state)
   end function manager_is_ready

   !> \brief Mark the state manager as configured
   subroutine manager_set_configured(this)
      class(StateManagerType), intent(inout) :: this

      this%is_configured = .true.
   end subroutine manager_set_configured

   !> \brief Get pointer to config manager for modification
   function manager_get_config_ptr(this) result(config_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(ConfigManagerType), pointer :: config_ptr

      if (associated(this%config)) then
         config_ptr => this%config
      else
         nullify(config_ptr)
      endif
   end function manager_get_config_ptr

   !> \brief Set the config manager to use an external instance
   subroutine manager_set_config(this, external_config, rc)
      class(StateManagerType), intent(inout) :: this
      type(ConfigManagerType), intent(in), target :: external_config
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Point to the external config (share the same instance)
      this%config => external_config
   end subroutine manager_set_config

   !> \brief Get pointer to met state for modification
   function manager_get_met_state_ptr(this) result(met_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(MetStateType), pointer :: met_ptr

      if (allocated(this%met_state)) then
         met_ptr => this%met_state
      else
         nullify(met_ptr)
      endif
   end function manager_get_met_state_ptr

   !> \brief Get pointer to time state for modification
   function manager_get_time_state_ptr(this) result(time_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(TimeStateType), pointer :: time_ptr

      if (allocated(this%time_state)) then
         time_ptr => this%time_state
      else
         nullify(time_ptr)
      endif
   end function manager_get_time_state_ptr

   !> \brief Get pointer to chem state for modification
   function manager_get_chem_state_ptr(this) result(chem_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(ChemStateType), pointer :: chem_ptr

      if (allocated(this%chem_state)) then
         chem_ptr => this%chem_state
      else
         nullify(chem_ptr)
      endif
   end function manager_get_chem_state_ptr

   !> \brief Get pointer to error manager
   function manager_get_error_manager(this) result(error_mgr_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(ErrorManagerType), pointer :: error_mgr_ptr

      error_mgr_ptr => this%error_mgr
   end function manager_get_error_manager

   !> \brief Get pointer to grid manager
   function manager_get_grid_manager(this) result(grid_mgr_ptr)
      class(StateManagerType), intent(in), target :: this
      type(GridManagerType), pointer :: grid_mgr_ptr

      if (associated(this%grid_mgr)) then
         grid_mgr_ptr => this%grid_mgr
      else
         nullify(grid_mgr_ptr)
      endif
   end function manager_get_grid_manager

   !> \brief Set grid manager pointer
   subroutine manager_set_grid_manager(this, grid_mgr_ptr, rc)
      class(StateManagerType), intent(inout) :: this
      type(GridManagerType), pointer, intent(in) :: grid_mgr_ptr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (associated(grid_mgr_ptr)) then
         this%grid_mgr => grid_mgr_ptr
      else
         rc = CC_FAILURE
      endif
   end subroutine manager_set_grid_manager

   !> \brief Get pointer to diagnostic manager
   function manager_get_diagnostic_manager(this) result(diag_mgr_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(DiagnosticManagerType), pointer :: diag_mgr_ptr

      if (associated(this%diag_mgr)) then
         diag_mgr_ptr => this%diag_mgr
      else
         nullify(diag_mgr_ptr)
      endif
   end function manager_get_diagnostic_manager

   !> \brief Set pointer to diagnostic manager (called by CATChemCore)
   subroutine manager_set_diagnostic_manager(this, diag_mgr_ptr)
      class(StateManagerType), intent(inout) :: this
      type(DiagnosticManagerType), intent(in), target :: diag_mgr_ptr

      this%diag_mgr => diag_mgr_ptr
   end subroutine manager_set_diagnostic_manager

   !> \brief Create virtual column for grid position (i,j)
   !! \details Creates a virtual column data container and populates it
   !! with data from the 3D grid at position (i,j)
   subroutine manager_create_virtual_column(this, i, j, virtual_col, rc)
      class(StateManagerType), intent(inout), target :: this
      integer, intent(in) :: i, j
      type(VirtualColumnType), intent(out) :: virtual_col
      integer, intent(out) :: rc

      integer :: nlev, nspec_chem, nspec_emis
      real(fp) :: lat, lon, area
      integer :: nx, ny, temp_nlev

      rc = CC_SUCCESS

      ! Get dimensions from MetState if available
      if (allocated(this%met_state)) then
         call this%met_state%get_dimensions(nx, ny, temp_nlev)
         nlev = temp_nlev
      else
         nlev = 50  ! Default fallback
      endif

      ! Get number of chemical species
      if (allocated(this%chem_state) .and. allocated(this%chem_state%ChemSpecies)) then
         nspec_chem = size(this%chem_state%ChemSpecies)
      else
         nspec_chem = 20  ! Default fallback
      endif

      ! For now, set emissions species to 0 (would come from EmissionState)
      nspec_emis = 0

      ! Get position metadata from MetState if available
      if (allocated(this%met_state)) then
         if (allocated(this%met_state%LAT) .and. allocated(this%met_state%LON)) then
            lat = this%met_state%LAT(i, j)
            lon = this%met_state%LON(i, j)
         else
            lat = 0.0_fp
            lon = 0.0_fp
         endif

         if (allocated(this%met_state%AREA_M2)) then
            area = this%met_state%AREA_M2(i, j)
         else
            area = 1.0_fp
         endif
      else
         lat = 0.0_fp
         lon = 0.0_fp
         area = 1.0_fp
      endif

      ! Initialize the virtual column data container
      call virtual_col%init(nlev, nspec_chem, nspec_emis, i, j, lat, lon, area, rc)
      if (rc /= CC_SUCCESS) return

      ! Populate with data from 3D grid
      call this%populate_virtual_column(virtual_col, rc)

   end subroutine manager_create_virtual_column

   !> \brief Apply virtual column changes back to 3D state
   !! \details Updates 3D state arrays with modified virtual column data
   subroutine manager_apply_virtual_column(this, virtual_col, rc)
      class(StateManagerType), intent(inout) :: this
      type(VirtualColumnType), intent(in) :: virtual_col
      integer, intent(out) :: rc

      integer :: grid_i, grid_j, k, ispec
      integer :: nlev, nspec_chem, nspec_emis
      real(fp) :: met_value, chem_value

      rc = CC_SUCCESS

      ! Check if states are allocated
      if (.not. allocated(this%met_state) .or. .not. allocated(this%chem_state)) then
         rc = CC_FAILURE
         return
      endif

      ! Get column position and dimensions
      call virtual_col%get_position(grid_i, grid_j)
      call virtual_col%get_dimensions(nlev, nspec_chem, nspec_emis)

      ! NOTE: Meteorological data copy-back is not needed because VirtualMetType
      ! uses pointers directly to the 3D arrays. Changes are automatically reflected.

      ! Apply chemical species data back to 3D arrays
      if (allocated(this%chem_state%ChemSpecies) .and. nspec_chem > 0) then
         do ispec = 1, min(nspec_chem, size(this%chem_state%ChemSpecies))
            if (associated(this%chem_state%ChemSpecies(ispec)%conc)) then
               do k = 1, nlev
                  ! Get modified concentration from virtual column
                  chem_value = virtual_col%get_chem_field(ispec, k)
                  ! Apply back to the 3D concentration array
                  this%chem_state%ChemSpecies(ispec)%conc(grid_i, grid_j, k) = chem_value
               end do
            endif
         end do
      endif

   end subroutine manager_apply_virtual_column

   !> \brief Populate virtual column from 3D state
   !! \details Uses VirtualMetType with direct pointers to MetState fields
   subroutine populate_virtual_column(this, virtual_col, rc)
      class(StateManagerType), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: virtual_col
      integer, intent(out) :: rc

      integer :: grid_i, grid_j, k, ispec
      integer :: nlev, nspec_chem, nspec_emis
      real(fp), pointer :: column_ptr(:)
      integer, pointer :: column_ptr_int(:)
      logical, pointer :: column_ptr_logical(:)
      real(fp) :: scalar_val
      integer :: scalar_val_int
      logical :: scalar_val_logical
      real(fp) :: chem_value
      integer :: field_rc

      rc = CC_SUCCESS

      ! Check if states are allocated
      if (.not. allocated(this%met_state) .or. .not. allocated(this%chem_state)) then
         rc = CC_FAILURE
         return
      endif

      ! Get column position and dimensions
      call virtual_col%get_position(grid_i, grid_j)
      call virtual_col%get_dimensions(nlev, nspec_chem, nspec_emis)

      ! Populate VirtualMetType using generated macro
#include "virtualmet_populate.inc"

      ! Extract chemical species data (unchanged)
      if (allocated(this%chem_state%ChemSpecies) .and. nspec_chem > 0) then
         do ispec = 1, min(nspec_chem, size(this%chem_state%ChemSpecies))
            if (associated(this%chem_state%ChemSpecies(ispec)%conc)) then
               do k = 1, nlev
                  chem_value = this%chem_state%ChemSpecies(ispec)%conc(grid_i, grid_j, k)
                  call virtual_col%set_chem_field(k, ispec, chem_value)
               end do
            else
               do k = 1, nlev
                  call virtual_col%set_chem_field(k, ispec, 0.0_fp)
               end do
            endif
         end do
      endif

   end subroutine populate_virtual_column

   !> \brief Set manager name
   subroutine manager_set_name(this, name)
      class(StateManagerType), intent(inout) :: this
      character(len=*), intent(in) :: name

      this%name = trim(name)
   end subroutine manager_set_name

   !> \brief Print manager information
   subroutine manager_print_info(this)
      class(StateManagerType), intent(in) :: this

      write(*,'(A)') '=== StateManager Information ==='
      write(*,'(A,A)') 'Name: ', trim(this%name)
      write(*,'(A,L1)') 'Initialized: ', this%is_initialized
      write(*,'(A,L1)') 'Config manager associated: ', associated(this%config)
      write(*,'(A,L1)') 'Met state allocated: ', allocated(this%met_state)
      write(*,'(A,L1)') 'Time state allocated: ', allocated(this%time_state)
      write(*,'(A,L1)') 'Chem state allocated: ', allocated(this%chem_state)
      write(*,'(A)') '================================='

   end subroutine manager_print_info

   !> \brief Finalize the state manager (alias for cleanup)
   subroutine manager_finalize(this, rc)
      class(StateManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%cleanup(rc)
   end subroutine manager_finalize

   !> \brief Get approximate memory usage in bytes
   function manager_get_memory_usage(this) result(memory_bytes)
      class(StateManagerType), intent(in) :: this
      integer(8) :: memory_bytes

      ! Simplified calculation - real implementation would query each state object
      memory_bytes = 0_8

      if (associated(this%config)) memory_bytes = memory_bytes + 1024_8
      if (allocated(this%met_state)) memory_bytes = memory_bytes + 102400_8
      if (allocated(this%time_state)) memory_bytes = memory_bytes + 32_8
      if (allocated(this%chem_state)) memory_bytes = memory_bytes + 1048576_8
   end function manager_get_memory_usage

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
      class(MetStateType), intent(inout) :: met_state
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc

      call met_state%allocate_field(field_name, rc)
   end subroutine allocate_met_field

end module StateManager_Mod
