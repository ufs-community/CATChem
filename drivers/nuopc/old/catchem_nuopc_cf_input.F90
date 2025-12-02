! \file catchem_nuopc_cf_input.F90
! \brief CF-compliant data input module for CATChem NUOPC cap
!>
! \details
! This module provides comprehensive CF-compliant data input capabilities for the
! CATChem NUOPC cap. It supports automatic regridding using ESMF, time interpolation,
! refresh scheduling, and quality control. Input configuration is controlled via
! YAML configuration file with flexible dataset and field management.
!>
! Features:
! - CF-compliant NetCDF data input
! - ESMF-based regridding with multiple methods
! - Time interpolation and extrapolation
! - Automatic refresh scheduling based on templates
! - Data validation and quality control
! - Parallel I/O capabilities
! - Caching for performance optimization
! - Unit conversion and coordinate mapping
!>
! \author CATChem Development Team
! \date 11/2024
! \ingroup catchem_nuopc_group

module catchem_nuopc_cf_input

   use ESMF
   use NUOPC
   use CATChem
   use catchem_types, only: catchem_container_type
   use netcdf
   use QfYaml_Mod
   use QfYaml_Mod, only: cc_yaml_get => QFYAML_Get

   implicit none

   private

   public :: cf_input_init
   public :: cf_input_update
   public :: cf_input_finalize
   public :: cf_input_config_type
   public :: cf_dataset_type
   public :: cf_input_field_type
   public :: refresh_template_type

   ! Maximum string lengths
   integer, parameter :: MAX_STRING_LEN = 256
   integer, parameter :: MAX_NAME_LEN = 128
   integer, parameter :: MAX_DATASETS = 50
   integer, parameter :: MAX_FIELDS = 200

   ! Refresh template type
   type :: refresh_template_type
      character(len=64) :: frequency        ! daily, hourly, 6hourly, simulation_start, etc.
      integer :: offset_hours = 0           ! Offset from base time
      integer :: lookahead_hours = 0        ! Hours to look ahead
      integer :: lookback_hours = 0         ! Hours to look back
      integer :: lookahead_days = 0         ! Days to look ahead
      integer :: lookback_days = 0          ! Days to look back
   end type refresh_template_type

   ! Dataset configuration type
   type :: cf_dataset_type
      character(len=MAX_NAME_LEN) :: name
      character(len=MAX_STRING_LEN) :: description
      character(len=MAX_STRING_LEN) :: file_template
      character(len=64) :: time_frequency
      type(refresh_template_type) :: refresh_template

      ! Coordinate mapping
      character(len=MAX_NAME_LEN) :: lon_name = "longitude"
      character(len=MAX_NAME_LEN) :: lat_name = "latitude"
      character(len=MAX_NAME_LEN) :: time_name = "time"
      character(len=MAX_NAME_LEN) :: level_name = "level"

      ! Dataset-level regridding options
      character(len=64) :: regrid_method = "bilinear"
      character(len=64) :: pole_method = "none"
      character(len=64) :: unmapped_action = "error"
      character(len=64) :: extrap_method = "nearest"
      logical :: cache_weights = .true.
      logical :: check_conservation = .false.
      real(ESMF_KIND_R8) :: conservation_tolerance = 1.0e-6_ESMF_KIND_R8
      integer :: conservative_order = 1
      logical :: conservative_frac_area = .true.
      logical :: enforce_positive = .false.
      logical :: apply_floor_values = .false.

      ! Weight file options
      character(len=MAX_STRING_LEN) :: weight_file = ""
      logical :: use_weight_file = .false.
      logical :: save_weights = .false.
      character(len=MAX_STRING_LEN) :: weight_file_template = ""

      ! Runtime data
      character(len=MAX_STRING_LEN) :: current_file
      type(ESMF_Time) :: last_refresh_time
      type(ESMF_Time) :: next_refresh_time
      logical :: needs_refresh = .true.
      integer :: ncid = -1                  ! NetCDF file ID
      logical :: file_open = .false.
   end type cf_dataset_type

   ! Input field configuration type
   type :: cf_input_field_type
      character(len=MAX_NAME_LEN) :: variable_name
      character(len=MAX_NAME_LEN) :: dataset_name
      character(len=MAX_STRING_LEN) :: cf_standard_name
      character(len=MAX_STRING_LEN) :: cf_long_name
      character(len=64) :: cf_units
      character(len=MAX_STRING_LEN) :: target_var
      integer :: dimensions
      logical :: required = .true.
      logical :: enabled = .true.

      ! Regridding options
      character(len=64) :: regrid_method = "bilinear"
      logical :: time_interpolation = .true.
      logical :: boundary_only = .false.
      logical :: initial_only = .false.

      ! Unit conversion
      character(len=64) :: from_units
      character(len=64) :: to_units
      real(ESMF_KIND_R8) :: scale_factor = 1.0_ESMF_KIND_R8
      real(ESMF_KIND_R8) :: offset = 0.0_ESMF_KIND_R8

      ! Runtime data
      type(ESMF_Field) :: source_field
      type(ESMF_Field) :: target_field
      type(ESMF_RouteHandle) :: regrid_route
      logical :: regrid_initialized = .false.
      integer :: varid = -1                 ! NetCDF variable ID
      type(ESMF_Time) :: last_update_time
   end type cf_input_field_type

   ! Main configuration type
   type :: cf_input_config_type
      ! General settings
      logical :: enable_cf_input = .true.
      character(len=MAX_STRING_LEN) :: input_directory = "./input_data"
      character(len=MAX_STRING_LEN) :: cache_directory = "./cache"
      logical :: strict_cf_compliance = .true.
      logical :: allow_missing_fields = .false.
      logical :: fill_missing_with_default = .true.

      ! Regridding settings
      character(len=64) :: regrid_method = "bilinear"
      logical :: cache_regridded_data = .true.
      character(len=64) :: regrid_extrap_method = "nearest"

      ! Time interpolation
      character(len=64) :: time_interp_method = "linear"
      character(len=64) :: time_extrap_method = "constant"

      ! Datasets and fields
      integer :: n_datasets = 0
      integer :: n_fields = 0
      type(cf_dataset_type), allocatable :: datasets(:)
      type(cf_input_field_type), allocatable :: fields(:)

      ! YAML parser handle
      type(QFYAML_t) :: yaml_parser
      type(QFYAML_t) :: yaml_anchored
      logical :: yaml_initialized = .false.
   end type cf_input_config_type

   ! Module-level variables
   type(cf_input_config_type), save :: cf_config
   logical, save :: module_initialized = .false.
   type(ESMF_VM), save :: vm
   type(ESMF_Grid), save :: model_grid
   integer, save :: local_pet, pet_count

contains

   ! Initialize CF data input system
   !!
   !! \param config_file Path to YAML configuration file
   !! \param grid ESMF grid for the model
   !! \param rc Return code
   subroutine cf_input_init(config_file, grid, rc)
      character(len=*), intent(in) :: config_file
      type(ESMF_Grid), intent(in) :: grid
      integer, intent(out) :: rc

      ! Local variables
      integer :: i, status
      character(len=ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS

      ! Get VM information
      call ESMF_VMGetGlobal(vm, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_VMGet(vm, localPet=local_pet, petCount=pet_count, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Store model grid
      model_grid = grid

      ! Load configuration
      call load_cf_input_config(config_file, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Initialize datasets
      do i = 1, cf_config%n_datasets
         call initialize_dataset(cf_config%datasets(i), rc)
         if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      end do

      ! Initialize fields
      do i = 1, cf_config%n_fields
         call initialize_field(cf_config%fields(i), rc)
         if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      end do

      module_initialized = .true.

      if (local_pet == 0) then
         write(msg, '(A,I0,A,I0,A)') 'CF input system initialized with ', &
            cf_config%n_datasets, ' datasets and ', cf_config%n_fields, ' fields'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=status)
      end if

   end subroutine cf_input_init

   ! Update CF input data
   !!
   !! \param current_time Current simulation time
   !! \param container CATChem container for data storage
   !! \param rc Return code
   subroutine cf_input_update(current_time, container, rc)
      type(ESMF_Time), intent(in) :: current_time
      type(catchem_container_type), intent(inout) :: container
      integer, intent(out) :: rc

      ! Local variables
      integer :: i, status
      logical :: needs_update
      character(len=ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS

      if (.not. module_initialized) then
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="CF input system not initialized", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if

      ! Check if datasets need refresh
      do i = 1, cf_config%n_datasets
         call check_dataset_refresh(cf_config%datasets(i), current_time, &
            needs_update, rc)
         if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

         if (needs_update) then
            call refresh_dataset(cf_config%datasets(i), current_time, rc)
            if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
         end if
      end do

      ! Update fields
      do i = 1, cf_config%n_fields
         if (cf_config%fields(i)%enabled) then
            call update_field(cf_config%fields(i), current_time, container, rc)
            if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
         end if
      end do

   end subroutine cf_input_update

   ! Finalize CF input system
   !!
   !! \param rc Return code
   subroutine cf_input_finalize(rc)
      integer, intent(out) :: rc

      ! Local variables
      integer :: i, status

      rc = ESMF_SUCCESS

      if (.not. module_initialized) return

      ! Clean up fields
      do i = 1, cf_config%n_fields
         call cleanup_field(cf_config%fields(i), rc)
         ! Continue cleanup even if error occurs
      end do

      ! Clean up datasets
      do i = 1, cf_config%n_datasets
         call cleanup_dataset(cf_config%datasets(i), rc)
         ! Continue cleanup even if error occurs
      end do

      ! Clean up YAML parser
      if (cf_config%yaml_initialized) then
         call QFYAML_CleanUp(cf_config%yaml_parser)
         call QFYAML_CleanUp(cf_config%yaml_anchored)
         cf_config%yaml_initialized = .false.
      end if

      ! Deallocate arrays
      if (allocated(cf_config%datasets)) deallocate(cf_config%datasets)
      if (allocated(cf_config%fields)) deallocate(cf_config%fields)

      module_initialized = .false.

   end subroutine cf_input_finalize

   ! Load CF input configuration from YAML file
   !!
   !! \param config_file Path to YAML configuration file
   !! \param rc Return code
   subroutine load_cf_input_config(config_file, rc)
      character(len=*), intent(in) :: config_file
      integer, intent(out) :: rc

      ! Local variables
      integer :: status, i
      character(len=MAX_STRING_LEN) :: key, value
      character(len=16) :: field_key, dataset_key

      rc = ESMF_SUCCESS

      ! Initialize YAML parser
      call QFYAML_Init(trim(config_file), cf_config%yaml_parser, &
         cf_config%yaml_anchored, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_FILE_READ, &
            msg="Failed to initialize YAML parser for CF input config", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if
      cf_config%yaml_initialized = .true.

      ! Load general settings
      call load_input_settings(rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Load datasets
      call QFYAML_Get(cf_config%yaml_parser, "datasets/n_datasets", &
         cf_config%n_datasets, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read n_datasets from config", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if

      if (cf_config%n_datasets > 0) then
         allocate(cf_config%datasets(cf_config%n_datasets))

         do i = 1, cf_config%n_datasets
            write(dataset_key, '(A,I0)') "dataset_", i
            call load_dataset_config(dataset_key, cf_config%datasets(i), rc)
            if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
         end do
      end if

      ! Load fields
      call QFYAML_Get(cf_config%yaml_parser, "input_fields/n_fields", &
         cf_config%n_fields, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read n_fields from config", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if

      if (cf_config%n_fields > 0) then
         allocate(cf_config%fields(cf_config%n_fields))

         do i = 1, cf_config%n_fields
            write(field_key, '(A,I0)') "field_", i
            call load_field_config(field_key, cf_config%fields(i), rc)
            if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
         end do
      end if

   end subroutine load_cf_input_config

   ! Load input settings from YAML
   !!
   !! \param rc Return code
   subroutine load_input_settings(rc)
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      character(len=MAX_STRING_LEN) :: string_val

      rc = ESMF_SUCCESS

      ! Load settings with defaults
      call QFYAML_Get(cf_config%yaml_parser, "input_settings/enable_cf_input", &
         cf_config%enable_cf_input, status)

      call QFYAML_Get(cf_config%yaml_parser, "input_settings/input_directory", &
         string_val, status)
      if (status == QFYAML_Success) cf_config%input_directory = trim(string_val)

      call QFYAML_Get(cf_config%yaml_parser, "input_settings/cache_directory", &
         string_val, status)
      if (status == QFYAML_Success) cf_config%cache_directory = trim(string_val)

      call QFYAML_Get(cf_config%yaml_parser, "input_settings/strict_cf_compliance", &
         cf_config%strict_cf_compliance, status)

      call QFYAML_Get(cf_config%yaml_parser, "input_settings/allow_missing_fields", &
         cf_config%allow_missing_fields, status)

      call QFYAML_Get(cf_config%yaml_parser, "input_settings/regrid_method", &
         string_val, status)
      if (status == QFYAML_Success) cf_config%regrid_method = trim(string_val)

      call QFYAML_Get(cf_config%yaml_parser, "input_settings/time_interp_method", &
         string_val, status)
      if (status == QFYAML_Success) cf_config%time_interp_method = trim(string_val)

   end subroutine load_input_settings

   ! Load dataset configuration from YAML
   !!
   !! \param dataset_key YAML key for the dataset
   !! \param dataset Dataset configuration object
   !! \param rc Return code
   subroutine load_dataset_config(dataset_key, dataset, rc)
      character(len=*), intent(in) :: dataset_key
      type(cf_dataset_type), intent(inout) :: dataset
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      character(len=MAX_STRING_LEN) :: base_key, key, string_val

      rc = ESMF_SUCCESS
      base_key = "datasets/" // trim(dataset_key)

      ! Load required fields
      key = trim(base_key) // "/name"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read dataset name", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if
      dataset%name = trim(string_val)

      key = trim(base_key) // "/file_template"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read dataset file_template", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if
      dataset%file_template = trim(string_val)

      ! Load optional fields
      key = trim(base_key) // "/description"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%description = trim(string_val)

      key = trim(base_key) // "/time_frequency"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%time_frequency = trim(string_val)

      ! Load coordinate mapping
      key = trim(base_key) // "/coordinate_mapping/longitude"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%lon_name = trim(string_val)

      key = trim(base_key) // "/coordinate_mapping/latitude"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%lat_name = trim(string_val)

      key = trim(base_key) // "/coordinate_mapping/time"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%time_name = trim(string_val)

      key = trim(base_key) // "/coordinate_mapping/level"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%level_name = trim(string_val)

      ! Load refresh template
      call load_refresh_template(trim(base_key) // "/refresh_template", &
         dataset%refresh_template, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Load regrid options
      call load_dataset_regrid_options(trim(base_key) // "/regrid_options", &
         dataset, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

   end subroutine load_dataset_config

   ! Load field configuration from YAML
   !!
   !! \param field_key YAML key for the field
   !! \param field Field configuration object
   !! \param rc Return code
   subroutine load_field_config(field_key, field, rc)
      character(len=*), intent(in) :: field_key
      type(cf_input_field_type), intent(inout) :: field
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      character(len=MAX_STRING_LEN) :: base_key, key, string_val
      real(ESMF_KIND_R8) :: real_val

      rc = ESMF_SUCCESS
      base_key = "input_fields/" // trim(field_key)

      ! Load required fields
      key = trim(base_key) // "/variable_name"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read field variable_name", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if
      field%variable_name = trim(string_val)

      key = trim(base_key) // "/dataset"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read field dataset", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if
      field%dataset_name = trim(string_val)

      key = trim(base_key) // "/target_var"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read field target_var", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if
      field%target_var = trim(string_val)

      key = trim(base_key) // "/dimensions"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), field%dimensions, status)
      if (status /= QFYAML_Success) then
         call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="Failed to read field dimensions", &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if

      ! Load optional fields with defaults
      key = trim(base_key) // "/cf_standard_name"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) field%cf_standard_name = trim(string_val)

      key = trim(base_key) // "/cf_long_name"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) field%cf_long_name = trim(string_val)

      key = trim(base_key) // "/cf_units"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) field%cf_units = trim(string_val)

      key = trim(base_key) // "/required"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), field%required, status)

      key = trim(base_key) // "/regrid_method"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) field%regrid_method = trim(string_val)

      key = trim(base_key) // "/time_interpolation"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), field%time_interpolation, status)

      key = trim(base_key) // "/boundary_only"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), field%boundary_only, status)

      key = trim(base_key) // "/initial_only"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), field%initial_only, status)

      ! Load unit conversion
      key = trim(base_key) // "/unit_conversion/from_units"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) field%from_units = trim(string_val)

      key = trim(base_key) // "/unit_conversion/to_units"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) field%to_units = trim(string_val)

      key = trim(base_key) // "/unit_conversion/scale_factor"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), real_val, status)
      if (status == QFYAML_Success) field%scale_factor = real_val

   end subroutine load_field_config

   ! Load refresh template from YAML
   !!
   !! \param base_key Base YAML key
   !! \param refresh_template Refresh template object
   !! \param rc Return code
   subroutine load_refresh_template(base_key, refresh_template, rc)
      character(len=*), intent(in) :: base_key
      type(refresh_template_type), intent(inout) :: refresh_template
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      character(len=MAX_STRING_LEN) :: key, string_val

      rc = ESMF_SUCCESS

      key = trim(base_key) // "/frequency"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) refresh_template%frequency = trim(string_val)

      key = trim(base_key) // "/offset_hours"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), refresh_template%offset_hours, status)

      key = trim(base_key) // "/lookahead_hours"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), refresh_template%lookahead_hours, status)

      key = trim(base_key) // "/lookback_hours"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), refresh_template%lookback_hours, status)

      key = trim(base_key) // "/lookahead_days"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), refresh_template%lookahead_days, status)

      key = trim(base_key) // "/lookback_days"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), refresh_template%lookback_days, status)

   end subroutine load_refresh_template

   ! Load dataset regrid options from YAML
   !!
   !! \param base_key YAML base key for regrid options
   !! \param dataset Dataset configuration object
   !! \param rc Return code
   subroutine load_dataset_regrid_options(base_key, dataset, rc)
      character(len=*), intent(in) :: base_key
      type(cf_dataset_type), intent(inout) :: dataset
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      character(len=MAX_STRING_LEN) :: key, string_val
      logical :: logical_val
      real(ESMF_KIND_R8) :: real_val
      integer :: int_val

      rc = ESMF_SUCCESS

      ! Load regridding method
      key = trim(base_key) // "/method"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%regrid_method = trim(string_val)

      ! Load pole method
      key = trim(base_key) // "/pole_method"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%pole_method = trim(string_val)

      ! Load unmapped action
      key = trim(base_key) // "/unmapped_action"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%unmapped_action = trim(string_val)

      ! Load extrapolation method
      key = trim(base_key) // "/extrap_method"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%extrap_method = trim(string_val)

      ! Load caching options
      key = trim(base_key) // "/cache_weights"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), logical_val, status)
      if (status == QFYAML_Success) dataset%cache_weights = logical_val

      ! Load conservation checking
      key = trim(base_key) // "/check_conservation"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), logical_val, status)
      if (status == QFYAML_Success) dataset%check_conservation = logical_val

      key = trim(base_key) // "/conservation_tolerance"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), real_val, status)
      if (status == QFYAML_Success) dataset%conservation_tolerance = real_val

      ! Load conservative regridding options
      key = trim(base_key) // "/conservative_order"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), int_val, status)
      if (status == QFYAML_Success) dataset%conservative_order = int_val

      key = trim(base_key) // "/conservative_frac_area"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), logical_val, status)
      if (status == QFYAML_Success) dataset%conservative_frac_area = logical_val

      ! Load special handling options
      key = trim(base_key) // "/special_handling/enforce_positive"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), logical_val, status)
      if (status == QFYAML_Success) dataset%enforce_positive = logical_val

      key = trim(base_key) // "/special_handling/apply_floor_values"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), logical_val, status)
      if (status == QFYAML_Success) dataset%apply_floor_values = logical_val

      ! Load weight file options
      key = trim(base_key) // "/weight_file"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) then
         dataset%weight_file = trim(string_val)
         dataset%use_weight_file = .true.
      end if

      key = trim(base_key) // "/weight_file_template"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), string_val, status)
      if (status == QFYAML_Success) dataset%weight_file_template = trim(string_val)

      key = trim(base_key) // "/save_weights"
      call cc_yaml_get(cf_config%yaml_parser, trim(key), logical_val, status)
      if (status == QFYAML_Success) dataset%save_weights = logical_val

   end subroutine load_dataset_regrid_options

   ! Initialize a dataset
   !!
   !! \param dataset Dataset to initialize
   !! \param rc Return code
   subroutine initialize_dataset(dataset, rc)
      type(cf_dataset_type), intent(inout) :: dataset
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      type(ESMF_TimeInterval) :: timeInterval

      rc = ESMF_SUCCESS

      ! Initialize times
      call ESMF_TimeSet(dataset%last_refresh_time, yy=1900, mm=1, dd=1, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_TimeSet(dataset%next_refresh_time, yy=1900, mm=1, dd=1, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      dataset%needs_refresh = .true.
      dataset%ncid = -1
      dataset%file_open = .false.

   end subroutine initialize_dataset

   ! Initialize a field
   !!
   !! \param field Field to initialize
   !! \param rc Return code
   subroutine initialize_field(field, rc)
      type(cf_input_field_type), intent(inout) :: field
      integer, intent(out) :: rc

      ! Local variables
      integer :: status

      rc = ESMF_SUCCESS

      ! Initialize field state
      field%regrid_initialized = .false.
      field%varid = -1

      call ESMF_TimeSet(field%last_update_time, yy=1900, mm=1, dd=1, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Create ESMF fields will be done during first update when we know the grid

   end subroutine initialize_field

   ! Read grid information from file using ESMF I/O
   !!
   !! \param io_obj ESMF I/O object for the file
   !! \param dataset Dataset configuration
   !! \param file_grid Output grid object
   !! \param rc Return code
   subroutine read_grid_from_file_esmf(io_obj, dataset, file_grid, rc)
      type(ESMF_IO), intent(in) :: io_obj
      type(cf_dataset_type), intent(in) :: dataset
      type(ESMF_Grid), intent(out) :: file_grid
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      integer :: lon_size, lat_size, lev_size
      real(ESMF_KIND_R8), allocatable :: lon_data(:), lat_data(:)
      real(ESMF_KIND_R8), pointer :: grid_lon(:,:), grid_lat(:,:)
      integer :: i, j
      character(len=ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS

      ! Read dimension sizes using ESMF I/O
      call ESMF_IOGetDimSize(io_obj, dimName=trim(dataset%lon_name), dimSize=lon_size, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_IOGetDimSize(io_obj, dimName=trim(dataset%lat_name), dimSize=lat_size, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Create grid
      file_grid = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/), &
         maxIndex=(/lon_size,lat_size/), &
         regDecomp=(/1,pet_count/), &
         coordSys=ESMF_COORDSYS_SPH_DEG, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Add coordinates to grid
      call ESMF_GridAddCoord(file_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Read coordinate data
      allocate(lon_data(lon_size), lat_data(lat_size))

      ! Read longitude and latitude using ESMF I/O
      call ESMF_IORead(io_obj, variableName=trim(dataset%lon_name), array=lon_data, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_IORead(io_obj, variableName=trim(dataset%lat_name), array=lat_data, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Get coordinate pointers and populate grid
      call ESMF_GridGetCoord(file_grid, coordDim=1, localDe=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=grid_lon, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_GridGetCoord(file_grid, coordDim=2, localDe=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=grid_lat, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Populate coordinate arrays
      do j = 1, size(grid_lat, 2)
         do i = 1, size(grid_lon, 1)
            if (i <= lon_size) grid_lon(i,j) = lon_data(i)
            if (j <= lat_size) grid_lat(i,j) = lat_data(j)
         end do
      end do

      deallocate(lon_data, lat_data)

      if (local_pet == 0) then
         write(msg, '(A,I0,A,I0)') 'Created source grid from file: ', &
            lon_size, ' x ', lat_size
         call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=status)
      end if

   end subroutine read_grid_from_file_esmf

   ! Get dataset regrid options for a field
   !!
   !! \param field Field configuration
   !! \param dataset_opts Output dataset regrid options
   !! \param rc Return code
   subroutine get_dataset_regrid_options(field, dataset_opts, rc)
      type(cf_input_field_type), intent(in) :: field
      type(cf_dataset_type), intent(out) :: dataset_opts
      integer, intent(out) :: rc

      ! Local variables
      integer :: dataset_idx

      rc = ESMF_SUCCESS

      ! Find the dataset that contains this field
      call find_dataset_for_field(field, dataset_idx, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Copy dataset regrid options
      dataset_opts = cf_config%datasets(dataset_idx)

   end subroutine get_dataset_regrid_options

   ! Find dataset that contains a given field
   !!
   !! \param field Field configuration
   !! \param dataset_idx Output dataset index
   !! \param rc Return code
   subroutine find_dataset_for_field(field, dataset_idx, rc)
      type(cf_input_field_type), intent(in) :: field
      integer, intent(out) :: dataset_idx
      integer, intent(out) :: rc

      ! Local variables
      integer :: i, j

      rc = ESMF_SUCCESS
      dataset_idx = -1

      ! Search through all datasets for this field
      do i = 1, cf_config%n_datasets
         do j = 1, cf_config%datasets(i)%n_fields
            if (trim(cf_config%datasets(i)%fields(j)%variable_name) == &
               trim(field%variable_name)) then
               dataset_idx = i
               return
            end if
         end do
      end do

      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
         msg="No dataset found for field: " // trim(field%variable_name), &
         ESMF_CONTEXT, rcToReturn=rc)

   end subroutine find_dataset_for_field

   ! Resolve regrid method for a field (field-specific overrides dataset default)
   !!
   !! \param field Field configuration
   !! \param resolved_method Output resolved regrid method string
   !! \param rc Return code
   subroutine resolve_field_regrid_method(field, resolved_method, rc)
      type(cf_input_field_type), intent(in) :: field
      character(len=*), intent(out) :: resolved_method
      integer, intent(out) :: rc

      ! Local variables
      integer :: dataset_idx

      rc = ESMF_SUCCESS

      ! Check if field has its own regrid method specified
      if (len_trim(field%regrid_method) > 0) then
         resolved_method = trim(field%regrid_method)
         return
      end if

      ! Use dataset default
      call find_dataset_for_field(field, dataset_idx, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      resolved_method = trim(cf_config%datasets(dataset_idx)%regrid_method)

   end subroutine resolve_field_regrid_method

   ! Open dataset file using ESMF I/O (for compatibility with existing logic)
   !!
   !! \param dataset Dataset configuration
   !! \param filename Filename to open
   !! \param rc Return code
   subroutine open_dataset_file_esmf(dataset, filename, rc)
      type(cf_dataset_type), intent(inout) :: dataset
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      ! Local variables
      logical :: file_exists
      character(len=ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS

      ! Check if file exists
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
            msg="Input file not found: " // trim(filename), &
            ESMF_CONTEXT, rcToReturn=rc)
         return
      end if

      ! With ESMF I/O, we don't need to keep files open - they are opened per operation
      ! Just update the dataset state
      dataset%current_file = trim(filename)
      dataset%file_open = .true.
      dataset%ncid = 1  ! Dummy value for compatibility

      if (local_pet == 0) then
         write(msg, '(A,A)') 'Prepared dataset file for ESMF I/O: ', trim(filename)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
      end if

   end subroutine open_dataset_file_esmf

   ! Check if dataset needs refresh based on time and frequency
   !!
   !! \param dataset Dataset to check
   !! \param current_time Current simulation time
   !! \param needs_refresh Output refresh status
   !! \param rc Return code
   subroutine check_dataset_refresh_needed(dataset, current_time, needs_refresh, rc)
      type(cf_dataset_type), intent(in) :: dataset
      type(ESMF_Time), intent(in) :: current_time
      logical, intent(out) :: needs_refresh
      integer, intent(out) :: rc

      ! Local variables
      integer :: status
      character(len=ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS
      needs_refresh = .false.

      ! Check if we've passed the next refresh time
      if (current_time >= dataset%next_refresh_time) then
         needs_refresh = .true.
      end if

      ! Always refresh on first call
      if (dataset%needs_refresh) then
         needs_refresh = .true.
      end if

      if (local_pet == 0 .and. needs_refresh) then
         write(msg, '(A,A)') 'Dataset refresh needed for: ', trim(dataset%name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=status)
      end if

   end subroutine check_dataset_refresh_needed

   ! Update dataset file path based on current time
   !!
   !! \param dataset Dataset to update
   !! \param current_time Current simulation time
   !! \param rc Return code
   subroutine update_dataset_file_path(dataset, current_time, rc)
      type(cf_dataset_type), intent(inout) :: dataset
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      ! Local variables
      character(len=MAX_STRING_LEN) :: new_filename
      character(len=ESMF_MAXSTR) :: msg
      integer :: status

      rc = ESMF_SUCCESS

      ! Generate filename from template
      call generate_filename(dataset%file_template, current_time, new_filename, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Update dataset if filename changed
      if (trim(new_filename) /= trim(dataset%current_file)) then
         if (local_pet == 0) then
            write(msg, '(A,A,A,A)') 'Updating dataset file from ', &
               trim(dataset%current_file), ' to ', trim(new_filename)
            call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=status)
         end if

         ! Close old file and open new one
         call cleanup_dataset(dataset, status)
         call open_dataset_file_esmf(dataset, new_filename, rc)
         if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      end if

   end subroutine update_dataset_file_path

   ! Check if dataset needs refresh (wrapper for check_dataset_refresh_needed)
   !!
   !! \param dataset Dataset to check
   !! \param current_time Current simulation time
   !! \param needs_update Output refresh status
   !! \param rc Return code
   subroutine check_dataset_refresh(dataset, current_time, needs_update, rc)
      type(cf_dataset_type), intent(in) :: dataset
      type(ESMF_Time), intent(in) :: current_time
      logical, intent(out) :: needs_update
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      ! Use the existing helper function
      call check_dataset_refresh_needed(dataset, current_time, needs_update, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

   end subroutine check_dataset_refresh

   ! Refresh dataset by updating file path and reopening if needed
   !!
   !! \param dataset Dataset to refresh
   !! \param current_time Current simulation time
   !! \param rc Return code
   subroutine refresh_dataset(dataset, current_time, rc)
      type(cf_dataset_type), intent(inout) :: dataset
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_TimeInterval) :: refresh_interval
      integer :: status
      character(len=ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS

      ! Update file path based on current time
      call update_dataset_file_path(dataset, current_time, rc)
      if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      ! Update refresh timing
      call ESMF_TimeIntervalSet(refresh_interval, &
         h=dataset%refresh_freq_hours, rc=status)
      if (ESMF_LogFoundError(status, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

      dataset%next_refresh_time = current_time + refresh_interval
      dataset%needs_refresh = .false.

      if (local_pet == 0) then
         write(msg, '(A,A)') 'Refreshed dataset: ', trim(dataset%name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=status)
      end if

   end subroutine refresh_dataset

end module catchem_nuopc_cf_input
