!> \file catchem_nuopc_interface.F90
!! \brief NUOPC interface for CATChem atmospheric chemistry model
!!
!! \details
!! This module provides the core interface routines for the CATChem NUOPC cap,
!! handling data transformations between ESMF/NUOPC fields and CATChem states.
!! It follows similar patterns as the CCPP interface but is adapted for NUOPC
!! framework requirements including ESMF grids, fields, and parallel operations.
!!
!! Key functionalities include:
!! - CATChem model initialization and finalization within NUOPC framework
!! - Data transformation between ESMF fields and CATChem state objects
!! - Field mapping configuration and management
!! - Integration with CF-compliant input and NetCDF output systems
!! - Support for flexible grid configurations and parallel decomposition
!! - Error handling and logging for NUOPC applications
!!
!! The interface supports both sequential and parallel execution modes
!! and provides standardized data exchange capabilities for coupling
!! with other Earth system model components.
!!
!! \author Barry Baker, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

module catchem_nuopc_interface

   use iso_c_binding, only: c_loc, c_null_char, c_null_ptr, c_char, c_double, c_ptr, c_int, c_long_long, c_associated, c_f_pointer, c_intptr_t
   use ESMF
   use NUOPC
   use MPI
   use CATChem_API, only: CATChem_Model
   use catchem_bridge_precision, only: fp, is_exact_zero
   use catchem_bridge_constants, only: g0, Rd, Re, AIRMW
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use catchem_bridge_error, only : CC_SUCCESS, CC_FAILURE
   use catchem_bridge_error, only: ErrorManagerType
   use catchem_nuopc_emis_data_mod, only: ExtEmisDataType, ExtEmisFieldType  ! External emissions data types
   use aqmio, only: AQMIO_Create, AQMIO_Destroy, AQMIO_Write, AQMIO_Close, AQMIO_Write1D, AQMIO_FMT_NETCDF, &
      AQMIO_LatlonInit, AQMIO_LatlonCleanup
   use catchem_latlon_output_mod, only: latlon_diag_set_time, latlon_diag_is_init
   use catchem_nuopc_emis_mod

   implicit none

   integer, parameter :: DIAG_REAL_SCALAR = 0, DIAG_REAL_1D = 1, DIAG_REAL_2D = 2, DIAG_REAL_3D = 3
   integer, parameter :: TRACER_HOST_OWNED = 0, TRACER_CHEMICAL = 1, TRACER_DIAGNOSTIC = 2

   interface
      integer(c_int) function catchem_state_get_pointer_3d_checked(state_ptr, name, ptr_out) &
         bind(C, name="catchem_state_get_pointer_3d_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
         type(c_ptr), intent(out) :: ptr_out
      end function

      integer(c_int) function catchem_state_get_species_count_checked(state_ptr, count_out) &
         bind(C, name="catchem_state_get_species_count_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), intent(out) :: count_out
      end function

      integer(c_int) function catchem_state_get_species_index_checked(state_ptr, name, index_out) &
         bind(C, name="catchem_state_get_species_index_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: index_out
      end function

      integer(c_int) function catchem_state_get_species_conc_pointer_checked(state_ptr, species_index, dim1, dim2, ptr_out) &
         bind(C, name="catchem_state_get_species_conc_pointer_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: species_index, dim1, dim2
         type(c_ptr), intent(out) :: ptr_out
      end function

      integer(c_int) function catchem_diag_get_pointer_checked(core_ptr, name, rank, dims, ptr_out) &
         bind(C, name="catchem_diag_get_pointer_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: rank
         integer(c_int), intent(in) :: dims(*)
         type(c_ptr), intent(out) :: ptr_out
      end function

      ! Diagnostic enumeration used by write_process_diagnostics to discover
      ! the fields each process registered in the C++ DiagnosticManager.
      integer(c_int) function catchem_diag_get_count_checked(core_ptr, count_out) &
         bind(C, name="catchem_diag_get_count_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), intent(out) :: count_out
      end function

      integer(c_int) function catchem_diag_get_name_at_checked(core_ptr, index, name_out, name_length) &
         bind(C, name="catchem_diag_get_name_at_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index, name_length
         character(kind=c_char), intent(out) :: name_out(*)
      end function

      integer(c_int) function catchem_diag_get_rank_checked(core_ptr, name, rank_out) &
         bind(C, name="catchem_diag_get_rank_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: rank_out
      end function

      integer(c_int) function catchem_diag_get_dims_checked(core_ptr, name, dims_out, dims_length) &
         bind(C, name="catchem_diag_get_dims_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: dims_out(*)
         integer(c_int), value :: dims_length
      end function

      integer(c_int) function catchem_diag_get_units_checked(core_ptr, name, units_out, units_length) &
         bind(C, name="catchem_diag_get_units_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         character(kind=c_char), intent(out) :: units_out(*)
         integer(c_int), value :: units_length
      end function

      integer(c_int) function catchem_diag_get_description_checked(core_ptr, name, desc_out, desc_length) &
         bind(C, name="catchem_diag_get_description_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         character(kind=c_char), intent(out) :: desc_out(*)
         integer(c_int), value :: desc_length
      end function

      subroutine catchem_diag_sync_to_host(core_ptr) &
         bind(C, name="catchem_diag_sync_to_host")
         import :: c_ptr
         type(c_ptr), value :: core_ptr
      end subroutine

      ! Species metadata predicates (macro-generated in catchem_api_config.cpp;
      ! they return the flag directly and yield 0 for an invalid index).
      integer(c_int) function catchem_state_is_species_dust(state_ptr, index) &
         bind(C, name="catchem_state_is_species_dust")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_is_species_seasalt(state_ptr, index) &
         bind(C, name="catchem_state_is_species_seasalt")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_name_at_checked(state_ptr, index, name_out, name_length) &
         bind(C, name="catchem_state_get_species_name_at_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index, name_length
         character(kind=c_char), intent(out) :: name_out(*)
      end function

      integer(c_int) function catchem_state_is_species_gas_checked(state_ptr, index, value_out) &
         bind(C, name="catchem_state_is_species_gas_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         integer(c_int), intent(out) :: value_out
      end function

      integer(c_int) function catchem_state_get_species_mw_checked(state_ptr, index, molecular_weight_out) &
         bind(C, name="catchem_state_get_species_mw_checked")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         real(c_double), intent(out) :: molecular_weight_out
      end function

      integer(c_int) function catchem_state_is_species_aerosol_checked(state_ptr, index, value_out) &
         bind(C, name="catchem_state_is_species_aerosol_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
         integer(c_int), intent(out) :: value_out
      end function

      integer(c_int) function catchem_state_set_time_checked(state_ptr, yr, mo, dy, hr, mn, sc, doy, tstep) &
         bind(C, name="catchem_state_set_time_checked")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: yr, mo, dy, hr, mn, sc, doy
         real(c_double), value :: tstep
      end function

      integer(c_int) function catchem_state_derive_airden_dry_checked(state_ptr) &
         bind(C, name="catchem_state_derive_airden_dry_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      integer(c_int) function catchem_state_derive_bxheight_checked(state_ptr) &
         bind(C, name="catchem_state_derive_bxheight_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      integer(c_int) function catchem_state_begin_import_generation(state_ptr) &
         bind(C, name="catchem_state_begin_import_generation")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      integer(c_int) function catchem_config_get_yaml_bool(core_ptr, yaml_path, default_val) &
         bind(C, name="catchem_config_get_yaml_bool")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: yaml_path(*)
         integer(c_int), value :: default_val
      end function catchem_config_get_yaml_bool

      subroutine catchem_set_config_echo_enabled(enabled) &
         bind(C, name="catchem_set_config_echo_enabled")
         import :: c_int
         integer(c_int), value :: enabled
      end subroutine catchem_set_config_echo_enabled

      integer(c_int) function catchem_core_get_timestep_outcome(core_ptr, status, timestep, duration, &
         import_generation, process_index, state_classification, process_name, process_name_len, cause, cause_len) &
         bind(C, name="catchem_core_get_timestep_outcome")
         import :: c_ptr, c_int, c_long_long, c_double, c_char
         type(c_ptr), value :: core_ptr
         integer(c_int), intent(out) :: status, process_index, state_classification
         integer(c_long_long), intent(out) :: timestep, import_generation
         real(c_double), intent(out) :: duration
         character(kind=c_char), intent(out) :: process_name(*), cause(*)
         integer(c_int), value :: process_name_len, cause_len
      end function catchem_core_get_timestep_outcome

      ! Emission-config query API used to resolve a MET_* map target back to the
      ! emission field that declares it (config is the source of truth).  These
      ! bind to the same global C symbols declared privately in
      ! catchem_nuopc_emis_mod; redeclaring here keeps this module self-contained.
      integer(c_int) function catchem_config_has_emission_mapping(core_ptr) &
         bind(C, name="catchem_config_has_emission_mapping")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function catchem_config_has_emission_mapping

      integer(c_int) function catchem_config_get_emission_category_count(core_ptr) &
         bind(C, name="catchem_config_get_emission_category_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function catchem_config_get_emission_category_count

      subroutine catchem_config_get_emission_category_name_at(core_ptr, index, name_out, max_len) &
         bind(C, name="catchem_config_get_emission_category_name_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
         integer(c_int), value :: max_len
      end subroutine catchem_config_get_emission_category_name_at

      integer(c_int) function catchem_config_get_emission_field_count(core_ptr, category_name) &
         bind(C, name="catchem_config_get_emission_field_count")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
      end function catchem_config_get_emission_field_count

      subroutine catchem_config_get_emission_field_name_at(core_ptr, category_name, field_idx, name_out, max_len) &
         bind(C, name="catchem_config_get_emission_field_name_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
         integer(c_int), value :: field_idx
         character(kind=c_char), intent(out) :: name_out(*)
         integer(c_int), value :: max_len
      end subroutine catchem_config_get_emission_field_name_at

      integer(c_int) function catchem_config_get_emission_species_map_count(core_ptr, category_name, field_name) &
         bind(C, name="catchem_config_get_emission_species_map_count")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
         character(kind=c_char), intent(in) :: field_name(*)
      end function catchem_config_get_emission_species_map_count

      subroutine catchem_config_get_emission_species_map_at(core_ptr, category_name, field_name, map_idx, &
         target_species_out, max_len, scale_out, species_idx_out) &
         bind(C, name="catchem_config_get_emission_species_map_at")
         import :: c_ptr, c_char, c_int, c_double
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
         character(kind=c_char), intent(in) :: field_name(*)
         integer(c_int), value :: map_idx
         character(kind=c_char), intent(out) :: target_species_out(*)
         integer(c_int), value :: max_len
         real(c_double), intent(out) :: scale_out
         integer(c_int), intent(out) :: species_idx_out
      end subroutine catchem_config_get_emission_species_map_at
   end interface

   private

   public :: catchem_nuopc_init
   public :: catchem_nuopc_run
   public :: catchem_nuopc_finalize
   public :: transform_nuopc_to_catchem
   public :: transform_catchem_to_nuopc
   public :: load_field_config
   public :: catchem_diagnostics_write
   public :: write_chem_diagnostics
   !public :: get_cc_wrap  ! Accessor for process-local wrapper
   public :: get_n_import_fields, get_import_field_info  ! Safe field_config access
   public :: get_n_export_fields, get_export_field_info  ! Safe field_config access
   public :: update_pm_diagnostics  ! Exposed for the NUOPC transform test harness
   public :: TRACER_HOST_OWNED, TRACER_CHEMICAL, TRACER_DIAGNOSTIC  ! Tracer-map contract values for the test harness
   public :: catchem_nuopc_get_physical_validation_report

   !> \brief Field mapping configuration structure
   !!
   !! Defines the mapping between NUOPC standard names and CATChem variables,
   !! including metadata for proper data transformation and validation.
   !! \{
   type :: field_mapping_type
      character(len=128) :: standard_name !< NUOPC/CF standard field name
      character(len=128) :: catchem_var   !< Corresponding CATChem variable path
      integer :: dimensions               !< Number of spatial dimensions (2D/3D)
      character(len=64) :: units          !< Physical units for conversion
      character(len=32) :: vertical_axis = 'level' !< Vertical coordinate semantic for 3D fields
      character(len=128) :: host_tracer_name = '' !< Host tracer to expose as a CATChem met field
      character(len=128) :: host_tracer_var = ''  !< CATChem met name for host_tracer_name
      logical :: optional = .false.       !< Whether field is required or optional
      logical :: advertise = .false.      !< Advertise an optional field for a host coupling contract
   end type field_mapping_type
   !! \}

   !> \brief Field configuration structure for NUOPC interface
   !!
   !! Contains the complete field mapping configuration including both
   !! import and export field definitions with associated metadata.
   !! \{
   type :: field_config_type
      integer :: n_import_fields = 0                         !< Number of import fields
      integer :: n_export_fields = 0                         !< Number of export fields
      type(field_mapping_type), allocatable :: import_fields(:) !< Import field mapping array
      type(field_mapping_type), allocatable :: export_fields(:) !< Export field mapping array
   end type field_config_type

   !> field_config can be shared among MPIs
   type(field_config_type), public :: field_config

   !> \brief tracer mapping between CATChem and NUOPC
   !!
   !! Defines the tracer mapping index between NUOPC and chem_state if CATChem
   !! \{
   type :: tracer_index_map
      integer, allocatable :: nuopc_to_cc(:)  !< mapping index from NUOPC to CATChem
      integer, allocatable :: entry_kind(:)   !< 1=CATChem species, 2=diagnostic pseudo-tracer, 0=host-owned
      character(len=128), allocatable :: names(:) !< NUOPC tracer name
      character(len=128), allocatable :: units(:) !< NUOPC tracer unit
      real(c_double), allocatable :: host_to_catchem(:) !< kg kg-1 to CATChem native unit
      real(c_double), allocatable :: catchem_to_host(:) !< CATChem native unit to kg kg-1
   end type tracer_index_map
   !! \}

   type :: buffer_2d_type
      real(c_double), allocatable :: data(:,:)
   end type buffer_2d_type

   type :: buffer_3d_type
      real(c_double), allocatable :: data(:,:,:)
   end type buffer_3d_type

   !> Container for process-private CATChem state to avoid MPI sharing
   type :: cc_wrap_type
      type(CATChem_Model) :: catchem_model
      type(ExtEmisDataType) :: ext_emis ! External emissions data object
      type(field_config_type) :: field_config  ! Moved from module level for MPI safety
      type(tracer_index_map) :: tracer_map
      type(ESMF_Grid) :: grid
      type(buffer_2d_type), allocatable :: met_buf_2d(:)
      type(buffer_3d_type), allocatable :: met_buf_3d(:)
      real(c_double), allocatable :: chem_buf_4d(:,:,:,:)
      real(c_double), allocatable :: host_tracer_buf_4d(:,:,:,:)
      real(c_double), allocatable :: lat(:,:)
      real(c_double), allocatable :: lon(:,:)
      real(c_double), allocatable :: area_m2(:,:)
      real(c_double), allocatable :: z0_m(:,:)
      real(c_double), allocatable :: dust_clayfrac(:,:)
      real(c_double), allocatable :: dust_sandfrac(:,:)
      real(c_double), allocatable :: dust_ssm(:,:)
      real(c_double), allocatable :: dust_rdrag(:,:)
      real(c_double), allocatable :: dust_ustar_threshold(:,:)
      logical :: initialized = .false.
      logical :: verbose_logging = .false. !< Runtime YAML switch: simulation/verbose/activate
      ! Diagnostic output variables (moved from module level for MPI safety)
      type(ESMF_Time) :: last_output_time
      type(ESMF_Time) :: startTime
      type(ESMF_Time) :: endTime
      type(ESMF_TimeInterval) :: output_interval
      type(ESMF_TimeInterval) :: timeStep
      logical :: output_timing_initialized = .false.
      integer :: timestep_counter = 0  !< Per-component run counter; never shared across CATChem instances.
      character(len=256) :: output_directory = './output'
      character(len=64) :: output_prefix = 'catchem_diag'
      integer :: output_frequency = 3600  ! Default: 1 hour in seconds
      integer :: compress_lev = 0         !< Compression level for output NC files (0-9)
      type(ESMF_GridComp) :: iocomp
      ! Time slice tracking for NetCDF output
      integer :: current_time_slice = 0
      logical :: pm_diag_registered = .false.  !< Track PM diagnostic registration per-instance
   end type cc_wrap_type

   type CATChem_InternalState
      type(cc_wrap_type), pointer :: wrap => null()
   end type CATChem_InternalState
   public :: cc_wrap_type, CATChem_InternalState


contains

   pure function lowercase(value) result(lower)
      character(len=*), intent(in) :: value
      character(len=len(value)) :: lower
      integer :: p, code
      lower = value
      do p = 1, len(value)
         code = iachar(lower(p:p))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(p:p) = achar(code + 32)
      end do
   end function lowercase

   pure logical function is_mass_mixing_ratio_unit(unit)
      character(len=*), intent(in) :: unit
      character(len=128) :: normalized
      normalized = trim(adjustl(lowercase(unit)))
      is_mass_mixing_ratio_unit = normalized == 'kg kg-1' .or. normalized == 'kg/kg' .or. &
         normalized == 'kg kg^-1' .or. normalized == 'kg kg**-1'
   end function is_mass_mixing_ratio_unit

   pure logical function is_ppm_unit(unit)
      character(len=*), intent(in) :: unit
      character(len=128) :: normalized
      normalized = trim(adjustl(lowercase(unit)))
      is_ppm_unit = normalized == 'ppm'
   end function is_ppm_unit

   pure logical function is_micro_mass_mixing_ratio_unit(unit)
      character(len=*), intent(in) :: unit
      character(len=128) :: normalized
      normalized = trim(adjustl(lowercase(unit)))
      is_micro_mass_mixing_ratio_unit = normalized == 'ug kg-1' .or. normalized == 'ug/kg' .or. &
         normalized == 'ug kg^-1' .or. normalized == 'ug kg**-1'
   end function is_micro_mass_mixing_ratio_unit

   ! PM2.5 and PM10 are diagnostic mass concentrations, unlike the
   ! transported aerosol tracers, which are mass mixing ratios.
   pure logical function is_micro_mass_concentration_unit(unit)
      character(len=*), intent(in) :: unit
      character(len=128) :: normalized
      normalized = trim(adjustl(lowercase(unit)))
      is_micro_mass_concentration_unit = normalized == 'ug m-3' .or. normalized == 'ug/m3' .or. &
         normalized == 'ug m^-3' .or. normalized == 'ug m**-3'
   end function is_micro_mass_concentration_unit

   pure logical function is_pm_diagnostic_name(name)
      character(len=*), intent(in) :: name
      character(len=128) :: normalized
      normalized = trim(adjustl(lowercase(name)))
      is_pm_diagnostic_name = normalized == 'pm25' .or. normalized == 'pm10'
   end function is_pm_diagnostic_name

   subroutine catchem_c_string_to_fortran(c_value, value)
      character(kind=c_char), intent(in) :: c_value(*)
      character(len=*), intent(out) :: value
      integer :: i

      value = ''
      do i = 1, len(value)
         if (c_value(i) == c_null_char) exit
         value(i:i) = c_value(i)
      end do
   end subroutine catchem_c_string_to_fortran

   subroutine catchem_append_timestep_failure(core_ptr, errmsg)
      type(c_ptr), intent(in) :: core_ptr
      character(len=*), intent(inout) :: errmsg
      character(kind=c_char) :: c_process(128), c_cause(512)
      character(len=128) :: process_name
      character(len=512) :: cause
      integer(c_int) :: status, process_index, state_classification, outcome_rc
      integer(c_long_long) :: timestep, import_generation
      real(c_double) :: duration

      if (.not. c_associated(core_ptr)) return
      outcome_rc = catchem_core_get_timestep_outcome(core_ptr, status, timestep, duration, import_generation, &
         process_index, state_classification, c_process, int(size(c_process), c_int), c_cause, int(size(c_cause), c_int))
      if (outcome_rc /= 0_c_int) return
      call catchem_c_string_to_fortran(c_process, process_name)
      call catchem_c_string_to_fortran(c_cause, cause)
      if (len_trim(process_name) > 0 .or. len_trim(cause) > 0) then
         errmsg = trim(errmsg) // ': process=' // trim(process_name) // ' cause=' // trim(cause)
      end if
   end subroutine catchem_append_timestep_failure

   !> Emit a low-volume, per-PET run-phase marker when enabled in runtime YAML.
   !! The final marker in an ESMF PET log identifies the phase that stalled.
   subroutine catchem_log_run_phase(cc_wrap, step, phase)
      type(cc_wrap_type), intent(in) :: cc_wrap
      integer, intent(in) :: step
      character(len=*), intent(in) :: phase
      character(len=256) :: message
      integer :: localrc

      if (.not. cc_wrap%verbose_logging) return
      write(message, '(A,I0,A,A)') 'CATChem run step=', step, ' phase=', trim(phase)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO, rc=localrc)
   end subroutine catchem_log_run_phase

   subroutine catchem_nuopc_get_physical_validation_report(cc_wrap, issue_count, detail, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(out) :: issue_count
      character(len=*), intent(out) :: detail
      integer, intent(out) :: rc

      call cc_wrap%catchem_model%get_physical_validation_report(issue_count, detail, rc)
   end subroutine catchem_nuopc_get_physical_validation_report

   !> Initialize CATChem model for NUOPC interface
   !!
   !! This routine performs comprehensive initialization of the CATChem model
   !! within the NUOPC framework, including configuration loading, memory
   !! allocation, and setup of I/O systems for external data and diagnostics.
   !!
   !! @param config CATChem configuration object to initialize
   !! @param catchem_states Main container for all CATChem state variables
   !! @param dustState Dust emission and transport state object
   !! @param seaSaltState Sea salt emission and transport state object
   !! @param dryDepState Dry deposition process state object
   !! @param im Number of horizontal grid points
   !! @param config_file Path to CATChem configuration file
   !! @param grid ESMF grid for regridding and I/O operations
   !! @param errflg Error flag (CC_SUCCESS on success)
   !! @param errmsg Error message string if errflg indicates failure
   !!
   !! This routine performs the following initialization steps:
   !! - Reads and validates the CATChem configuration file
   !! - Initializes all CATChem state objects and process modules
   !! - Allocates memory for MetState, EmisState, ChemState, and DiagState arrays
   !! - Sets up field mapping configuration for NUOPC data exchange
   !! - Initializes CF-compliant input system for external data
   !! - Initializes NetCDF output system for diagnostic data
   !! - Validates grid compatibility and spatial dimensions
   !!
   !! The routine is similar to catchem_init in the CCPP interface but includes
   !! additional NUOPC-specific functionality such as ESMF grid integration
   !! and YAML-based field configuration management.
   !!
   !! @note This routine must be called before any CATChem calculations
   !!       and requires valid ESMF grid and configuration files
   !!
   !! @warning Proper error checking should be performed on errflg after calling
   subroutine catchem_nuopc_init(model, config_file, lat, lon, nlev, tracerinfo, input_grid, startTime,stopTime, timeStep, clock, nsoil, nsoiltype, nsurftype, rc)

      type(ESMF_GridComp)  :: model
      character(len=*), intent(in) :: config_file
      real(ESMF_KIND_R8), dimension(:,:), intent(in) :: lat
      real(ESMF_KIND_R8), dimension(:,:), intent(in) :: lon
      integer, intent(in) :: nlev
      ! Absent on the standalone path (no coupling partner advertises a tracer
      ! field); required in practice for coupled runs, enforced below.
      type(ESMF_Info), intent(in), optional :: tracerinfo
      type(ESMF_Grid), intent(in) :: input_grid
      type(ESMF_Time), intent(in), optional :: startTime,stopTime
      type(ESMF_TimeInterval), intent(in), optional :: timeStep
      type(ESMF_Clock), intent(in), optional :: clock
      integer, intent(in), optional :: nsoil, nsoiltype, nsurftype
      integer, intent(out) :: rc

      ! Local variables
      integer :: nx, ny, num_processes, stat, i, j
      integer(c_int) :: catchem_status, species_index, is_gas, is_aerosol
      real(c_double) :: molecular_weight, conversion_factor
      integer(ESMF_KIND_I8) :: tstep_seconds
      character(len=128), allocatable :: tracer_names(:) !< NUOPC tracer name
      character(len=128), allocatable :: tracer_units(:) !< NUOPC tracer unit
      type(CATChem_InternalState) :: is
      type(CATChem_InternalState) :: verify_is
      type(cc_wrap_type), pointer:: cc_wrap
      integer :: verify_rc
      type(ESMF_VM) :: vm
      integer :: localPet, vmrc

      ! Initialize
      rc = CC_SUCCESS

      ! -- allocate memory for the internal state and store it into component
      allocate(is%wrap, stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return !bail out

      cc_wrap => is%wrap

      !get nx, ny
      nx = size(lat, 1)
      ny = size(lat, 2)

      ! Every PET parses the same YAML, so echo the effective configuration to
      ! stdout only on the root PET; the C++ core prints it during the
      ! initialize() call below (see catchem_set_config_echo_enabled).
      call ESMF_VMGetCurrent(vm, rc=vmrc)
      if (vmrc == ESMF_SUCCESS) then
         call ESMF_VMGet(vm, localPet=localPet, rc=vmrc)
      end if
      if (vmrc /= ESMF_SUCCESS) localPet = 0  ! single-PET or query failed: print
      if (localPet /= 0) call catchem_set_config_echo_enabled(0_c_int)

      ! Initialize catchem using process-local variable
      if (present(nsoil) .and. present(nsoiltype) .and. present(nsurftype)) then
         call cc_wrap%catchem_model%initialize(config_file, nx, ny, nlev, nsoil, nsoiltype, nsurftype, rc)
      else
         call cc_wrap%catchem_model%initialize(config_file, nx, ny, nlev, rc = rc)
      end if

      if (rc /= CC_SUCCESS) then
         call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="CATChem initialization failed: "//trim(cc_wrap%catchem_model%last_error), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return  ! bail out
      end if

      cc_wrap%verbose_logging = catchem_config_get_yaml_bool( &
         cc_wrap%catchem_model%cpp_core_ptr, 'simulation/verbose/activate' // c_null_char, 0_c_int) /= 0

      !assign lat and lon directly to C++ StateManager persistently in cc_wrap
      allocate(cc_wrap%lat(nx, ny))
      allocate(cc_wrap%lon(nx, ny))
      cc_wrap%lat = real(lat, c_double)
      cc_wrap%lon = real(lon, c_double)
      where (cc_wrap%lon > 180.0_c_double)
         cc_wrap%lon = cc_wrap%lon - 360.0_c_double
      end where
      call cc_wrap%catchem_model%bind_met_2d("LAT", cc_wrap%lat, rc)
      if (rc /= CC_SUCCESS) return
      call cc_wrap%catchem_model%bind_met_2d("LON", cc_wrap%lon, rc)
      if (rc /= CC_SUCCESS) return

      ! Populate grid-cell areas [m2] used for point-source emissions
      allocate(cc_wrap%area_m2(nx, ny))
      cc_wrap%area_m2 = 0.0_c_double
      block
         type(ESMF_Field) :: areaField
         real(ESMF_KIND_R8), pointer :: areaPtr(:,:)
         integer :: arc
         logical :: areaIsPresent
         nullify(areaPtr)
         areaIsPresent = .false.
         call ESMF_GridGetItem(input_grid, itemflag=ESMF_GRIDITEM_AREA, &
            staggerloc=ESMF_STAGGERLOC_CENTER, isPresent=areaIsPresent, rc=arc)
         if (arc == ESMF_SUCCESS .and. areaIsPresent) then
            call ESMF_GridGetItem(input_grid, itemflag=ESMF_GRIDITEM_AREA, &
               staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=areaPtr, rc=arc)
         else
            arc = ESMF_RC_NOT_FOUND
         end if
         if (arc == ESMF_SUCCESS .and. associated(areaPtr)) then
            if (size(areaPtr,1) == nx .and. size(areaPtr,2) == ny) then
               cc_wrap%area_m2 = real(areaPtr, c_double)
            end if
         else
            ! Fall back to computing cell areas from the grid geometry.
            nullify(areaPtr)
            areaField = ESMF_FieldCreate(input_grid, typekind=ESMF_TYPEKIND_R8, &
               staggerloc=ESMF_STAGGERLOC_CENTER, rc=arc)
            if (arc == ESMF_SUCCESS) call ESMF_FieldRegridGetArea(areaField, rc=arc)
            if (arc == ESMF_SUCCESS) call ESMF_FieldGet(areaField, farrayPtr=areaPtr, rc=arc)
            if (arc == ESMF_SUCCESS .and. associated(areaPtr)) then
               if (size(areaPtr,1) == nx .and. size(areaPtr,2) == ny) then
                  cc_wrap%area_m2 = real(areaPtr, c_double) * real(Re * Re, c_double)
               end if
            else
               call ESMF_LogWrite('catchem_nuopc_init: could not determine grid-cell '// &
                  'areas; AREA_M2 left unset (point emissions will be skipped)', &
                  ESMF_LOGMSG_WARNING, rc=arc)
            end if
            call ESMF_FieldDestroy(areaField, rc=arc)
         end if
         call cc_wrap%catchem_model%bind_met_2d("AREA_M2", cc_wrap%area_m2, rc)
         if (rc /= CC_SUCCESS) return
      end block

      !initialize extemission data here
      call catchem_emis_init(cc_wrap%ext_emis, cc_wrap%catchem_model%cpp_core_ptr, nx, ny, nlev, clock, rc)

      !get output information from config
      cc_wrap%output_frequency = cc_wrap%catchem_model%get_output_frequency()
      cc_wrap%compress_lev = cc_wrap%catchem_model%get_compress_level()
      call cc_wrap%catchem_model%get_output_directory(cc_wrap%output_directory)
      call cc_wrap%catchem_model%get_output_prefix(cc_wrap%output_prefix)

      !populate tracer mapping using process-local tracer_map
      if (present(tracerinfo)) then
         call TracerInfoGet(tracerinfo, 'tracerNames', tracer_names, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__)) return  ! bail out

         if (.not.allocated(tracer_names)) then
            call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
               msg="CATChem requires tracerNames metadata on its rank-4 host tracer field", &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            return
         end if

         ! - import tracer units if available
         call TracerInfoGet(tracerinfo, 'tracerUnits', tracer_units, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__)) return  ! bail out

         if (.not.allocated(tracer_units)) then
            allocate(tracer_units(size(tracer_names)), stat=stat)
            if (ESMF_LogFoundAllocError(statusToCheck=stat, &
               msg="Unable to allocate internal workspace", &
               line=__LINE__,  file=__FILE__)) return  ! bail out
            tracer_units = 'n/a'
         else if (size(tracer_units) /= size(tracer_names)) then
            call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
               msg='CATChem tracerUnits length does not match tracerNames length', &
               line=__LINE__, file=__FILE__, rcToReturn=rc)
            return
         end if
      else
         ! Standalone (no coupling partner): there is no imported tracer list to
         ! map onto. Build empty name/unit lists so the (unused) NUOPC<->CATChem
         ! tracer map is well-defined. A single CATChem model owns its species
         ! through its YAML configuration rather than through coupling.
         allocate(tracer_names(0), stat=stat)
         if (ESMF_LogFoundAllocError(statusToCheck=stat, &
            msg="Unable to allocate empty tracer name list", &
            line=__LINE__,  file=__FILE__)) return  ! bail out
         allocate(tracer_units(0), stat=stat)
         if (ESMF_LogFoundAllocError(statusToCheck=stat, &
            msg="Unable to allocate empty tracer unit list", &
            line=__LINE__,  file=__FILE__)) return  ! bail out
      end if

      !copy to cc_wrap
      allocate(cc_wrap%tracer_map%names, stat=stat, source=tracer_names)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Unable to allocate tracers name", &
         line=__LINE__,  file=__FILE__, rcToReturn=rc)) return  ! bail out

      allocate(cc_wrap%tracer_map%units, stat=stat, source=tracer_units)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Unable to allocate tracers unit", &
         line=__LINE__,  file=__FILE__, rcToReturn=rc)) return  ! bail out

      allocate(cc_wrap%tracer_map%nuopc_to_cc(size(tracer_names)), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Unable to allocate nuopc_to_cc mapping", &
         line=__LINE__,  file=__FILE__, rcToReturn=rc)) return  ! bail out
      allocate(cc_wrap%tracer_map%entry_kind(size(tracer_names)), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Unable to allocate tracer classifications", &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      allocate(cc_wrap%tracer_map%host_to_catchem(size(tracer_names)), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Unable to allocate tracer import factors", &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      allocate(cc_wrap%tracer_map%catchem_to_host(size(tracer_names)), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Unable to allocate tracer export factors", &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      cc_wrap%tracer_map%host_to_catchem = 1.0_c_double
      cc_wrap%tracer_map%catchem_to_host = 1.0_c_double

      ! assign mapping index directly from C++ StateManager
      do i = 1, size(cc_wrap%tracer_map%names)
         catchem_status = catchem_state_get_species_index_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, trim(cc_wrap%tracer_map%names(i)) // c_null_char, species_index)
         if (catchem_status == 0_c_int) then
            cc_wrap%tracer_map%nuopc_to_cc(i) = int(species_index)
         else
            cc_wrap%tracer_map%nuopc_to_cc(i) = 0
         end if
         if (cc_wrap%tracer_map%nuopc_to_cc(i) > 0) then
            cc_wrap%tracer_map%entry_kind(i) = TRACER_CHEMICAL
            is_gas = 0_c_int
            is_aerosol = 0_c_int
            catchem_status = catchem_state_is_species_gas_checked(cc_wrap%catchem_model%state_mgr_ptr, &
               species_index, is_gas)
            if (catchem_status == 0_c_int) catchem_status = catchem_state_is_species_aerosol_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, species_index, is_aerosol)
            if (catchem_status /= 0_c_int) then
               call ESMF_LogWrite('Unable to classify CATChem tracer: ' // trim(cc_wrap%tracer_map%names(i)), &
                  ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
            if ((is_gas == 0_c_int .and. is_aerosol == 0_c_int) .or. &
               (is_gas /= 0_c_int .and. is_aerosol /= 0_c_int)) then
               call ESMF_LogWrite('CATChem tracer must be exactly one of gas or aerosol: ' // &
                  trim(cc_wrap%tracer_map%names(i)), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
            if (is_gas /= 0_c_int) then
               if (.not. is_mass_mixing_ratio_unit(cc_wrap%tracer_map%units(i)) .and. &
                  .not. is_ppm_unit(cc_wrap%tracer_map%units(i))) then
                  call ESMF_LogWrite('Unsupported units for CATChem gas tracer ' // trim(cc_wrap%tracer_map%names(i)) // &
                     ': "' // trim(cc_wrap%tracer_map%units(i)) // '"; expected kg kg-1 or ppm', ESMF_LOGMSG_ERROR, rc=rc)
                  rc = ESMF_FAILURE
                  return
               end if
               molecular_weight = 0.0_c_double
               catchem_status = catchem_state_get_species_mw_checked(cc_wrap%catchem_model%state_mgr_ptr, &
                  species_index, molecular_weight)
               if (catchem_status /= 0_c_int .or. .not. ieee_is_finite(molecular_weight) .or. molecular_weight <= 0.0_c_double) then
                  call ESMF_LogWrite('Invalid molecular weight for CATChem gas tracer: ' // &
                     trim(cc_wrap%tracer_map%names(i)), ESMF_LOGMSG_ERROR, rc=rc)
                  rc = ESMF_FAILURE
                  return
               end if
               if (is_ppm_unit(cc_wrap%tracer_map%units(i))) then
                  ! UFS/GOCART supplies these gas tracers in ppmv.
                  conversion_factor = 1.0_c_double
               else
                  ! Convert a host kg/kg mass mixing ratio to CATChem ppmv.
                  conversion_factor = real(AIRMW, c_double) / molecular_weight * 1.0e6_c_double
               end if
            else
               if (.not. is_mass_mixing_ratio_unit(cc_wrap%tracer_map%units(i)) .and. &
                  .not. is_micro_mass_mixing_ratio_unit(cc_wrap%tracer_map%units(i)) .and. &
                  .not. is_ppm_unit(cc_wrap%tracer_map%units(i))) then
                  call ESMF_LogWrite('Unsupported units for CATChem aerosol tracer ' // trim(cc_wrap%tracer_map%names(i)) // &
                     ': "' // trim(cc_wrap%tracer_map%units(i)) // '"; expected kg kg-1, ug kg-1, or ppm', ESMF_LOGMSG_ERROR, rc=rc)
                  rc = ESMF_FAILURE
                  return
               end if
               if (is_micro_mass_mixing_ratio_unit(cc_wrap%tracer_map%units(i))) then
                  ! UFS/GOCART supplies these aerosol tracers in ug/kg.
                  conversion_factor = 1.0_c_double
               else if (is_ppm_unit(cc_wrap%tracer_map%units(i))) then
                  ! The UFS/GOCART tracer convention defines ppm as a mass
                  ! scale (not a species-dependent molar mixing ratio):
                  ! 1 ppm = 1000 ug/kg.  This matches the legacy
                  ! AerosolTracerGetUnitsConv table.
                  conversion_factor = 1.0e3_c_double
               else
                  ! Convert a host kg/kg mass mixing ratio to CATChem ug/kg.
                  conversion_factor = 1.0e9_c_double
               end if
            end if
            if (.not. ieee_is_finite(conversion_factor) .or. conversion_factor <= 0.0_c_double) then
               call ESMF_LogWrite('Invalid conversion factor for CATChem tracer: ' // &
                  trim(cc_wrap%tracer_map%names(i)), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
            cc_wrap%tracer_map%host_to_catchem(i) = conversion_factor
            cc_wrap%tracer_map%catchem_to_host(i) = 1.0_c_double / conversion_factor
            do j = 1, i - 1
               if (cc_wrap%tracer_map%nuopc_to_cc(j) == cc_wrap%tracer_map%nuopc_to_cc(i)) then
                  call ESMF_LogWrite("Duplicate tracer mapping for species: " // &
                     trim(cc_wrap%tracer_map%names(i)), ESMF_LOGMSG_ERROR, rc=rc)
                  rc = ESMF_FAILURE
                  return
               end if
            end do
         else if (is_pm_diagnostic_name(cc_wrap%tracer_map%names(i))) then
            if (.not. is_micro_mass_concentration_unit(cc_wrap%tracer_map%units(i))) then
               call ESMF_LogWrite('Unsupported units for PM diagnostic ' // trim(cc_wrap%tracer_map%names(i)) // &
                  ': "' // trim(cc_wrap%tracer_map%units(i)) // '"; expected ug m-3', ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
            cc_wrap%tracer_map%entry_kind(i) = TRACER_DIAGNOSTIC
         else
            ! NUOPC tracer metadata can include host prognostics that are not
            ! part of the active chemistry mechanism (for example sphum).
            ! Keep those slots in the host array, but leave them untouched by
            ! CATChem.  Chemical membership remains entirely mechanism- and
            ! configuration-driven; no host tracer names are hardcoded here.
            cc_wrap%tracer_map%entry_kind(i) = TRACER_HOST_OWNED
            call ESMF_LogWrite("Ignoring host-owned tracer not present in active mechanism: " // &
               trim(cc_wrap%tracer_map%names(i)), ESMF_LOGMSG_INFO, rc=rc)
         end if
      end do

      !copy fields to cc_wrap
      cc_wrap%field_config = field_config
      ! Set the process-local grid variable
      cc_wrap%grid = input_grid

      ! Initialize AQMIO component if not done
      if (.not. ESMF_GridCompIsCreated(cc_wrap%iocomp)) then
         cc_wrap%iocomp = AQMIO_Create(cc_wrap%grid, rc =rc)
         if (rc /= CC_SUCCESS) return
      end if

      ! Initialize lat/lon stitched output if configured (multi-tile only)
      if (cc_wrap%catchem_model%is_latlon_output_enabled()) then
         call AQMIO_LatlonInit(cc_wrap%grid, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) then
            call ESMF_LogWrite('AQMIO_LatlonInit failed, lat/lon output disabled', &
               ESMF_LOGMSG_WARNING, rc=rc)
            rc = ESMF_SUCCESS  ! Non-fatal
         end if
      end if

      ! Set time information if provided
      if (present(stopTime)) then
         cc_wrap%endTime = stopTime
      end if
      if (present(startTime)) then
         cc_wrap%startTime = startTime
      end if

      if (present(timeStep)) then
         cc_wrap%timeStep = timeStep
         call ESMF_TimeIntervalGet(timeStep, s_i8=tstep_seconds, rc=rc)
         catchem_status = catchem_state_set_time_checked(cc_wrap%catchem_model%state_mgr_ptr, &
            0_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int, real(tstep_seconds, c_double))
         if (catchem_status /= 0_c_int) then
            call ESMF_LogWrite("Failed to set CATChem time state: "// &
               trim(cc_wrap%catchem_model%last_error), ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if
      end if

      ! Verify process registration.  The processes activated by the runtime
      ! YAML were already registered and created during model_initialize;
      ! add_process is a compatibility no-op that checks the Core handle.
      call cc_wrap%catchem_model%add_process(rc)
      num_processes = cc_wrap%catchem_model%get_num_processes()
      ! A process-free configuration is a valid lifecycle/contract run.  It
      ! permits an embedding application to initialize, exchange fields, and
      ! finalize CATChem before runtime YAML activates any processes.  Only a
      ! failed process-registration call is an initialization error.
      if (rc /= CC_SUCCESS) then
         call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="CATChem initialization failed", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return  ! bail out
      end if
      call catchem_log_run_phase(cc_wrap, 0, 'initialize: processes registered')

      ! Mark this process as initialized
      cc_wrap%initialized = .true.

      call ESMF_GridCompSetInternalState(model, is, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__)) return  ! bail out
      nullify(verify_is%wrap)
      call ESMF_GridCompGetInternalState(model, verify_is, verify_rc)
      if (verify_rc /= ESMF_SUCCESS .or. .not. associated(verify_is%wrap, cc_wrap)) then
         call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
            msg="CATChem internal state registration did not round-trip", &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end if

      !deallocate variable
      if (allocated(tracer_names)) deallocate(tracer_names)
      if (allocated(tracer_units)) deallocate(tracer_units)
      if (allocated(field_config%import_fields)) deallocate(field_config%import_fields)
      if (allocated(field_config%export_fields)) deallocate(field_config%export_fields)

      ! ! Initialize CF input system
      ! call cf_input_init('catchem_input_config.yml', grid, errflg)
      ! if (errflg /= ESMF_SUCCESS) then
      !   errmsg = 'Error initializing CF input system'
      !   errflg = CC_FAILURE
      !   return
      ! end if

      ! ! Initialize NetCDF output system
      ! call output_diagnostics_init('catchem_output_config.yml', grid, errflg)
      ! if (errflg /= ESMF_SUCCESS) then
      !   errmsg = 'Error initializing NetCDF output system'
      !   errflg = CC_FAILURE
      !   return
      ! end if

   end subroutine catchem_nuopc_init

   !> Get process-local CATChem wrapper (guaranteed thread/process safe)
   !!
   !! This function provides access to the process-local CATChem state
   !! using static local variables which are guaranteed to be process-private
   !! in MPI environments.
   !!
   !! \return Pointer to process-local CATChem wrapper
   ! function get_cc_wrap() result(wrap_ptr)
   !   type(ESMF_GridComp)         :: model
   !   type(cc_wrap_type), pointer :: wrap_ptr

   !   type(CATChem_InternalState), target :: is
   !   integer                    :: verbosity, localrc
   !   character(len=128) :: name

   !   ! -- get component's information
   !   call NUOPC_CompGet(model, rc=localrc)
   !   if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
   !     line=__LINE__,  file=__FILE__)) return  ! bail out

   !   ! -- get component's internal state
   !   call ESMF_GridCompGetInternalState(model, is, localrc)
   !   if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
   !     line=__LINE__,  file=__FILE__))  return  ! bail out

   !   wrap_ptr => is%wrap
   ! end function get_cc_wrap

   ! Run CATChem processes for NUOPC interface
   !!
   !! \param config         CATChem configuration
   !! \param catchem_states CATChem container
   !! \param dustState      Dust process state
   !! \param seaSaltState   Sea salt process state
   !! \param dryDepState    Dry deposition process state
   !! \param    dt             Time step (seconds)
   !! \param    current_time   Current model time
   !! \param   errflg         Error flag
   !! \param   errmsg         Error message
   !!
   subroutine catchem_nuopc_run( cc_wrap, dt, current_time, errmsg, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      real(ESMF_KIND_R8), intent(in) :: dt
      type(ESMF_Time), intent(in) :: current_time
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      errmsg = ''

      cc_wrap%timestep_counter = cc_wrap%timestep_counter + 1
      call catchem_log_run_phase(cc_wrap, cc_wrap%timestep_counter, 'enter emissions update')

      ! Update extemission data first
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("catchem_emis_update", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call catchem_emis_update(cc_wrap%ext_emis, cc_wrap%catchem_model%cpp_core_ptr, current_time, cc_wrap%catchem_model%nz, cc_wrap%iocomp, cc_wrap%grid, real(dt, fp), rc)

      if (rc == CC_SUCCESS) then
         call catchem_log_run_phase(cc_wrap, cc_wrap%timestep_counter, 'enter static-met bind')
         call bind_static_met_from_aqmio(cc_wrap, rc)
      end if
      if (rc /= CC_SUCCESS) then
         errmsg = 'Error binding AQMIO static met fields'
         return
      end if
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("catchem_emis_update", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      !Run CATChem processes
      call catchem_log_run_phase(cc_wrap, cc_wrap%timestep_counter, 'enter core process dispatch')
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("cc_wrap%catchem_model%run_timestep", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call cc_wrap%catchem_model%run_timestep(cc_wrap%timestep_counter, real(dt, fp), rc)
      if (rc /= CC_SUCCESS) then
         write(errmsg, '(A,I0)') 'Error in run_timestep at timestep = ', cc_wrap%timestep_counter
         call catchem_append_timestep_failure(cc_wrap%catchem_model%cpp_core_ptr, errmsg)
         return
      end if
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("cc_wrap%catchem_model%run_timestep", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      ! Update PM2.5/PM10 aerosol diagnostics. These are stored in the
      ! DiagnosticManager so they are available both for NetCDF output and for
      ! NUOPC export, and must be computed after run_timestep (so concentrations
      ! are current) and before the export transform.
      call catchem_log_run_phase(cc_wrap, cc_wrap%timestep_counter, 'enter PM diagnostics')
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("update_pm_diagnostics", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call update_pm_diagnostics(cc_wrap, rc)
      if (rc /= CC_SUCCESS) then
         errmsg = 'Error updating PM2.5/PM10 aerosol diagnostics'
         return
      end if
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("update_pm_diagnostics", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      ! Write NetCDF output diagnostics if needed
      call catchem_log_run_phase(cc_wrap, cc_wrap%timestep_counter, 'enter diagnostics output')
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionEnter("catchem_diagnostics_write", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif
      call catchem_diagnostics_write(cc_wrap, current_time, rc)
      if (rc /= ESMF_SUCCESS) then
         errmsg = 'Error writing NetCDF output diagnostics'
         return
      end if
#ifdef CATCHEM_TRACE_NUOPC
      call ESMF_TraceRegionExit("catchem_diagnostics_write", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
#endif

      call catchem_log_run_phase(cc_wrap, cc_wrap%timestep_counter, 'complete')

   end subroutine catchem_nuopc_run

   !> \brief Bind static dust-support met fields from AQMIO-loaded data.
   !!
   !! Some dust inputs (e.g., clay/sand fractions and threshold friction
   !! velocity) come from static files handled via the external-emissions AQMIO
   !! path. This routine maps those fields into CATChem met state when present.
   subroutine bind_static_met_from_aqmio(cc_wrap, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      if (cc_wrap%ext_emis%n_categories == 0) then
#ifdef CATCHEM_TRACE_NUOPC
         write(*,'(A)') '[CATCHEM DEBUG] bind_static_met_from_aqmio: no emission categories available'
         call flush(6)
#endif
         return
      end if

      ! Bind each FENGSHA dust input by the EXACT canonical field name that the
      ! Bind each FENGSHA dust input by its canonical emission-map TARGET
      ! (MET_CLAYFRAC, MET_SSM, MET_RDRAG, ...).  bind_static_field resolves the
      ! target back to the emission field that declares it in `map:` using the
      ! emission config itself, so the config is the single source of truth for
      ! the field-name-to-target mapping (dust.sep -> MET_SSM,
      ! fengsha.PC/albedo_drag -> MET_RDRAG, ...).  No hardcoded field-name
      ! aliases and no substring/fuzzy matching.
      call bind_static_field(cc_wrap, 'MET_CLAYFRAC', &
         'CLAYFRAC', cc_wrap%dust_clayfrac, 1.0_c_double, rc)
      if (rc /= CC_SUCCESS) return

      call bind_static_field(cc_wrap, 'MET_SANDFRAC', &
         'SNDFRC', cc_wrap%dust_sandfrac, 1.0_c_double, rc)
      if (rc /= CC_SUCCESS) return

      call bind_static_field(cc_wrap, 'MET_SSM', &
         'SSM', cc_wrap%dust_ssm, 1.0_c_double, rc)
      if (rc /= CC_SUCCESS) return

      call bind_static_field(cc_wrap, 'MET_RDRAG', &
         'RDRAG', cc_wrap%dust_rdrag, 1.0_c_double, rc)
      if (rc /= CC_SUCCESS) return

      call bind_static_field(cc_wrap, 'MET_USTAR_THRESHOLD', &
         'USTAR_THRESHOLD', cc_wrap%dust_ustar_threshold, 1.0_c_double, rc)

   end subroutine bind_static_met_from_aqmio

   !> \brief Resolve the emission field NAME that maps to a canonical target.
   !!
   !! Consults the emission configuration (the single source of truth): scans
   !! every category/field and returns the first field whose `map:` list
   !! contains `target_name` (e.g. 'MET_SSM').  Returns '' if none maps to it.
   !! This removes any hardcoded field-name alias from the binding code -- the
   !! config alone decides which field (dust.sep, fengsha.PC/albedo_drag, ...)
   !! supplies each MET_* target.
   function resolve_emission_field_for_target(cc_wrap, target_name) result(field_out)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: target_name
      character(len=64) :: field_out

      type(c_ptr) :: core_ptr
      character(len=64) :: category_name, field_name, mapped_species
      integer(c_int) :: n_categories, n_fields, n_maps
      integer :: icat, ifield, imap
      real(c_double) :: scale_factor
      integer(c_int) :: species_index

      field_out = ''
      core_ptr = cc_wrap%catchem_model%cpp_core_ptr
      if (.not. c_associated(core_ptr)) return
      if (catchem_config_has_emission_mapping(core_ptr) == 0_c_int) return

      n_categories = catchem_config_get_emission_category_count(core_ptr)
      do icat = 0, n_categories - 1
         call catchem_config_get_emission_category_name_at(core_ptr, icat, category_name, 64_c_int)
         call trim_at_null(category_name)
         n_fields = catchem_config_get_emission_field_count(core_ptr, trim(category_name) // c_null_char)
         do ifield = 0, n_fields - 1
            call catchem_config_get_emission_field_name_at(core_ptr, trim(category_name) // c_null_char, &
               ifield, field_name, 64_c_int)
            call trim_at_null(field_name)
            n_maps = catchem_config_get_emission_species_map_count(core_ptr, &
               trim(category_name) // c_null_char, trim(field_name) // c_null_char)
            do imap = 0, n_maps - 1
               call catchem_config_get_emission_species_map_at(core_ptr, &
                  trim(category_name) // c_null_char, trim(field_name) // c_null_char, imap, &
                  mapped_species, 64_c_int, scale_factor, species_index)
               call trim_at_null(mapped_species)
               if (trim(mapped_species) == trim(target_name)) then
                  field_out = trim(field_name)
                  return
               end if
            end do
         end do
      end do

   contains
      !> Truncate a C-filled buffer at its NUL terminator and blank-pad.
      subroutine trim_at_null(str)
         character(len=*), intent(inout) :: str
         integer :: idx
         idx = index(str, c_null_char)
         if (idx > 0) str(idx:) = ' '
      end subroutine trim_at_null

   end function resolve_emission_field_for_target

   !> \brief Helper: bind one static AQMIO field into CATChem met state.
   subroutine bind_static_field(cc_wrap, target_name, met_name, met_buffer, scale, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: target_name  !< canonical emission-map target, e.g. 'MET_SSM'
      character(len=*), intent(in) :: met_name
      real(c_double), allocatable, intent(inout) :: met_buffer(:,:)
      real(c_double), intent(in) :: scale
      integer, intent(out) :: rc

      type(ExtEmisFieldType), pointer :: src_field
      character(len=64) :: resolved_field

      rc = CC_SUCCESS

      ! Resolve the emission field NAME that declares target_name in its `map:`
      ! by consulting the emission CONFIG (the single source of truth), then
      ! bind that field by its exact name.  This keeps the field-name <-> target
      ! mapping entirely in configuration (dust.sep -> MET_SSM,
      ! fengsha.PC/albedo_drag -> MET_RDRAG, ...), with no hardcoded field-name
      ! aliases in code and no substring/fuzzy matching.
      resolved_field = resolve_emission_field_for_target(cc_wrap, target_name)

      src_field => null()
      if (len_trim(resolved_field) > 0) &
         src_field => cc_wrap%ext_emis%find_emission_field(trim(resolved_field))
      if (associated(src_field)) then
         call bind_static_field_data(cc_wrap, src_field, met_name, met_buffer, scale, rc)
         return
      end if

      ! Not found: leave the met buffer unbound (the dust process validates the
      ! required inputs and will error clearly if a mandatory field is missing).
#ifdef CATCHEM_TRACE_NUOPC
      write(*,'(A,A,A,I0)') '[CATCHEM DEBUG] bind_static_field missing met=', trim(met_name), &
         ' n_categories=', cc_wrap%ext_emis%n_categories
      call flush(6)
      do i = 1, cc_wrap%ext_emis%n_categories
         write(*,'(A,A,A,I0)') '[CATCHEM DEBUG]   category=', trim(cc_wrap%ext_emis%categories(i)%category_name), &
            ' n_fields=', cc_wrap%ext_emis%categories(i)%n_fields
         call flush(6)
         if (.not. allocated(cc_wrap%ext_emis%categories(i)%fields)) cycle
         do j = 1, cc_wrap%ext_emis%categories(i)%n_fields
            write(*,'(A,A,A,L1,A,L1)') '[CATCHEM DEBUG]     field=', &
               trim(cc_wrap%ext_emis%categories(i)%fields(j)%field_name), &
               ' emission_data=', allocated(cc_wrap%ext_emis%categories(i)%fields(j)%emission_data), &
               ' interp_t1=', allocated(cc_wrap%ext_emis%categories(i)%fields(j)%interp_data_t1)
            call flush(6)
         end do
      end do
#endif

   end subroutine bind_static_field

   !> \brief Copy one loaded static AQMIO field into a persistent CATChem met buffer.
   subroutine bind_static_field_data(cc_wrap, src_field, met_name, met_buffer, scale, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ExtEmisFieldType), intent(in) :: src_field
      character(len=*), intent(in) :: met_name
      real(c_double), allocatable, intent(inout) :: met_buffer(:,:)
      real(c_double), intent(in) :: scale
      integer, intent(out) :: rc

      integer :: nx, ny

      rc = CC_SUCCESS

      if (allocated(src_field%interp_data_t1)) then
         nx = size(src_field%interp_data_t1, 1)
         ny = size(src_field%interp_data_t1, 2)
         if (size(src_field%interp_data_t1, 3) < 1 .or. size(src_field%interp_data_t1, 4) < 1) return
         if (.not. allocated(met_buffer) .or. size(met_buffer, 1) /= nx .or. size(met_buffer, 2) /= ny) then
            if (allocated(met_buffer)) deallocate(met_buffer)
            allocate(met_buffer(nx, ny))
         end if
         met_buffer = real(src_field%interp_data_t1(:,:,1,1), c_double) * scale
      else if (allocated(src_field%emission_data)) then
         nx = size(src_field%emission_data, 1)
         ny = size(src_field%emission_data, 2)
         if (size(src_field%emission_data, 3) < 1 .or. size(src_field%emission_data, 4) < 1) return
         if (.not. allocated(met_buffer) .or. size(met_buffer, 1) /= nx .or. size(met_buffer, 2) /= ny) then
            if (allocated(met_buffer)) deallocate(met_buffer)
            allocate(met_buffer(nx, ny))
         end if
         met_buffer = real(src_field%emission_data(:,:,1,1), c_double) * scale
      else
         return
      end if

      call cc_wrap%catchem_model%bind_met_2d(trim(met_name), met_buffer)

   end subroutine bind_static_field_data

   ! Finalize CATChem for NUOPC interface
   !!
   !! \param catchem_states CATChem container
   !! \param dustState      Dust process state
   !! \param seaSaltState   Sea salt process state
   !! \param dryDepState    Dry deposition process state
   !! \param   errflg         Error flag
   !! \param   errmsg         Error message
   !!
   subroutine catchem_nuopc_finalize(cc_wrap, rc, errmsg)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(out) :: rc
      character(len=*), intent(out) :: errmsg

      ! Get process-local state
      !type(cc_wrap_type), pointer :: cc_wrap
      !cc_wrap => get_cc_wrap()

      rc = CC_SUCCESS
      errmsg = ''

      ! Finalize CF input system
      !call cf_input_finalize()

      ! Finalize NetCDF output system
      !call output_diagnostics_finalize()

      ! Finalize CATChem model
      if (cc_wrap%initialized) then
         !finalize extemission data
         call catchem_emis_finalize(cc_wrap%ext_emis, rc)
         !finalize catchem model
         call cc_wrap%catchem_model%finalize(rc)
         if (rc /= CC_SUCCESS) then
            errmsg = 'Error in calling cc_wrap%catchem_model%finalize!'
            return
         end if
         cc_wrap%initialized = .false.
      end if

      ! Deallocate field mappings
      if (allocated(cc_wrap%field_config%import_fields)) deallocate(cc_wrap%field_config%import_fields)
      if (allocated(cc_wrap%field_config%export_fields)) deallocate(cc_wrap%field_config%export_fields)

      ! Deallocate tracer mapping
      if (allocated(cc_wrap%tracer_map%nuopc_to_cc)) deallocate(cc_wrap%tracer_map%nuopc_to_cc)
      if (allocated(cc_wrap%tracer_map%entry_kind)) deallocate(cc_wrap%tracer_map%entry_kind)
      if (allocated(cc_wrap%tracer_map%names)) deallocate(cc_wrap%tracer_map%names)
      if (allocated(cc_wrap%tracer_map%units)) deallocate(cc_wrap%tracer_map%units)
      if (allocated(cc_wrap%tracer_map%host_to_catchem)) deallocate(cc_wrap%tracer_map%host_to_catchem)
      if (allocated(cc_wrap%tracer_map%catchem_to_host)) deallocate(cc_wrap%tracer_map%catchem_to_host)

      if (allocated(cc_wrap%lat)) deallocate(cc_wrap%lat)
      if (allocated(cc_wrap%lon)) deallocate(cc_wrap%lon)
      if (allocated(cc_wrap%area_m2)) deallocate(cc_wrap%area_m2)
      if (allocated(cc_wrap%z0_m)) deallocate(cc_wrap%z0_m)
      if (allocated(cc_wrap%dust_clayfrac)) deallocate(cc_wrap%dust_clayfrac)
      if (allocated(cc_wrap%dust_sandfrac)) deallocate(cc_wrap%dust_sandfrac)
      if (allocated(cc_wrap%dust_ssm)) deallocate(cc_wrap%dust_ssm)
      if (allocated(cc_wrap%dust_rdrag)) deallocate(cc_wrap%dust_rdrag)
      if (allocated(cc_wrap%dust_ustar_threshold)) deallocate(cc_wrap%dust_ustar_threshold)
      if (allocated(cc_wrap%met_buf_2d)) deallocate(cc_wrap%met_buf_2d)
      if (allocated(cc_wrap%met_buf_3d)) deallocate(cc_wrap%met_buf_3d)
      if (allocated(cc_wrap%chem_buf_4d)) deallocate(cc_wrap%chem_buf_4d)
      if (allocated(cc_wrap%host_tracer_buf_4d)) deallocate(cc_wrap%host_tracer_buf_4d)

      ! Clean up lat/lon stitched output resources
      call AQMIO_LatlonCleanup(rc=rc)

      ! Destroy the IO component and its per-tile taskComps
      ! Must happen before the parent component's VM is torn down
      call AQMIO_Destroy(cc_wrap%iocomp, rc=rc)

   end subroutine catchem_nuopc_finalize

   ! Transform NUOPC import fields to CATChem states
   !!
   !! \param    importState    NUOPC import state
   !! \param catchem_states CATChem container
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_nuopc_to_catchem(cc_wrap, importState, currTime, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_State), intent(in) :: importState
      type(ESMF_Time), intent(in) :: currTime
      integer, intent(out) :: rc

      type(ESMF_Field) :: field
      logical, allocatable :: set_required_met(:)
      integer(ESMF_KIND_I8) :: timestep_seconds
      integer :: year, month, day, hour, minute, second
      integer :: i, n, n_met
      logical :: required_from_import
      integer(c_int) :: catchem_status

      rc = ESMF_SUCCESS

      if (catchem_state_begin_import_generation(cc_wrap%catchem_model%state_mgr_ptr) /= 0_c_int) then
         call ESMF_LogWrite("Failed to begin CATChem import generation", ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE
         return
      end if

      ! assign time to catchem model's time state
      call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
         h=hour, m=minute, s=second, rc=rc)
      call ESMF_TimeIntervalGet(cc_wrap%timeStep, s_i8=timestep_seconds, rc=rc)
      catchem_status = catchem_state_set_time_checked(cc_wrap%catchem_model%state_mgr_ptr, &
         int(year, c_int), int(month, c_int), int(day, c_int), &
         int(hour, c_int), int(minute, c_int), int(second, c_int), &
         0_c_int, real(timestep_seconds, c_double))
      if (catchem_status /= 0_c_int) then
         call ESMF_LogWrite("Failed to set CATChem time state", ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE
         return
      end if

      ! This is to check if all required met fields in CATChem are set
      if (allocated(cc_wrap%catchem_model%required_fields)) then
         n_met = size(cc_wrap%catchem_model%required_fields)
         if (n_met > 0) then
            allocate(set_required_met(n_met))
            set_required_met = .false.
         end if
      end if

      ! Loop through all import fields and transform to CATChem states
      if (cc_wrap%field_config%n_import_fields > 0) then
         if (.not. allocated(cc_wrap%met_buf_2d)) then
            allocate(cc_wrap%met_buf_2d(cc_wrap%field_config%n_import_fields))
         end if
         if (.not. allocated(cc_wrap%met_buf_3d)) then
            allocate(cc_wrap%met_buf_3d(cc_wrap%field_config%n_import_fields))
         end if
      end if

      do n = 1, cc_wrap%field_config%n_import_fields

         ! Keep an unadvertised optional input out of ESMF_StateGet.  Calling
         ! StateGet for a field the YAML explicitly chose not to advertise
         ! produces an ESMF error log even though absence is valid.
         if (cc_wrap%field_config%import_fields(n)%optional .and. &
            .not. cc_wrap%field_config%import_fields(n)%advertise) cycle

         ! Try to get field from import state (will fail if not present)
         call ESMF_StateGet(importState, trim(cc_wrap%field_config%import_fields(n)%standard_name), field, rc=rc)

         if (rc /= ESMF_SUCCESS) then
            if (.not. cc_wrap%field_config%import_fields(n)%optional) then
               call ESMF_LogWrite("Required field not found: "// &
                  trim(cc_wrap%field_config%import_fields(n)%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            else
               cycle  ! Skip optional fields that are not present
            end if
         end if

         ! Transform based on field type and dimensions
         call transform_field_to_catchem(cc_wrap, field, cc_wrap%field_config%import_fields(n), &
            cc_wrap%field_config%import_fields(n)%optional, set_required_met, rc, n)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

      end do

      ! derive required met fields in C++ StateManager
      catchem_status = catchem_state_derive_airden_dry_checked(cc_wrap%catchem_model%state_mgr_ptr)
      if (catchem_status == 0_c_int) &
         catchem_status = catchem_state_derive_bxheight_checked(cc_wrap%catchem_model%state_mgr_ptr)
      if (catchem_status /= 0_c_int) then
         call ESMF_LogWrite("CATChem physical derivation failed", ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE
         return
      end if

      !check if all require met fields are set
      if (allocated(cc_wrap%catchem_model%required_fields) .and. allocated(set_required_met)) then
         do i = 1, n_met
            ! Only validate requirements owned by the NUOPC import contract
            ! here.  Derived fields and static AQMIO inputs are validated by
            ! the core execution plan after their respective preparation
            ! stages; they must never be replaced with fallback values.
            required_from_import = .false.
            do n = 1, cc_wrap%field_config%n_import_fields
               if (.not. cc_wrap%field_config%import_fields(n)%optional .and. &
                  (trim(cc_wrap%field_config%import_fields(n)%catchem_var) == &
                  trim(cc_wrap%catchem_model%required_fields(i)) .or. &
                  trim(cc_wrap%field_config%import_fields(n)%host_tracer_var) == &
                  trim(cc_wrap%catchem_model%required_fields(i)))) then
                  required_from_import = .true.
                  exit
               end if
            end do
            if (.not. required_from_import) cycle
            if (.not. set_required_met(i)) then
               !write(*,*) 'Wait. A required field is not set: ' // trim(cc_wrap%catchem_model%required_fields(i))
               call ESMF_LogWrite("Required met field not set yet: "// &
                  trim(cc_wrap%catchem_model%required_fields(i)), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
         end do
         !deallocate array
         deallocate(set_required_met)
      end if

   end subroutine transform_nuopc_to_catchem

   ! Transform CATChem states to NUOPC export fields
   !!
   !! \param exportState    NUOPC export state
   !! \param    catchem_states CATChem container
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_catchem_to_nuopc(cc_wrap, exportState, rc, currTime)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_State), intent(inout) :: exportState
      integer, intent(out) :: rc
      type(ESMF_Time), intent(in), optional :: currTime

      type(ESMF_Field) :: field
      integer :: n
      !type(cc_wrap_type), pointer :: cc_wrap

      rc = ESMF_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Loop through all export fields and transform from CATChem states
      do n = 1, cc_wrap%field_config%n_export_fields

         ! Try to get field from export state (will fail if not present)
         call ESMF_StateGet(exportState, trim(cc_wrap%field_config%export_fields(n)%standard_name), field, rc=rc)

         if (rc /= ESMF_SUCCESS) then
            if (.not. cc_wrap%field_config%export_fields(n)%optional) then
               call ESMF_LogWrite("Required export field not found: "// &
                  trim(cc_wrap%field_config%export_fields(n)%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            else
               cycle  ! Skip optional fields that are not present
            end if
         end if

         ! Transform from CATChem to field
         call transform_catchem_to_field(cc_wrap, field, cc_wrap%field_config%export_fields(n), rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         if (present(currTime)) then
            ! NUOPC validates each shared Field's timestamp, not just the
            ! enclosing State.  Stamp the completed member at the current
            ! coupling time so the consumer observes a valid export.
            call NUOPC_SetTimestamp(field, currTime, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
         end if

      end do

   end subroutine transform_catchem_to_nuopc

   ! Transform individual field to CATChem state
   !!
   !! \param    field          ESMF field
   !! \param    field_map      Field mapping information
   !! \param catchem_states CATChem container
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_field_to_catchem(cc_wrap, field, field_map, required, is_met_set, rc, fidx)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Field), intent(in) :: field
      type(field_mapping_type), intent(in) :: field_map
      logical, allocatable, intent(inout) :: is_met_set(:)
      logical, intent(in) :: required
      integer, intent(out) :: rc
      integer, intent(in) :: fidx

      !local vars
      real(ESMF_KIND_R8), pointer :: fptr4d(:,:,:,:), fptr3d(:,:,:), fptr2d(:,:)
      integer :: met_index, v_cc, v, found_index, expected_levels, expected_tracers, localrc
      integer(c_int) :: catchem_status, species_count
      logical :: tracer_shape_valid
      type(ESMF_Info) :: field_info
      character(len=64) :: observed_units

      ! `required` is carried in the signature for the import loop but the
      ! missing-field policy is enforced by the caller; reference it so the
      ! interface stays as documented without an unused-argument warning.
      associate(unused_required => required); end associate

      rc = ESMF_SUCCESS

      observed_units = ''
      localrc = ESMF_SUCCESS
      call ESMF_InfoGetFromHost(field, field_info, rc=localrc)
      if (localrc == ESMF_SUCCESS) then
         call ESMF_InfoGet(field_info, key="units", value=observed_units, default="", rc=localrc)
         if (localrc == ESMF_SUCCESS .and. len_trim(observed_units) > 0 .and. len_trim(field_map%units) > 0) then
            if (trim(observed_units) /= trim(field_map%units)) then
               call ESMF_LogWrite("Unit mismatch for import field " // trim(field_map%standard_name) // &
                  ": expected " // trim(field_map%units) // ", observed " // trim(observed_units), &
                  ESMF_LOGMSG_ERROR, rc=rc)
               rc = ESMF_FAILURE
               return
            end if
         end if
      end if

      ! Transform based on field mapping
      select case (field_map%dimensions)

         ! 2D meteorological fields
       case (2)
         nullify(fptr2d)
         call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (.not. associated(fptr2d)) then
#ifdef CATCHEM_TRACE_NUOPC
            write(*, '(A,A,A,A)') '[CATCHEM DEBUG] transform_field_to_catchem 2D UNASSOCIATED: ', &
               trim(field_map%standard_name), ' -> ', trim(field_map%catchem_var)
            call flush(6)
#endif
            call ESMF_LogWrite("transform_field_to_catchem: fptr2d NOT associated for "// &
               trim(field_map%standard_name)//" -> "//trim(field_map%catchem_var), &
               ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         if (size(fptr2d, 1) /= cc_wrap%catchem_model%nx .or. &
            size(fptr2d, 2) /= cc_wrap%catchem_model%ny) then
            call ESMF_LogWrite("Shape mismatch for 2D import field: " // trim(field_map%standard_name) // &
               " -> " // trim(field_map%catchem_var), ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         if (.not. allocated(cc_wrap%met_buf_2d(fidx)%data)) then
            allocate(cc_wrap%met_buf_2d(fidx)%data(size(fptr2d, 1), size(fptr2d, 2)))
         end if
         cc_wrap%met_buf_2d(fidx)%data = real(fptr2d, c_double)

         ! Standard pointer mapping to C++ core via persistent contiguous buffer
         if (trim(field_map%catchem_var) == 'Z0') then
            if (.not. allocated(cc_wrap%z0_m)) then
               allocate(cc_wrap%z0_m(size(fptr2d, 1), size(fptr2d, 2)))
            end if
            cc_wrap%z0_m = cc_wrap%met_buf_2d(fidx)%data * 0.01_c_double
            call cc_wrap%catchem_model%bind_met_2d("Z0", cc_wrap%z0_m, rc)
         else
            call cc_wrap%catchem_model%bind_met_2d(trim(field_map%catchem_var), cc_wrap%met_buf_2d(fidx)%data, rc)
         end if
         if (rc /= CC_SUCCESS) return

         if (allocated(cc_wrap%catchem_model%required_fields) .and. allocated(is_met_set)) then
            met_index = cc_wrap%catchem_model%get_required_met_index( trim(field_map%catchem_var) )
            if (met_index > 0 .and. met_index <= size(is_met_set)) then
               is_met_set(met_index) = .true.
            end if
         end if

         ! 3D meteorological fields
       case (3)
         nullify(fptr3d)
         call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (.not. associated(fptr3d)) then
#ifdef CATCHEM_TRACE_NUOPC
            write(*, '(A,A,A,A)') '[CATCHEM DEBUG] transform_field_to_catchem 3D UNASSOCIATED: ', &
               trim(field_map%standard_name), ' -> ', trim(field_map%catchem_var)
            call flush(6)
#endif
            call ESMF_LogWrite("transform_field_to_catchem: fptr3d NOT associated for "// &
               trim(field_map%standard_name)//" -> "//trim(field_map%catchem_var), &
               ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         select case (trim(field_map%vertical_axis))
          case ('level')
            expected_levels = cc_wrap%catchem_model%nz
          case ('interface')
            expected_levels = cc_wrap%catchem_model%nz + 1
          case ('level_to_interface')
            ! UFS provides the physical interfaces 1:nz.  CATChem retains
            ! the legacy 0:nz storage convention, so the cap supplies the
            ! zero-valued index-0 surface slot while importing nlev values.
            expected_levels = cc_wrap%catchem_model%nz
          case ('soil_layer')
            expected_levels = size(fptr3d, 3)
          case default
            call ESMF_LogWrite("Unsupported vertical_axis for 3D import field: " // &
               trim(field_map%standard_name) // " (" // trim(field_map%vertical_axis) // ")", &
               ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end select
         if (size(fptr3d,1) /= cc_wrap%catchem_model%nx .or. &
            size(fptr3d,2) /= cc_wrap%catchem_model%ny .or. size(fptr3d,3) /= expected_levels) then
            call ESMF_LogWrite("Shape mismatch for 3D import field: " // trim(field_map%standard_name), &
               ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         if (trim(field_map%vertical_axis) == 'level_to_interface') then
            if (allocated(cc_wrap%met_buf_3d(fidx)%data)) then
               if (size(cc_wrap%met_buf_3d(fidx)%data, 1) /= size(fptr3d, 1) .or. &
                  size(cc_wrap%met_buf_3d(fidx)%data, 2) /= size(fptr3d, 2) .or. &
                  size(cc_wrap%met_buf_3d(fidx)%data, 3) /= size(fptr3d, 3) + 1) then
                  deallocate(cc_wrap%met_buf_3d(fidx)%data)
               end if
            end if
            if (.not. allocated(cc_wrap%met_buf_3d(fidx)%data)) then
               allocate(cc_wrap%met_buf_3d(fidx)%data(size(fptr3d, 1), size(fptr3d, 2), size(fptr3d, 3) + 1))
            end if
            ! Match the upstream cap's interface reconstruction: the UFS
            ! nlev values occupy CATChem entries 1:nlev and the final value
            ! is repeated at nlev+1.
            cc_wrap%met_buf_3d(fidx)%data(:,:,1:size(fptr3d,3)) = real(fptr3d, c_double)
            cc_wrap%met_buf_3d(fidx)%data(:,:,size(fptr3d,3)+1) = real(fptr3d(:,:,size(fptr3d,3)), c_double)
         else
            if (.not. allocated(cc_wrap%met_buf_3d(fidx)%data)) then
               allocate(cc_wrap%met_buf_3d(fidx)%data(size(fptr3d, 1), size(fptr3d, 2), size(fptr3d, 3)))
            end if
            cc_wrap%met_buf_3d(fidx)%data = real(fptr3d, c_double)
         end if

         ! UFS provides geopotential for Z and ZMID; CATChem process
         ! kernels consume geometric height in metres, as in upstream.
         if (trim(field_map%catchem_var) == 'Z' .or. trim(field_map%catchem_var) == 'ZMID') then
            cc_wrap%met_buf_3d(fidx)%data = cc_wrap%met_buf_3d(fidx)%data / real(g0, c_double)
         end if

         ! Bind through the checked semantic contract.  The field mapping, not
         ! a variable-name special case, defines whether vertical extent is
         ! atmospheric levels, interfaces, or host-defined soil layers.
         select case (trim(field_map%vertical_axis))
          case ('level')
            call cc_wrap%catchem_model%bind_met_3d_axis(trim(field_map%catchem_var), &
               cc_wrap%met_buf_3d(fidx)%data, 0, rc)
          case ('interface')
            call cc_wrap%catchem_model%bind_met_3d_axis(trim(field_map%catchem_var), &
               cc_wrap%met_buf_3d(fidx)%data, 1, rc)
          case ('level_to_interface')
            call cc_wrap%catchem_model%bind_met_3d_axis(trim(field_map%catchem_var), &
               cc_wrap%met_buf_3d(fidx)%data, 1, rc)
          case ('soil_layer')
            call cc_wrap%catchem_model%bind_met_3d_axis(trim(field_map%catchem_var), &
               cc_wrap%met_buf_3d(fidx)%data, 2, rc)
         end select
         if (rc /= CC_SUCCESS) return

         if (allocated(cc_wrap%catchem_model%required_fields) .and. allocated(is_met_set)) then
            met_index = cc_wrap%catchem_model%get_required_met_index( trim(field_map%catchem_var) )
            if (met_index > 0 .and. met_index <= size(is_met_set)) then
               is_met_set(met_index) = .true.
            end if
         end if

         ! 4D tracer concentrations
       case (4)
         nullify(fptr4d)
         call ESMF_FieldGet(field, farrayPtr=fptr4d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (.not. associated(fptr4d)) then
#ifdef CATCHEM_TRACE_NUOPC
            write(*, '(A,A,A,A)') '[CATCHEM DEBUG] transform_field_to_catchem 4D UNASSOCIATED: ', &
               trim(field_map%standard_name), ' -> ', trim(field_map%catchem_var)
            call flush(6)
#endif
            call ESMF_LogWrite("transform_field_to_catchem: fptr4d NOT associated for "// &
               trim(field_map%standard_name)//" -> "//trim(field_map%catchem_var), &
               ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         ! Preserve every host tracer, including host-owned tracers that are
         ! not in the active CATChem mechanism (for example moisture).  The
         ! export state is not guaranteed to alias the import state, so an
         ! unowned slot must never be left undefined for the host component.
         ! Fortran does not short-circuit .or. expressions, so do not query
         ! SIZE on the unallocated buffer in the same logical expression.
         if (.not. allocated(cc_wrap%host_tracer_buf_4d)) then
            allocate(cc_wrap%host_tracer_buf_4d(size(fptr4d,1), size(fptr4d,2), &
               size(fptr4d,3), size(fptr4d,4)))
         else if (size(cc_wrap%host_tracer_buf_4d,1) /= size(fptr4d,1) .or. &
            size(cc_wrap%host_tracer_buf_4d,2) /= size(fptr4d,2) .or. &
            size(cc_wrap%host_tracer_buf_4d,3) /= size(fptr4d,3) .or. &
            size(cc_wrap%host_tracer_buf_4d,4) /= size(fptr4d,4)) then
            deallocate(cc_wrap%host_tracer_buf_4d)
            allocate(cc_wrap%host_tracer_buf_4d(size(fptr4d,1), size(fptr4d,2), &
               size(fptr4d,3), size(fptr4d,4)))
         end if
         cc_wrap%host_tracer_buf_4d = real(fptr4d, c_double)

         catchem_status = catchem_state_get_species_count_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, species_count)
         if (catchem_status /= 0_c_int) then
            rc = ESMF_FAILURE
            return
         end if
         v_cc = int(species_count)
         if (v_cc <= 0) v_cc = size(fptr4d, 4)
         expected_tracers = v_cc
         tracer_shape_valid = size(fptr4d,4) <= v_cc
         if (allocated(cc_wrap%tracer_map%nuopc_to_cc)) then
            expected_tracers = size(cc_wrap%tracer_map%nuopc_to_cc)
            tracer_shape_valid = size(fptr4d,4) == expected_tracers
         end if

         if (size(fptr4d,1) /= cc_wrap%catchem_model%nx .or. &
            size(fptr4d,2) /= cc_wrap%catchem_model%ny .or. &
            size(fptr4d,3) /= cc_wrap%catchem_model%nz .or. &
            .not. tracer_shape_valid) then
            call ESMF_LogWrite("Shape mismatch for 4D chemistry import field: " // &
               trim(field_map%standard_name), ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         if (.not. allocated(cc_wrap%chem_buf_4d)) then
            allocate(cc_wrap%chem_buf_4d(size(fptr4d, 1), size(fptr4d, 2), size(fptr4d, 3), v_cc))
            cc_wrap%chem_buf_4d = 0.0_c_double
         else if (size(cc_wrap%chem_buf_4d, 1) /= size(fptr4d, 1) .or. &
            size(cc_wrap%chem_buf_4d, 2) /= size(fptr4d, 2) .or. &
            size(cc_wrap%chem_buf_4d, 3) /= size(fptr4d, 3) .or. &
            size(cc_wrap%chem_buf_4d, 4) /= v_cc) then
            deallocate(cc_wrap%chem_buf_4d)
            allocate(cc_wrap%chem_buf_4d(size(fptr4d, 1), size(fptr4d, 2), size(fptr4d, 3), v_cc))
            cc_wrap%chem_buf_4d = 0.0_c_double
         end if

         if (allocated(cc_wrap%tracer_map%nuopc_to_cc)) then
            do v = 1, min(size(fptr4d, 4), size(cc_wrap%tracer_map%nuopc_to_cc))
               if (cc_wrap%tracer_map%entry_kind(v) /= TRACER_CHEMICAL) cycle
               found_index = cc_wrap%tracer_map%nuopc_to_cc(v)
               if (found_index > 0 .and. found_index <= v_cc) then
                  cc_wrap%chem_buf_4d(:,:,:, found_index) = real(fptr4d(:,:,:, v), c_double) * &
                     cc_wrap%tracer_map%host_to_catchem(v)
               end if
            end do
         else
            call ESMF_LogWrite('Missing validated tracer mapping for rank-4 chemistry import', ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         ! Direct pointer mapping to C++ core StateManager via persistent contiguous buffer
         call cc_wrap%catchem_model%bind_unified_chemistry(cc_wrap%chem_buf_4d, rc)
         if (rc /= CC_SUCCESS) return

         ! A host-owned tracer may be a meteorological prerequisite (e.g.,
         ! specific humidity) rather than an active chemical species.  The
         ! YAML mapping declares this relationship, keeping host conventions
         ! out of the mechanism and core process code.
         if (len_trim(field_map%host_tracer_name) > 0 .and. len_trim(field_map%host_tracer_var) > 0) then
            found_index = 0
            if (allocated(cc_wrap%tracer_map%names)) then
               do v = 1, min(size(fptr4d, 4), size(cc_wrap%tracer_map%names))
                  if (trim(cc_wrap%tracer_map%names(v)) == trim(field_map%host_tracer_name)) then
                     found_index = v
                     exit
                  end if
               end do
            end if
            if (found_index <= 0) then
               ! The relationship is declarative, but a host can legitimately
               ! omit an optional prognostic (as the transform test does).
               ! Do not manufacture a meteorological profile: processes that
               ! truly require this field will fail their own contract.
               call ESMF_LogWrite("Configured host tracer is unavailable; leaving met field unbound: " // &
                  trim(field_map%host_tracer_name), ESMF_LOGMSG_INFO, rc=rc)
               rc = ESMF_SUCCESS
               return
            end if
            if (.not. allocated(cc_wrap%met_buf_3d(fidx)%data)) then
               allocate(cc_wrap%met_buf_3d(fidx)%data(size(fptr4d, 1), size(fptr4d, 2), size(fptr4d, 3)))
            end if
            cc_wrap%met_buf_3d(fidx)%data = real(fptr4d(:,:,:,found_index), c_double)
            call cc_wrap%catchem_model%bind_met_3d_axis(trim(field_map%host_tracer_var), &
               cc_wrap%met_buf_3d(fidx)%data, 0, rc)
            if (rc /= CC_SUCCESS) return
            if (allocated(cc_wrap%catchem_model%required_fields) .and. allocated(is_met_set)) then
               met_index = cc_wrap%catchem_model%get_required_met_index(trim(field_map%host_tracer_var))
               if (met_index > 0 .and. met_index <= size(is_met_set)) is_met_set(met_index) = .true.
            end if
         end if

       case default
         call ESMF_LogWrite("Unknown field mapping dimension for: " // trim(field_map%catchem_var), &
            ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE

      end select

   end subroutine transform_field_to_catchem

   ! Transform CATChem state to ESMF field
   !!
   !! \param    catchem_states CATChem container
   !! \param    field_map      Field mapping information
   !! \param field          ESMF field
   !! \param    im             Horizontal dimension
   !! \param    kme            Vertical dimension
   !! \param   rc             ESMF return code
   !!
   subroutine transform_catchem_to_field(cc_wrap, field, field_map, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(field_mapping_type), intent(in) :: field_map
      type(ESMF_Field), intent(inout) :: field
      integer, intent(out) :: rc

      real(ESMF_KIND_R8), pointer :: fptr4d(:,:,:,:), fptr3d(:,:,:), fptr2d(:,:)
      real(c_double), pointer :: cc_species_conc(:,:)
      real(fp), allocatable :: cc_diag_data(:,:,:)
      type(c_ptr) :: raw_species_ptr
      character(len=128), allocatable :: diagnostic_names(:)
      integer :: i, j, k, v, col, ni, nj, nk, kk, nv, found_index
      character(len=256) :: export_msg

      rc = ESMF_SUCCESS

      !TODO: we assume all the export fields are from DiagManager
      call cc_wrap%catchem_model%get_diagnostic_names(diagnostic_names, rc = rc)

      ! Transform based on field mapping
      select case (field_map%dimensions)

         ! 2D  fields
       case (2)
         nullify(fptr2d)
         if(allocated(cc_diag_data)) deallocate(cc_diag_data)
         ! get field pointer
         call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (.not. associated(fptr2d)) return
         !find index of diagnostic name in the format of process_name.field_name
         found_index = cc_wrap%catchem_model%get_diag_index_from_field(field_map%catchem_var)
         if (found_index > 0) then
            call cc_wrap%catchem_model%get_diagnostic(diagnostic_names(found_index), cc_diag_data, rc)
            if (rc == ESMF_SUCCESS .and. allocated(cc_diag_data)) then
               if (size(cc_diag_data, 1) /= size(fptr2d, 1) .or. size(cc_diag_data, 2) /= size(fptr2d, 2)) then
                  call ESMF_LogWrite("Shape mismatch for 2D export field: " // trim(field_map%catchem_var) // "; zeroing field", &
                     ESMF_LOGMSG_WARNING, rc=rc)
                  fptr2d = 0.0_ESMF_KIND_R8
                  rc = ESMF_SUCCESS
               else
                  fptr2d = cc_diag_data(:,:,1)
               end if
            else
               call ESMF_LogWrite("Could not retrieve diagnostic data for: " // trim(diagnostic_names(found_index)) // "; zeroing field", &
                  ESMF_LOGMSG_WARNING, rc=rc)
               fptr2d = 0.0_ESMF_KIND_R8
               rc = ESMF_SUCCESS
            end if

         else
            fptr2d = 0.0_ESMF_KIND_R8  ! Species not found
         end if

         ! 3D  fields
       case (3)
         nullify(fptr3d)
         if(allocated(cc_diag_data)) deallocate(cc_diag_data)
         ! get field pointer
         call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (.not. associated(fptr3d)) return
         !find index of diagnostic name in the format of process_name.field_name
         found_index = cc_wrap%catchem_model%get_diag_index_from_field(field_map%catchem_var)
         if (found_index > 0) then
            call cc_wrap%catchem_model%get_diagnostic(diagnostic_names(found_index), cc_diag_data, rc)
            if (rc == ESMF_SUCCESS .and. allocated(cc_diag_data)) then
               ni = size(fptr3d, 1)
               nj = size(fptr3d, 2)
               nk = size(fptr3d, 3)
               if (size(cc_diag_data, 1) /= ni .or. size(cc_diag_data, 2) /= nj .or. size(cc_diag_data, 3) /= nk) then
                  call ESMF_LogWrite("Shape mismatch for 3D export field: " // trim(field_map%catchem_var) // "; zeroing field", &
                     ESMF_LOGMSG_WARNING, rc=rc)
                  fptr3d = 0.0_ESMF_KIND_R8
                  rc = ESMF_SUCCESS
               else
                  do k = 1, nk
                     kk = k
                     do j = 1, nj
                        do i = 1, ni
                           fptr3d(i,j,kk) = cc_diag_data(i,j,k)
                        end do
                     end do
                  end do
               end if
            else
               call ESMF_LogWrite("Could not retrieve diagnostic data for: " // trim(diagnostic_names(found_index)) // "; zeroing field", &
                  ESMF_LOGMSG_WARNING, rc=rc)
               fptr3d = 0.0_ESMF_KIND_R8
               rc = ESMF_SUCCESS
            end if

         else
            fptr3d = 0.0_ESMF_KIND_R8  ! Species not found
         end if

         ! 4D  fields for chemical species
       case (4)
         nullify(fptr4d)
         if(allocated(cc_diag_data)) deallocate(cc_diag_data)
         ! get field pointer
         call ESMF_FieldGet(field, farrayPtr=fptr4d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (.not. associated(fptr4d)) return

         write(export_msg, '(A,A,A,I0,A,I0,A,I0,A,I0)') &
            'CATChem export field ', trim(field_map%standard_name), ' shape=', &
            size(fptr4d,1), 'x', size(fptr4d,2), 'x', size(fptr4d,3), 'x', size(fptr4d,4)
         call ESMF_LogWrite(trim(export_msg), ESMF_LOGMSG_INFO, rc=rc)
         if (rc /= ESMF_SUCCESS) return

         ! As above, avoid relying on logical short-circuiting when checking
         ! an allocatable buffer before querying its shape.
         if (allocated(cc_wrap%host_tracer_buf_4d)) then
            if (all(shape(cc_wrap%host_tracer_buf_4d) == shape(fptr4d))) then
               fptr4d = real(cc_wrap%host_tracer_buf_4d, ESMF_KIND_R8)
            end if
         end if

         ni = size(fptr4d, 1)
         nj = size(fptr4d, 2)
         nk = size(fptr4d, 3)
         nv = size(fptr4d, 4)

         ! Copy updated concentrations from the live C++ ChemState buffer back to ESMF tracer array.
         ! The imported 4D buffer is only a staging source for bind_unified_chemistry; process kernels update
         ! ChemState in-place, so export must read the C++ state rather than replaying the original staging buffer.
         if (allocated(cc_wrap%tracer_map%nuopc_to_cc)) then
            do v = 1, min(nv, size(cc_wrap%tracer_map%nuopc_to_cc))
               if (cc_wrap%tracer_map%entry_kind(v) /= TRACER_CHEMICAL) cycle
               found_index = cc_wrap%tracer_map%nuopc_to_cc(v)
               if (found_index <= 0) cycle
               if (catchem_state_get_species_conc_pointer_checked(cc_wrap%catchem_model%state_mgr_ptr, &
                  int(found_index, c_int), int(ni * nj, c_int), int(nk, c_int), raw_species_ptr) /= 0_c_int) cycle
               call c_f_pointer(raw_species_ptr, cc_species_conc, [ni * nj, nk])
               do k = 1, nk
                  do j = 1, nj
                     do i = 1, ni
                        col = i + (j - 1) * ni
                        fptr4d(i,j,k,v) = cc_species_conc(col,k) * cc_wrap%tracer_map%catchem_to_host(v)
                     end do
                  end do
               end do
               nullify(cc_species_conc)
            end do
         else
            call ESMF_LogWrite('Missing validated tracer mapping for rank-4 chemistry export', ESMF_LOGMSG_ERROR, rc=rc)
            rc = ESMF_FAILURE
            return
         end if

         if (allocated(cc_wrap%tracer_map%names)) then
            do v = 1, min(nv, size(cc_wrap%tracer_map%names))
               if (is_pm_diagnostic_name(cc_wrap%tracer_map%names(v))) then
                  found_index = cc_wrap%catchem_model%get_diag_index_from_field( &
                     trim(lowercase(cc_wrap%tracer_map%names(v))))
                  if (found_index > 0) then
                     if (allocated(cc_diag_data)) deallocate(cc_diag_data)
                     call cc_wrap%catchem_model%get_diagnostic(diagnostic_names(found_index), cc_diag_data, rc)
                     if (rc /= ESMF_SUCCESS) then
                        call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, &
                           msg="Failed to get diagnostic data for: " // trim(diagnostic_names(found_index)), &
                           line=__LINE__, file=__FILE__, rcToReturn=rc)
                        return
                     end if
                     fptr4d(:,:,:,v) = cc_diag_data(:,:,:)
                  end if
               end if
            end do   !nv
         end if

       case default
         call ESMF_LogWrite("Unknown export field dimension for: "//trim(field_map%catchem_var), &
            ESMF_LOGMSG_WARNING, rc=rc)

      end select

   end subroutine transform_catchem_to_field

   !> \brief Write diagnostic output using AQMIO
   !!
   !! This subroutine writes CATChem diagnostic fields to NetCDF files using the
   !! working AQMIO module. It retrieves field data, metadata (description, units),
   !! and configuration from the DiagnosticManager and creates properly documented
   !! NetCDF output files.
   !!
   !! \param current_time Current simulation time
   !! \param rc Return code
   subroutine catchem_diagnostics_write(cc_wrap, current_time, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap
      type(ESMF_Time) :: time_on_file
      logical :: time_to_write
      character(len=256) :: filename

      rc = CC_SUCCESS

      ! Check top-level diagnostics/output/enabled switch before doing anything
      if (.not. cc_wrap%catchem_model%is_diag_enabled()) then
         return
      end if

      ! Initialize output timing if not done
      if (.not. cc_wrap%output_timing_initialized) then
         call initialize_output_timing(cc_wrap, current_time, rc)
         if (rc /= CC_SUCCESS) return
      end if

      ! Check if it's time to write output
      call check_diagnostic_output_time(cc_wrap, current_time, time_to_write, time_on_file, rc)
      if (rc /= CC_SUCCESS) return
      if (.not. time_to_write) return

      ! Use grid (must be set during initialization)
      if (.not. ESMF_GridIsCreated(cc_wrap%grid)) then
         rc = CC_FAILURE
         write(*,'(A)') 'Error: grid not initialized.'
         return
      end if

      ! Generate filename for current time
      call generate_diagnostic_filename(cc_wrap, time_on_file, filename, rc)
      if (rc /= CC_SUCCESS) return

      ! Update time variable in NetCDF file and get the time slice to use
      call update_time_variable(cc_wrap, filename, time_on_file, cc_wrap%current_time_slice, rc)
      if (rc /= CC_SUCCESS) return

      !write extemission fields if needed
      call catchem_emis_write_diagnostics(cc_wrap%ext_emis, cc_wrap%current_time_slice, cc_wrap%iocomp, cc_wrap%grid, filename, rc)
      if (rc /= CC_SUCCESS) then
         write(*,'(A)') 'Error: Failed to write external emission diagnostics.'
         return
      end if

      ! Write chemical species concentration diagnostics
      call write_chem_diagnostics(cc_wrap, filename, rc)
      if (rc /= CC_SUCCESS) then
         write(*,'(A)') 'Error: Failed to write chemical species diagnostics.'
         return
      end if

      ! Write per-process diagnostic variables (dust/seasalt emissions,
      ! fluxes, thresholds) when diagnostics.output/process_diagnostics is on.
      call write_process_diagnostics(cc_wrap, 'all', filename, rc)
      if (rc /= CC_SUCCESS) then
         write(*,'(A)') 'Warning: Failed to write process diagnostics.'
         rc = CC_SUCCESS
      end if

      ! Update last output time
      cc_wrap%last_output_time = current_time

      write(*,'(A,A)') 'CATChem: Wrote diagnostic output to ', trim(filename)

   end subroutine catchem_diagnostics_write

   !> \brief Write diagnostics for a specific process
   !!
   !! Discovers the fields the C++ process layer registered in the
   !! DiagnosticManager and writes the dust_ and seasalt_ prefixed ones to
   !! the NetCDF diagnostic file.  Per-column fields (registered shape
   !! [ncols, 1]) are written as 2D (nx, ny) variables; per-bin fields
   !! (shape [ncols, nbin]) are written as a single 3D (nx, ny, nbin)
   !! variable whose third axis follows the mechanism's dust/seasalt bin
   !! order (the bin species names are appended to the description).
   !!
   !! The diagnostic storage is column-major with the same flattened column
   !! index (col = i + (j-1)*nx) the science bridges use, so the [ncols, n]
   !! buffer reinterprets directly as (nx, ny[, nbin]) with no reshaping.
   !!
   !! Gated by diagnostics/output/enabled (checked by the caller) and
   !! diagnostics/output/process_diagnostics.
   !!
   !! \param cc_wrap CATChem wrapper containing model state and configuration
   !! \param process_name Name of the process ('all' selects every process)
   !! \param filename Output filename
   !! \param rc Return code
   subroutine write_process_diagnostics(cc_wrap, process_name, filename, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      ! Local variables
      integer(c_int) :: c_status, c_count, i, rank
      integer(c_int) :: dims(3), n_dust, n_seasalt, n_species
      integer :: nx, ny, k
      character(kind=c_char) :: c_name(64), c_species_name(64)
      character(kind=c_char) :: c_units(32), c_desc(256)
      character(len=64) :: field_name, species_name
      character(len=32) :: units_str
      character(len=256) :: desc_str, bin_list
      type(c_ptr) :: raw_ptr
      real(fp), pointer :: ptr_2d(:,:) => null()
      real(fp), pointer :: ptr_3d(:,:,:) => null()
      logical :: is_dust_field, is_seasalt_field

      rc = CC_SUCCESS

      ! Gate: per-process diagnostic output is opt-in.
      if (.not. cc_wrap%catchem_model%is_process_diag_enabled()) return

      nx = cc_wrap%catchem_model%nx
      ny = cc_wrap%catchem_model%ny

      ! Bring every device-resident diagnostic view back to the host once
      ! before reading raw pointers (no-op in host-only builds).
      call catchem_diag_sync_to_host(cc_wrap%catchem_model%cpp_core_ptr)

      c_status = catchem_diag_get_count_checked(cc_wrap%catchem_model%cpp_core_ptr, c_count)
      if (c_status /= 0_c_int) then
         rc = CC_FAILURE
         return
      end if

      ! Bin counts used to validate per-bin field extents, from the loaded
      ! species mechanism metadata.
      n_dust = 0
      n_seasalt = 0
      n_species = 0
      c_status = catchem_state_get_species_count_checked(cc_wrap%catchem_model%state_mgr_ptr, n_species)
      if (c_status == 0_c_int) then
         do i = 1, n_species
            if (catchem_state_is_species_dust(cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int)) /= 0) then
               n_dust = n_dust + 1
            end if
            if (catchem_state_is_species_seasalt(cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int)) /= 0) then
               n_seasalt = n_seasalt + 1
            end if
         end do
      end if

      do i = 0, c_count - 1
         c_status = catchem_diag_get_name_at_checked(cc_wrap%catchem_model%cpp_core_ptr, &
            int(i, c_int), c_name, 64_c_int)
         if (c_status /= 0_c_int) cycle
         call catchem_c_string_to_fortran(c_name, field_name)

         is_dust_field = (index(field_name, 'dust_') == 1)
         is_seasalt_field = (index(field_name, 'seasalt_') == 1)
         if (.not. (is_dust_field .or. is_seasalt_field)) cycle
         if (trim(process_name) /= 'all') then
            if (.not. ((trim(process_name) == 'dust' .and. is_dust_field) .or. &
               (trim(process_name) == 'seasalt' .and. is_seasalt_field))) cycle
         end if

         rank = 0
         c_status = catchem_diag_get_rank_checked(cc_wrap%catchem_model%cpp_core_ptr, &
            trim(field_name) // c_null_char, rank)
         if (c_status /= 0_c_int .or. rank /= 2_c_int) cycle

         dims = 0
         c_status = catchem_diag_get_dims_checked(cc_wrap%catchem_model%cpp_core_ptr, &
            trim(field_name) // c_null_char, dims, 3_c_int)
         if (c_status /= 0_c_int) cycle
         ! Process diagnostics are flattened over columns: [ncols, 1] totals
         ! or [ncols, nbin] per-bin arrays.  Anything else (e.g. a 3D field
         ! a host registered under the same prefix) is not ours to write.
         if (dims(1) /= nx * ny .or. dims(2) < 1) then
            write(*,'(A,A)') 'Warning: Skipping process diagnostic with unexpected shape: ', trim(field_name)
            cycle
         end if
         if (is_dust_field .and. dims(2) > 1 .and. dims(2) /= n_dust) then
            write(*,'(A,A)') 'Warning: Skipping process diagnostic, bin count mismatch: ', trim(field_name)
            cycle
         end if
         if (is_seasalt_field .and. dims(2) > 1 .and. dims(2) /= n_seasalt) then
            write(*,'(A,A)') 'Warning: Skipping process diagnostic, bin count mismatch: ', trim(field_name)
            cycle
         end if

         units_str = ''
         c_status = catchem_diag_get_units_checked(cc_wrap%catchem_model%cpp_core_ptr, &
            trim(field_name) // c_null_char, c_units, 32_c_int)
         if (c_status == 0_c_int) call catchem_c_string_to_fortran(c_units, units_str)
         desc_str = ''
         c_status = catchem_diag_get_description_checked(cc_wrap%catchem_model%cpp_core_ptr, &
            trim(field_name) // c_null_char, c_desc, 256_c_int)
         if (c_status == 0_c_int) call catchem_c_string_to_fortran(c_desc, desc_str)

         raw_ptr = c_null_ptr
         c_status = catchem_diag_get_pointer_checked(cc_wrap%catchem_model%cpp_core_ptr, &
            trim(field_name) // c_null_char, 2_c_int, dims, raw_ptr)
         if (c_status /= 0_c_int .or. .not. c_associated(raw_ptr)) then
            write(*,'(A,A)') 'Warning: Could not map process diagnostic storage: ', trim(field_name)
            cycle
         end if

         if (dims(2) == 1) then
            ! Per-column field: reinterpret [ncols,1] as (nx,ny) 2D output.
            call c_f_pointer(raw_ptr, ptr_2d, [nx, ny])
            call write_diagnostic_field(cc_wrap, trim(field_name), DIAG_REAL_2D, 0.0_fp, &
               array_2d_ptr=ptr_2d, description=trim(desc_str), &
               units=trim(units_str), filename=filename, rc=rc)
            nullify(ptr_2d)
         else
            ! Per-bin field: reinterpret [ncols,nbin] as (nx,ny,nbin) 3D
            ! output.  The third axis follows mechanism declaration order, so
            ! list the bin species names in the description for consumers.
            bin_list = ''
            do k = 1, n_species
               c_status = catchem_state_get_species_name_at_checked( &
                  cc_wrap%catchem_model%state_mgr_ptr, int(k, c_int), c_species_name, 64_c_int)
               if (c_status /= 0_c_int) cycle
               call catchem_c_string_to_fortran(c_species_name, species_name)
               if (is_dust_field) then
                  if (catchem_state_is_species_dust(cc_wrap%catchem_model%state_mgr_ptr, &
                     int(k, c_int)) == 0) cycle
               else
                  if (catchem_state_is_species_seasalt(cc_wrap%catchem_model%state_mgr_ptr, &
                     int(k, c_int)) == 0) cycle
               end if
               if (len_trim(bin_list) > 0) bin_list = trim(bin_list) // ','
               bin_list = trim(bin_list) // trim(species_name)
            end do
            if (len_trim(bin_list) > 0) then
               desc_str = trim(desc_str) // ' [bins: ' // trim(bin_list) // ']'
            end if
            call c_f_pointer(raw_ptr, ptr_3d, [nx, ny, int(dims(2))])
            call write_diagnostic_field(cc_wrap, trim(field_name), DIAG_REAL_3D, 0.0_fp, &
               array_3d_ptr=ptr_3d, description=trim(desc_str), &
               units=trim(units_str), filename=filename, rc=rc)
            nullify(ptr_3d)
         end if

         if (rc /= CC_SUCCESS) then
            write(*,'(A,A)') 'Warning: Failed to write process diagnostic: ', trim(field_name)
            rc = CC_SUCCESS
         end if
      end do

   end subroutine write_process_diagnostics

   !> \brief Write individual diagnostic field to NetCDF
   !!
   !! \param field_name Field name for NetCDF variable
   !! \param data_type Type of diagnostic data
   !! \param scalar_value Scalar value (if applicable)
   !! \param array_1d_ptr 1D array pointer (if applicable)
   !! \param array_2d_ptr 2D array pointer (if applicable)
   !! \param array_3d_ptr 3D array pointer (if applicable)
   !! \param description Field description for metadata
   !! \param units Field units for metadata
   !! \param filename Output filename
   !! \param rc Return code
   subroutine write_diagnostic_field(cc_wrap, field_name, data_type, scalar_value, &
      array_1d_ptr, array_2d_ptr, array_3d_ptr, &
      description, units, filename, rc)

      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: data_type
      real(fp), intent(in) :: scalar_value
      real(fp), pointer, optional, intent(in) :: array_1d_ptr(:)
      real(fp), pointer, optional, intent(in) :: array_2d_ptr(:,:)
      real(fp), pointer, optional, intent(in) :: array_3d_ptr(:,:,:)
      character(len=*), intent(in) :: description
      character(len=*), intent(in) :: units
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      type(ESMF_Field) :: esmf_field
      type(ESMF_Info) :: info
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      real(ESMF_KIND_R4), pointer :: field_data_3d(:,:,:) => null()
      integer :: i, j, k, time_slice

      ! Only the 2D/3D array kinds are written today; the scalar and 1D inputs
      ! are part of the generic writer signature and intentionally unused.
      ! (array_1d_ptr is an optional pointer, so it may only be referenced via
      ! present(); scalar_value is a plain intent(in) value.)
      associate(unused_scalar => scalar_value); end associate
      if (present(array_1d_ptr)) continue

      rc = CC_SUCCESS

      ! Get current time slice - this will be the same for all fields in this diagnostic write
      time_slice = cc_wrap%current_time_slice

      ! Create appropriate ESMF field based on data type
      select case (data_type)
       case (DIAG_REAL_2D)
         if (.not. present(array_2d_ptr)) then
            rc = CC_FAILURE
            return
         end if
         if (.not. associated(array_2d_ptr)) then
            rc = CC_FAILURE
            return
         end if
         esmf_field = ESMF_FieldCreate(cc_wrap%grid, &
            name=trim(field_name), &
            typekind=ESMF_TYPEKIND_R4, &
            rc=rc)
         if (rc /= ESMF_SUCCESS) return

         !set some info for the field
         call ESMF_InfoGetFromHost(esmf_field, info, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out

         call ESMF_InfoSet(info, "units", trim(units), rc=rc)
         call ESMF_InfoSet(info, "description", trim(description), rc=rc)

         !set values
         call ESMF_FieldGet(esmf_field, farrayPtr=field_data_2d, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         do j = 1, size(array_2d_ptr, 2)
            do i = 1, size(array_2d_ptr, 1)
               field_data_2d(i, j) = real(array_2d_ptr(i, j), ESMF_KIND_R4)
            end do
         end do
         call AQMIO_Write(cc_wrap%iocomp, (/esmf_field/), timeSlice=time_slice, compressLev=cc_wrap%compress_lev, &
            fileName=trim(filename), iofmt=AQMIO_FMT_NETCDF, rc=rc)

       case (DIAG_REAL_3D)
         if (.not. present(array_3d_ptr)) then
            rc = CC_FAILURE
            return
         end if
         if (.not. associated(array_3d_ptr)) then
            rc = CC_FAILURE
            return
         end if
         esmf_field = ESMF_FieldCreate(cc_wrap%grid, &
            name=trim(field_name), &
            typekind=ESMF_TYPEKIND_R4, &
            ungriddedLBound=(/1/), &
            ungriddedUBound=(/size(array_3d_ptr, 3)/), &
            rc=rc)
         if (rc /= ESMF_SUCCESS) return

         !set some info for the field
         call ESMF_InfoGetFromHost(esmf_field, info, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out

         call ESMF_InfoSet(info, "units", trim(units), rc=rc)
         call ESMF_InfoSet(info, "description", trim(description), rc=rc)

         !set values
         call ESMF_FieldGet(esmf_field, farrayPtr=field_data_3d, rc=rc)
         if (rc /= ESMF_SUCCESS) return
         do k = 1, size(array_3d_ptr, 3)
            do j = 1, size(array_3d_ptr, 2)
               do i = 1, size(array_3d_ptr, 1)
                  field_data_3d(i, j, k) = real(array_3d_ptr(i, j, k), ESMF_KIND_R4)
               end do
            end do
         end do
         call AQMIO_Write(cc_wrap%iocomp, (/esmf_field/), timeSlice=time_slice, compressLev=cc_wrap%compress_lev, &
            fileName=trim(filename), iofmt=AQMIO_FMT_NETCDF, rc=rc)

       case default
         rc = CC_FAILURE
         return
      end select

      ! TODO: Add NetCDF attributes for description and units
      ! This would require extending AQMIO or using NetCDF directly
      ! For now, we rely on the working AQMIO functionality

      ! Clean up
      if (ESMF_FieldIsCreated(esmf_field)) then
         call ESMF_FieldDestroy(esmf_field, rc=rc)
      end if

   end subroutine write_diagnostic_field

   !> \brief Write chemical species diagnostics to NetCDF file
   !!
   !! This function saves chemical species concentrations as diagnostic output
   !! based on the diag_species configuration. It handles both individual species
   !! and the 'All' option to save all available species. Units are properly
   !! converted - aerosols from ug/kg to ug/m3 using air density, and gases
   !! are output in ppm.
   !!
   !! \param cc_wrap CATChem wrapper containing model state and configuration
   !! \param filename Output NetCDF filename
   !! \param rc Return code
   subroutine write_chem_diagnostics(cc_wrap, filename, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      ! Local variables
      character(len=64), allocatable :: diag_species(:)
      integer :: num_diag_species, i, species_idx, num_total_species, dims(3)
      character(len=64) :: species_name, field_name, units_str
      character(kind=c_char) :: c_species_name(64)
      character(len=128) :: description
      logical :: save_all_species, is_gas, is_aerosol
      real(fp), pointer :: conc_data(:,:,:) => null()
      real(fp), pointer :: converted_conc(:,:,:) => null()
      real(fp), pointer :: air_density(:,:,:) => null()
      type(c_ptr) :: raw_airden_ptr
      integer(c_int) :: catchem_status
      integer(c_int) :: species_count, species_index_out, gas_value, aerosol_value

      ! Initialize return code
      rc = CC_SUCCESS

      if (.not. cc_wrap%catchem_model%is_diag_enabled()) then
         ! Chemistry diagnostics not enabled, skip
         return
      end if

      dims = [cc_wrap%catchem_model%nx, cc_wrap%catchem_model%ny, cc_wrap%catchem_model%nz]

      ! Get diagnostic species configuration
      num_diag_species = cc_wrap%catchem_model%get_diag_species_count()
      if (num_diag_species > 0) then
         allocate(diag_species(num_diag_species))
         do i = 1, num_diag_species
            call cc_wrap%catchem_model%get_diag_species_at(i, diag_species(i))
         end do
      else
         ! No species configured for diagnostics
         return
      end if

      ! Check if 'All' species should be saved
      save_all_species = .false.
      if (num_diag_species > 0) then
         if (trim(diag_species(1)) == 'All' .or. trim(diag_species(1)) == 'ALL') then
            save_all_species = .true.
         end if
      end if

      ! Get air density field for unit conversion
      catchem_status = catchem_state_get_pointer_3d_checked( &
         cc_wrap%catchem_model%state_mgr_ptr, "AIRDEN" // c_null_char, raw_airden_ptr)
      if (catchem_status /= 0_c_int) then
         catchem_status = catchem_state_get_pointer_3d_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, "AIRDEN_DRY" // c_null_char, raw_airden_ptr)
      end if
      if (c_associated(raw_airden_ptr)) then
         call c_f_pointer(raw_airden_ptr, air_density, dims)
      end if

      catchem_status = catchem_state_get_species_count_checked( &
         cc_wrap%catchem_model%state_mgr_ptr, species_count)
      if (catchem_status /= 0_c_int) then
         rc = CC_FAILURE
         return
      end if
      num_total_species = int(species_count)

      if (save_all_species) then
         ! Save all available chemical species
         do i = 1, num_total_species
            catchem_status = catchem_state_get_species_name_at_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int), c_species_name, 64_c_int)
            if (catchem_status == 0_c_int) catchem_status = catchem_state_is_species_gas_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int), gas_value)
            if (catchem_status == 0_c_int) catchem_status = catchem_state_is_species_aerosol_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int), aerosol_value)
            if (catchem_status /= 0_c_int) cycle
            call catchem_c_string_to_fortran(c_species_name, species_name)
            field_name = 'conc_' // trim(species_name)

            is_gas = (gas_value /= 0_c_int)
            is_aerosol = (aerosol_value /= 0_c_int)

            if (is_gas) then
               units_str = 'ppm'
               description = 'Gas phase concentration of ' // trim(species_name)
            else if (is_aerosol) then
               units_str = 'ug/m3'
               description = 'Aerosol mass concentration of ' // trim(species_name)
            else
               cycle
            end if

            call cc_wrap%catchem_model%get_species_conc_ptr(i, conc_data, dims, rc)
            if (rc /= 0 .or. .not. associated(conc_data)) cycle

            if (is_aerosol .and. associated(air_density)) then
               allocate(converted_conc(dims(1), dims(2), dims(3)))
               converted_conc = conc_data * air_density
               call write_diagnostic_field(cc_wrap, field_name, DIAG_REAL_3D, 0.0_fp, &
                  array_3d_ptr=converted_conc, description=trim(description), &
                  units=trim(units_str), filename=filename, rc=rc)
               deallocate(converted_conc)
            else
               call write_diagnostic_field(cc_wrap, field_name, DIAG_REAL_3D, 0.0_fp, &
                  array_3d_ptr=conc_data, description=trim(description), &
                  units=trim(units_str), filename=filename, rc=rc)
            end if

            if (rc /= CC_SUCCESS) then
               write(*,'(A,A)') 'Warning: Failed to write diagnostics for species: ', trim(species_name)
               rc = CC_SUCCESS
            end if
         end do

      else
         ! Save only specified species
         do i = 1, num_diag_species
            species_name = trim(diag_species(i))
            catchem_status = catchem_state_get_species_index_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, trim(species_name) // c_null_char, species_index_out)
            species_idx = int(species_index_out)

            if (catchem_status /= 0_c_int .or. species_idx <= 0) then
               write(*,'(A,A)') 'Warning: Requested diagnostic species not found: ', trim(species_name)
               cycle
            end if

            field_name = 'conc_' // trim(species_name)

            catchem_status = catchem_state_is_species_gas_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, int(species_idx, c_int), gas_value)
            if (catchem_status == 0_c_int) catchem_status = catchem_state_is_species_aerosol_checked( &
               cc_wrap%catchem_model%state_mgr_ptr, int(species_idx, c_int), aerosol_value)
            if (catchem_status /= 0_c_int) cycle
            is_gas = (gas_value /= 0_c_int)
            is_aerosol = (aerosol_value /= 0_c_int)

            if (is_gas) then
               units_str = 'ppm'
               description = 'Gas phase concentration of ' // trim(species_name)
            else if (is_aerosol) then
               units_str = 'ug/m3'
               description = 'Aerosol mass concentration of ' // trim(species_name)
            else
               write(*,'(A,A)') 'Warning: Species is neither gas nor aerosol: ', trim(species_name)
               cycle
            end if

            call cc_wrap%catchem_model%get_species_conc_ptr(species_idx, conc_data, dims, rc)
            if (rc /= 0 .or. .not. associated(conc_data)) then
               write(*,'(A,A)') 'Warning: Concentration data not available for species: ', trim(species_name)
               cycle
            end if

            if (is_aerosol .and. associated(air_density)) then
               allocate(converted_conc(dims(1), dims(2), dims(3)))
               converted_conc = conc_data * air_density
               call write_diagnostic_field(cc_wrap, field_name, DIAG_REAL_3D, 0.0_fp, &
                  array_3d_ptr=converted_conc, description=trim(description), &
                  units=trim(units_str), filename=filename, rc=rc)
               deallocate(converted_conc)
            else
               call write_diagnostic_field(cc_wrap, field_name, DIAG_REAL_3D, 0.0_fp, &
                  array_3d_ptr=conc_data, description=trim(description), &
                  units=trim(units_str), filename=filename, rc=rc)
            end if

            if (rc /= CC_SUCCESS) then
               write(*,'(A,A)') 'Warning: Failed to write diagnostics for species: ', trim(species_name)
               rc = CC_SUCCESS
            end if
         end do
      end if

   end subroutine write_chem_diagnostics

   !> \brief PM mass weight for a single aerosol species/bin
   !!
   !! \details
   !! Returns the fractional contribution of an aerosol species (or size bin)
   !! to a given particulate-matter size class (PM2.5 or PM10). This mirrors
   !! the weighted-sum approach of the GOCART UFS Aerosol_Diag_Mod ComputePM /
   !! PMGetTracerWeight routine, but is keyed on CATChem per-bin species
   !! short_names (dust1..dust5, seas1..seas5, so4, bc1/bc2, oc1/oc2, and the
   !! NO3an* nitrate aerosols) instead of contiguous tracer indices.
   !!
   !! A weight of 0 means the species does not contribute to that size class.
   !!
   !! \param name    Aerosol species short_name (e.g. 'dust2', 'seas3', 'so4')
   !! \param pm_size Size class string: 'PM25' or 'PM10'
   !! \return w      Mass weight (dimensionless multiplier)
   function pm_tracer_weight(name, pm_size) result(w)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: pm_size
      real(fp) :: w

      ! Partial-bin mass fractions (log-ratio of size cutoff to bin upper edge),
      ! taken directly from the GOCART Aerosol_Diag_Mod PMGetTracerWeight routine.
      real(fp), parameter :: one        = 1.0_fp
      real(fp), parameter :: w25_du2    = log(1.250_fp) / log(1.8_fp)
      real(fp), parameter :: w_du4      = log(1.667_fp) / log(2.0_fp)
      real(fp), parameter :: w25_ss3    = log(2.50_fp)  / log(3.0_fp)
      real(fp), parameter :: w_so4      = 132.14_fp / 96.06_fp
      real(fp), parameter :: w_no3      = 80.043_fp / 62.0_fp
      real(fp), parameter :: w10_no3an2 = 0.808_fp * w_no3
      real(fp), parameter :: w25_no3an2 = 0.138_fp * w_no3
      real(fp), parameter :: w10_no3an3 = 0.164_fp * w_no3

      logical :: is25

      w = 0.0_fp
      is25 = (trim(pm_size) == 'PM25')

      select case (trim(name))
         ! --- Mineral dust (5 bins) ---
       case ('dust1', 'DUST1')
         w = one                                   ! fully in PM2.5 and PM10
       case ('dust2', 'DUST2')
         if (is25) then
            w = w25_du2                            ! partial in PM2.5
         else
            w = one
         end if
       case ('dust3', 'DUST3')
         if (.not. is25) w = one                   ! PM10 only
       case ('dust4', 'DUST4')
         if (.not. is25) w = w_du4                 ! partial in PM10
       case ('dust5', 'DUST5')
         w = 0.0_fp                                ! coarser than PM10

         ! --- Sea salt (5 bins) ---
       case ('seas1', 'SEAS1', 'seas2', 'SEAS2')
         w = one
       case ('seas3', 'SEAS3')
         if (is25) then
            w = w25_ss3                            ! partial in PM2.5
         else
            w = one
         end if
       case ('seas4', 'SEAS4')
         if (.not. is25) w = one                   ! PM10 only
       case ('seas5', 'SEAS5')
         w = 0.0_fp                                ! coarser than PM10

         ! --- Sulfate ---
       case ('so4', 'SO4')
         w = w_so4                                 ! (NH4)2SO4 mass scaling

         ! --- Nitrate aerosols (present in extended mechanisms) ---
       case ('NO3an1', 'no3an1','NO3AN1')
         w = w_no3
       case ('NO3an2', 'no3an2', 'NO3AN2')
         if (is25) then
            w = w25_no3an2
         else
            w = w10_no3an2
         end if
       case ('NO3an3', 'no3an3', 'NO3AN3')
         w = w10_no3an3                            ! same weight for PM2.5/PM10

         ! --- Carbonaceous aerosols (BC/OC, all fine mode) ---
       case ('bc1', 'bc2', 'oc1', 'oc2', 'BC1', 'BC2', 'OC1', 'OC2')
         w = one

       case default
         w = 0.0_fp
      end select

   end function pm_tracer_weight

   !> \brief Compute 3D PM2.5 and PM10 aerosol mass concentrations
   !!
   !! \details
   !! Computes particulate-matter mass concentrations (ug m-3) as a weighted
   !! sum over aerosol species:  PM = sum_s w_s * conc_s * air_density, where
   !! conc_s is the aerosol mixing ratio (ug kg-1) and air_density is the dry
   !! air density (kg m-3). Weights are obtained from pm_tracer_weight and the
   !! sum runs over every aerosol species in the chemistry state. No vertical
   !! flip is applied (CATChem concentrations and AIRDEN share orientation).
   !!
   !! \param cc_wrap CATChem wrapper containing the model state
   !! \param pm25    (out) allocatable 3D PM2.5 mass concentration (ug m-3)
   !! \param pm10    (out) allocatable 3D PM10 mass concentration (ug m-3)
   !! \param rc      Return code
   subroutine compute_pm_diagnostics(cc_wrap, pm25, pm10, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      real(fp), allocatable, intent(out) :: pm25(:,:,:)
      real(fp), allocatable, intent(out) :: pm10(:,:,:)
      integer, intent(out) :: rc

      ! StateManager owns its met and chemistry buffers as C++ double arrays.
      ! Do not use CATChem_Model%get_species_conc_ptr here: that legacy wrapper
      ! presents the raw C++ storage as real(fp), which is 4-byte by default.
      real(c_double), pointer :: air_density(:,:,:) => null()
      real(c_double), pointer :: conc_data(:,:,:) => null()
      integer :: i, num_total_species, dims(3)
      real(fp) :: w25, w10
      character(len=64) :: species_name
      character(kind=c_char) :: c_species_name(64)
      type(c_ptr) :: raw_airden_ptr, raw_conc_ptr
      integer(c_int) :: catchem_status, species_count, aerosol_value

      rc = CC_SUCCESS

      dims = [cc_wrap%catchem_model%nx, cc_wrap%catchem_model%ny, cc_wrap%catchem_model%nz]

      catchem_status = catchem_state_get_pointer_3d_checked( &
         cc_wrap%catchem_model%state_mgr_ptr, "AIRDEN" // c_null_char, raw_airden_ptr)
      if (catchem_status /= 0_c_int) then
         catchem_status = catchem_state_get_pointer_3d_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, "AIRDEN_DRY" // c_null_char, raw_airden_ptr)
      end if
      if (.not. c_associated(raw_airden_ptr)) then
         write(*,'(A)') 'Error: AIRDEN not available for PM diagnostics'
         rc = CC_FAILURE
         return
      end if
      call c_f_pointer(raw_airden_ptr, air_density, dims)

      allocate(pm25(dims(1), dims(2), dims(3)))
      allocate(pm10(dims(1), dims(2), dims(3)))
      pm25 = 0.0_fp
      pm10 = 0.0_fp

      catchem_status = catchem_state_get_species_count_checked( &
         cc_wrap%catchem_model%state_mgr_ptr, species_count)
      if (catchem_status /= 0_c_int) then
         rc = CC_FAILURE
         return
      end if
      num_total_species = int(species_count)

      do i = 1, num_total_species
         catchem_status = catchem_state_is_species_aerosol_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int), aerosol_value)
         if (catchem_status /= 0_c_int .or. aerosol_value == 0_c_int) cycle

         catchem_status = catchem_state_get_species_name_at_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int), c_species_name, 64_c_int)
         if (catchem_status /= 0_c_int) cycle
         call catchem_c_string_to_fortran(c_species_name, species_name)
         w25 = pm_tracer_weight(trim(species_name), 'PM25')
         w10 = pm_tracer_weight(trim(species_name), 'PM10')
         if (is_exact_zero(w25) .and. is_exact_zero(w10)) cycle

         catchem_status = catchem_state_get_species_conc_pointer_checked( &
            cc_wrap%catchem_model%state_mgr_ptr, int(i, c_int), &
            int(dims(1) * dims(2), c_int), int(dims(3), c_int), raw_conc_ptr)
         if (catchem_status /= 0_c_int .or. .not. c_associated(raw_conc_ptr)) cycle
         call c_f_pointer(raw_conc_ptr, conc_data, dims)

         if (.not. is_exact_zero(w25)) pm25 = pm25 + real(real(w25, c_double) * conc_data * air_density, fp)
         if (.not. is_exact_zero(w10)) pm10 = pm10 + real(real(w10, c_double) * conc_data * air_density, fp)

         nullify(conc_data)
      end do

   end subroutine compute_pm_diagnostics

   !> \brief Register (lazily) and update PM2.5/PM10 diagnostics each timestep
   !!
   !! \details
   !! On first invocation this registers an 'aerosol' diagnostic process in the
   !! DiagnosticManager and creates two 3D fields, 'pm25' and 'pm10'. On every
   !! invocation it recomputes the PM mass concentrations and stores them in the
   !! DiagnosticManager. Storing the fields here makes them available both for
   !! NetCDF file output (via the standard process-diagnostics writer) and for
   !! NUOPC export (via transform_catchem_to_field). Must be called after the
   !! chemistry timestep has run and before the export transform.
   !!
   !! \param cc_wrap CATChem wrapper containing the model state
   !! \param rc      Return code
   subroutine update_pm_diagnostics(cc_wrap, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(out) :: rc

      real(fp), allocatable :: pm25(:,:,:), pm10(:,:,:)
      real(c_double), pointer :: diag_ptr(:,:,:) => null()
      integer :: ni, nj, nk
      integer(c_int) :: diag_dims(3), catchem_status
      type(c_ptr) :: raw_diag_ptr

      rc = CC_SUCCESS

      ! Compute current PM mass concentrations
      call compute_pm_diagnostics(cc_wrap, pm25, pm10, rc)
      if (rc /= CC_SUCCESS) return

      ni = size(pm25, 1)
      nj = size(pm25, 2)
      nk = size(pm25, 3)

      ! Lazily register the PM fields in the C++ diagnostic manager.  Both
      ! fields are rewritten in full (diag_ptr = pm25/pm10 below) on every
      ! step, so they are registered Persistent to skip the blanket
      ! per-step reset that Instantaneous fields pay in begin_timestep().
      if (.not. cc_wrap%pm_diag_registered) then
         call cc_wrap%catchem_model%register_diagnostic('pm25', &
            'PM2.5 aerosol mass concentration', 'ug m-3', (/ni, nj, nk/), rc, persistent=.true.)
         if (rc /= 0) then
            write(*,'(A)') 'Error: could not register pm25 diagnostic'
            rc = CC_FAILURE
            return
         end if
         call cc_wrap%catchem_model%register_diagnostic('pm10', &
            'PM10 aerosol mass concentration', 'ug m-3', (/ni, nj, nk/), rc, persistent=.true.)
         if (rc /= 0) then
            write(*,'(A)') 'Error: could not register pm10 diagnostic'
            rc = CC_FAILURE
            return
         end if
         cc_wrap%pm_diag_registered = .true.
      end if

      ! Write current values into the C++-owned diagnostic storage
      diag_dims = int([ni, nj, nk], c_int)
      catchem_status = catchem_diag_get_pointer_checked(cc_wrap%catchem_model%cpp_core_ptr, &
         'pm25' // c_null_char, 3_c_int, diag_dims, raw_diag_ptr)
      if (catchem_status /= 0_c_int .or. .not. c_associated(raw_diag_ptr)) then
         write(*,'(A)') 'Error: could not map pm25 diagnostic storage'
         rc = CC_FAILURE
         return
      end if
      call c_f_pointer(raw_diag_ptr, diag_ptr, [ni, nj, nk])
      diag_ptr = pm25
      nullify(diag_ptr)

      catchem_status = catchem_diag_get_pointer_checked(cc_wrap%catchem_model%cpp_core_ptr, &
         'pm10' // c_null_char, 3_c_int, diag_dims, raw_diag_ptr)
      if (catchem_status /= 0_c_int .or. .not. c_associated(raw_diag_ptr)) then
         write(*,'(A)') 'Error: could not map pm10 diagnostic storage'
         rc = CC_FAILURE
         return
      end if
      call c_f_pointer(raw_diag_ptr, diag_ptr, [ni, nj, nk])
      diag_ptr = pm10
      nullify(diag_ptr)

      if (allocated(pm25)) deallocate(pm25)
      if (allocated(pm10)) deallocate(pm10)

   end subroutine update_pm_diagnostics

   !> \brief Update time variable in NetCDF file
   !!
   !! This function handles proper time series management by either:
   !! - Reading existing time values and appending a new time step, or
   !! - Creating the time variable for the first time
   !! Uses 1D LocStream-based fields and new AQMIO_Write1D function for proper
   !! 1D field handling without MPI operations.
   !!
   !! \param cc_wrap CATChem wrapper containing ESMF components
   !! \param filename NetCDF filename
   !! \param current_time Current simulation time to append
   !! \param time_slice Returns the time slice number that was written
   !! \param rc Return code
   subroutine update_time_variable(cc_wrap, filename, current_time, time_slice, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      character(len=*), intent(in) :: filename
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: time_slice
      integer, intent(out) :: rc

      ! Local variables for direct data I/O approach
      integer(ESMF_KIND_I4) :: new_time_data(1)
      type(ESMF_Time) :: reference_time
      type(ESMF_TimeInterval) :: time_diff
      integer(ESMF_KIND_I8) :: time_seconds
      type(ESMF_VM) :: vm
      type(ESMF_Grid) :: grid
      integer :: ibuf(1)  ! Buffer for MPI broadcast
      integer :: tileCount, tile
      character(len=256) :: tileFilename
      character(len=16) :: tileSuffix
      integer :: dotpos

      rc = CC_SUCCESS

      ! Create reference time (Unix epoch: 1970-01-01 00:00:00)
      call ESMF_TimeSet(reference_time, yy=1970, mm=1, dd=1, h=0, m=0, s=0, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Convert current time to epoch seconds
      time_diff = current_time - reference_time
      call ESMF_TimeIntervalGet(time_diff, s_i8=time_seconds, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      new_time_data(1) = int(time_seconds, ESMF_KIND_I4)

      ! Store time value for lat/lon coordinate output
      if (latlon_diag_is_init()) call latlon_diag_set_time(new_time_data(1))

      ! Determine tile count to match AQMIO's per-tile file naming
      call ESMF_GridCompGet(cc_wrap%iocomp, grid=grid, rc=rc)
      if (rc /= ESMF_SUCCESS) return
      call ESMF_GridGet(grid, tileCount=tileCount, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      if (tileCount > 1 .and. index(filename, '<tile>') == 0) then
         ! Multi-tile without <tile> placeholder: write time to each per-tile file
         ! Must match AQMIO_FileNameGet auto-tile naming: "file.nc" -> "file.tileN.nc"
         do tile = 1, tileCount
            write(tileSuffix, '(".tile",I0)') tile
            dotpos = index(filename, '.', back=.true.)
            if (dotpos > 1) then
               tileFilename = filename(1:dotpos-1) // trim(tileSuffix) // trim(filename(dotpos:))
            else
               tileFilename = trim(filename) // trim(tileSuffix)
            end if
            call AQMIO_Write1D(tileFilename, "time", append=.true., del_old_file=.true., rc=rc, &
               data_i4=new_time_data, current_size=time_slice, &
               iocomp=cc_wrap%iocomp)
            if (rc /= ESMF_SUCCESS) return
         end do
      else
         ! Single tile or filename has <tile> placeholder
         call AQMIO_Write1D(filename, "time", append=.true., del_old_file=.true., rc=rc, &
            data_i4=new_time_data, current_size=time_slice, &
            iocomp=cc_wrap%iocomp)
         if (rc /= ESMF_SUCCESS) return
      end if

      ! Broadcast time_slice from I/O PET to all other PETs so they have the correct value
      ! Get VM from the IOComp for broadcasting
      call ESMF_GridCompGet(cc_wrap%iocomp, vm=vm, rc=rc)
      if (rc == ESMF_SUCCESS) then
         ! Use buffer for broadcast (ESMF_VMBroadcast expects arrays)
         ibuf(1) = time_slice
         call ESMF_VMBroadcast(vm, ibuf, 1, 0, rc=rc)
         time_slice = ibuf(1)
      end if

   end subroutine update_time_variable

   !> \brief Initialize output timing
   !!
   !! \param start_time Simulation start time
   !! \param rc Return code
   subroutine initialize_output_timing(cc_wrap, start_time, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: start_time
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap

      rc = CC_SUCCESS

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Create output interval from frequency (in seconds)
      call ESMF_TimeIntervalSet(cc_wrap%output_interval, s=cc_wrap%output_frequency, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Set initial output time (start time minus interval so first check will trigger output)
      cc_wrap%last_output_time = start_time - cc_wrap%output_interval

      cc_wrap%output_timing_initialized = .true.

   end subroutine initialize_output_timing

   !> \brief Check if it's time to write diagnostic output
   !!
   !! \param current_time Current simulation time
   !! \param time_to_write True if it's time to write
   !! \param rc Return code
   subroutine check_diagnostic_output_time(cc_wrap, current_time, time_to_write, time_on_file, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: current_time
      logical, intent(out) :: time_to_write
      type(ESMF_Time), intent(out) :: time_on_file
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap
      type(ESMF_Time) :: next_output_time

      rc = CC_SUCCESS
      time_to_write = .false.

      ! Get process-local state
      !cc_wrap => get_cc_wrap()

      ! Calculate next output time
      next_output_time = cc_wrap%last_output_time + cc_wrap%output_interval

      ! Check if current time is at or past next output time
      if (current_time >= next_output_time) then
         time_to_write = .true.
         time_on_file = current_time
      end if

      !Note that the last run cycle is one time step before the end time, so we do not have the last hour saved out.
      !To ensure the last time step is written, we can force output if we are within one time step of the end time.
      ! next_time = current_time + cc_wrap%timeStep
      ! if ( (cc_wrap%endTime == next_time) ) then
      !   time_to_write = .true.
      !   time_on_file = cc_wrap%endTime
      ! end if


      ! !!!!!log to debug
      ! block
      !   integer :: year, month, day, hour, minute, second
      !   character(len=64) :: time_str_current, time_str_next, time_str_end

      !   call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, &
      !                    h=hour, m=minute, s=second, rc=rc)
      !   write(time_str_current, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      !         year, month, day, hour, minute, second

      !   call ESMF_TimeGet(cc_wrap%endTime, yy=year, mm=month, dd=day, &
      !                    h=hour, m=minute, s=second, rc=rc)
      !   write(time_str_next, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
      !         year, month, day, hour, minute, second

      !   write(*,'(A,A,A)') 'Debug: Current time: ', trim(time_str_current), ', End time: ', trim(time_str_next)

      ! end block

   end subroutine check_diagnostic_output_time

   !> \brief Generate filename for diagnostic output
   !!
   !! \param current_time Current simulation time
   !! \param filename Generated filename
   !! \param rc Return code
   subroutine generate_diagnostic_filename(cc_wrap, time_on_file, filename, rc)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      type(ESMF_Time), intent(in) :: time_on_file
      character(len=*), intent(out) :: filename
      integer, intent(out) :: rc

      !type(cc_wrap_type), pointer :: cc_wrap
      type(ESMF_VM) :: vm
      integer :: year, month, day, hour, minute, second, localPet
      character(len=256) :: time_string
      logical :: dir_exists

      rc = CC_SUCCESS

      ! Only have PET 0 create directories to avoid race conditions
      call ESMF_GridCompGet(cc_wrap%iocomp, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_VMGet(vm, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (localPet == 0) then
         ! check if output directory exists, create if not
         inquire(file=trim(cc_wrap%output_directory), exist=dir_exists)
         if (.not. dir_exists) then
            ! Create directory (execute_command_line is the F2008 standard form;
            ! the `system` extension has no explicit interface and ifx warns)
            call execute_command_line('mkdir -p ' // trim(cc_wrap%output_directory))
         end if
      end if

      ! Barrier to ensure directory is created before any PET tries to write
      call ESMF_VMBarrier(vm, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Get time components
      call ESMF_TimeGet(time_on_file, yy=year, mm=month, dd=day, &
         h=hour, m=minute, s=second, rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Create filename: output_directory/prefix_YYYYMMDD_HHMMSS.nc
      write(time_string, '(I4.4,I2.2,I2.2,A,I2.2,I2.2,I2.2)') &
         year, month, day, '_', hour, minute, second

      filename = trim(cc_wrap%output_directory) // '/' // trim(cc_wrap%output_prefix) // '_' // trim(time_string) // '.nc'

      !This will save all hours to the same file
      !filename = trim(cc_wrap%output_directory) // '/' // trim(cc_wrap%output_prefix) // '.nc'

   end subroutine generate_diagnostic_filename


   ! Load field configuration from YAML file
   !!
   !! \param info  tracerinfo from NUOPC
   !! \param key   infor intended to get
   !! \param values  infor values returned
   !! \param rc  return status
   !!
   subroutine TracerInfoGet(info, key, values, rc)
      ! -- interface variables
      type(ESMF_Info),               intent(in)  :: info
      character(len=*),              intent(in)  :: key
      character(len=*), allocatable, intent(out) :: values(:)
      integer,          optional,    intent(out) :: rc

      ! -- local variables
      integer :: localrc
      logical :: isKeyFound

      ! -- begin
      if (present(rc)) rc = ESMF_SUCCESS

      ! -- check if metadata key is present in ESMF_Info object
      isKeyFound = ESMF_InfoIsPresent(info, trim(key), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__, &
         rcToReturn=rc)) return

      if (isKeyFound) then
         isKeyFound = ESMF_InfoIsSet(info, trim(key), rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
      end if

      if (isKeyFound) then
         call ESMF_InfoGetAlloc(info, trim(key), values, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)) return
      end if

   end subroutine TracerInfoGet

   ! Load field configuration from YAML file
   !!
   !! \param  config_file Configuration file path
   !! \param errflg      Error flag
   !! \param errmsg      Error message
   !!
   subroutine load_field_config(config_file, errflg, errmsg)

      character(len=*), intent(in) :: config_file
      integer, intent(out) :: errflg
      character(len=*), intent(out) :: errmsg

      !type(cc_wrap_type), pointer :: cc_wrap

      errflg = CC_SUCCESS
      errmsg = ''

      ! Parse import fields
      call parse_field_section(config_file, 'import_fields', field_config%import_fields, &
         field_config%n_import_fields, errflg, errmsg)
      if (errflg /= CC_SUCCESS) return

      ! Parse export fields
      call parse_field_section(config_file, 'export_fields', field_config%export_fields, &
         field_config%n_export_fields, errflg, errmsg)
      if (errflg /= CC_SUCCESS) return

   end subroutine load_field_config

   !> Parse a field section (import_fields or export_fields) from YAML file
   !!
   !! \param filename YAML configuration file
   !! \param section_name Section name ('import_fields' or 'export_fields')
   !! \param fields Array to store parsed fields
   !! \param n_fields Number of fields found
   !! \param errflg Error flag
   !! \param errmsg Error message
   !!
   subroutine parse_field_section(filename, section_name, fields, n_fields, errflg, errmsg)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: section_name
      type(field_mapping_type), allocatable, intent(out) :: fields(:)
      integer, intent(out) :: n_fields
      integer, intent(out) :: errflg
      character(len=*), intent(out) :: errmsg

      integer :: unit_num, io_stat, indent_level, section_indent
      character(len=256) :: line, trimmed_line, field_name, field_value
      logical :: in_section, found_section
      integer :: line_number, colon_pos
      type(field_mapping_type), allocatable :: temp_fields(:)
      type(field_mapping_type) :: current_field
      logical :: in_field_item, field_already_saved

      errflg = CC_SUCCESS
      errmsg = ''
      n_fields = 0
      in_section = .false.
      found_section = .false.
      section_indent = -1
      line_number = 0
      in_field_item = .false.
      field_already_saved = .false.

      ! Initialize current field
      current_field%standard_name = ''
      current_field%catchem_var = ''
      current_field%dimensions = 0
      current_field%units = ''
      current_field%vertical_axis = 'level'
      current_field%host_tracer_name = ''
      current_field%host_tracer_var = ''
      current_field%optional = .false.
      current_field%advertise = .false.

      ! Open file for reading
      open(newunit=unit_num, file=trim(filename), status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         write(errmsg, '(A,A)') 'Cannot open configuration file: ', trim(filename)
         errflg = CC_FAILURE
         return
      endif

      ! Allocate temporary storage and grow as needed for large mapping files
      allocate(temp_fields(64))

      ! Read file line by line
      do
         read(unit_num, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit  ! End of file or error

         line_number = line_number + 1

         ! Remove inline comments - everything after '#' character
         if (index(line, '#') > 0) then
            line = line(1:index(line, '#')-1)
         endif

         trimmed_line = trim(adjustl(line))

         ! Skip empty lines (comments have already been stripped)
         if (len_trim(trimmed_line) == 0) cycle

         ! Calculate indentation level
         do indent_level = 1, len_trim(line)
            if (line(indent_level:indent_level) /= ' ') exit
         end do
         indent_level = indent_level - 1

         ! Look for section header
         if (index(trimmed_line, ':') > 0) then
            colon_pos = index(trimmed_line, ':')
            field_name = trimmed_line(1:colon_pos-1)
            field_name = trim(adjustl(field_name))

            ! Check if we found our target section
            if (trim(field_name) == trim(section_name) .and. indent_level == 0) then
               in_section = .true.
               found_section = .true.
               section_indent = indent_level
               cycle
            endif

            ! If we're already in a section and encounter another top-level section, exit
            if (in_section .and. indent_level == 0 .and. trim(field_name) /= trim(section_name)) then
               ! Save the last field if we're still processing one
               if (in_field_item .and. current_field%standard_name /= '') then
                  call append_parsed_field(temp_fields, n_fields, current_field, errflg, errmsg)
                  if (errflg /= CC_SUCCESS) then
                     close(unit_num)
                     return
                  end if
                  field_already_saved = .true.  ! Mark that we've saved the field
               endif
               exit
            endif

            ! Process items within the section
            if (in_section .and. indent_level > section_indent) then

               ! Look for array items (lines starting with "- ")
               if (index(trimmed_line, '- ') == 1) then
                  ! Save previous field if we have one
                  if (in_field_item .and. current_field%standard_name /= '') then
                     call append_parsed_field(temp_fields, n_fields, current_field, errflg, errmsg)
                     if (errflg /= CC_SUCCESS) then
                        close(unit_num)
                        return
                     end if
                  endif

                  ! Start new field item
                  in_field_item = .true.
                  current_field%standard_name = ''
                  current_field%catchem_var = ''
                  current_field%dimensions = 0
                  current_field%units = ''
                  current_field%vertical_axis = 'level'
                  current_field%host_tracer_name = ''
                  current_field%host_tracer_var = ''
                  current_field%optional = .false.
                  current_field%advertise = .false.

                  ! Parse the first property if it's on the same line as the dash
                  if (len_trim(trimmed_line) > 2) then
                     trimmed_line = trim(adjustl(trimmed_line(3:)))  ! Remove "- "
                     if (index(trimmed_line, ':') > 0) then
                        colon_pos = index(trimmed_line, ':')
                        field_name = trim(adjustl(trimmed_line(1:colon_pos-1)))
                        field_value = trim(adjustl(trimmed_line(colon_pos+1:)))
                        call parse_field_property(field_name, field_value, current_field)
                     endif
                  endif

               elseif (in_field_item .and. indent_level > section_indent + 2) then
                  ! Parse field properties
                  if (index(trimmed_line, ':') > 0) then
                     colon_pos = index(trimmed_line, ':')
                     field_name = trim(adjustl(trimmed_line(1:colon_pos-1)))
                     field_value = trim(adjustl(trimmed_line(colon_pos+1:)))
                     call parse_field_property(field_name, field_value, current_field)
                  endif
               endif

            elseif (in_section .and. indent_level <= section_indent) then
               ! We've left our section
               if (in_field_item .and. current_field%standard_name /= '') then
                  call append_parsed_field(temp_fields, n_fields, current_field, errflg, errmsg)
                  if (errflg /= CC_SUCCESS) then
                     close(unit_num)
                     return
                  end if
                  field_already_saved = .true.  ! Mark that we've saved the field
               endif
               exit
            endif
         endif
      end do

      ! Save the last field if we're still processing one AND it hasn't been saved yet
      if (in_field_item .and. current_field%standard_name /= '' .and. .not. field_already_saved) then
         call append_parsed_field(temp_fields, n_fields, current_field, errflg, errmsg)
         if (errflg /= CC_SUCCESS) then
            close(unit_num)
            return
         end if
      endif

      close(unit_num)

      ! Check if we found the section
      if (.not. found_section) then
         write(errmsg, '(A,A,A)') 'Section "', trim(section_name), '" not found in configuration file'
         errflg = CC_FAILURE
         deallocate(temp_fields)
         return
      endif

      ! Allocate final array and copy data
      if (n_fields > 0) then
         allocate(fields(n_fields))
         fields(1:n_fields) = temp_fields(1:n_fields)
      endif

      deallocate(temp_fields)

   end subroutine parse_field_section

   !> \brief Append one parsed field to a dynamically growing temporary array.
   subroutine append_parsed_field(temp_fields, n_fields, current_field, errflg, errmsg)
      type(field_mapping_type), allocatable, intent(inout) :: temp_fields(:)
      integer, intent(inout) :: n_fields
      type(field_mapping_type), intent(in) :: current_field
      integer, intent(out) :: errflg
      character(len=*), intent(out) :: errmsg

      type(field_mapping_type), allocatable :: expanded_fields(:)
      integer :: old_size, new_size, stat

      errflg = CC_SUCCESS
      errmsg = ''

      n_fields = n_fields + 1
      if (n_fields > size(temp_fields)) then
         old_size = size(temp_fields)
         new_size = max(old_size * 2, n_fields)
         allocate(expanded_fields(new_size), stat=stat)
         if (stat /= 0) then
            errflg = CC_FAILURE
            errmsg = 'Unable to grow temporary field mapping storage'
            return
         end if
         expanded_fields(1:old_size) = temp_fields(1:old_size)
         call move_alloc(expanded_fields, temp_fields)
      end if

      temp_fields(n_fields) = current_field
   end subroutine append_parsed_field

   !> Parse a field property and set it in the field structure
   !!
   !! \param property_name Name of the property
   !! \param property_value Value of the property
   !! \param field Field structure to update
   !!
   subroutine parse_field_property(property_name, property_value, field)
      character(len=*), intent(in) :: property_name
      character(len=*), intent(in) :: property_value
      type(field_mapping_type), intent(inout) :: field

      character(len=256) :: clean_value
      integer :: read_stat

      ! Remove quotes from string values
      clean_value = property_value
      if (len_trim(clean_value) >= 2) then
         if ((clean_value(1:1) == '"' .and. clean_value(len_trim(clean_value):len_trim(clean_value)) == '"') .or. &
            (clean_value(1:1) == "'" .and. clean_value(len_trim(clean_value):len_trim(clean_value)) == "'")) then
            clean_value = clean_value(2:len_trim(clean_value)-1)
         endif
      endif

      select case (trim(property_name))
       case ('standard_name')
         field%standard_name = trim(clean_value)
       case ('catchem_var')
         field%catchem_var = trim(clean_value)
       case ('dimensions')
         read(clean_value, *, iostat=read_stat) field%dimensions
         if (read_stat /= 0) field%dimensions = 0
       case ('units')
         field%units = trim(clean_value)
       case ('vertical_axis')
         field%vertical_axis = trim(clean_value)
       case ('host_tracer_name')
         field%host_tracer_name = trim(clean_value)
       case ('host_tracer_var')
         field%host_tracer_var = trim(clean_value)
       case ('optional')
         select case (trim(clean_value))
          case ('true', 'True', 'TRUE', '.true.')
            field%optional = .true.
          case ('false', 'False', 'FALSE', '.false.')
            field%optional = .false.
          case default
            field%optional = .false.
         end select
       case ('advertise')
         select case (trim(clean_value))
          case ('true', 'True', 'TRUE', '.true.')
            field%advertise = .true.
          case ('false', 'False', 'FALSE', '.false.')
            field%advertise = .false.
          case default
            field%advertise = .false.
         end select
      end select

   end subroutine parse_field_property

   !> \brief Get number of import fields (MPI-safe accessor)
   !!
   !! \return Number of import fields configured
   function get_n_import_fields(cc_wrap) result(n_fields)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer :: n_fields
      !type(cc_wrap_type), pointer :: cc_wrap

      !cc_wrap => get_cc_wrap()
      n_fields = cc_wrap%field_config%n_import_fields
   end function get_n_import_fields

   !> \brief Get import field information (MPI-safe accessor)
   !!
   !! \param field_index Index of the field (1-based)
   !! \param standard_name NUOPC standard name
   !! \param optional Whether field is optional
   function get_import_field_info(cc_wrap, field_index, standard_name, optional) result(success)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(in) :: field_index
      character(len=*), intent(out) :: standard_name
      logical, intent(out) :: optional
      logical :: success
      !type(cc_wrap_type), pointer :: cc_wrap

      success = .false.
      !cc_wrap => get_cc_wrap()

      if (field_index > 0 .and. field_index <= cc_wrap%field_config%n_import_fields) then
         standard_name = cc_wrap%field_config%import_fields(field_index)%standard_name
         optional = cc_wrap%field_config%import_fields(field_index)%optional
         success = .true.
      end if
   end function get_import_field_info

   !> \brief Get number of export fields (MPI-safe accessor)
   !!
   !! \return Number of export fields configured
   function get_n_export_fields(cc_wrap) result(n_fields)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer :: n_fields
      !type(cc_wrap_type), pointer :: cc_wrap

      !cc_wrap => get_cc_wrap()
      n_fields = cc_wrap%field_config%n_export_fields
   end function get_n_export_fields

   !> \brief Get export field information (MPI-safe accessor)
   !!
   !! \param field_index Index of the field (1-based)
   !! \param standard_name NUOPC standard name
   !! \param optional Whether field is optional
   function get_export_field_info(cc_wrap, field_index, standard_name, optional) result(success)
      type(cc_wrap_type), intent(inout) :: cc_wrap
      integer, intent(in) :: field_index
      character(len=*), intent(out) :: standard_name
      logical, intent(out) :: optional
      logical :: success
      !type(cc_wrap_type), pointer :: cc_wrap

      success = .false.
      !cc_wrap => get_cc_wrap()

      if (field_index > 0 .and. field_index <= cc_wrap%field_config%n_export_fields) then
         standard_name = cc_wrap%field_config%export_fields(field_index)%standard_name
         optional = cc_wrap%field_config%export_fields(field_index)%optional
         success = .true.
      end if
   end function get_export_field_info

end module catchem_nuopc_interface
