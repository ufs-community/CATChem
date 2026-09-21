!> \file CATChem_API.F90
!! \brief Modernized high-level BIND(C) API for host model integration
!!
!! This module provides the standard, backward-compatible CATChem_Model
!! derived type used by Earth System drivers (like NUOPC and UFS).
!! It delegates all memory management, synchronization, and process scheduling
!! directly to the modernized top-down C++ Core via standard BIND(C) bridges.
!!
module CATChem_API
   use iso_c_binding
   use catchem_bridge_precision, only: fp
   use catchem_bridge_error, only: CC_SUCCESS, CC_FAILURE

   implicit none
   private

   public :: CATChem_Model, catchem_diag_register_contract_checked, catchem_diag_get_contract

   !=========================================================================
   ! C-API Interfaces to catchem_api.cpp
   !=========================================================================
   interface
      integer(c_int) function catchem_core_create_from_config_with_grid_checked(config_file, ncols, nlevels, core_out) &
         bind(C, name="catchem_core_create_from_config_with_grid_checked")
         import :: c_char, c_ptr, c_int
         character(kind=c_char), intent(in) :: config_file(*)
         integer(c_int), value :: ncols, nlevels
         type(c_ptr), intent(out) :: core_out
      end function

      integer(c_int) function catchem_core_destroy_checked(core_ptr) bind(C, name="catchem_core_destroy_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_core_get_state_manager_checked(core_ptr, state_out) &
         bind(C, name="catchem_core_get_state_manager_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
         type(c_ptr), intent(out) :: state_out
      end function

      integer(c_int) function catchem_get_last_error(buffer, max_len) bind(C, name="catchem_get_last_error")
         import :: c_char, c_int
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end function

      type(c_ptr) function catchem_core_create_from_config(config_file) bind(C, name="catchem_core_create_from_config")
         import :: c_char, c_ptr
         character(kind=c_char), intent(in) :: config_file(*)
      end function

      type(c_ptr) function catchem_core_create_from_config_with_grid(config_file, ncols, nlevels) &
         bind(C, name="catchem_core_create_from_config_with_grid")
         import :: c_char, c_ptr, c_int
         character(kind=c_char), intent(in) :: config_file(*)
         integer(c_int), value :: ncols, nlevels
      end function

      subroutine catchem_core_destroy(core_ptr) bind(C, name="catchem_core_destroy")
         import :: c_ptr
         type(c_ptr), value :: core_ptr
      end subroutine

      subroutine catchem_register_carbchem_cpp() bind(C, name="catchem_register_carbchem_cpp")
      end subroutine

      subroutine catchem_register_drydep_cpp() bind(C, name="catchem_register_drydep_cpp")
      end subroutine

      subroutine catchem_register_dust_cpp() bind(C, name="catchem_register_dust_cpp")
      end subroutine

      subroutine catchem_register_seasalt_cpp() bind(C, name="catchem_register_seasalt_cpp")
      end subroutine

      subroutine catchem_register_settling_cpp() bind(C, name="catchem_register_settling_cpp")
      end subroutine

      subroutine catchem_register_so4chem_cpp() bind(C, name="catchem_register_so4chem_cpp")
      end subroutine

      subroutine catchem_register_wetdep_cpp() bind(C, name="catchem_register_wetdep_cpp")
      end subroutine

      type(c_ptr) function catchem_core_get_state_manager(core_ptr) bind(C, name="catchem_core_get_state_manager")
         import :: c_ptr
         type(c_ptr), value :: core_ptr
      end function

      subroutine catchem_core_add_process_by_name(core_ptr, name) bind(C, name="catchem_core_add_process_by_name")
         import :: c_ptr, c_char
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
      end subroutine

      integer(c_int) function catchem_core_get_num_processes(core_ptr) bind(C, name="catchem_core_get_num_processes")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_core_get_num_processes_checked(core_ptr, count_out) &
         bind(C, name="catchem_core_get_num_processes_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), intent(out) :: count_out
      end function

      integer(c_int) function catchem_core_get_required_host_field_count_checked(core_ptr, count_out) &
         bind(C, name="catchem_core_get_required_host_field_count_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), intent(out) :: count_out
      end function

      integer(c_int) function catchem_core_get_required_host_field_name_checked(core_ptr, index, name_out, name_out_len) &
         bind(C, name="catchem_core_get_required_host_field_name_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index, name_out_len
         character(kind=c_char), intent(out) :: name_out(*)
      end function

      integer(c_int) function catchem_core_run_timestep(core_ptr, dt) bind(C, name="catchem_core_run_timestep")
         import :: c_ptr, c_double, c_int
         type(c_ptr), value :: core_ptr
         real(c_double), value :: dt
      end function

      subroutine catchem_state_bind_met_3d(state_ptr, name, ptr) bind(C, name="catchem_state_bind_met_3d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
         type(c_ptr), value :: ptr
      end subroutine

      integer(c_int) function catchem_state_bind_met_3d_checked(state_ptr, name, ptr, dim1, dim2, dim3) &
         bind(C, name="catchem_state_bind_met_3d_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr, ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: dim1, dim2, dim3
      end function

      integer(c_int) function catchem_state_bind_met_3d_axis_checked(state_ptr, name, ptr, dim1, dim2, dim3, axis) &
         bind(C, name="catchem_state_bind_met_3d_axis_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr, ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: dim1, dim2, dim3, axis
      end function

      subroutine catchem_state_bind_met_2d(state_ptr, name, ptr) bind(C, name="catchem_state_bind_met_2d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
         type(c_ptr), value :: ptr
      end subroutine

      integer(c_int) function catchem_state_bind_met_2d_checked(state_ptr, name, ptr, dim1, dim2) &
         bind(C, name="catchem_state_bind_met_2d_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr, ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: dim1, dim2
      end function

      subroutine catchem_state_bind_unified_chemistry(state_ptr, ptr) bind(C, name="catchem_state_bind_unified_chemistry")
         import :: c_ptr
         type(c_ptr), value :: state_ptr
         type(c_ptr), value :: ptr
      end subroutine

      integer(c_int) function catchem_state_bind_unified_chemistry_checked(state_ptr, ptr, dim1, dim2, dim3) &
         bind(C, name="catchem_state_bind_unified_chemistry_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr, ptr
         integer(c_int), value :: dim1, dim2, dim3
      end function

      subroutine catchem_state_sync_to_device(state_ptr) bind(C, name="catchem_state_sync_to_device")
         import :: c_ptr
         type(c_ptr), value :: state_ptr
      end subroutine

      integer(c_int) function catchem_state_sync_to_device_checked(state_ptr) &
         bind(C, name="catchem_state_sync_to_device_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      subroutine catchem_state_sync_to_host(state_ptr) bind(C, name="catchem_state_sync_to_host")
         import :: c_ptr
         type(c_ptr), value :: state_ptr
      end subroutine

      integer(c_int) function catchem_state_sync_to_host_checked(state_ptr) &
         bind(C, name="catchem_state_sync_to_host_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      type(c_ptr) function catchem_state_get_species_conc_pointer(state_ptr, index) &
         bind(C, name="catchem_state_get_species_conc_pointer")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_conc_pointer_checked(state_ptr, index, dim1, dim2, ptr_out) &
         bind(C, name="catchem_state_get_species_conc_pointer_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index, dim1, dim2
         type(c_ptr), intent(out) :: ptr_out
      end function

      subroutine catchem_get_grid_dimensions(core_ptr, nx, ny, nz) bind(C, name="catchem_get_grid_dimensions")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), intent(out) :: nx, ny, nz
      end subroutine

      real(c_double) function catchem_get_config_timestep(core_ptr) bind(C, name="catchem_get_config_timestep")
         import :: c_ptr, c_double
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_config_get_output_frequency(core_ptr) &
         bind(C, name="catchem_config_get_output_frequency")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_config_get_compress_level(core_ptr) &
         bind(C, name="catchem_config_get_compress_level")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      subroutine catchem_config_get_output_directory(core_ptr, buffer, max_len) &
         bind(C, name="catchem_config_get_output_directory")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      subroutine catchem_config_get_output_prefix(core_ptr, buffer, max_len) &
         bind(C, name="catchem_config_get_output_prefix")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      integer(c_int) function catchem_config_get_latlon_output(core_ptr) &
         bind(C, name="catchem_config_get_latlon_output")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_config_get_diag_enabled(core_ptr) &
         bind(C, name="catchem_config_get_diag_enabled")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_config_get_process_diagnostics_enabled(core_ptr) &
         bind(C, name="catchem_config_get_process_diagnostics_enabled")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_config_get_diag_species_count(core_ptr) &
         bind(C, name="catchem_config_get_diag_species_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      subroutine catchem_config_get_diag_species_at(core_ptr, index, buffer, max_len) &
         bind(C, name="catchem_config_get_diag_species_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      integer(c_int) function catchem_config_get_output_attribute_count(core_ptr) &
         bind(C, name="catchem_config_get_output_attribute_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      subroutine catchem_config_get_output_attribute_key_at(core_ptr, index, buffer, max_len) &
         bind(C, name="catchem_config_get_output_attribute_key_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      subroutine catchem_config_get_output_attribute_value_at(core_ptr, index, buffer, max_len) &
         bind(C, name="catchem_config_get_output_attribute_value_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      subroutine catchem_config_get_config_file_path(core_ptr, buffer, max_len) &
         bind(C, name="catchem_config_get_config_file_path")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      subroutine catchem_get_build_version(buffer, max_len) &
         bind(C, name="catchem_get_build_version")
         import :: c_char, c_int
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      subroutine catchem_get_build_commit(buffer, max_len) &
         bind(C, name="catchem_get_build_commit")
         import :: c_char, c_int
         character(kind=c_char), intent(out) :: buffer(*)
         integer(c_int), value :: max_len
      end subroutine

      integer(c_int) function catchem_config_get_process_active(core_ptr, process_name) &
         bind(C, name="catchem_config_get_process_active")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: process_name(*)
      end function

      integer(c_int) function catchem_config_has_emission_mapping(core_ptr) &
         bind(C, name="catchem_config_has_emission_mapping")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      type(c_ptr) function catchem_diag_get_pointer(core_ptr, name) bind(C, name="catchem_diag_get_pointer")
         import :: c_ptr, c_char
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
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

      integer(c_int) function catchem_diag_get_rank(core_ptr, name) bind(C, name="catchem_diag_get_rank")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function

      integer(c_int) function catchem_diag_get_rank_checked(core_ptr, name, rank_out) &
         bind(C, name="catchem_diag_get_rank_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: rank_out
      end function

      subroutine catchem_diag_get_dims(core_ptr, name, dims_out) bind(C, name="catchem_diag_get_dims")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: dims_out(*)
      end subroutine

      integer(c_int) function catchem_diag_get_dims_checked(core_ptr, name, dims_out, dims_length) &
         bind(C, name="catchem_diag_get_dims_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: dims_out(*)
         integer(c_int), value :: dims_length
      end function

      subroutine catchem_diag_register(core_ptr, name, desc, units, rank, dim1, dim2, dim3) &
         bind(C, name="catchem_diag_register")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         character(kind=c_char), intent(in) :: desc(*)
         character(kind=c_char), intent(in) :: units(*)
         integer(c_int), value :: rank, dim1, dim2, dim3
      end subroutine

      integer(c_int) function catchem_diag_register_checked(core_ptr, name, desc, units, rank, dim1, dim2, dim3) &
         bind(C, name="catchem_diag_register_checked")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*), desc(*), units(*)
         integer(c_int), value :: rank, dim1, dim2, dim3
      end function

      integer(c_int) function catchem_diag_register_contract_checked(core_ptr, name, desc, units, rank, dims, axes, &
         policy, reset_value) &
         bind(C, name="catchem_diag_register_contract_checked")
         import :: c_ptr, c_char, c_int, c_double
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*), desc(*), units(*)
         integer(c_int), value :: rank, policy
         integer(c_int), intent(in) :: dims(*), axes(*)
         real(c_double), value :: reset_value
      end function

      integer(c_int) function catchem_diag_get_contract(core_ptr, name, generation, availability, latest_writer, &
         policy) bind(C, name="catchem_diag_get_contract")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: generation, availability, latest_writer, policy
      end function

      subroutine catchem_diag_sync_to_host(core_ptr) bind(C, name="catchem_diag_sync_to_host")
         import :: c_ptr
         type(c_ptr), value :: core_ptr
      end subroutine

      integer(c_int) function catchem_diag_get_count(core_ptr) bind(C, name="catchem_diag_get_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_diag_get_count_checked(core_ptr, count_out) &
         bind(C, name="catchem_diag_get_count_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), intent(out) :: count_out
      end function

      subroutine catchem_diag_get_name_at(core_ptr, index, name_out) bind(C, name="catchem_diag_get_name_at")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
      end subroutine

      integer(c_int) function catchem_diag_get_name_at_checked(core_ptr, index, name_out, name_length) &
         bind(C, name="catchem_diag_get_name_at_checked")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index, name_length
         character(kind=c_char), intent(out) :: name_out(*)
      end function

      integer(c_int) function catchem_state_get_species_count(state_ptr) bind(C, name="catchem_state_get_species_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      integer(c_int) function catchem_state_set_physical_validation_policy_checked(state_ptr, policy) &
         bind(C, name="catchem_state_set_physical_validation_policy_checked")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: policy
      end function

      integer(c_int) function catchem_state_get_physical_validation_report_checked( &
         state_ptr, issue_count, detail, detail_length) &
         bind(C, name="catchem_state_get_physical_validation_report_checked")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: state_ptr
         integer(c_int), intent(out) :: issue_count
         character(kind=c_char), intent(out) :: detail(*)
         integer(c_int), value :: detail_length
      end function

      real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function
   end interface

   !=========================================================================
   ! CATChem_Model Derived Type
   !=========================================================================
   type :: CATChem_Model
      type(c_ptr) :: cpp_core_ptr = c_null_ptr
      type(c_ptr) :: state_mgr_ptr = c_null_ptr
      character(len=64), allocatable, public :: required_fields(:)
      integer :: nx = 0
      integer :: ny = 0
      integer :: nz = 0
      logical :: initialized = .false.
      character(len=512), public :: last_error = ''
   contains
      procedure :: initialize => model_initialize
      procedure :: finalize => model_finalize
      procedure :: add_process => model_add_process
      procedure :: get_num_processes => model_get_num_processes
      procedure :: load_required_fields => model_load_required_fields
      procedure :: run_timestep => model_run_timestep
      procedure :: get_diagnostic_names => model_get_diagnostic_names
      procedure :: get_diagnostic => model_get_diagnostic
      procedure :: register_diagnostic => model_register_diagnostic
      procedure :: get_diagnostic_ptr => model_get_diagnostic_ptr
      procedure :: get_species_conc_ptr => model_get_species_conc_ptr
      procedure :: get_diag_index_from_field => model_get_diag_index_from_field
      procedure :: get_required_met_index => model_get_required_met_index
      procedure :: get_grid_dimensions => model_get_grid_dimensions
      procedure :: is_initialized => model_is_initialized
      procedure :: bind_met_3d => model_bind_met_3d
      procedure :: bind_met_3d_axis => model_bind_met_3d_axis
      procedure :: bind_met_2d => model_bind_met_2d
      procedure :: bind_unified_chemistry_3d => model_bind_unified_chemistry_3d
      procedure :: bind_unified_chemistry_4d => model_bind_unified_chemistry_4d
      generic :: bind_unified_chemistry => bind_unified_chemistry_3d, bind_unified_chemistry_4d
      procedure :: get_output_frequency => model_get_output_frequency
      procedure :: get_compress_level => model_get_compress_level
      procedure :: get_output_directory => model_get_output_directory
      procedure :: get_output_prefix => model_get_output_prefix
      procedure :: is_latlon_output_enabled => model_is_latlon_output_enabled
      procedure :: is_diag_enabled => model_is_diag_enabled
      procedure :: is_process_diag_enabled => model_is_process_diag_enabled
      procedure :: get_diag_species_count => model_get_diag_species_count
      procedure :: get_diag_species_at => model_get_diag_species_at
      procedure :: get_output_attribute_count => model_get_output_attribute_count
      procedure :: get_output_attribute_at => model_get_output_attribute_at
      procedure :: get_config_file_path => model_get_config_file_path
      procedure :: get_build_version => model_get_build_version
      procedure :: get_build_commit => model_get_build_commit
      procedure :: is_process_active => model_is_process_active
      procedure :: has_emission_mapping => model_has_emission_mapping
      procedure :: set_physical_validation_policy => model_set_physical_validation_policy
      procedure :: get_physical_validation_report => model_get_physical_validation_report
   end type CATChem_Model

contains

   subroutine capture_boundary_error(this)
      class(CATChem_Model), intent(inout) :: this
      character(kind=c_char) :: buffer(512)
      integer :: i, ignored_status
      this%last_error = ''
      ignored_status = catchem_get_last_error(buffer, int(size(buffer), c_int))
      do i = 1, min(len(this%last_error), size(buffer))
         if (buffer(i) == c_null_char) exit
         this%last_error(i:i) = buffer(i)
      end do
   end subroutine capture_boundary_error

   ! Helper to convert standard Fortran string to null-terminated C char array
   subroutine to_c_string(f_str, c_arr)
      character(len=*), intent(in) :: f_str
      character(kind=c_char), intent(out) :: c_arr(*)
      integer :: i, f_len

      f_len = len_trim(f_str)
      ! Strip any trailing c_null_char if present in f_str
      if (f_len > 0) then
         if (f_str(f_len:f_len) == c_null_char) f_len = f_len - 1
      end if
      do i = 1, f_len
         c_arr(i) = f_str(i:i)
      end do
      c_arr(f_len+1) = c_null_char
   end subroutine to_c_string

   ! Initialize CATChem model and load configuration
   subroutine model_initialize(this, config_file, nx, ny, nz, nsoil, nsoiltype, nsurftype, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: config_file
      integer, intent(in) :: nx, ny, nz
      integer, intent(in), optional :: nsoil, nsoiltype, nsurftype
      integer, intent(out) :: rc

      character(kind=c_char) :: c_filename(512)
      integer(c_int) :: c_status

      call to_c_string(config_file, c_filename)

      call catchem_register_carbchem_cpp()
      call catchem_register_drydep_cpp()
      call catchem_register_dust_cpp()
      call catchem_register_seasalt_cpp()
      call catchem_register_settling_cpp()
      call catchem_register_so4chem_cpp()
      call catchem_register_wetdep_cpp()

      ! Retain the C++ configuration error if construction fails.  This is a
      ! plain C pointer output, so no ESMF descriptor crosses the ABI.
      c_status = catchem_core_create_from_config_with_grid_checked( &
         c_filename, int(nx*ny, c_int), int(nz, c_int), this%cpp_core_ptr)
      if (c_status /= 0_c_int) then
         call capture_boundary_error(this)
         rc = CC_FAILURE
         return
      end if
      if (.not. c_associated(this%cpp_core_ptr)) then
         rc = CC_FAILURE
         this%last_error = 'core_create_from_config_with_grid returned a null handle without a boundary error'
         return
      end if

      this%state_mgr_ptr = catchem_core_get_state_manager(this%cpp_core_ptr)
      if (.not. c_associated(this%state_mgr_ptr)) then
         call catchem_core_destroy(this%cpp_core_ptr)
         this%cpp_core_ptr = c_null_ptr
         rc = CC_FAILURE
         this%last_error = 'core_get_state_manager returned a null handle'
         return
      end if

      ! Host-local dimensions are authoritative
      this%nx = nx
      this%ny = ny
      this%nz = nz

      call this%load_required_fields(rc)
      if (rc /= CC_SUCCESS) then
         call catchem_core_destroy(this%cpp_core_ptr)
         this%cpp_core_ptr = c_null_ptr
         this%state_mgr_ptr = c_null_ptr
         return
      end if

      this%initialized = .true.
      rc = CC_SUCCESS
   end subroutine model_initialize

   ! Finalize model and release memory
   subroutine model_finalize(this, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: rc

      if (c_associated(this%cpp_core_ptr)) then
         call catchem_core_destroy(this%cpp_core_ptr)
         rc = CC_SUCCESS
         this%cpp_core_ptr = c_null_ptr
         this%state_mgr_ptr = c_null_ptr
      else
         rc = CC_SUCCESS
      end if
      this%initialized = .false.
      if (rc == CC_SUCCESS) this%last_error = ''
   end subroutine model_finalize

   ! Validate process registration state.
   !! NOTE: This routine does not add processes.  All science processes are
   !! registered with the C++ ProcessRegistry during model_initialize and are
   !! instantiated from the runtime YAML (processes/<name>/activate) while the
   !! Core is constructed.  This entry point is retained for API compatibility
   !! and only verifies that the underlying Core handle is present; process
   !! selection is configuration-driven, not call-driven.
   subroutine model_add_process(this, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: rc

      if (c_associated(this%cpp_core_ptr)) then
         rc = CC_SUCCESS
      else
         rc = CC_FAILURE
         this%last_error = 'model_add_process called before model_initialize: no Core handle'
      end if
   end subroutine model_add_process

   ! Get number of active processes
   function model_get_num_processes(this) result(num_processes)
      class(CATChem_Model), intent(inout) :: this
      integer :: num_processes
      num_processes = int(catchem_core_get_num_processes(this%cpp_core_ptr))
   end function model_get_num_processes

   ! Load active process-contract requirements for host-owned fields.
   subroutine model_load_required_fields(this, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: rc

      character(kind=c_char) :: c_name(64)
      integer(c_int) :: c_count, c_status
      integer :: i, j

      if (allocated(this%required_fields)) deallocate(this%required_fields)
      c_count = 0_c_int
      c_status = catchem_core_get_required_host_field_count_checked(this%cpp_core_ptr, c_count)
      if (c_status /= 0_c_int) then
         call capture_boundary_error(this)
         rc = CC_FAILURE
         return
      end if

      allocate(this%required_fields(int(c_count)))
      this%required_fields = ''
      do i = 1, int(c_count)
         c_name = c_null_char
         c_status = catchem_core_get_required_host_field_name_checked(this%cpp_core_ptr, int(i - 1, c_int), &
            c_name, int(size(c_name), c_int))
         if (c_status /= 0_c_int) then
            call capture_boundary_error(this)
            deallocate(this%required_fields)
            rc = CC_FAILURE
            return
         end if
         do j = 1, size(c_name)
            if (c_name(j) == c_null_char) exit
            this%required_fields(i)(j:j) = c_name(j)
         end do
      end do
      rc = CC_SUCCESS
   end subroutine model_load_required_fields

   ! Execute standard timestep
   subroutine model_run_timestep(this, timestep, dt, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(in) :: timestep
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      if (c_associated(this%cpp_core_ptr)) then
         call catchem_state_sync_to_device(this%state_mgr_ptr)
         rc = catchem_core_run_timestep(this%cpp_core_ptr, real(dt, c_double))
         if (rc == CC_SUCCESS) then
            call catchem_state_sync_to_host(this%state_mgr_ptr)
         end if
      else
         rc = CC_FAILURE
         this%last_error = 'run_timestep: model is not initialized'
      end if
   end subroutine model_run_timestep

   ! Get available diagnostic field names
   subroutine model_get_diagnostic_names(this, diagnostic_names, diagnostic_fields, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), allocatable, intent(out) :: diagnostic_names(:)
      character(len=*), allocatable, optional, intent(out) :: diagnostic_fields(:)
      integer, intent(out) :: rc

      integer :: i, count
      integer(c_int) :: c_count, status
      character(kind=c_char) :: c_name(64)
      character(len=64) :: f_name

      status = catchem_diag_get_count_checked(this%cpp_core_ptr, c_count)
      if (status /= 0_c_int) then
         allocate(diagnostic_names(0))
         if (present(diagnostic_fields)) allocate(diagnostic_fields(0))
         rc = int(status)
         call capture_boundary_error(this)
         return
      end if
      count = int(c_count)
      allocate(diagnostic_names(count))
      if (present(diagnostic_fields)) allocate(diagnostic_fields(count))

      do i = 1, count
         status = catchem_diag_get_name_at_checked(this%cpp_core_ptr, int(i - 1, c_int), c_name, 64_c_int)
         if (status /= 0_c_int) then
            rc = int(status)
            call capture_boundary_error(this)
            return
         end if
         f_name = ""
         block
            integer :: j
            do j = 1, 64
               if (c_name(j) == c_null_char) exit
               f_name(j:j) = c_name(j)
            end do
         end block
         diagnostic_names(i) = trim(f_name)
         if (present(diagnostic_fields)) diagnostic_fields(i) = trim(f_name)
      end do

      rc = CC_SUCCESS
   end subroutine model_get_diagnostic_names

   ! Retrieve individual diagnostic data mapped directly from C++ heap
   subroutine model_get_diagnostic(this, diagnostic_name, diagnostic_data, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: diagnostic_name
      real(fp), allocatable, intent(out) :: diagnostic_data(:,:,:)
      integer, intent(out) :: rc

      character(kind=c_char) :: c_name(64)
      type(c_ptr) :: raw_ptr
      integer(c_int) :: rank, dims(3)
      real(c_double), pointer :: f_ptr_2d(:,:) => null()
      real(c_double), pointer :: f_ptr_3d(:,:,:) => null()

      call to_c_string(diagnostic_name, c_name)
      raw_ptr = c_null_ptr
      rc = int(catchem_diag_get_rank_checked(this%cpp_core_ptr, c_name, rank))
      if (rc /= CC_SUCCESS) then
         call capture_boundary_error(this)
         return
      end if
      dims = 0
      rc = int(catchem_diag_get_dims_checked(this%cpp_core_ptr, c_name, dims, 3_c_int))
      if (rc /= CC_SUCCESS) then
         call capture_boundary_error(this)
         return
      end if
      rc = int(catchem_diag_get_pointer_checked(this%cpp_core_ptr, c_name, rank, dims, raw_ptr))
      if (rc /= CC_SUCCESS .or. .not. c_associated(raw_ptr)) then
         call capture_boundary_error(this)
         return
      end if

      allocate(diagnostic_data(this%nx, this%ny, this%nz))
      diagnostic_data = 0.0_fp

      if (rank == 2) then
         call c_f_pointer(raw_ptr, f_ptr_2d, [dims(1), dims(2)])
         if (dims(2) == this%nz .and. dims(1) == this%nx * this%ny) then
            diagnostic_data = real(reshape(f_ptr_2d, [this%nx, this%ny, this%nz]), fp)
         else if (dims(1) == this%nx * this%ny) then
            diagnostic_data(:,:,1) = real(reshape(f_ptr_2d(:,1), [this%nx, this%ny]), fp)
         end if
      else if (rank == 3) then
         call c_f_pointer(raw_ptr, f_ptr_3d, [dims(1), dims(2), dims(3)])
         if (dims(1) == this%nx .and. dims(2) == this%ny .and. dims(3) == this%nz) then
            diagnostic_data = real(f_ptr_3d, fp)
         else
            deallocate(diagnostic_data)
            rc = CC_FAILURE
            return
         end if
      end if

      rc = CC_SUCCESS
   end subroutine model_get_diagnostic

   ! Find index of registered diagnostic name
   function model_get_diag_index_from_field(this, field_name) result(found_index)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer :: found_index

      character(len=128), allocatable :: diagnostic_names(:)
      integer :: rc, i

      found_index = 0
      call this%get_diagnostic_names(diagnostic_names, rc=rc)
      if (allocated(diagnostic_names)) then
         do i = 1, size(diagnostic_names)
            if (trim(field_name) == trim(diagnostic_names(i))) then
               found_index = i
               exit
            end if
         end do
      end if
   end function model_get_diag_index_from_field

   ! Required met field utility
   function model_get_required_met_index(this, var_name) result(found_index)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer :: found_index, i

      found_index = 0
      if (.not. allocated(this%required_fields)) return
      do i = 1, size(this%required_fields)
         if (trim(var_name) == trim(this%required_fields(i))) then
            found_index = i
            return
         end if
      end do
   end function model_get_required_met_index

   ! Grid dimension query
   subroutine model_get_grid_dimensions(this, nx, ny, nz)
      class(CATChem_Model), intent(in) :: this
      integer, intent(out) :: nx, ny, nz

      nx = this%nx
      ny = this%ny
      nz = this%nz
   end subroutine model_get_grid_dimensions

   ! Status checking
   function model_is_initialized(this) result(is_initialized)
      class(CATChem_Model), intent(in) :: this
      logical :: is_initialized

      is_initialized = this%initialized
   end function model_is_initialized

   subroutine model_set_physical_validation_policy(this, policy, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(in) :: policy
      integer, intent(out) :: rc

      rc = int(catchem_state_set_physical_validation_policy_checked( &
         this%state_mgr_ptr, int(policy, c_int)))
      if (rc /= CC_SUCCESS) call capture_boundary_error(this)
   end subroutine model_set_physical_validation_policy

   subroutine model_get_physical_validation_report(this, issue_count, detail, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(out) :: issue_count
      character(len=*), intent(out) :: detail
      integer, intent(out) :: rc
      character(kind=c_char) :: c_detail(1024)
      integer(c_int) :: c_count, status
      integer :: i

      issue_count = 0
      detail = ''
      status = catchem_state_get_physical_validation_report_checked( &
         this%state_mgr_ptr, c_count, c_detail, int(size(c_detail), c_int))
      rc = int(status)
      issue_count = int(c_count)
      if (status /= 0_c_int) then
         call capture_boundary_error(this)
         return
      end if
      do i = 1, min(len(detail), size(c_detail))
         if (c_detail(i) == c_null_char) exit
         detail(i:i) = c_detail(i)
      end do
   end subroutine model_get_physical_validation_report

   ! Bind a 3D meteorological field
   subroutine model_bind_met_3d(this, name, arr, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(c_double), target, contiguous, intent(in) :: arr(:,:,:)
      integer, optional, intent(out) :: rc

      character(kind=c_char) :: c_name(64)
      integer(c_int) :: status

      call to_c_string(name, c_name)
      ! Route through the checked binder so an array whose extents do not
      ! match the state contract is rejected with the exact boundary status
      ! (e.g. extent mismatch) instead of being silently accepted.  CATChem
      ! flattens the two horizontal Fortran dimensions into one column axis.
      status = catchem_state_bind_met_3d_checked(this%state_mgr_ptr, c_name, c_loc(arr(1,1,1)), &
         int(size(arr, 1) * size(arr, 2), c_int), int(size(arr, 3), c_int), 1_c_int)
      if (present(rc)) then
         rc = int(status)
         if (status /= 0_c_int) call capture_boundary_error(this)
      end if
   end subroutine model_bind_met_3d

   ! Bind a 3D meteorological field with an explicit vertical semantic axis.
   ! CATChem stores horizontal Fortran dimensions as one flattened column axis.
   subroutine model_bind_met_3d_axis(this, name, arr, semantic_axis, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(c_double), target, contiguous, intent(in) :: arr(:,:,:)
      integer, intent(in) :: semantic_axis
      integer, optional, intent(out) :: rc

      character(kind=c_char) :: c_name(64)
      integer(c_int) :: status

      call to_c_string(name, c_name)
      status = catchem_state_bind_met_3d_axis_checked(this%state_mgr_ptr, c_name, c_loc(arr(1,1,1)), &
         int(size(arr, 1) * size(arr, 2), c_int), int(size(arr, 3), c_int), 1_c_int, int(semantic_axis, c_int))
      if (present(rc)) then
         if (status == 0_c_int) then
            rc = CC_SUCCESS
         else
            rc = CC_FAILURE
            call capture_boundary_error(this)
         end if
      end if
   end subroutine model_bind_met_3d_axis

   ! Bind a 2D meteorological field
   subroutine model_bind_met_2d(this, name, arr, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(c_double), target, contiguous, intent(in) :: arr(:,:)
      integer, optional, intent(out) :: rc

      character(kind=c_char) :: c_name(64)

      call to_c_string(name, c_name)
      call catchem_state_bind_met_2d(this%state_mgr_ptr, c_name, c_loc(arr(1,1)))
      if (present(rc)) rc = CC_SUCCESS
   end subroutine model_bind_met_2d

   ! Bind unified chemical concentrations 3D array
   subroutine model_bind_unified_chemistry_3d(this, arr, rc)
      class(CATChem_Model), intent(inout) :: this
      real(c_double), target, contiguous, intent(in) :: arr(:,:,:)
      integer, optional, intent(out) :: rc

      call catchem_state_bind_unified_chemistry(this%state_mgr_ptr, c_loc(arr(1,1,1)))
      if (present(rc)) rc = CC_SUCCESS
   end subroutine model_bind_unified_chemistry_3d

   ! Bind unified chemical concentrations 4D array
   subroutine model_bind_unified_chemistry_4d(this, arr, rc)
      class(CATChem_Model), intent(inout) :: this
      real(c_double), target, contiguous, intent(in) :: arr(:,:,:,:)
      integer, optional, intent(out) :: rc

      call catchem_state_bind_unified_chemistry(this%state_mgr_ptr, c_loc(arr(1,1,1,1)))
      if (present(rc)) rc = CC_SUCCESS
   end subroutine model_bind_unified_chemistry_4d

   ! Register a 3D diagnostic field in the C++ DiagnosticManager.
   !
   ! By default the field is DiagnosticPolicy::Instantaneous, which the C++
   ! manager resets (full-array memset) at the start of every timestep.  Hosts
   ! that fully overwrite the field each step (e.g. the PM2.5/PM10 diagnostics)
   ! can pass persistent=.true. to select DiagnosticPolicy::Persistent and skip
   ! that per-step reset; the writer must then guarantee every element is
   ! rewritten before the field is read.
   subroutine model_register_diagnostic(this, name, desc, units, dims, rc, persistent)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: name, desc, units
      integer, intent(in) :: dims(3)
      integer, intent(out) :: rc
      logical, optional, intent(in) :: persistent

      character(kind=c_char) :: c_name(64), c_desc(128), c_units(64)
      integer(c_int) :: c_dims(3), c_axes(3)
      integer, parameter :: POLICY_INSTANTANEOUS = 0, POLICY_PERSISTENT = 2
      integer, parameter :: AXIS_COLUMN = 0, AXIS_LEVEL = 1, AXIS_SPECIES = 4

      rc = CC_FAILURE
      if (.not. c_associated(this%cpp_core_ptr)) return
      call to_c_string(name, c_name)
      call to_c_string(desc, c_desc)
      call to_c_string(units, c_units)
      if (present(persistent) .and. persistent) then
         ! Match the axes DiagnosticManager::register_field() would derive for a
         ! 3D field (Column, Level, Species) so re-registration stays consistent.
         c_dims = int(dims, c_int)
         c_axes = [AXIS_COLUMN, AXIS_LEVEL, AXIS_SPECIES]
         rc = int(catchem_diag_register_contract_checked(this%cpp_core_ptr, c_name, c_desc, c_units, &
            3_c_int, c_dims, c_axes, int(POLICY_PERSISTENT, c_int), 0.0_c_double))
      else
         rc = int(catchem_diag_register_checked(this%cpp_core_ptr, c_name, c_desc, c_units, &
            3_c_int, int(dims(1), c_int), int(dims(2), c_int), int(dims(3), c_int)))
      end if
      if (rc /= CC_SUCCESS) call capture_boundary_error(this)
   end subroutine model_register_diagnostic

   subroutine model_get_species_conc_ptr(this, species_index, ptr3d, dims, rc)
      class(CATChem_Model), intent(inout) :: this
      integer, intent(in) :: species_index
      real(fp), pointer, intent(out) :: ptr3d(:,:,:)
      integer, intent(in) :: dims(3)
      integer, intent(out) :: rc
      type(c_ptr) :: raw_ptr
      integer(c_int) :: status

      rc = CC_FAILURE
      nullify(ptr3d)
      if (.not. c_associated(this%state_mgr_ptr)) return
      status = catchem_state_get_species_conc_pointer_checked(this%state_mgr_ptr, int(species_index, c_int), &
         int(dims(1) * dims(2), c_int), int(dims(3), c_int), raw_ptr)
      if (status /= 0_c_int .or. .not. c_associated(raw_ptr)) then
         rc = int(status)
         call capture_boundary_error(this)
         return
      end if
      call c_f_pointer(raw_ptr, ptr3d, dims)
      rc = CC_SUCCESS
   end subroutine model_get_species_conc_ptr

   ! Map the C++-owned storage of a registered 3D diagnostic for in-place
   ! writes (zero-copy; the same memory NUOPC export and NetCDF output read)
   subroutine model_get_diagnostic_ptr(this, name, ptr3d, dims, rc)
      class(CATChem_Model), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(fp), pointer, intent(out) :: ptr3d(:, :, :)
      integer, intent(in) :: dims(3)
      integer, intent(out) :: rc

      character(kind=c_char) :: c_name(64)
      type(c_ptr) :: raw_ptr
      integer(c_int) :: c_dims(3), status

      rc = CC_FAILURE
      nullify(ptr3d)
      if (.not. c_associated(this%cpp_core_ptr)) return
      call to_c_string(name, c_name)
      c_dims = int(dims, c_int)
      status = catchem_diag_get_pointer_checked(this%cpp_core_ptr, c_name, 3_c_int, c_dims, raw_ptr)
      if (status /= 0_c_int .or. .not. c_associated(raw_ptr)) then
         rc = int(status)
         call capture_boundary_error(this)
         return
      end if
      call c_f_pointer(raw_ptr, ptr3d, dims)
      rc = CC_SUCCESS
   end subroutine model_get_diagnostic_ptr

   function model_get_output_frequency(this) result(freq)
      class(CATChem_Model), intent(in) :: this
      integer :: freq
      freq = int(catchem_config_get_output_frequency(this%cpp_core_ptr))
   end function model_get_output_frequency

   function model_get_compress_level(this) result(clev)
      class(CATChem_Model), intent(in) :: this
      integer :: clev
      clev = int(catchem_config_get_compress_level(this%cpp_core_ptr))
   end function model_get_compress_level

   subroutine model_get_output_directory(this, dir_out)
      class(CATChem_Model), intent(in) :: this
      character(len=*), intent(out) :: dir_out
      character(kind=c_char) :: c_buf(256)
      integer :: i
      call catchem_config_get_output_directory(this%cpp_core_ptr, c_buf, 256_c_int)
      dir_out = ""
      do i = 1, 256
         if (c_buf(i) == c_null_char) exit
         dir_out(i:i) = c_buf(i)
      end do
   end subroutine model_get_output_directory

   subroutine model_get_output_prefix(this, prefix_out)
      class(CATChem_Model), intent(in) :: this
      character(len=*), intent(out) :: prefix_out
      character(kind=c_char) :: c_buf(256)
      integer :: i
      call catchem_config_get_output_prefix(this%cpp_core_ptr, c_buf, 256_c_int)
      prefix_out = ""
      do i = 1, 256
         if (c_buf(i) == c_null_char) exit
         prefix_out(i:i) = c_buf(i)
      end do
   end subroutine model_get_output_prefix

   function model_is_latlon_output_enabled(this) result(enabled)
      class(CATChem_Model), intent(in) :: this
      logical :: enabled
      enabled = (catchem_config_get_latlon_output(this%cpp_core_ptr) /= 0_c_int)
   end function model_is_latlon_output_enabled

   function model_is_diag_enabled(this) result(enabled)
      class(CATChem_Model), intent(in) :: this
      logical :: enabled
      enabled = (catchem_config_get_diag_enabled(this%cpp_core_ptr) /= 0_c_int)
   end function model_is_diag_enabled

   !> \brief True when diagnostics.output/process_diagnostics is enabled in the
   !! runtime YAML.  DEPRECATED: the NUOPC driver no longer consults this key
   !! -- per-process diagnostics are written whenever they are registered and
   !! runtime diagnostics are enabled (parity with the legacy Fortran core).
   !! The accessor is retained only so existing configs keep parsing.
   function model_is_process_diag_enabled(this) result(enabled)
      class(CATChem_Model), intent(in) :: this
      logical :: enabled
      enabled = (catchem_config_get_process_diagnostics_enabled(this%cpp_core_ptr) /= 0_c_int)
   end function model_is_process_diag_enabled

   function model_get_diag_species_count(this) result(count)
      class(CATChem_Model), intent(in) :: this
      integer :: count
      count = int(catchem_config_get_diag_species_count(this%cpp_core_ptr))
   end function model_get_diag_species_count

   subroutine model_get_diag_species_at(this, index, species_name)
      class(CATChem_Model), intent(in) :: this
      integer, intent(in) :: index
      character(len=*), intent(out) :: species_name
      character(kind=c_char) :: c_buf(128)
      integer :: i
      call catchem_config_get_diag_species_at(this%cpp_core_ptr, int(index - 1, c_int), c_buf, 128_c_int)
      species_name = ""
      do i = 1, 128
         if (c_buf(i) == c_null_char) exit
         species_name(i:i) = c_buf(i)
      end do
   end subroutine model_get_diag_species_at

   function model_get_output_attribute_count(this) result(count)
      class(CATChem_Model), intent(in) :: this
      integer :: count
      count = int(catchem_config_get_output_attribute_count(this%cpp_core_ptr))
   end function model_get_output_attribute_count

   subroutine model_get_output_attribute_at(this, index, attr_name, attr_value)
      class(CATChem_Model), intent(in) :: this
      integer, intent(in) :: index
      character(len=*), intent(out) :: attr_name, attr_value
      character(kind=c_char) :: c_buf(256)
      integer :: i

      call catchem_config_get_output_attribute_key_at(this%cpp_core_ptr, int(index - 1, c_int), c_buf, 256_c_int)
      attr_name = ""
      do i = 1, 256
         if (c_buf(i) == c_null_char) exit
         attr_name(i:i) = c_buf(i)
      end do

      call catchem_config_get_output_attribute_value_at(this%cpp_core_ptr, int(index - 1, c_int), c_buf, 256_c_int)
      attr_value = ""
      do i = 1, 256
         if (c_buf(i) == c_null_char) exit
         attr_value(i:i) = c_buf(i)
      end do
   end subroutine model_get_output_attribute_at

   subroutine model_get_config_file_path(this, path_out)
      class(CATChem_Model), intent(in) :: this
      character(len=*), intent(out) :: path_out
      character(kind=c_char) :: c_buf(512)
      integer :: i
      call catchem_config_get_config_file_path(this%cpp_core_ptr, c_buf, 512_c_int)
      path_out = ""
      do i = 1, 512
         if (c_buf(i) == c_null_char) exit
         path_out(i:i) = c_buf(i)
      end do
   end subroutine model_get_config_file_path

   subroutine model_get_build_version(this, version_out)
      class(CATChem_Model), intent(in) :: this
      character(len=*), intent(out) :: version_out
      character(kind=c_char) :: c_buf(64)
      integer :: i
      call catchem_get_build_version(c_buf, 64_c_int)
      version_out = ""
      do i = 1, 64
         if (c_buf(i) == c_null_char) exit
         version_out(i:i) = c_buf(i)
      end do
   end subroutine model_get_build_version

   subroutine model_get_build_commit(this, commit_out)
      class(CATChem_Model), intent(in) :: this
      character(len=*), intent(out) :: commit_out
      character(kind=c_char) :: c_buf(64)
      integer :: i
      call catchem_get_build_commit(c_buf, 64_c_int)
      commit_out = ""
      do i = 1, 64
         if (c_buf(i) == c_null_char) exit
         commit_out(i:i) = c_buf(i)
      end do
   end subroutine model_get_build_commit

   function model_is_process_active(this, process_name) result(active)
      class(CATChem_Model), intent(in) :: this
      character(len=*), intent(in) :: process_name
      logical :: active
      character(kind=c_char) :: c_name(128)
      call to_c_string(process_name, c_name)
      active = (catchem_config_get_process_active(this%cpp_core_ptr, c_name) /= 0_c_int)
   end function model_is_process_active

   function model_has_emission_mapping(this) result(has_mapping)
      class(CATChem_Model), intent(in) :: this
      logical :: has_mapping
      has_mapping = (catchem_config_has_emission_mapping(this%cpp_core_ptr) /= 0_c_int)
   end function model_has_emission_mapping

end module CATChem_API
