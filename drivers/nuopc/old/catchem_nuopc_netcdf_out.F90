! \file catchem_nuopc_netcdf_out.F90
! \brief NetCDF diagnostic output module for CATChem NUOPC cap
!>
! \details
! This module provides parallel NetCDF diagnostic output capabilities for the
! CATChem NUOPC cap. It leverages ESMF I/O capabilities for efficient parallel
! writing and ensures CF convention compliance. Output configuration is
! controlled via YAML configuration file.
!>
! Features:
! - ESMF-based parallel NetCDF I/O
! - CF convention compliance
! - YAML-based configuration
! - Time averaging and aggregation
! - Flexible field selection
! - Proper coordinate variable handling
!>
! \author Barry Baker
! \date 11/2024
! \ingroup catchem_nuopc_group

module catchem_nuopc_netcdf_out

   use ESMF
   use NUOPC
   use CATChem
   use catchem_types, only: catchem_container_type

   implicit none

   private

   public :: output_diagnostics_init
   public :: output_diagnostics_write
   public :: output_diagnostics_finalize
   public :: output_config_type
   public :: diagnostic_field_type

   ! Diagnostic field configuration type
   type :: diagnostic_field_type
      character(len=128) :: variable_name
      character(len=128) :: standard_name
      character(len=256) :: long_name
      character(len=128) :: catchem_var
      integer :: dimensions
      character(len=64) :: units
      character(len=128) :: cell_methods
      character(len=256) :: coordinates
      logical :: optional = .false.
      logical :: enabled = .true.
      type(ESMF_Field) :: esmf_field
   end type diagnostic_field_type

   ! Output configuration type
   type :: output_config_type
      ! Timing and frequency
      character(len=64) :: output_frequency
      character(len=128) :: output_start_time

      ! File naming
      character(len=256) :: filename_template
      character(len=256) :: output_directory
      logical :: split_by_time
      integer :: max_records_per_file

      ! NetCDF format
      character(len=32) :: netcdf_format
      integer :: compression_level
      logical :: shuffle

      ! CF conventions
      character(len=32) :: cf_conventions
      character(len=128) :: institution
      character(len=256) :: source
      character(len=512) :: references
      character(len=512) :: comment

      ! Fields
      integer :: n_fields
      type(diagnostic_field_type), allocatable :: fields(:)

      ! Aggregation
      logical :: time_averaging
      character(len=64) :: averaging_period
      logical :: output_instantaneous

      ! Parallel I/O
      logical :: use_esmf_io
      logical :: collective_io
      integer :: io_root_rank
      logical :: output_on_native_grid
   end type output_config_type

   ! Module variables
   type(output_config_type), save :: output_config
   type(QFYAML_t), save :: config_yaml
   type(ESMF_FieldBundle), save :: output_bundle
   type(ESMF_Grid), save :: output_grid
   type(ESMF_Time), save :: last_output_time
   type(ESMF_TimeInterval), save :: output_interval
   character(len=256), save :: current_filename
   logical, save :: output_initialized = .false.

contains

   ! Initialize diagnostic output system
   !!
   !! \param  grid         ESMF grid for output
   !! \param  config_file  Output configuration file
   !! \param  start_time   Simulation start time
   !! \param rc           ESMF return code
   !!
   subroutine output_diagnostics_init(grid, config_file, start_time, rc)

      type(ESMF_Grid), intent(in) :: grid
      character(len=*), intent(in) :: config_file
      type(ESMF_Time), intent(in) :: start_time
      integer, intent(out) :: rc

      character(len=*), parameter :: routine = 'output_diagnostics_init'
      integer :: i, yaml_rc
      character(len=512) :: errmsg

      rc = ESMF_SUCCESS

      ! Load configuration from YAML file
      call cc_yaml_init(config_yaml, trim(config_file), yaml_rc)
      if (yaml_rc /= QFYAML_Success) then
         call ESMF_LogWrite("Failed to load output config: " // trim(config_file), &
            ESMF_LOGMSG_ERROR, rc=rc)
         rc = ESMF_FAILURE
         return
      end if

      ! Parse configuration
      call parse_output_config(rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Set up output grid
      output_grid = grid

      ! Create output directory if needed
      call create_output_directory(rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Set up output timing
      call setup_output_timing(start_time, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Create output fields and bundle
      call create_output_fields(rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      output_initialized = .true.

      call ESMF_LogWrite("CATChem: Initialized diagnostic output", &
         ESMF_LOGMSG_INFO, rc=rc)

   end subroutine output_diagnostics_init

   ! Write diagnostic output
   !!
   !! \param  catchem_states  CATChem state container
   !! \param  current_time    Current simulation time
   !! \param rc              ESMF return code
   !!
   subroutine output_diagnostics_write(catchem_states, current_time, rc)

      type(catchem_container_type), intent(in) :: catchem_states
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      logical :: time_to_write
      character(len=*), parameter :: routine = 'output_diagnostics_write'

      rc = ESMF_SUCCESS

      if (.not. output_initialized) then
         call ESMF_LogWrite("Output system not initialized", ESMF_LOGMSG_WARNING, rc=rc)
         return
      end if

      ! Check if it's time to write output
      call check_output_time(current_time, time_to_write, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (.not. time_to_write) return

      ! Copy data to output fields
      call populate_output_fields(catchem_states, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Generate filename for this time
      call generate_filename(current_time, current_filename, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Write the data using ESMF I/O
      call write_netcdf_file(current_time, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Update last output time
      last_output_time = current_time

      call ESMF_LogWrite("CATChem: Wrote diagnostic output to " // trim(current_filename), &
         ESMF_LOGMSG_INFO, rc=rc)

   end subroutine output_diagnostics_write

   ! Finalize diagnostic output system
   !!
   !! \param rc  ESMF return code
   !!
   subroutine output_diagnostics_finalize(rc)

      integer, intent(out) :: rc
      integer :: i

      rc = ESMF_SUCCESS

      if (.not. output_initialized) return

      ! Clean up ESMF objects
      if (ESMF_FieldBundleIsCreated(output_bundle)) then
         call ESMF_FieldBundleDestroy(output_bundle, rc=rc)
      end if

      ! Clean up individual fields
      if (allocated(output_config%fields)) then
         do i = 1, output_config%n_fields
            if (ESMF_FieldIsCreated(output_config%fields(i)%esmf_field)) then
               call ESMF_FieldDestroy(output_config%fields(i)%esmf_field, rc=rc)
            end if
         end do
         deallocate(output_config%fields)
      end if

      ! Clean up YAML parser
      call cc_yaml_cleanup(config_yaml)

      output_initialized = .false.

      call ESMF_LogWrite("CATChem: Finalized diagnostic output", &
         ESMF_LOGMSG_INFO, rc=rc)

   end subroutine output_diagnostics_finalize

   ! Parse output configuration from YAML
   !!
   !! \param rc  Return code
   !!
   subroutine parse_output_config(rc)

      integer, intent(out) :: rc
      integer :: i, yaml_rc
      character(len=128) :: field_key, temp_str

      rc = CC_SUCCESS

      ! Parse output settings
      call cc_yaml_get(config_yaml, 'output_settings%output_frequency', &
         output_config%output_frequency, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%output_frequency = "1h"

      call cc_yaml_get(config_yaml, 'output_settings%filename_template', &
         output_config%filename_template, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%filename_template = "CATChem_diag_%Y%m%d_%H%M%S.nc"

      call cc_yaml_get(config_yaml, 'output_settings%output_directory', &
         output_config%output_directory, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%output_directory = "./output"

      call cc_yaml_get(config_yaml, 'output_settings%netcdf_format', &
         output_config%netcdf_format, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%netcdf_format = "NETCDF4_CLASSIC"

      call cc_yaml_get(config_yaml, 'output_settings%compression_level', &
         output_config%compression_level, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%compression_level = 4

      call cc_yaml_get(config_yaml, 'output_settings%cf_conventions', &
         output_config%cf_conventions, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%cf_conventions = "CF-1.8"

      call cc_yaml_get(config_yaml, 'output_settings%institution', &
         output_config%institution, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%institution = "NOAA"

      call cc_yaml_get(config_yaml, 'output_settings%source', &
         output_config%source, yaml_rc)
      if (yaml_rc /= QFYAML_Success) output_config%source = "CATChem atmospheric chemistry model"

      ! Parse parallel I/O settings
      call cc_yaml_get(config_yaml, 'parallel_io%use_esmf_io', temp_str, yaml_rc)
      if (yaml_rc == QFYAML_Success) then
         output_config%use_esmf_io = (trim(temp_str) == 'true')
      else
         output_config%use_esmf_io = .true.
      end if

      call cc_yaml_get(config_yaml, 'parallel_io%collective_io', temp_str, yaml_rc)
      if (yaml_rc == QFYAML_Success) then
         output_config%collective_io = (trim(temp_str) == 'true')
      else
         output_config%collective_io = .true.
      end if

      ! Parse diagnostic fields
      call cc_yaml_get(config_yaml, 'diagnostic_fields%n_fields', &
         output_config%n_fields, yaml_rc)
      if (yaml_rc /= QFYAML_Success) then
         rc = CC_FAILURE
         return
      end if

      allocate(output_config%fields(output_config%n_fields))

      do i = 1, output_config%n_fields
         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%variable_name'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%variable_name, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%standard_name'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%standard_name, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%long_name'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%long_name, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%catchem_var'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%catchem_var, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%dimensions'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%dimensions, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%units'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%units, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%coordinates'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%coordinates, yaml_rc)

         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%cell_methods'
         call cc_yaml_get(config_yaml, trim(field_key), &
            output_config%fields(i)%cell_methods, yaml_rc)

         ! Optional field
         write(field_key, '(A,I0,A)') 'diagnostic_fields%field_', i, '%optional'
         call cc_yaml_get(config_yaml, trim(field_key), temp_str, yaml_rc)
         if (yaml_rc == QFYAML_Success) then
            output_config%fields(i)%optional = (trim(temp_str) == 'true')
         end if
      end do

   end subroutine parse_output_config

   ! Create output directory
   !!
   !! \param rc  ESMF return code
   !!
   subroutine create_output_directory(rc)

      integer, intent(out) :: rc
      type(ESMF_VM) :: vm
      integer :: localPet

      rc = ESMF_SUCCESS

      ! Get VM info
      call ESMF_VMGetGlobal(vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_VMGet(vm, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Only root PET creates directory
      if (localPet == 0) then
         call system('mkdir -p ' // trim(output_config%output_directory))
      end if

      ! Barrier to ensure directory is created before other PETs proceed
      call ESMF_VMBarrier(vm, rc=rc)

   end subroutine create_output_directory

   ! Set up output timing
   !!
   !! \param  start_time  Simulation start time
   !! \param rc          ESMF return code
   !!
   subroutine setup_output_timing(start_time, rc)

      type(ESMF_Time), intent(in) :: start_time
      integer, intent(out) :: rc

      integer :: frequency_seconds

      rc = ESMF_SUCCESS

      ! Parse output frequency
      call parse_frequency_string(output_config%output_frequency, frequency_seconds, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Create output interval
      call ESMF_TimeIntervalSet(output_interval, s=frequency_seconds, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Set initial output time
      last_output_time = start_time - output_interval

   end subroutine setup_output_timing

   ! Create output fields and bundle
   !!
   !! \param rc  ESMF return code
   !!
   subroutine create_output_fields(rc)

      integer, intent(out) :: rc
      integer :: i

      rc = ESMF_SUCCESS

      ! Create field bundle
      output_bundle = ESMF_FieldBundleCreate(name="CATChem_diagnostics", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Create individual fields
      do i = 1, output_config%n_fields
         call create_diagnostic_field(output_config%fields(i), rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         ! Add field to bundle
         call ESMF_FieldBundleAdd(output_bundle, (/output_config%fields(i)%esmf_field/), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
      end do

   end subroutine create_output_fields

   ! Create a single diagnostic field
   !!
   !! \param field_config  Field configuration
   !! \param   rc            ESMF return code
   !!
   subroutine create_diagnostic_field(field_config, rc)

      type(diagnostic_field_type), intent(inout) :: field_config
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      ! Create field based on dimensions
      if (field_config%dimensions == 3) then
         field_config%esmf_field = ESMF_FieldCreate(output_grid, &
            typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1/), ungriddedUBound=(/64/), & ! Assuming 64 levels
            name=trim(field_config%variable_name), rc=rc)
      else
         field_config%esmf_field = ESMF_FieldCreate(output_grid, &
            typekind=ESMF_TYPEKIND_R8, &
            name=trim(field_config%variable_name), rc=rc)
      end if

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Set field attributes for CF compliance
      call set_field_attributes(field_config, rc)

   end subroutine create_diagnostic_field

   ! Set CF-compliant field attributes
   !!
   !! \param  field_config  Field configuration
   !! \param rc            ESMF return code
   !!
   subroutine set_field_attributes(field_config, rc)

      type(diagnostic_field_type), intent(in) :: field_config
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      ! Set CF convention attributes
      call ESMF_AttributeSet(field_config%esmf_field, &
         name="standard_name", value=trim(field_config%standard_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field_config%esmf_field, &
         name="long_name", value=trim(field_config%long_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field_config%esmf_field, &
         name="units", value=trim(field_config%units), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      if (len_trim(field_config%coordinates) > 0) then
         call ESMF_AttributeSet(field_config%esmf_field, &
            name="coordinates", value=trim(field_config%coordinates), rc=rc)
      end if

      if (len_trim(field_config%cell_methods) > 0) then
         call ESMF_AttributeSet(field_config%esmf_field, &
            name="cell_methods", value=trim(field_config%cell_methods), rc=rc)
      end if

   end subroutine set_field_attributes

   ! Check if it's time to write output
   !!
   !! \param  current_time   Current simulation time
   !! \param time_to_write  True if it's time to write
   !! \param rc             ESMF return code
   !!
   subroutine check_output_time(current_time, time_to_write, rc)

      type(ESMF_Time), intent(in) :: current_time
      logical, intent(out) :: time_to_write
      integer, intent(out) :: rc

      type(ESMF_Time) :: next_output_time

      rc = ESMF_SUCCESS
      time_to_write = .false.

      ! Calculate next output time
      next_output_time = last_output_time + output_interval

      ! Check if current time is at or past next output time
      if (current_time >= next_output_time) then
         time_to_write = .true.
      end if

   end subroutine check_output_time

   ! Generate filename for current time
   !!
   !! \param  current_time  Current simulation time
   !! \param filename      Generated filename
   !! \param rc            ESMF return code
   !!
   subroutine generate_filename(current_time, filename, rc)

      type(ESMF_Time), intent(in) :: current_time
      character(len=*), intent(out) :: filename
      integer, intent(out) :: rc

      integer :: year, month, day, hour, minute, second
      character(len=256) :: template_str, time_str

      rc = ESMF_SUCCESS

      ! Get time components
      call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, &
         h=hour, m=minute, s=second, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Create filename from template
      template_str = output_config%filename_template

      ! Replace time format specifiers
      write(time_str, '(I4.4)') year
      call replace_string(template_str, '%Y', trim(time_str))

      write(time_str, '(I2.2)') month
      call replace_string(template_str, '%m', trim(time_str))

      write(time_str, '(I2.2)') day
      call replace_string(template_str, '%d', trim(time_str))

      write(time_str, '(I2.2)') hour
      call replace_string(template_str, '%H', trim(time_str))

      write(time_str, '(I2.2)') minute
      call replace_string(template_str, '%M', trim(time_str))

      write(time_str, '(I2.2)') second
      call replace_string(template_str, '%S', trim(time_str))

      ! Combine with output directory
      filename = trim(output_config%output_directory) // '/' // trim(template_str)

   end subroutine generate_filename

   ! Populate output fields with data from CATChem states
   !!
   !! \param  catchem_states  CATChem state container
   !! \param rc              ESMF return code
   !!
   subroutine populate_output_fields(catchem_states, rc)

      type(catchem_container_type), intent(in) :: catchem_states
      integer, intent(out) :: rc

      integer :: i

      rc = ESMF_SUCCESS

      ! Copy data for each diagnostic field
      do i = 1, output_config%n_fields
         if (output_config%fields(i)%enabled) then
            call copy_field_data(output_config%fields(i), catchem_states, rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
         end if
      end do

   end subroutine populate_output_fields

   ! Copy data from CATChem state to output field
   !!
   !! \param  field_config    Field configuration
   !! \param  catchem_states  CATChem state container
   !! \param rc              ESMF return code
   !!
   subroutine copy_field_data(field_config, catchem_states, rc)

      type(diagnostic_field_type), intent(in) :: field_config
      type(catchem_container_type), intent(in) :: catchem_states
      integer, intent(out) :: rc

      real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:), fptr2d(:,:)
      integer :: species_id

      rc = ESMF_SUCCESS

      ! Copy data based on field mapping (similar to transform routines)
      select case (trim(field_config%catchem_var))

         ! 3D chemical species fields
       case ("ChemState%Species%O3")
         call ESMF_FieldGet(field_config%esmf_field, farrayPtr=fptr3d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         call cc_find_species_by_name(catchem_states%ChemState(1), "O3", species_id)
         if (species_id > 0) then
            ! Copy data from all columns
            call copy_3d_species_data(catchem_states, species_id, fptr3d)
         else
            fptr3d = 0.0_ESMF_KIND_R8
         end if

         ! Add more cases for other species and diagnostic fields...
         ! (Similar pattern to the transform routines)

       case default
         call ESMF_LogWrite("Unknown diagnostic field: " // trim(field_config%catchem_var), &
            ESMF_LOGMSG_WARNING, rc=rc)

      end select

   end subroutine copy_field_data

   ! Write NetCDF file using ESMF I/O
   !!
   !! \param  current_time  Current simulation time
   !! \param rc            ESMF return code
   !!
   subroutine write_netcdf_file(current_time, rc)

      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      type(ESMF_IO) :: io_obj

      rc = ESMF_SUCCESS

      ! Create ESMF I/O object
      io_obj = ESMF_IOCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Write the field bundle to NetCDF
      call ESMF_IOWrite(io_obj, filename=trim(current_filename), &
         fieldbundle=output_bundle, &
         convention="CF", purpose="history", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Add global attributes for CF compliance
      call add_global_attributes(io_obj, current_time, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Clean up I/O object
      call ESMF_IODestroy(io_obj, rc=rc)

   end subroutine write_netcdf_file

   ! Add global attributes to NetCDF file
   !!
   !! \param  io_obj        ESMF I/O object
   !! \param  current_time  Current simulation time
   !! \param rc            ESMF return code
   !!
   subroutine add_global_attributes(io_obj, current_time, rc)

      type(ESMF_IO), intent(in) :: io_obj
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      character(len=256) :: time_string, history_string

      rc = ESMF_SUCCESS

      ! Set CF convention attributes
      call ESMF_AttributeSet(io_obj, name="Conventions", &
         value=trim(output_config%cf_conventions), rc=rc)

      call ESMF_AttributeSet(io_obj, name="institution", &
         value=trim(output_config%institution), rc=rc)

      call ESMF_AttributeSet(io_obj, name="source", &
         value=trim(output_config%source), rc=rc)

      call ESMF_AttributeSet(io_obj, name="references", &
         value=trim(output_config%references), rc=rc)

      call ESMF_AttributeSet(io_obj, name="comment", &
         value=trim(output_config%comment), rc=rc)

      ! Add creation time
      call ESMF_TimeGet(current_time, timeStringISOFrac=time_string, rc=rc)
      write(history_string, '(A,A)') "Created on ", trim(time_string)
      call ESMF_AttributeSet(io_obj, name="history", value=trim(history_string), rc=rc)

   end subroutine add_global_attributes

   ! Helper routines

   ! Copy 3D species data from CATChem states
   subroutine copy_3d_species_data(catchem_states, species_id, output_array)
      type(catchem_container_type), intent(in) :: catchem_states
      integer, intent(in) :: species_id
      real(ESMF_KIND_R8), intent(out) :: output_array(:,:,:)

      integer :: i, j, k

      ! Copy data from all columns to output array
      do i = 1, catchem_states%im
         do k = 1, size(output_array, 3)
            do j = 1, size(output_array, 2)
               output_array(i,j,k) = catchem_states%ChemState(i)%Species(species_id)%Conc(j,k)
            end do
         end do
      end do
   end subroutine copy_3d_species_data

   ! Parse frequency string (e.g., "1h", "30m")
   subroutine parse_frequency_string(freq_str, seconds, rc)
      character(len=*), intent(in) :: freq_str
      integer, intent(out) :: seconds
      integer, intent(out) :: rc

      integer :: value, pos
      character(len=8) :: unit

      rc = ESMF_SUCCESS

      ! Find position of unit character
      pos = len_trim(freq_str)
      unit = freq_str(pos:pos)
      read(freq_str(1:pos-1), *) value

      select case (unit)
       case ('s')
         seconds = value
       case ('m')
         seconds = value * 60
       case ('h')
         seconds = value * 3600
       case ('d')
         seconds = value * 86400
       case default
         seconds = 3600  ! Default to 1 hour
         rc = ESMF_FAILURE
      end select
   end subroutine parse_frequency_string

   ! Replace substring in string
   subroutine replace_string(string, old_substr, new_substr)
      character(len=*), intent(inout) :: string
      character(len=*), intent(in) :: old_substr, new_substr

      integer :: pos

      pos = index(string, old_substr)
      if (pos > 0) then
         string = string(1:pos-1) // new_substr // string(pos+len(old_substr):)
      end if
   end subroutine replace_string

end module catchem_nuopc_netcdf_out
