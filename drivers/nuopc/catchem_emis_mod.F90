!> \file catchem_emis_mod.F90
!! \brief External emission data management for CATChem NUOPC interface
!! \ingroup catchem_nuopc_group
!!
!! \details
!! This module provides comprehensive external emission data management
!! for the CATChem NUOPC interface, following patterns from the AQM emission
!! module. It handles:
!! - Initialization of emission data structures from YAML configuration
!! - Reading and parsing emission files using AQMIO module
!! - Time interpolation and data updates for temporal emission data
!! - Population of ExtEmisDataType structures for CATChem processes
!! - Cleanup and finalization of emission resources
!!
!! The module integrates with CATChem's ConfigManager and ExtEmisData
!! infrastructure to provide consistent emission handling across different
!! emission categories (anthropogenic, point sources, fires, etc.)
!!
!! \author CATChem Development Team
!! \date January 2026
!! \version 1.0
!!
!! \section catchem_emis_usage Usage Example
!! \code{.f90}
!! use catchem_emis_mod
!! integer :: rc
!! call catchem_emis_init(cc_wrap, config_file, rc)
!! call catchem_emis_update(cc_wrap, current_time, rc)
!! call catchem_emis_finalize(cc_wrap, rc)
!! \endcode

module catchem_emis_mod

   use ESMF
   use NUOPC
   use aqmio
   use netcdf
   use catchem_regrid_mod, only: RegridCache, catchem_regrid_field, catchem_regrid_cleanup
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use ConfigManager_Mod, only: ConfigManagerType, ConfigDataType, EmissionCategoryMapping, &
      EmisSpeciesMappingEntry, EmissionMappingConfig
   use StateManager_Mod, only: StateManagerType
   use ChemState_Mod, only: ChemStateType
   use MetState_Mod, only: MetStateType
   use ExtEmisData_Mod, only: ExtEmisDataType, ExtEmisCategoryType, ExtEmisFieldType
   use Constants, only: AIRMW, AVO

   implicit none
   private

   ! Public interfaces
   public :: catchem_emis_init
   public :: catchem_emis_update
   public :: catchem_emis_finalize
   public :: catchem_emis_write_diagnostics
   public :: catchem_map_points_to_grid


   !> \brief Parameters for emission handling
   integer, parameter :: EMIS_MAXSTR = 256
   integer, parameter :: EMIS_MAXFIELDS = 100
   real(fp), parameter :: EMIS_MISSING = -999.0_fp

   !> Module-level regrid cache (weights computed once, reused)
   type(RegridCache), save :: emis_regrid_cache
   real(fp), parameter :: EMIS_ACCEPT = 1.e+15_fp ! Same as MAPL library "undefval"

   !> \brief Emission timing and alarm information
contains

   !> \brief Initialize emission data from configuration
   !!
   !! This subroutine initializes the external emission data system by:
   !! - Using already-loaded emission configuration from ConfigManagerType
   !! - Setting up emission categories and timing alarms
   !! - Populating ExtEmisDataType structures
   !!
   !! \param[inout] ext_emis_data External emission data container
   !! \param[in] config_manager Already loaded CATChem configuration manager
   !! \param[in] grid ESMF grid for I/O operations
   !! \param[out] rc Return code
   subroutine catchem_emis_init(ext_emis_data, config_manager, nx, ny, nlev, clock, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      type(ConfigManagerType), pointer, intent(in) :: config_manager
      integer, intent(in) :: nx, ny, nlev
      type(ESMF_Clock), intent(in) :: clock
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc,  icat
      logical :: extemis_activate
      character(len=EMIS_MAXSTR) :: msg
      character(len=*), parameter :: pName = 'catchem_emis_init'

      ! Initialize
      rc = CC_SUCCESS

      ! Check top-level processes/extemis/activate switch
      call config_manager%get_logical('processes/extemis/activate', extemis_activate, localrc, .true.)
      if (.not. extemis_activate) then
         call ESMF_LogWrite(trim(pName)//': External emissions disabled (processes/extemis/activate=false)', &
            ESMF_LOGMSG_INFO, rc=localrc)
         return
      end if

      ! Check if emission mapping is loaded
      if (.not. config_manager%config_data%emission_mapping%is_loaded) then
         write(msg, '(A,A)') trim(pName), ': Emission mapping not loaded in ConfigManager'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ! Initialize ExtEmisDataType with 0 to allow push-back population
      ! We start with 0 and let add_category grow the array incrementally
      call ext_emis_data%init(0, 'CATChem NUOPC Emission Data', localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to initialize ExtEmisDataType'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ! Enable global emission diagnostics - read from configuration or default to true
      call config_manager%get_logical('processes/extemis/global_diagnostics', ext_emis_data%diagnostic, localrc, .true.)

      ! Populate emission categories from already-loaded configuration
      do icat = 1, config_manager%config_data%emission_mapping%n_categories

         if (config_manager%config_data%emission_mapping%categories(icat)%is_active) then
            call catchem_emis_populate_category(ext_emis_data, &
               config_manager%config_data%emission_mapping%categories(icat), &
               config_manager, nx, ny, nlev, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': Failed to populate category ', &
                  trim(config_manager%config_data%emission_mapping%categories(icat)%category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
               rc = CC_FAILURE
               return
            end if

            call catchem_emis_setup_timing(ext_emis_data%categories(icat), clock, localrc)

         end if
      end do

      call ESMF_LogWrite(trim(pName)//': Emission initialization completed', &
         ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_init

   !> \brief Update emission data for current time
   !!
   !! Checks emission alarms and reads new emission data when needed.
   !! Handles time interpolation for temporal emission data.
   !!
   !! \param[inout] ext_emis_data External emission data container
   !! \param[in] current_time Current model time
   !! \param[out] rc Return code
   subroutine catchem_emis_update(ext_emis_data, current_time, state_manager, IO, grid, dt, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      type(ESMF_Time), intent(in) :: current_time
      type(StateManagerType), intent(inout) :: state_manager
      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      ! Local variables
      type(ConfigManagerType),pointer :: config_manager
      type(ErrorManagerType), pointer :: error_manager
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      integer :: localrc, i, period_key
      integer :: blo_year, blo_month
      real(fp) :: bfrac
      character(len=EMIS_MAXSTR) :: msg, timeString
      character(len=*), parameter :: pName = 'catchem_emis_update'

      rc = CC_SUCCESS

      ! Skip if no emission categories were initialized (e.g. extemis disabled)
      if (ext_emis_data%n_categories == 0) return

      ! Get managers from state manager
      config_manager => state_manager%get_config_ptr()
      error_manager => state_manager%get_error_manager()
      met_state => state_manager%get_met_state_ptr()
      chem_state => state_manager%get_chem_state_ptr()

      ! Loop through all emission categories and check if updates are needed
      do i = 1, ext_emis_data%n_categories
         if (.not. ext_emis_data%categories(i)%is_active) cycle

         ! Determine the current calendar period key for this category's frequency.
         ! A period key encodes the calendar unit that triggers a new read:
         !   daily -> yyyymmdd, monthly -> yyyymm, hourly -> yyyymmddhh, static -> 0.
         ! When the key differs from last_period_key (including the sentinel -1 on
         ! the first call), the emission data must be re-read.  This approach is
         ! immune to the alarm-drift problem that occurs when simulations do not
         ! start at a "natural" boundary (e.g. 06:00 start with daily data).
         call catchem_emis_period_key(ext_emis_data%categories(i)%frequency, &
            current_time, period_key, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return

         ! For a single multi-record monthly climatology with linear time
         ! interpolation, the interpolation bracket changes at mid-month (not at
         ! the month start), so re-read the two bracketing slices at mid-month.
         if (trim(ext_emis_data%categories(i)%frequency) == 'monthly' .and. &
            trim(ext_emis_data%categories(i)%time_interpolation) == 'linear' .and. &
            ext_emis_data%categories(i)%n_times >= 2) then
            call catchem_emis_month_bracket(current_time, blo_year, blo_month, bfrac, localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return
            period_key = blo_year*100 + blo_month
         end if

         if (period_key /= ext_emis_data%categories(i)%last_period_key) then

            call ESMF_TimeGet(current_time, timeString=timeString, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return

            call ESMF_LogWrite(trim(pName)//': reading emission for '// &
               trim(ext_emis_data%categories(i)%category_name)// &
               " @ "//trim(timeString), ESMF_LOGMSG_INFO, rc=localrc)

            ! For files without time-coordinate matching and no filename template,
            ! advance irec sequentially (one slice per period).
            if (ext_emis_data%categories(i)%n_times == 0 .and. &
               index(trim(ext_emis_data%categories(i)%source_file), '%') == 0) then
               ext_emis_data%categories(i)%irec = ext_emis_data%categories(i)%irec + 1
            end if

            call catchem_emis_read(ext_emis_data%categories(i), IO, grid, &
               met_state%NLEVS, current_time, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': Failed to read data for category: ', &
                  trim(ext_emis_data%categories(i)%category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
            end if

            ext_emis_data%categories(i)%last_period_key = period_key
         end if

         ! Recompute temporal blend weights every timestep for time-interpolated categories
         if (ext_emis_data%categories(i)%needs_time_blend) then
            call catchem_emis_blend_time(ext_emis_data%categories(i), current_time, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': Failed to blend time for category: ', &
                  trim(ext_emis_data%categories(i)%category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
            end if
         end if

         ! Apply emissions to chemical state every timestep
         ! (data is read only when the period changes, but applied every step)
         call catchem_emis_apply(ext_emis_data%categories(i), i, ext_emis_data%global_scale, config_manager, error_manager, chem_state, met_state, dt, current_time, localrc)
         if (localrc /= CC_SUCCESS) then
            write(msg, '(A,A,A)') trim(pName), ': Failed to apply emissions for category: ', &
               trim(ext_emis_data%categories(i)%category_name)
            call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
         end if
      end do
      nullify(config_manager, met_state, chem_state) ! Clean up pointers

      call ESMF_LogWrite(trim(pName)//': Emission data updated', &
         ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_update

   !> \brief Read emission data from files
   !!
   !! Reads emission data from NetCDF files using AQMIO module
   !! and populates the ExtEmisFieldType structures.
   !!
   !! \param[inout] ext_emis_data External emission data container
   !! \param[in] category_name Name of emission category to read
   !! \param[out] rc Return code
   subroutine catchem_emis_read(category, IO, grid, nlev, curr_time, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      integer, intent(in) :: nlev
      type(ESMF_Time), intent(in) :: curr_time
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc,   ifield
      character(len=EMIS_MAXSTR) :: msg, filename
      character(len=64) :: category_name
      type(ESMF_Field) :: esmf_field
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      real(ESMF_KIND_R4), pointer :: field_data_3d(:,:,:) => null()
      character(len=*), parameter :: pName = 'catchem_emis_read'
      logical :: use_regrid
      logical :: didRegrid
      logical :: file_exists

      rc = CC_SUCCESS

      category_name = trim(category%category_name)

      ! Point/volcanic categories are read from an ASCII point table (.rc), not a
      ! gridded NetCDF file, and are injected directly into the 3D column at apply
      ! time.  Dispatch to the dedicated reader and skip the gridded I/O path.
      if (is_point_category(category)) then
         call catchem_emis_read_points(category, curr_time, rc)
         return
      end if

      ! Resolve filename: substitute date tokens if the template contains '%'
      if (index(trim(category%source_file), '%') > 0) then
         call resolve_filename_template(category%source_file, curr_time, filename, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      else
         filename = trim(category%source_file)
      end if

      if (len_trim(filename) == 0) then
         write(msg, '(A,A,A)') trim(pName), ': No source file specified for category: ', trim(category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ! For template-resolved filenames, check that the file exists before attempting I/O.
      ! If missing, log a warning and keep the last loaded data unchanged.
      if (index(trim(category%source_file), '%') > 0) then
         inquire(file=trim(filename), exist=file_exists)
         if (.not. file_exists) then
            write(msg, '(A,A,A)') trim(pName), ': File not found (holding last data): ', trim(filename)
            call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_WARNING, rc=localrc)
            return
         end if
      end if

      ! Populate time-coordinate cache if the file has changed (or first call)
      if (trim(filename) /= trim(category%last_resolved_file) .or. category%n_times == 0) then
         if (trim(category%frequency) /= 'static') then
            call catchem_emis_read_time_coord(filename, category, localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) return
         end if
         category%last_resolved_file = trim(filename)
      end if

      ! Compute the correct time-slice index from cached time coordinates
      if (category%n_times > 0) then
         call catchem_emis_find_time_index(category, curr_time, category%frequency, &
            category%irec, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      end if

      ! Determine if this category needs runtime regridding.
      ! When regrid_method is set to anything other than 'none' (e.g.
      ! bilinear, neareststod, conserve, ...) the file is assumed to be
      ! on a different grid and will be regridded to the model grid.
      use_regrid = (trim(category%regrid_method) /= 'none' .and. &
         trim(category%regrid_method) /= 'NONE' .and. &
         len_trim(category%regrid_method) > 0)

      if (use_regrid) then
         if (len_trim(category%latname) == 0 .or. len_trim(category%lonname) == 0) then
            write(msg, '(A,A,A)') trim(pName), &
               ': regrid_method set but lat_name/lon_name missing for category: ', &
               trim(category_name)
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
            rc = CC_FAILURE
            return
         end if
         call catchem_emis_read_regrid(category, grid, nlev, filename, curr_time, rc)
         return
      end if

      !open file (Note: although AQMIO_Read can open file in its source code, it gives zeros for some reason.
      !           So we have to open it here first.)
      call AQMIO_Open(IO, filename, iomode="read", iofmt=AQMIO_FMT_NETCDF, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

      do ifield = 1, category%n_fields
         !create field to receive data
         if (category%is_2d) then
            esmf_field = ESMF_FieldCreate(grid, name=trim(category%fields(ifield)%field_name), &
               typekind=ESMF_TYPEKIND_R4, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
         else !3D field
            esmf_field = ESMF_FieldCreate(grid, name=trim(category%fields(ifield)%field_name), &
               typekind=ESMF_TYPEKIND_R4, ungriddedLBound=(/1/), ungriddedUBound=(/nlev/), rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
         end if

         !read data into field
         call AQMIO_Read(IO, (/ esmf_field /), fieldNameList=(/ trim(category%fields(ifield)%field_name) /), &
            timeSlice=category % irec, iofmt=AQMIO_FMT_NETCDF, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) then
            ! Clean up field before returning
            call ESMF_FieldDestroy(esmf_field, rc=localrc)
            return  ! bail out
         end if

         if (category%is_2d) then
            !get data pointer and assign to emission field array
            call ESMF_FieldGet(esmf_field, farrayPtr=field_data_2d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) then
               ! Clean up field before returning
               call ESMF_FieldDestroy(esmf_field, rc=localrc)
               return  ! bail out
            end if
            !!TODO: We should check unit conversion in the future. Here we make sure the gridded emission is in kg/m2/s already
            category%fields(ifield)%emission_data(:,:,1,1) = real(field_data_2d(:,:), fp)  !assuming 2D data for now
         else !3D field
            !get data pointer and assign to emission field array
            call ESMF_FieldGet(esmf_field, farrayPtr=field_data_3d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) then
               ! Clean up field before returning
               call ESMF_FieldDestroy(esmf_field, rc=localrc)
               return  ! bail out
            end if
            !!TODO: We should check unit conversion in the future. Here we make sure the gridded emission is in kg/m2/s already
            if (category%reverse_vertical) then
               category%fields(ifield)%emission_data(:,:,:,1) = real(field_data_3d(:,:,nlev:1:-1), fp)  !reverse vertical level
            else
               category%fields(ifield)%emission_data(:,:,:,1) = real(field_data_3d(:,:,:), fp)
            end if
         end if

         category%fields(ifield)%is_loaded = .true.   !set to true; otherwise diagnostics will not be saved.

         ! Clean up ESMF field after data transfer
         call ESMF_FieldDestroy(esmf_field, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

         ! Nullify pointer for safety
         field_data_2d => null()
         field_data_3d => null()
      end do

      !!not sure why this write will crash the model
      write(msg, '(A,A,A)') trim(pName), ': Successfully read emission data for category ', &
         trim(category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_read

   !> \brief Read emission data with runtime regridding from lat-lon to model grid
   !!
   !! Reads global lat-lon emission data and regrids it onto the model
   !! grid using ESMF bilinear regridding.  Route handles are cached in
   !! the module-level emis_regrid_cache so weights are computed only once.
   subroutine catchem_emis_read_regrid(category, grid, nlev, filename, curr_time, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_Grid),          intent(in)    :: grid
      integer,                  intent(in)    :: nlev
      character(len=*),         intent(in)    :: filename
      type(ESMF_Time),          intent(in)    :: curr_time
      integer,                  intent(out)   :: rc

      ! Local variables
      integer :: localrc, ifield, klev
      character(len=EMIS_MAXSTR) :: msg
      character(len=64) :: category_name
      type(ESMF_Field) :: esmf_field
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      logical :: didRegrid
      character(len=*), parameter :: pName = 'catchem_emis_read_regrid'

      ! Temporal interpolation variables
      logical :: do_time_interp
      logical :: multi_file_interp   ! t2 comes from a different file
      integer :: irec_next
      integer :: nx, ny
      character(len=EMIS_MAXSTR) :: filename_next
      type(ESMF_Time) :: next_time
      type(ESMF_TimeInterval) :: period_step
      logical :: next_file_exists

      rc = CC_SUCCESS
      category_name = trim(category%category_name)

      ! Determine if temporal interpolation is needed.
      ! Two modes:
      !   (a) single multi-record file: n_times >= 2
      !   (b) separate files per period (template with %): n_times <= 1
      multi_file_interp = .false.
      do_time_interp = (trim(category%time_interpolation) == 'linear')
      if (do_time_interp) then
         if (category%n_times >= 2) then
            ! Single file with multiple time records — use next record in same file
            multi_file_interp = .false.
         else if (index(trim(category%source_file), '%') > 0) then
            ! Template-based: each file has 1 record, read t2 from next-period file
            multi_file_interp = .true.
         else
            ! Single-record file, no template — cannot interpolate
            do_time_interp = .false.
         end if
      end if

      ! Compute next record index or resolve next-period filename
      if (do_time_interp) then
         if (multi_file_interp) then
            ! Advance curr_time by one period to resolve the next file
            irec_next = 1  ! next file's first (only) record
            select case (trim(category%frequency))
             case ('monthly')
               call ESMF_TimeIntervalSet(period_step, mm=1, rc=localrc)
             case ('daily')
               call ESMF_TimeIntervalSet(period_step, d=1, rc=localrc)
             case ('hourly')
               call ESMF_TimeIntervalSet(period_step, h=1, rc=localrc)
             case default
               call ESMF_TimeIntervalSet(period_step, mm=1, rc=localrc)
            end select
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) return

            next_time = curr_time + period_step
            call resolve_filename_template(category%source_file, next_time, filename_next, localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) return

            ! Check next file exists; fall back to no interpolation if missing
            inquire(file=trim(filename_next), exist=next_file_exists)
            if (.not. next_file_exists) then
               write(msg, '(A,A,A)') trim(pName), &
                  ': next-period file not found, disabling time_interp: ', trim(filename_next)
               call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_WARNING, rc=localrc)
               do_time_interp = .false.
               multi_file_interp = .false.
            else
               write(msg, '(A,A,A,A)') trim(pName), &
                  ': multi-file time_interp for ', trim(category_name), &
                  ' next_file='//trim(filename_next)
               call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=localrc)
            end if
         else
            ! Same-file interpolation
            if (category%irec < category%n_times) then
               irec_next = category%irec + 1
            else
               irec_next = 1  ! wrap around for climatological data (Dec -> Jan)
            end if

            write(msg, '(A,A,A,I3,A,I3)') trim(pName), &
               ': time_interp read for ', trim(category_name), &
               category%irec, '  and next=', irec_next
            call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=localrc)
         end if
      end if

      do ifield = 1, category%n_fields
         ! Create 2D destination field on the model grid
         esmf_field = ESMF_FieldCreate(grid, &
            name=trim(category%fields(ifield)%field_name), &
            typekind=ESMF_TYPEKIND_R4, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         if (category%is_2d) then
            ! --- 2D field ---
            call catchem_regrid_field( &
               cache     = emis_regrid_cache, &
               filename  = trim(filename), &
               varname   = trim(category%fields(ifield)%field_name), &
               dstField  = esmf_field, &
               latname   = trim(category%latname), &
               lonname   = trim(category%lonname), &
               regrid_method_name = trim(category%regrid_method), &
               timeSlice = category%irec, &
               didRegrid = didRegrid, &
               rc        = localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) then
               call ESMF_FieldDestroy(esmf_field, rc=localrc)
               return
            end if

            call ESMF_FieldGet(esmf_field, farrayPtr=field_data_2d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) then
               call ESMF_FieldDestroy(esmf_field, rc=localrc)
               return
            end if

            if (do_time_interp) then
               ! Store current time slice in interp_data_t1
               nx = size(field_data_2d, 1)
               ny = size(field_data_2d, 2)
               if (.not. allocated(category%fields(ifield)%interp_data_t1)) then
                  allocate(category%fields(ifield)%interp_data_t1(nx, ny, 1, 1))
               end if
               category%fields(ifield)%interp_data_t1(:,:,1,1) = real(field_data_2d(:,:), fp)

               ! Regrid the next time slice (from same file or next-period file)
               if (multi_file_interp) then
                  call catchem_regrid_field( &
                     cache     = emis_regrid_cache, &
                     filename  = trim(filename_next), &
                     varname   = trim(category%fields(ifield)%field_name), &
                     dstField  = esmf_field, &
                     latname   = trim(category%latname), &
                     lonname   = trim(category%lonname), &
                     regrid_method_name = trim(category%regrid_method), &
                     timeSlice = irec_next, &
                     didRegrid = didRegrid, &
                     rc        = localrc)
               else
                  call catchem_regrid_field( &
                     cache     = emis_regrid_cache, &
                     filename  = trim(filename), &
                     varname   = trim(category%fields(ifield)%field_name), &
                     dstField  = esmf_field, &
                     latname   = trim(category%latname), &
                     lonname   = trim(category%lonname), &
                     regrid_method_name = trim(category%regrid_method), &
                     timeSlice = irec_next, &
                     didRegrid = didRegrid, &
                     rc        = localrc)
               end if
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)) then
                  call ESMF_FieldDestroy(esmf_field, rc=localrc)
                  return
               end if

               ! Store next time slice in interp_data_t2
               if (.not. allocated(category%fields(ifield)%interp_data_t2)) then
                  allocate(category%fields(ifield)%interp_data_t2(nx, ny, 1, 1))
               end if
               category%fields(ifield)%interp_data_t2(:,:,1,1) = real(field_data_2d(:,:), fp)

               ! Mark category for per-timestep blending
               category%needs_time_blend = .true.

               ! Initialize emission_data with a placeholder (will be overwritten by blend)
               category%fields(ifield)%emission_data(:,:,1,1) = &
                  category%fields(ifield)%interp_data_t1(:,:,1,1)
            else
               category%fields(ifield)%emission_data(:,:,1,1) = real(field_data_2d(:,:), fp)
            end if
         else
            ! --- 3D field: regrid each vertical level as a 2D slab ---
            do klev = 1, nlev
               call catchem_regrid_field( &
                  cache     = emis_regrid_cache, &
                  filename  = trim(filename), &
                  varname   = trim(category%fields(ifield)%field_name), &
                  dstField  = esmf_field, &
                  latname   = trim(category%latname), &
                  lonname   = trim(category%lonname), &
                  regrid_method_name = trim(category%regrid_method), &
                  timeSlice = category%irec, &
                  levelSlice = klev, &
                  didRegrid = didRegrid, &
                  rc        = localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)) then
                  call ESMF_FieldDestroy(esmf_field, rc=localrc)
                  return
               end if

               call ESMF_FieldGet(esmf_field, farrayPtr=field_data_2d, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__, rcToReturn=rc)) then
                  call ESMF_FieldDestroy(esmf_field, rc=localrc)
                  return
               end if

               if (do_time_interp) then
                  nx = size(field_data_2d, 1)
                  ny = size(field_data_2d, 2)
                  if (.not. allocated(category%fields(ifield)%interp_data_t1)) then
                     allocate(category%fields(ifield)%interp_data_t1(nx, ny, nlev, 1))
                  end if
                  category%fields(ifield)%interp_data_t1(:,:,klev,1) = real(field_data_2d(:,:), fp)

                  if (multi_file_interp) then
                     call catchem_regrid_field( &
                        cache     = emis_regrid_cache, &
                        filename  = trim(filename_next), &
                        varname   = trim(category%fields(ifield)%field_name), &
                        dstField  = esmf_field, &
                        latname   = trim(category%latname), &
                        lonname   = trim(category%lonname), &
                        regrid_method_name = trim(category%regrid_method), &
                        timeSlice = irec_next, &
                        levelSlice = klev, &
                        didRegrid = didRegrid, &
                        rc        = localrc)
                  else
                     call catchem_regrid_field( &
                        cache     = emis_regrid_cache, &
                        filename  = trim(filename), &
                        varname   = trim(category%fields(ifield)%field_name), &
                        dstField  = esmf_field, &
                        latname   = trim(category%latname), &
                        lonname   = trim(category%lonname), &
                        regrid_method_name = trim(category%regrid_method), &
                        timeSlice = irec_next, &
                        levelSlice = klev, &
                        didRegrid = didRegrid, &
                        rc        = localrc)
                  end if
                  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                     line=__LINE__, file=__FILE__, rcToReturn=rc)) then
                     call ESMF_FieldDestroy(esmf_field, rc=localrc)
                     return
                  end if

                  if (.not. allocated(category%fields(ifield)%interp_data_t2)) then
                     allocate(category%fields(ifield)%interp_data_t2(nx, ny, nlev, 1))
                  end if
                  category%fields(ifield)%interp_data_t2(:,:,klev,1) = real(field_data_2d(:,:), fp)

                  category%needs_time_blend = .true.
                  category%fields(ifield)%emission_data(:,:,klev,1) = &
                     category%fields(ifield)%interp_data_t1(:,:,klev,1)
               else
                  category%fields(ifield)%emission_data(:,:,klev,1) = real(field_data_2d(:,:), fp)
               end if
            end do

            ! Reverse vertical levels if configured (apply to both stored slices)
            if (category%reverse_vertical) then
               category%fields(ifield)%emission_data(:,:,:,1) = &
                  category%fields(ifield)%emission_data(:,:,nlev:1:-1,1)
               if (do_time_interp) then
                  category%fields(ifield)%interp_data_t1(:,:,:,1) = &
                     category%fields(ifield)%interp_data_t1(:,:,nlev:1:-1,1)
                  category%fields(ifield)%interp_data_t2(:,:,:,1) = &
                     category%fields(ifield)%interp_data_t2(:,:,nlev:1:-1,1)
               end if
            end if
         end if

         category%fields(ifield)%is_loaded = .true.

         call ESMF_FieldDestroy(esmf_field, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         field_data_2d => null()
      end do

      write(msg, '(A,A,A)') trim(pName), &
         ': Successfully read & regridded emission data for category ', trim(category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_read_regrid

   !> \brief Get emission data for a specific field and location
   !!
   !! Returns emission rates for specified field at grid location.
   !! Provides interface similar to aqm_emis_get.
   !!
   !! \param[in] ext_emis_data External emission data container
   !! \param[in] category_name Name of emission category
   !! \param[in] field_name Name of emission field
   !! \param[in] i Longitude index
   !! \param[in] j Latitude index
   !! \param[in] k Vertical index (optional)
   !! \return Emission rate [kg/m2/s]
   function catchem_emis_get(ext_emis_data, field_name, i, j, k) result(emission_rate)
      implicit none

      type(ExtEmisDataType), intent(in) :: ext_emis_data
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: i, j
      integer, intent(in), optional :: k
      real(fp) :: emission_rate

      ! Local variables
      integer :: kk
      real(fp) :: rate

      kk = 1
      if (present(k)) kk = k

      ! Get emission rate from ExtEmisDataType
      rate = ext_emis_data%get_emission_rate(field_name, i, j, kk)

      ! Apply any additional scaling or processing
      emission_rate = rate

   end function catchem_emis_get

   !> \brief Apply emission data to chemical state
   !!
   !> \brief Distribute 2D surface emissions vertically based on specified method
   !!
   !! Based on GOCART2G SulfateDistributeEmissions and distribute_aviation_emissions.
   !! Supports four methods:
   !! - P100: Distribute uniformly from surface to 100m altitude
   !! - P500: Distribute uniformly from 100m to 500m altitude
   !! - Ppbl: Distribute uniformly from surface to PBL height (min 100m)
   !! - aviation: Distribute using hardcoded aviation layers
   !!   [0.0, 100.0, 9000.0, 10000.0] m covering LTO/CDS/CRS ranges
   !! - aviation_lto: LTO only (0 - 100m)
   !! - aviation_cds: CDS only (100 - 9000m)
   !! - aviation_crs: CRS only (9000 - 10000m)
   !!
   !! Uses pressure-based vertical fractions following the GOCART2G approach:
   !! compute pressure at target altitudes by walking from surface upward,
   !! then assign layer fractions proportional to the pressure overlap
   !! with the target range.
   !!
   !! UFS/FV3 convention: k=1 = surface, k=nz = top of atmosphere.
   !! Loops walk from k=1 (surface) to k=nz (top) for altitude calculations.
   !!
   !! \param[inout] emission_flux 3D emission flux array (nx,ny,nz); surface data in k=1
   !! \param[in] met_state Meteorological state (DELP, AIRDEN, PBLH)
   !! \param[in] vertical_dist Distribution method name
   !! \param[in] nx,ny,nz Grid dimensions
   subroutine distribute_emissions_vertical(emission_flux, met_state, vertical_dist, nx, ny, nz)
      use Constants, only: g0
      implicit none

      real(fp), intent(inout) :: emission_flux(:,:,:)
      type(MetStateType), intent(in) :: met_state
      character(len=*), intent(in) :: vertical_dist
      integer, intent(in) :: nx, ny, nz

      ! Local variables
      integer :: i, j, k
      real(fp) :: ps, p0, p1, z0_col, z1_col, dz, deltaz, deltap
      real(fp) :: p100, p500, pPBL, p9000, p10000, zpbl
      real(fp) :: f_dist, emis_sfc
      real(fp) :: p_top, p_bot  ! pressure range for distribution

      ! Hardcoded aviation emission layers [m] following GOCART2G convention:
      !   aviation_layers = [LTO_bot, LTO_top/CDS_bot, CDS_top/CRS_bot, CRS_top]
      !   LTO (Landing/Take-Off):     0 -   100 m
      !   CDS (Climb/Descent):      100 -  9000 m
      !   CRS (Cruise):            9000 - 10000 m
      real(fp), parameter :: AVN_LTO_BOT =     0.0_fp
      real(fp), parameter :: AVN_LTO_TOP =   100.0_fp
      real(fp), parameter :: AVN_CDS_TOP =  9.0e3_fp
      real(fp), parameter :: AVN_CRS_TOP = 10.0e3_fp

      ! UFS/FV3 convention: k=1 = surface, k=nz = top of atmosphere
      ! Return immediately if no distribution needed
      select case (trim(vertical_dist))
       case ('none', 'NONE', 'None', '')
         return
       case ('P100', 'p100', 'P500', 'p500', 'Ppbl', 'ppbl', 'PBL', 'pbl', &
          'aviation', 'AVIATION', &
          'aviation_lto', 'AVIATION_LTO', &
          'aviation_cds', 'AVIATION_CDS', &
          'aviation_crs', 'AVIATION_CRS')
         ! proceed
       case default
         call ESMF_LogWrite('distribute_emissions_vertical: Unrecognized vertical_dist option: '// &
            trim(vertical_dist)//'; skipping vertical distribution', &
            ESMF_LOGMSG_WARNING)
         return
      end select

      do j = 1, ny
         do i = 1, nx
            ! Save surface emission value (2D data is stored in k=1 slot)
            emis_sfc = emission_flux(i, j, 1)
            if (emis_sfc == 0.0_fp) cycle

            ! Compute surface pressure by summing all layer thicknesses
            ps = 0.0_fp
            do k = 1, nz
               ps = ps + met_state%DELP(i, j, k)
            end do

            ! Find pressure at target altitudes by walking from surface (k=1) upward (k=nz)
            p0 = ps
            z0_col = 0.0_fp
            p100   = 0.0_fp
            p500   = 0.0_fp
            pPBL   = 0.0_fp
            p9000  = 0.0_fp
            p10000 = 0.0_fp

            do k = 1, nz
               p1 = p0 - met_state%DELP(i, j, k)
               dz = met_state%DELP(i, j, k) / (met_state%AIRDEN(i, j, k) * g0)
               z1_col = z0_col + dz

               if (p100 == 0.0_fp .and. z0_col < 100.0_fp .and. z1_col >= 100.0_fp) then
                  deltaz = z1_col - 100.0_fp
                  deltap = deltaz * met_state%AIRDEN(i, j, k) * g0
                  p100 = p1 + deltap
               end if

               if (p500 == 0.0_fp .and. z0_col < 500.0_fp .and. z1_col >= 500.0_fp) then
                  deltaz = z1_col - 500.0_fp
                  deltap = deltaz * met_state%AIRDEN(i, j, k) * g0
                  p500 = p1 + deltap
               end if

               zpbl = max(met_state%PBLH(i, j), 100.0_fp)
               if (pPBL == 0.0_fp .and. z0_col < zpbl .and. z1_col >= zpbl) then
                  deltaz = z1_col - zpbl
                  deltap = deltaz * met_state%AIRDEN(i, j, k) * g0
                  pPBL = p1 + deltap
               end if

               if (p9000 == 0.0_fp .and. z0_col < AVN_CDS_TOP .and. z1_col >= AVN_CDS_TOP) then
                  deltaz = z1_col - AVN_CDS_TOP
                  deltap = deltaz * met_state%AIRDEN(i, j, k) * g0
                  p9000 = p1 + deltap
               end if

               if (p10000 == 0.0_fp .and. z0_col < AVN_CRS_TOP .and. z1_col >= AVN_CRS_TOP) then
                  deltaz = z1_col - AVN_CRS_TOP
                  deltap = deltaz * met_state%AIRDEN(i, j, k) * g0
                  p10000 = p1 + deltap
               end if

               p0 = p1
               z0_col = z1_col
            end do

            ! Fallback: if target height was never reached, use top-of-atmosphere pressure
            if (p100   == 0.0_fp) p100   = p0
            if (p500   == 0.0_fp) p500   = p0
            if (pPBL   == 0.0_fp) pPBL   = p0
            if (p9000  == 0.0_fp) p9000  = p0
            if (p10000 == 0.0_fp) p10000 = p0

            ! Determine pressure range for this distribution type
            ! p_bot = higher pressure (lower altitude), p_top = lower pressure (higher altitude)
            select case (trim(vertical_dist))
             case ('P100', 'p100')
               ! Surface to 100m
               p_bot = ps
               p_top = p100
             case ('P500', 'p500')
               ! 100m to 500m
               p_bot = p100
               p_top = p500
             case ('Ppbl', 'ppbl', 'PBL', 'pbl')
               ! Surface to PBL height
               p_bot = ps
               p_top = pPBL
             case ('aviation', 'AVIATION')
               ! Full aviation range: surface to CRS top (0 - 10000 m)
               ! Covers LTO (0-100m) + CDS (100-9000m) + CRS (9000-10000m)
               p_bot = ps
               p_top = p10000
             case ('aviation_lto', 'AVIATION_LTO')
               ! LTO only: surface to 100m
               p_bot = ps
               p_top = p100
             case ('aviation_cds', 'AVIATION_CDS')
               ! CDS only: 100m to 9000m
               p_bot = p100
               p_top = p9000
             case ('aviation_crs', 'AVIATION_CRS')
               ! CRS only: 9000m to 10000m
               p_bot = p9000
               p_top = p10000
             case default
               cycle
            end select

            ! Guard against zero or negative pressure range
            if (p_bot - p_top <= 0.0_fp) cycle

            ! Zero out all levels, then distribute using pressure fractions
            ! Walk from surface (k=1) to top (k=nz)
            emission_flux(i, j, :) = 0.0_fp

            p0 = ps
            do k = 1, nz
               p1 = p0 - met_state%DELP(i, j, k)

               ! Compute fractional overlap of this model layer with the target pressure range
               ! p0 = pressure at layer bottom (higher pressure, lower altitude)
               ! p1 = pressure at layer top (lower pressure, higher altitude)
               f_dist = 0.0_fp

               if (p0 <= p_bot .and. p1 >= p_top) then
                  ! Layer fully within target range
                  f_dist = met_state%DELP(i, j, k) / (p_bot - p_top)
               else if (p0 > p_bot .and. p1 >= p_top .and. p1 < p_bot) then
                  ! Layer straddles bottom boundary (extends below target)
                  f_dist = (p_bot - max(p1, p_top)) / (p_bot - p_top)
               else if (p0 <= p_bot .and. p0 > p_top .and. p1 < p_top) then
                  ! Layer straddles top boundary (extends above target)
                  f_dist = (min(p0, p_bot) - p_top) / (p_bot - p_top)
               else if (p0 > p_bot .and. p1 < p_top) then
                  ! Layer fully encompasses the target range
                  f_dist = 1.0_fp
               end if
               ! Otherwise: layer entirely outside range, f_dist = 0

               emission_flux(i, j, k) = emis_sfc * f_dist

               p0 = p1
            end do

         end do
      end do

   end subroutine distribute_emissions_vertical

   !> \brief Compute biomass burning emission scaling factor using Mie data
   !!
   !! Prevents unrealistically high aerosol optical thickness from biomass
   !! burning emissions by computing the extinction AOT that the emission
   !! would produce and scaling down if it exceeds max_bb_exttau (30.0).
   !! Follows GOCART2G CAEmission pattern.
   !!
   !! \param[in]  emission_flux  3D emission flux after vertical distribution [kg/m2/s]
   !! \param[in]  scale_factor   Species-specific scale factor from mapping
   !! \param[in]  dt             Time step [s]
   !! \param[in]  met_state      Meteorological state (for RH)
   !! \param[in]  chem_state     Chemical state (for MieData)
   !! \param[in]  species_idx    Species index in chem_state
   !! \param[out] f_bb           2D scaling factor [0..1] per column
   !! \param[out] rc             Return code
   subroutine compute_bb_emission_factor(emission_flux, scale_factor, dt, &
      met_state, chem_state, species_idx, &
      f_bb, rc)
      use Constants, only: g0
      implicit none

      real(fp), intent(in)    :: emission_flux(:,:,:)
      real(fp), intent(in)    :: scale_factor
      real(fp), intent(in)    :: dt
      type(MetStateType), intent(in)  :: met_state
      type(ChemStateType), intent(in) :: chem_state
      integer, intent(in)    :: species_idx
      real(fp), intent(out)   :: f_bb(:,:)
      integer, intent(out)   :: rc

      ! Local variables
      integer :: nx, ny, nz, i, j, k, mie_idx, ibin
      real, allocatable :: q_mass(:,:,:), rh_r4(:,:,:), tau(:,:,:)
      real(fp) :: exttau_bb, cutoff_bb_exttau
      integer :: localrc
      character(len=*), parameter :: pName = 'compute_bb_emission_factor'
      character(len=EMIS_MAXSTR) :: msg

      ! Parameters following GOCART2G CAEmission
      real(fp), parameter :: max_bb_exttau = 30.0_fp  ! daily maximum AOT from BB
      integer, parameter  :: nbin = 2  ! hardcoded for carbonaceous aerosols

      rc = CC_SUCCESS
      f_bb = 1.0_fp

      ! Scale daily max AOT to per-timestep cutoff (GOCART2G: cdt / (24*3600) * max_bb_exttau)
      cutoff_bb_exttau = (dt / 86400.0_fp) * max_bb_exttau

      ! Check species has Mie data
      if (.not. allocated(chem_state%SpcMieMap)) return
      if (species_idx < 1 .or. species_idx > size(chem_state%SpcMieMap)) return
      mie_idx = chem_state%SpcMieMap(species_idx)
      if (mie_idx <= 0) return

      nx = size(emission_flux, 1)
      ny = size(emission_flux, 2)
      nz = size(emission_flux, 3)

      ! Allocate working arrays as default real (GOCART2G_Mie uses default real)
      allocate(q_mass(nx, ny, nz), rh_r4(nx, ny, nz), tau(nx, ny, nz))

      ! Relative humidity clamped to [0, 0.99] for Mie table lookup
      rh_r4 = real(min(max(met_state%RH, 0.0_fp), 0.99_fp))

      ! Column mass from emission [kg/m2]: flux [kg/m2/s] * scale * dt [s]
      q_mass = real(emission_flux * scale_factor * dt)

      ! Sum extinction optical depth over all Mie bins
      do j = 1, ny
         do i = 1, nx
            exttau_bb = 0.0_fp
            do ibin = 1, min(nbin, chem_state%MieData(mie_idx)%nbin)
               call chem_state%MieData(mie_idx)%Query( &
                  550.0e-9, ibin, q_mass(i:i,j:j,:), rh_r4(i:i,j:j,:), &
                  tau=tau(i:i,j:j,:), rc=localrc)
               if (localrc /= CC_SUCCESS) then
                  write(msg, '(A,A,I0,A,I0)') trim(pName), &
                     ': Mie Query failed for species ', species_idx, ' bin ', ibin
                  call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
                  cycle
               end if
               do k = 1, nz
                  exttau_bb = exttau_bb + real(tau(i,j,k), fp)
               end do
            end do
            if (exttau_bb > cutoff_bb_exttau) then
               f_bb(i,j) = cutoff_bb_exttau / exttau_bb
            end if
         end do
      end do

      deallocate(q_mass, rh_r4, tau)

   end subroutine compute_bb_emission_factor

   !> \brief Apply diurnal cycle to biomass burning emissions
   !!
   !! Ported from GOCART2G Chem_BiomassDiurnal. Modulates daily-mean fire
   !! emissions using a GOES-12 derived diurnal profile (2003-2007).
   !! - NonBoreal (lat < 30): strong afternoon peak (~2x), near-zero at night
   !! - Boreal (lat >= 50): flat cycle (effectively 1.0)
   !! - Transition (30-50): linear blend
   !! Normalization ensures daily total emission is preserved.
   !!
   !! \param[inout] emission_2d  2D surface emission field [kg/m2/s], modified in place
   !! \param[in]    lons         2D longitude array [degrees]
   !! \param[in]    lats         2D latitude array [degrees]
   !! \param[in]    current_time ESMF_Time for current model time
   !! \param[in]    nx, ny       Grid dimensions
   !! \param[out]   rc           Return code
   subroutine apply_biomass_diurnal(emission_2d, lons, lats, current_time, nx, ny, rc)
      implicit none

      real(fp), intent(inout) :: emission_2d(:,:)
      real(fp), intent(in)    :: lons(:,:)
      real(fp), intent(in)    :: lats(:,:)
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(in)     :: nx, ny
      integer, intent(out)    :: rc

      ! Parameters: N=240 time bins per day, DT=360 seconds
      integer, parameter :: N = 240
      real(fp), parameter :: DT_DIURNAL = 86400.0_fp / N

      ! Boreal: flat diurnal cycle (no modulation for lat >= 50)
      real(fp), parameter :: Boreal(N) = 1.0_fp

      ! NonBoreal: GOES-12 derived diurnal profile (2003-2007)
      real(fp), parameter :: NonBoreal(N) = (/ &
         0.0121_fp, 0.0150_fp, 0.0172_fp, 0.0185_fp, 0.0189_fp, 0.0184_fp, &
         0.0174_fp, 0.0162_fp, 0.0151_fp, 0.0141_fp, 0.0133_fp, 0.0126_fp, &
         0.0121_fp, 0.0117_fp, 0.0115_fp, 0.0114_fp, 0.0114_fp, 0.0116_fp, &
         0.0120_fp, 0.0126_fp, 0.0133_fp, 0.0142_fp, 0.0151_fp, 0.0159_fp, &
         0.0167_fp, 0.0174_fp, 0.0180_fp, 0.0184_fp, 0.0187_fp, 0.0189_fp, &
         0.0190_fp, 0.0190_fp, 0.0191_fp, 0.0192_fp, 0.0192_fp, 0.0193_fp, &
         0.0194_fp, 0.0194_fp, 0.0193_fp, 0.0192_fp, 0.0190_fp, 0.0187_fp, &
         0.0185_fp, 0.0182_fp, 0.0180_fp, 0.0178_fp, 0.0177_fp, 0.0176_fp, &
         0.0174_fp, 0.0172_fp, 0.0169_fp, 0.0166_fp, 0.0162_fp, 0.0158_fp, &
         0.0153_fp, 0.0149_fp, 0.0144_fp, 0.0138_fp, 0.0132_fp, 0.0126_fp, &
         0.0118_fp, 0.0109_fp, 0.0101_fp, 0.0092_fp, 0.0085_fp, 0.0081_fp, &
         0.0080_fp, 0.0083_fp, 0.0091_fp, 0.0102_fp, 0.0117_fp, 0.0135_fp, &
         0.0157_fp, 0.0182_fp, 0.0210_fp, 0.0240_fp, 0.0273_fp, 0.0308_fp, &
         0.0345_fp, 0.0387_fp, 0.0432_fp, 0.0483_fp, 0.0540_fp, 0.0606_fp, &
         0.0683_fp, 0.0775_fp, 0.0886_fp, 0.1022_fp, 0.1188_fp, 0.1388_fp, &
         0.1625_fp, 0.1905_fp, 0.2229_fp, 0.2602_fp, 0.3025_fp, 0.3500_fp, &
         0.4031_fp, 0.4623_fp, 0.5283_fp, 0.6016_fp, 0.6824_fp, 0.7705_fp, &
         0.8650_fp, 0.9646_fp, 1.0676_fp, 1.1713_fp, 1.2722_fp, 1.3662_fp, &
         1.4491_fp, 1.5174_fp, 1.5685_fp, 1.6014_fp, 1.6173_fp, 1.6200_fp, &
         1.6150_fp, 1.6082_fp, 1.6040_fp, 1.6058_fp, 1.6157_fp, 1.6353_fp, &
         1.6651_fp, 1.7045_fp, 1.7513_fp, 1.8024_fp, 1.8541_fp, 1.9022_fp, &
         1.9429_fp, 1.9738_fp, 1.9947_fp, 2.0072_fp, 2.0132_fp, 2.0141_fp, &
         2.0096_fp, 1.9994_fp, 1.9829_fp, 1.9604_fp, 1.9321_fp, 1.8977_fp, &
         1.8562_fp, 1.8052_fp, 1.7419_fp, 1.6646_fp, 1.5738_fp, 1.4734_fp, &
         1.3693_fp, 1.2676_fp, 1.1724_fp, 1.0851_fp, 1.0052_fp, 0.9317_fp, &
         0.8637_fp, 0.8004_fp, 0.7414_fp, 0.6862_fp, 0.6348_fp, 0.5871_fp, &
         0.5434_fp, 0.5037_fp, 0.4682_fp, 0.4368_fp, 0.4097_fp, 0.3864_fp, &
         0.3667_fp, 0.3499_fp, 0.3355_fp, 0.3231_fp, 0.3123_fp, 0.3029_fp, &
         0.2944_fp, 0.2862_fp, 0.2773_fp, 0.2670_fp, 0.2547_fp, 0.2402_fp, &
         0.2238_fp, 0.2061_fp, 0.1882_fp, 0.1712_fp, 0.1562_fp, 0.1434_fp, &
         0.1332_fp, 0.1251_fp, 0.1189_fp, 0.1141_fp, 0.1103_fp, 0.1071_fp, &
         0.1043_fp, 0.1018_fp, 0.0996_fp, 0.0979_fp, 0.0968_fp, 0.0964_fp, &
         0.0966_fp, 0.0970_fp, 0.0973_fp, 0.0970_fp, 0.0959_fp, 0.0938_fp, &
         0.0909_fp, 0.0873_fp, 0.0831_fp, 0.0784_fp, 0.0732_fp, 0.0676_fp, &
         0.0618_fp, 0.0565_fp, 0.0521_fp, 0.0491_fp, 0.0475_fp, 0.0473_fp, &
         0.0480_fp, 0.0492_fp, 0.0504_fp, 0.0514_fp, 0.0519_fp, 0.0521_fp, &
         0.0520_fp, 0.0517_fp, 0.0513_fp, 0.0510_fp, 0.0507_fp, 0.0507_fp, &
         0.0508_fp, 0.0512_fp, 0.0515_fp, 0.0518_fp, 0.0519_fp, 0.0518_fp, &
         0.0513_fp, 0.0506_fp, 0.0496_fp, 0.0482_fp, 0.0465_fp, 0.0443_fp, &
         0.0418_fp, 0.0387_fp, 0.0351_fp, 0.0310_fp, 0.0263_fp, 0.0214_fp /)

      ! Local variables
      integer :: i, j, k, localrc, hh, mm, ss, ndt, NN, kk
      real(fp) :: secs, secs_local, aBoreal, aNonBoreal, alpha
      real(fp) :: fBoreal, fNonBoreal
      integer :: nhms
      character(len=*), parameter :: pName = 'apply_biomass_diurnal'

      rc = CC_SUCCESS

      ! Get HHMMSS from current ESMF time
      call ESMF_TimeGet(current_time, h=hh, m=mm, s=ss, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      nhms = hh * 10000 + mm * 100 + ss

      ! Compute normalization factors (depend on model timestep via ndt=1 for 360s bins)
      ! Use ndt=1 since we sample one bin per call (consistent with GOCART default)
      ndt = 1
      fBoreal = 0.0_fp
      fNonBoreal = 0.0_fp
      NN = 0
      do kk = 1, N, ndt
         NN = NN + 1
         fBoreal    = fBoreal    + Boreal(kk)
         fNonBoreal = fNonBoreal + NonBoreal(kk)
      end do
      fBoreal    = fBoreal / real(NN, fp)
      fNonBoreal = fNonBoreal / real(NN, fp)

      ! Find number of seconds since beginning of the day (GMT)
      secs = 3600.0_fp * hh + 60.0_fp * mm + ss

      ! Apply diurnal factors depending on latitude
      do j = 1, ny
         do i = 1, nx
            if (emission_2d(i,j) == 0.0_fp) cycle

            ! Find corresponding index in diurnal cycle array
            ! 240 = 24*60*60 / 360 (seconds per degree of longitude)
            secs_local = secs + 240.0_fp * lons(i,j)
            k = 1 + mod(nint(secs_local / DT_DIURNAL), N)
            if (k < 1) k = N + k

            ! Compute scaling factors normalized to preserve daily mean
            aBoreal    = Boreal(k) / fBoreal
            aNonBoreal = NonBoreal(k) / fNonBoreal

            ! Apply based on latitude band
            if (lats(i,j) >= 50.0_fp) then
               emission_2d(i,j) = aBoreal * emission_2d(i,j)
            else if (lats(i,j) >= 30.0_fp) then
               alpha = (lats(i,j) - 30.0_fp) / 20.0_fp
               emission_2d(i,j) = (1.0_fp - alpha) * aNonBoreal * emission_2d(i,j) + &
                  alpha * aBoreal * emission_2d(i,j)
            else
               emission_2d(i,j) = aNonBoreal * emission_2d(i,j)
            end if
         end do
      end do

   end subroutine apply_biomass_diurnal

   !! Applies emission data from ExtEmisDataType to the chemical state
   !! using species mapping from emission configuration. Processes entire
   !! arrays at once for efficiency and handles proper unit conversion.
   !!
   !! \param[in] ext_emis_data External emission data container
   !! \param[in] config_manager Configuration manager with emission mapping
   !! \param[inout] chem_state Chemical state to apply emissions to
   !! \param[in] met_state Meteorological state for unit conversion
   !! \param[in] dt Time step [s]
   !! \param[out] rc Return code
   subroutine catchem_emis_apply(category, icat, global_scale, config_manager, error_manager, chem_state, met_state, dt, current_time, rc)
      use Constants, only: g0, AIRMW  ! Gravitational acceleration and air molecular weight
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      integer, intent(in) :: icat !category index in the ext_emis_data
      real(fp), intent(in) :: global_scale
      type(ConfigManagerType), intent(in) :: config_manager
      type(ErrorManagerType), pointer, intent(inout) :: error_manager
      type(ChemStateType), intent(inout) :: chem_state
      type(MetStateType), intent(inout) :: met_state
      real(fp), intent(in) :: dt
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, ifield, ispec, n_mapped_species, species_idx
      integer :: nx, ny, nz, n_species, i, j, k
      character(len=EMIS_MAXSTR) :: msg, field_name, category_name
      character(len=64) :: mapped_species_name  ! Single species name
      real(fp) :: scale_factor  ! Single scale factor
      integer :: species_index  ! Single species index in chem_state
      character(len=*), parameter :: pName = 'catchem_emis_apply'

      ! Arrays for full domain processing
      real(fp), allocatable :: concentrations(:,:,:,:)  ! (nx,ny,nz,n_species)
      real(fp), allocatable :: emission_flux(:,:,:)       ! (nx,ny,nz) - emission rate [kg/m2/s]
      real(fp), allocatable :: species_tendency(:,:,:)  ! (nx,ny,nz) - species tendency [mol/mol/s]
      real(fp), allocatable :: f_bb(:,:)                 ! (nx,ny) - BB emission scaling factor
      real(fp) :: converter

      rc = CC_SUCCESS

      ! Point/volcanic categories inject directly into the 3D column at their
      ! mapped grid cells and plume altitude; dispatch to the dedicated handler.
      if (is_point_category(category)) then
         call catchem_emis_apply_points(category, icat, global_scale, config_manager, &
            chem_state, met_state, dt, rc)
         return
      end if

      ! Get dimensions
      nx = size(met_state%DELP, 1)
      ny = size(met_state%DELP, 2)
      nz = size(met_state%DELP, 3)
      n_species = chem_state%nSpecies

      ! Get current concentrations for all species
      allocate(concentrations(nx, ny, nz, n_species))
      call chem_state%get_all_concentrations(concentrations, localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to get concentrations from chem_state'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         deallocate(concentrations)
         return
      end if

      ! Allocate working arrays
      allocate(emission_flux(nx, ny, nz))
      allocate(species_tendency(nx, ny, nz))

      ! Get category name
      category_name = trim(category%category_name)

      ! Loop through all fields in this category
      do ifield = 1, category%n_fields
         if (.not. category%fields(ifield)%is_loaded .or. .not. allocated(category%fields(ifield)%emission_data)) cycle

         field_name = trim(category%fields(ifield)%field_name)

         ! Get emission data for entire domain [kg/m2/s]
         ! Assuming surface emissions (k=1, t=1) for now
         emission_flux(:,:,:) = category%fields(ifield)%emission_data(:,:,:,1)

         ! Apply category and global scaling factors
         emission_flux = emission_flux * category%global_scale * global_scale

         ! Apply diurnal biomass burning cycle if enabled (before vertical distribution)
         if (category%diurnal_bb) then
            call apply_biomass_diurnal(emission_flux(:,:,1), met_state%LON, met_state%LAT, &
               current_time, nx, ny, localrc)
         end if

         ! Apply vertical distribution if configured (redistributes 2D surface emission to 3D)
         if (trim(category%vertical_dist) /= 'none' .and. trim(category%vertical_dist) /= '') then
            call distribute_emissions_vertical(emission_flux, met_state, category%vertical_dist, nx, ny, nz)
         end if

         ! Direct mapping access using same indices (one-to-one correspondence)
         ! Add sanity checks to ensure category and field names match
         if (icat > config_manager%config_data%emission_mapping%n_categories) then
            write(msg, '(A,A,I0,A,I0)') trim(pName), ': Category index out of bounds: ', &
               icat, ' > ', config_manager%config_data%emission_mapping%n_categories
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
            cycle
         end if

         if (ifield > config_manager%config_data%emission_mapping%categories(icat)%n_emission_species) then
            write(msg, '(A,A,I0,A,I0)') trim(pName), ': Field index out of bounds: ', &
               ifield, ' > ', config_manager%config_data%emission_mapping%categories(icat)%n_emission_species
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
            cycle
         end if

         ! Sanity check: verify category names match
         if (trim(category_name) /= trim(config_manager%config_data%emission_mapping%categories(icat)%category_name)) then
            write(msg, '(A,A,A,A,A)') trim(pName), ': Category name mismatch: ', &
               trim(category_name), ' != ', &
               trim(config_manager%config_data%emission_mapping%categories(icat)%category_name)
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
            cycle
         end if

         ! Sanity check: verify field names match
         if (trim(field_name) /= trim(config_manager%config_data%emission_mapping%categories(icat)%species_mappings(ifield)%emission_field)) then
            write(msg, '(A,A,A,A,A)') trim(pName), ': Field name mismatch: ', &
               trim(field_name), ' != ', &
               trim(config_manager%config_data%emission_mapping%categories(icat)%species_mappings(ifield)%emission_field)
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
            cycle
         end if

         ! Direct access to species mapping data (no search needed)
         n_mapped_species = config_manager%config_data%emission_mapping%categories(icat)%species_mappings(ifield)%n_mappings

         ! Apply emissions to each mapped species
         do ispec = 1, n_mapped_species
            ! Get mapping data directly for this species
            mapped_species_name = config_manager%config_data%emission_mapping%categories(icat)%species_mappings(ifield)%map(ispec)
            scale_factor = config_manager%config_data%emission_mapping%categories(icat)%species_mappings(ifield)%scale(ispec)
            species_index = config_manager%config_data%emission_mapping%categories(icat)%species_mappings(ifield)%index(ispec)

            if (len_trim(mapped_species_name) == 0) cycle

            ! Get species index from mapping (or lookup if fallback was used)
            species_idx = species_index
            if (species_idx <= 0) then
               !check if this is to map to metstate variable since we read in some met variables from emissin reading too.
               !In the emission map yaml file, if the mapped_species_name starts with "MET_" or "met_", we will treat it as a met variable
               !and set the met state instead of chem state. The rest of the name after "MET_" should match the field name in met state.
               if (len_trim(mapped_species_name) > 4 .and. (trim(mapped_species_name(1:4)) == 'MET_' .or. trim(mapped_species_name(1:4)) == 'met_')) then
                  ! This is a mapping to a meteorological variable, not a chemical species. Skip applying to chem_state.
                  if (category%is_2d) then
                     call met_state%set_field(trim(mapped_species_name(5:)), emission_flux(:,:,1) * scale_factor, error_manager, localrc)
                  else
                     call met_state%set_field(trim(mapped_species_name(5:)), emission_flux * scale_factor, error_manager, localrc)
                  end if
                  if (localrc /= CC_SUCCESS) then
                     write(msg, '(A,A)') trim(pName), ': Failed to set met_state'
                     call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
                     rc = CC_FAILURE
                  end if
                  cycle !do not move to chemstate below
               end if
               ! Fallback case - need to lookup species index
               species_idx = chem_state%find_species(trim(mapped_species_name))
               if (species_idx <= 0) then
                  write(msg, '(A,A,A)') trim(pName), ': Species not found in chem_state: ', &
                     trim(mapped_species_name)
                  call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
                  cycle
               end if
            end if

            ! Unit conversion factor: emission flux (kg/m2/s) -> mass mixing ratio (kg/kg) -> model concentration units
            !   Step 1 (in loop below): kg/m2/s * dt[s] * g0[m/s2] / DELP[Pa] = kg/kg  (mass mixing ratio)
            !   Step 2 (here):           kg/kg * converter = final model units
            ! For gas species:    converter = AIRMW/MW_species * 1e6 => kg/kg -> ppmv (parts per million by volume)
            ! For aerosol species: converter = 1e9                  => kg/kg -> ug/kg (micrograms per kilogram)
            if (chem_state%ChemSpecies(species_idx)%is_gas) then
               converter = AIRMW / chem_state%ChemSpecies(species_idx)%mw_g * 1.0e6_fp
            else
               converter = 1.0e9_fp
            end if
            species_tendency = 0.0_fp

            do j = 1, ny
               do i = 1, nx
                  do k = 1, nz
                     if (emission_flux(i,j,k) > 0.0_fp) then
                        select case (trim(category%fields(ifield)%units))
                         case('nmol/l', 'nmol/L', 'NMOL/L')
                           ! Special case for DMS read in with nmol/L unit (Note: this is in water)
                           species_tendency(i,j,k) = emission_flux(i,j,k) * scale_factor
                         case ('1/cm3', '1/cm^3', '#/cm3', 'molec/cm3')
                           ! Special case for GMI oxidants OH which is in #/cm3 in the file (TODO:make sure the input file unit).
                           ! convert from #/cm3 to ppm to keep consistent with other species units
                           species_tendency(i,j,k) = emission_flux(i,j,k) * scale_factor / AVO * AIRMW / met_state%AIRDEN(i,j,k) * 1.e3
                         case ('mol/mol', 'MOL/MOL')
                           ! GMI NO3 and H2O2 are in mol/mol volume mixing ratio. Change to ppm
                           species_tendency(i,j,k) = emission_flux(i,j,k) * scale_factor * 1.e6_fp
                         case ('kg/m2/s', 'KG/M2/S')
                           ! Unit chain: [kg/m2/s] * scale * dt[s] * g0[m/s2] / DELP[Pa] * converter
                           !           = [kg/m2/s] * [s] * [m/s2] / [kg/m/s2 / m2] * converter
                           !           = [kg/kg] * converter
                           !           = [ug/kg] for aerosols (converter=1e9)
                           !           = [ppmv]  for gases    (converter=AIRMW/MW*1e6)

                           !safety check following GOCART
                           if (1.01_fp * emission_flux(i,j,k) / category%global_scale / global_scale > EMIS_ACCEPT) cycle
                           species_tendency(i,j,k) = emission_flux(i,j,k) * scale_factor *dt * g0 / met_state%DELP(i,j,k) * converter
                         case default
                           write(msg, '(A,A,A)') trim(pName), ': Unrecognized emission field units: ', &
                              trim(category%fields(ifield)%units)
                           call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
                        end select
                     end if
                  end do
               end do
            end do

            ! Add tendency to concentrations
            ! Apply Mie-based BB emission scaling factor if enabled
            ! Only for OC and BrC species (matching GOCART: prefix=='OC' or 'BR')
            if (category%use_oc_fbb .and. &
               .not. chem_state%ChemSpecies(species_idx)%is_gas .and. &
               (mapped_species_name(1:2) == 'oc' .or. mapped_species_name(1:2) == 'OC' .or. &
               mapped_species_name(1:2) == 'br' .or. mapped_species_name(1:2) == 'BR')) then
               if (.not. allocated(f_bb)) allocate(f_bb(nx, ny))
               call compute_bb_emission_factor(emission_flux, scale_factor, dt, &
                  met_state, chem_state, species_idx, f_bb, localrc)
               if (localrc == CC_SUCCESS) then
                  do k = 1, nz
                     species_tendency(:,:,k) = species_tendency(:,:,k) * f_bb(:,:)
                  end do
               end if
            end if
            ! Apply tendency: 'add' accumulates, 'replace' overwrites concentration
            if (trim(category%apply_method) == 'replace') then
               concentrations(:,:,:,species_idx) = species_tendency(:,:,:)
            else
               concentrations(:,:,:,species_idx) = concentrations(:,:,:,species_idx) + species_tendency(:,:,:)
            end if

         end do !end of mapped species loop

      end do ! end of field loop


      ! Set updated concentrations back to chemical state
      call chem_state%set_all_concentrations(concentrations, localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to set concentrations in chem_state'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
      end if

      ! Clean up
      deallocate(concentrations, emission_flux, species_tendency)
      if (allocated(f_bb)) deallocate(f_bb)
   end subroutine catchem_emis_apply

   !> \brief Lowercase a string (ASCII only)
   pure function emis_lower(str) result(low)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: low
      integer :: i, ic
      low = str
      do i = 1, len(str)
         ic = iachar(str(i:i))
         if (ic >= iachar('A') .and. ic <= iachar('Z')) low(i:i) = achar(ic + 32)
      end do
   end function emis_lower

   !> \brief Test whether an emission category is a point/volcanic source
   !!
   !! Point categories are read from an ASCII point table (.rc) rather than a
   !! gridded NetCDF file and are injected directly into the 3D column at the
   !! mapped grid cell(s).  A category is treated as a point source when its
   !! `format` is volcano / point / point_rc / volcano_rc (case-insensitive).
   pure logical function is_point_category(category) result(is_point)
      type(ExtEmisCategoryType), intent(in) :: category
      character(len=len(category%format)) :: fmt
      fmt = trim(emis_lower(category%format))
      is_point = (fmt == 'volcano'    .or. &
         fmt == 'point'      .or. &
         fmt == 'point_rc'   .or. &
         fmt == 'volcano_rc')
   end function is_point_category

   !> \brief Find the model layer (1..nz) whose edge interval contains altitude h
   !!
   !! Orientation-agnostic: each layer k spans [min(zedge(k),zedge(k+1)),
   !! max(zedge(k),zedge(k+1))].  Altitudes below the lowest edge clamp to the
   !! surface layer (k=1, UFS/FV3 convention) and altitudes above the highest
   !! edge clamp to the model top (k=nz).
   pure integer function find_point_layer(zedge, h, nz) result(ksel)
      real(fp), intent(in) :: zedge(:)   ! geopotential height at layer edges [m], size nz+1
      real(fp), intent(in) :: h          ! target altitude [m above sea level]
      integer,  intent(in) :: nz
      integer :: k
      real(fp) :: zb, zt, zmin, zmax
      ksel = 1
      zmin = min(zedge(1), zedge(nz+1))
      zmax = max(zedge(1), zedge(nz+1))
      if (h <= zmin) then
         ksel = 1                 ! at/below surface
         return
      else if (h >= zmax) then
         ksel = nz                ! at/above model top
         return
      end if
      do k = 1, nz
         zb = min(zedge(k), zedge(k+1))
         zt = max(zedge(k), zedge(k+1))
         if (h >= zb .and. h < zt) then
            ksel = k
            return
         end if
      end do
   end function find_point_layer

   !> \brief Read a point-source (.rc) emission table for a point/volcanic category
   !!
   !! Parses the ASCII `label::` ... `::` table block of a GOCART-style point
   !! emission resource file.  Each data row is
   !!   LAT  LON  EMIS  BASE_ELEVATION  TOP_ELEVATION
   !! where EMIS is the per-point source rate in the file's native mass units
   !! (kg S/s for the CARN volcanic degassing file), BASE/TOP are altitudes in
   !! metres above sea level (BASE==TOP for degassing; TOP>BASE for an explosive
   !! plume column).  All PEs read the full list; ownership is resolved later in
   !! catchem_map_points_to_grid.  The block label defaults to 'volcano' and may
   !! be overridden via the category's plume_rise key.
   subroutine catchem_emis_read_points(category, curr_time, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_Time), intent(in) :: curr_time
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, iounit, ios, npts, n, ifield
      logical :: file_exists, in_table
      character(len=EMIS_MAXSTR) :: filename, msg, line, label
      real(fp) :: vlat, vlon, vemis, vbot, vtop
      real(fp), allocatable :: tlat(:), tlon(:), temis(:), tbot(:), ttop(:)
      character(len=*), parameter :: pName = 'catchem_emis_read_points'

      rc = CC_SUCCESS

      ! Resolve filename: substitute date tokens if the template contains '%'
      if (index(trim(category%source_file), '%') > 0) then
         call resolve_filename_template(category%source_file, curr_time, filename, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      else
         filename = trim(category%source_file)
      end if

      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         call ESMF_LogWrite(trim(pName)//': point file not found (holding last data): '// &
            trim(filename), ESMF_LOGMSG_WARNING, rc=localrc)
         return
      end if

      ! Block label to read inside the .rc file (GOCART convention: 'volcano')
      label = 'volcano'

      open(newunit=iounit, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call ESMF_LogWrite(trim(pName)//': cannot open point file: '//trim(filename), &
            ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ! ---- First pass: count data rows inside the label:: ... :: block ----
      in_table = .false.
      npts = 0
      do
         read(iounit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         if (.not. in_table) then
            if (index(line, trim(label)//'::') > 0) in_table = .true.
            cycle
         else
            if (trim(line) == '::') exit
            npts = npts + 1
         end if
      end do

      if (npts <= 0) then
         call ESMF_LogWrite(trim(pName)//': no points found in block "'//trim(label)// &
            '" of '//trim(filename), ESMF_LOGMSG_WARNING, rc=localrc)
         close(iounit)
         return
      end if

      allocate(tlat(npts), tlon(npts), temis(npts), tbot(npts), ttop(npts))

      ! ---- Second pass: parse the data rows ----
      rewind(iounit)
      in_table = .false.
      n = 0
      do
         read(iounit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         if (.not. in_table) then
            if (index(line, trim(label)//'::') > 0) in_table = .true.
            cycle
         else
            if (trim(line) == '::') exit
            read(line, *, iostat=ios) vlat, vlon, vemis, vbot, vtop
            if (ios /= 0) then
               call ESMF_LogWrite(trim(pName)//': skipping malformed row: '//trim(line), &
                  ESMF_LOGMSG_WARNING, rc=localrc)
               cycle
            end if
            n = n + 1
            tlat(n)  = real(vlat,  fp)
            tlon(n)  = real(vlon,  fp)
            temis(n) = real(vemis, fp)
            tbot(n)  = real(vbot,  fp)
            ttop(n)  = real(vtop,  fp)
         end if
      end do
      close(iounit)
      npts = n

      ! ---- Store the point geometry on every field of this category ----
      ! (per-species partitioning happens at apply time via the species map scale)
      do ifield = 1, category%n_fields
         if (allocated(category%fields(ifield)%lat))   deallocate(category%fields(ifield)%lat)
         if (allocated(category%fields(ifield)%lon))   deallocate(category%fields(ifield)%lon)
         if (allocated(category%fields(ifield)%pemis)) deallocate(category%fields(ifield)%pemis)
         if (allocated(category%fields(ifield)%pbot))  deallocate(category%fields(ifield)%pbot)
         if (allocated(category%fields(ifield)%ptop))  deallocate(category%fields(ifield)%ptop)
         ! Drop any stale grid mapping so apply re-maps for the new point set
         if (allocated(category%fields(ifield)%ip))    deallocate(category%fields(ifield)%ip)
         if (allocated(category%fields(ifield)%jp))    deallocate(category%fields(ifield)%jp)

         allocate(category%fields(ifield)%lat(npts),   source=tlat(1:npts))
         allocate(category%fields(ifield)%lon(npts),   source=tlon(1:npts))
         allocate(category%fields(ifield)%pemis(npts), source=temis(1:npts))
         allocate(category%fields(ifield)%pbot(npts),  source=tbot(1:npts))
         allocate(category%fields(ifield)%ptop(npts),  source=ttop(1:npts))
         category%fields(ifield)%npts = npts
         category%fields(ifield)%is_loaded = .true.
      end do

      deallocate(tlat, tlon, temis, tbot, ttop)

      write(msg, '(A,I0,A,A)') trim(pName)//': read ', npts, ' point sources from ', trim(filename)
      call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_read_points

   !> \brief Map point (lat,lon) locations to local grid (i,j) indices
   !!
   !! General-purpose, distributed-memory point-to-grid locator (the CATChem
   !! analogue of MAPL_GetHorzIJIndex).  For each point it finds the nearest
   !! local cell centre using chordal distance on the unit sphere (robust to
   !! 0..360 vs -180..180 longitude conventions and to dateline/pole wrapping),
   !! then resolves global ownership with an ESMF VM all-reduce so that every
   !! point is owned by exactly one PE.  Points not owned by this PE return
   !! ip = jp = -1.
   !!
   !! \param[in]  plat,plon  Point coordinates [degrees], size npts
   !! \param[in]  npts       Number of points (identical on all PEs)
   !! \param[in]  gridlat,gridlon  Local cell-centre coordinates [degrees] (nx,ny)
   !! \param[out] ip,jp      Local i,j indices for owned points, else -1 (allocated here)
   !! \param[out] rc         Return code
   subroutine catchem_map_points_to_grid(plat, plon, npts, gridlat, gridlon, ip, jp, rc)
      implicit none

      real(fp), intent(in) :: plat(:), plon(:)
      integer,  intent(in) :: npts
      real(fp), intent(in) :: gridlat(:,:), gridlon(:,:)
      integer, allocatable, intent(out) :: ip(:), jp(:)
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_VM) :: vm
      integer :: localrc, i, j, it, nx, ny, localPet
      real(fp) :: dtor, px, py, pz, gx, gy, gz, d, dmin
      real(ESMF_KIND_R8), allocatable :: locmind(:), glomind(:), locpet(:), glopet(:)
      integer, allocatable :: lmi(:), lmj(:)
      real(fp), parameter :: dtol = 1.0e-9_fp
      character(len=*), parameter :: pName = 'catchem_map_points_to_grid'

      rc = CC_SUCCESS

      nx = size(gridlat, 1)
      ny = size(gridlat, 2)

      allocate(ip(max(npts,1)), jp(max(npts,1)))
      ip = -1
      jp = -1
      if (npts <= 0) return

      allocate(locmind(npts), glomind(npts), locpet(npts), glopet(npts), lmi(npts), lmj(npts))
      locmind = huge(1.0_ESMF_KIND_R8)
      lmi = -1
      lmj = -1
      dtor = acos(-1.0_fp) / 180.0_fp

      ! Local nearest cell-centre search (chordal distance on the unit sphere)
      do it = 1, npts
         px = cos(plat(it)*dtor) * cos(plon(it)*dtor)
         py = cos(plat(it)*dtor) * sin(plon(it)*dtor)
         pz = sin(plat(it)*dtor)
         dmin = huge(1.0_fp)
         do j = 1, ny
            do i = 1, nx
               gx = cos(gridlat(i,j)*dtor) * cos(gridlon(i,j)*dtor)
               gy = cos(gridlat(i,j)*dtor) * sin(gridlon(i,j)*dtor)
               gz = sin(gridlat(i,j)*dtor)
               d = (px-gx)**2 + (py-gy)**2 + (pz-gz)**2
               if (d < dmin) then
                  dmin = d
                  lmi(it) = i
                  lmj(it) = j
               end if
            end do
         end do
         locmind(it) = real(dmin, ESMF_KIND_R8)
      end do

      ! Global minimum distance per point across all PEs
      call ESMF_VMGetCurrent(vm, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call ESMF_VMGet(vm, localPet=localPet, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      call ESMF_VMAllReduce(vm, locmind, glomind, npts, ESMF_REDUCE_MIN, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! PET tiebreak: among PEs holding the global-min cell, the lowest PET owns
      do it = 1, npts
         if (locmind(it) <= glomind(it) + real(dtol, ESMF_KIND_R8)) then
            locpet(it) = real(localPet, ESMF_KIND_R8)
         else
            locpet(it) = huge(1.0_ESMF_KIND_R8)
         end if
      end do
      call ESMF_VMAllReduce(vm, locpet, glopet, npts, ESMF_REDUCE_MIN, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      do it = 1, npts
         if (locmind(it) <= glomind(it) + real(dtol, ESMF_KIND_R8) .and. &
            nint(glopet(it)) == localPet) then
            ip(it) = lmi(it)
            jp(it) = lmj(it)
         end if
      end do

      deallocate(locmind, glomind, locpet, glopet, lmi, lmj)

   end subroutine catchem_map_points_to_grid

   !> \brief Inject point/volcanic emissions into the 3D chemical state
   !!
   !! For each point owned by this PE (ip>0), converts the per-point source rate
   !! to a column-integrated flux [kg/m2/s] (rate / cell area), distributes it in
   !! the vertical, and adds the resulting tendency to each mapped species.
   !! Degassing points (TOP==BASE) deposit all mass in the layer containing the
   !! vent elevation; explosive points (TOP>BASE) spread mass over the top third
   !! of the cloud column, following GOCART2G's SUvolcanicEmissions.  The per-point
   !! rate (pemis) is taken in the file's native units [kg/s]; any mass conversion
   !! to the target species (e.g. kg S/s -> kg SO2/s, scale=2.0) is supplied through
   !! the species-map scale factor, exactly as for gridded emissions.
   subroutine catchem_emis_apply_points(category, icat, global_scale, config_manager, &
      chem_state, met_state, dt, rc)
      use Constants, only: g0, AIRMW
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      integer, intent(in) :: icat
      real(fp), intent(in) :: global_scale
      type(ConfigManagerType), intent(in) :: config_manager
      type(ChemStateType), intent(inout) :: chem_state
      type(MetStateType), intent(inout) :: met_state
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, nx, ny, nz, n_species, ifield, ispec, it, k, i, j
      integer :: npts, species_idx, ksel, n_mapped
      real(fp), allocatable :: concentrations(:,:,:,:)
      real(fp) :: area, fluxcol, hlow, hup, dzv, zb, zt, ovlp, frac
      real(fp) :: converter, scale_factor, dmr
      character(len=64) :: mapped_species_name
      character(len=EMIS_MAXSTR) :: msg
      character(len=*), parameter :: pName = 'catchem_emis_apply_points'

      rc = CC_SUCCESS

      nx = size(met_state%DELP, 1)
      ny = size(met_state%DELP, 2)
      nz = size(met_state%DELP, 3)
      n_species = chem_state%nSpecies

      ! Map points to the local grid once per read (ip is dropped on each re-read)
      do ifield = 1, category%n_fields
         if (category%fields(ifield)%npts <= 0) cycle
         if (.not. allocated(category%fields(ifield)%ip)) then
            call catchem_map_points_to_grid(category%fields(ifield)%lat, &
               category%fields(ifield)%lon, category%fields(ifield)%npts, &
               met_state%LAT, met_state%LON, category%fields(ifield)%ip, &
               category%fields(ifield)%jp, localrc)
            if (localrc /= CC_SUCCESS) then
               call ESMF_LogWrite(trim(pName)//': point-to-grid mapping failed', &
                  ESMF_LOGMSG_ERROR, rc=localrc)
               rc = CC_FAILURE
               return
            end if
         end if
      end do

      allocate(concentrations(nx, ny, nz, n_species))
      call chem_state%get_all_concentrations(concentrations, localrc)
      if (localrc /= CC_SUCCESS) then
         call ESMF_LogWrite(trim(pName)//': failed to get concentrations', &
            ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         deallocate(concentrations)
         return
      end if

      do ifield = 1, category%n_fields
         if (.not. category%fields(ifield)%is_loaded) cycle
         npts = category%fields(ifield)%npts
         if (npts <= 0) cycle

         n_mapped = config_manager%config_data%emission_mapping% &
            categories(icat)%species_mappings(ifield)%n_mappings

         do ispec = 1, n_mapped
            mapped_species_name = config_manager%config_data%emission_mapping% &
               categories(icat)%species_mappings(ifield)%map(ispec)
            scale_factor = config_manager%config_data%emission_mapping% &
               categories(icat)%species_mappings(ifield)%scale(ispec)
            species_idx = config_manager%config_data%emission_mapping% &
               categories(icat)%species_mappings(ifield)%index(ispec)

            if (len_trim(mapped_species_name) == 0) cycle
            if (species_idx <= 0) species_idx = chem_state%find_species(trim(mapped_species_name))
            if (species_idx <= 0) then
               call ESMF_LogWrite(trim(pName)//': species not found: '// &
                  trim(mapped_species_name), ESMF_LOGMSG_WARNING, rc=localrc)
               cycle
            end if

            ! kg/kg -> model units (ppmv for gases, ug/kg for aerosols)
            if (chem_state%ChemSpecies(species_idx)%is_gas) then
               converter = AIRMW / chem_state%ChemSpecies(species_idx)%mw_g * 1.0e6_fp
            else
               converter = 1.0e9_fp
            end if

            do it = 1, npts
               i = category%fields(ifield)%ip(it)
               j = category%fields(ifield)%jp(it)
               if (i < 1 .or. j < 1) cycle      ! not owned by this PE

               area = met_state%AREA_M2(i,j)
               if (area <= 1.0_fp) cycle

               ! Column-integrated flux [kg species/m2/s]: raw per-point rate
               ! divided by cell area, then category + global + species-map scaling
               ! (the map scale converts file units to the target species mass).
               fluxcol = category%fields(ifield)%pemis(it) / area * &
                  scale_factor * category%global_scale * global_scale
               if (fluxcol <= 0.0_fp) cycle

               hlow = category%fields(ifield)%pbot(it)
               hup  = category%fields(ifield)%ptop(it)

               if (hup > hlow) then
                  ! Explosive plume: emit in the top third of the cloud column
                  hlow = hup - (hup - hlow) / 3.0_fp
                  dzv  = max(hup - hlow, tiny(1.0_fp))
                  do k = 1, nz
                     zb = min(met_state%Z(i,j,k), met_state%Z(i,j,k+1))
                     zt = max(met_state%Z(i,j,k), met_state%Z(i,j,k+1))
                     ovlp = min(zt, hup) - max(zb, hlow)
                     if (ovlp <= 0.0_fp) cycle
                     frac = ovlp / dzv
                     dmr = fluxcol * frac * dt * g0 / met_state%DELP(i,j,k)
                     concentrations(i,j,k,species_idx) = &
                        concentrations(i,j,k,species_idx) + dmr * converter
                  end do
               else
                  ! Degassing: deposit all mass in the layer containing the vent
                  ksel = find_point_layer(met_state%Z(i,j,:), hlow, nz)
                  dmr = fluxcol * dt * g0 / met_state%DELP(i,j,ksel)
                  concentrations(i,j,ksel,species_idx) = &
                     concentrations(i,j,ksel,species_idx) + dmr * converter
               end if
            end do
         end do
      end do

      call chem_state%set_all_concentrations(concentrations, localrc)
      if (localrc /= CC_SUCCESS) then
         call ESMF_LogWrite(trim(pName)//': failed to set concentrations', &
            ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
      end if

      deallocate(concentrations)

   end subroutine catchem_emis_apply_points
   !!
   !! Loops through all emission categories and fields, writing diagnostic
   !! output for fields where diagnostics are enabled. Uses AQMIO for NetCDF output.
   !!
   !! \param[in] ext_emis_data External emission data container
   !! \param[inout] IO ESMF GridComp for I/O operations
   !! \param[in] grid ESMF grid for field creation
   !! \param[in] filename Output filename for diagnostics
   !! \param[out] rc Return code
   subroutine catchem_emis_write_diagnostics(ext_emis_data, time_slice, IO, grid, filename, rc)
      implicit none

      type(ExtEmisDataType), intent(in) :: ext_emis_data
      integer, intent(in) :: time_slice
      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      character(len=*), intent(in) :: filename
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, icat, ifield
      character(len=EMIS_MAXSTR) :: msg
      character(len=64) :: field_name, category_name
      character(len=128) :: description
      character(len=32) :: units
      character(len=*), parameter :: pName = 'catchem_emis_write_diagnostics'

      rc = CC_SUCCESS

      ! Check if diagnostics are enabled globally
      if (.not. ext_emis_data%diagnostic) then
         call ESMF_LogWrite(trim(pName)//': Global emission diagnostics disabled', &
            ESMF_LOGMSG_INFO, rc=localrc)
         return
      end if


      ! Loop through all emission categories
      do icat = 1, ext_emis_data%n_categories
         if (.not. ext_emis_data%categories(icat)%is_active) cycle
         if (.not. ext_emis_data%categories(icat)%diagnostic) cycle

         category_name = trim(ext_emis_data%categories(icat)%category_name)

         ! Loop through all fields in this category
         do ifield = 1, ext_emis_data%categories(icat)%n_fields
            field_name = trim(ext_emis_data%categories(icat)%fields(ifield)%field_name)

            if (.not. ext_emis_data%categories(icat)%fields(ifield)%diagnostic) cycle
            if (.not. ext_emis_data%categories(icat)%fields(ifield)%is_loaded) cycle
            if (.not. allocated(ext_emis_data%categories(icat)%fields(ifield)%emission_data)) cycle

            field_name = trim(ext_emis_data%categories(icat)%fields(ifield)%field_name)
            field_name = "emis_" // trim(category_name) // "_" // trim(field_name)  ! Prefix for diagnostics
            description = trim(ext_emis_data%categories(icat)%fields(ifield)%long_name)
            units = trim(ext_emis_data%categories(icat)%fields(ifield)%units)

            ! Write field based on whether it's gridded (2D) or not (3D)
            if (ext_emis_data%categories(icat)%gridded .and. ext_emis_data%categories(icat)%is_2d) then
               ! 2D gridded emission field
               call write_emission_field_2d(IO, grid, field_name, &
                  ext_emis_data%categories(icat)%fields(ifield)%emission_data(:,:,1,1), &
                  description, units, filename, time_slice, localrc)
            else
               ! 3D point source or vertical emission field
               call write_emission_field_3d(IO, grid, field_name, &
                  ext_emis_data%categories(icat)%fields(ifield)%emission_data(:,:,:,1), &
                  description, units, filename, time_slice, localrc)
            end if

            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A,A,A)') trim(pName), ': Failed to write emission field ', &
                  trim(field_name), ' from category ', trim(category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
               ! Continue with other fields
            else
               write(msg, '(A,A,A,A,A)') trim(pName), ': Wrote emission field ', &
                  trim(field_name), ' from category ', trim(category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)
            end if
         end do
      end do

      call ESMF_LogWrite(trim(pName)//': Emission diagnostics written to '//trim(filename), &
         ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_write_diagnostics

   !> \brief Write 2D emission field to NetCDF
   !!
   !! Helper subroutine to write 2D emission data using AQMIO.
   !!
   !! \param[inout] IO ESMF GridComp for I/O operations
   !! \param[in] grid ESMF grid for field creation
   !! \param[in] field_name Name of the emission field
   !! \param[in] emission_data 2D emission data array
   !! \param[in] description Field description for metadata
   !! \param[in] units Field units for metadata
   !! \param[in] filename Output filename
   !! \param[in] time_slice Time slice for NetCDF output
   !! \param[out] rc Return code
   subroutine write_emission_field_2d(IO, grid, field_name, emission_data, &
      description, units, filename, time_slice, rc)
      implicit none

      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: emission_data(:,:)
      character(len=*), intent(in) :: description
      character(len=*), intent(in) :: units
      character(len=*), intent(in) :: filename
      integer, intent(in) :: time_slice
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_Field) :: esmf_field
      type(ESMF_Info) :: info
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      integer :: i, j
      !character(len=*), parameter :: pName = 'write_emission_field_2d'

      rc = CC_SUCCESS

      ! Create 2D ESMF field
      esmf_field = ESMF_FieldCreate(grid, &
         name=trim(field_name), &
         typekind=ESMF_TYPEKIND_R4, &
         rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Set field metadata
      call ESMF_InfoGetFromHost(esmf_field, info, rc=rc)
      if (rc == ESMF_SUCCESS) then
         call ESMF_InfoSet(info, "units", trim(units), rc=rc)
         call ESMF_InfoSet(info, "description", trim(description), rc=rc)
      end if

      ! Get field data pointer and copy emission data
      call ESMF_FieldGet(esmf_field, farrayPtr=field_data_2d, rc=rc)
      if (rc /= ESMF_SUCCESS) then
         call ESMF_FieldDestroy(esmf_field, rc=rc)
         return
      end if

      ! Copy data (convert from fp to ESMF_KIND_R4)
      do j = 1, size(emission_data, 2)
         do i = 1, size(emission_data, 1)
            field_data_2d(i, j) = real(emission_data(i, j), ESMF_KIND_R4)
         end do
      end do

      ! Write to NetCDF using AQMIO
      call AQMIO_Write(IO, (/esmf_field/), timeSlice=time_slice, fileName=trim(filename), &
         iofmt=AQMIO_FMT_NETCDF, rc=rc)

      ! Clean up
      call ESMF_FieldDestroy(esmf_field, rc=rc)

   end subroutine write_emission_field_2d

   !> \brief Write 3D emission field to NetCDF
   !!
   !! Helper subroutine to write 3D emission data using AQMIO.
   !!
   !! \param[inout] IO ESMF GridComp for I/O operations
   !! \param[in] grid ESMF grid for field creation
   !! \param[in] field_name Name of the emission field
   !! \param[in] emission_data 3D emission data array
   !! \param[in] description Field description for metadata
   !! \param[in] units Field units for metadata
   !! \param[in] filename Output filename
   !! \param[in] time_slice Time slice for NetCDF output
   !! \param[out] rc Return code
   subroutine write_emission_field_3d(IO, grid, field_name, emission_data, &
      description, units, filename, time_slice, rc)
      implicit none

      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: emission_data(:,:,:)
      character(len=*), intent(in) :: description
      character(len=*), intent(in) :: units
      character(len=*), intent(in) :: filename
      integer, intent(in) :: time_slice
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_Field) :: esmf_field
      type(ESMF_Info) :: info
      real(ESMF_KIND_R4), pointer :: field_data_3d(:,:,:) => null()
      integer :: i, j, k
      !character(len=*), parameter :: pName = 'write_emission_field_3d'

      rc = CC_SUCCESS

      ! Create 3D ESMF field
      esmf_field = ESMF_FieldCreate(grid, &
         name=trim(field_name), &
         typekind=ESMF_TYPEKIND_R4, &
         ungriddedLBound=(/1/), &
         ungriddedUBound=(/size(emission_data, 3)/), &
         rc=rc)
      if (rc /= ESMF_SUCCESS) return

      ! Set field metadata
      call ESMF_InfoGetFromHost(esmf_field, info, rc=rc)
      if (rc == ESMF_SUCCESS) then
         call ESMF_InfoSet(info, "units", trim(units), rc=rc)
         call ESMF_InfoSet(info, "description", trim(description), rc=rc)
      end if

      ! Get field data pointer and copy emission data
      call ESMF_FieldGet(esmf_field, farrayPtr=field_data_3d, rc=rc)
      if (rc /= ESMF_SUCCESS) then
         call ESMF_FieldDestroy(esmf_field, rc=rc)
         return
      end if

      ! Copy data (convert from fp to ESMF_KIND_R4)
      do k = 1, size(emission_data, 3)
         do j = 1, size(emission_data, 2)
            do i = 1, size(emission_data, 1)
               field_data_3d(i, j, k) = real(emission_data(i, j, k), ESMF_KIND_R4)
            end do
         end do
      end do

      ! Write to NetCDF using AQMIO
      call AQMIO_Write(IO, (/esmf_field/), timeSlice=time_slice, fileName=trim(filename), &
         iofmt=AQMIO_FMT_NETCDF, rc=rc)

      ! Clean up
      call ESMF_FieldDestroy(esmf_field, rc=rc)

   end subroutine write_emission_field_3d


   !> \brief Finalize emission data and clean up resources
   !!
   !! Deallocates emission data structures and destroys ESMF objects.
   !! Should be called during model finalization.
   !!
   !! \param[inout] ext_emis_data External emission data container
   !! \param[out] rc Return code
   subroutine catchem_emis_finalize(ext_emis_data, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc
      character(len=*), parameter :: pName = 'catchem_emis_finalize'

      rc = CC_SUCCESS

      ! Clean up ExtEmisDataType
      call ext_emis_data%cleanup(localrc)
      if (localrc /= CC_SUCCESS) then
         call ESMF_LogWrite(trim(pName)//': Warning - ExtEmisDataType cleanup failed', &
            ESMF_LOGMSG_WARNING, rc=localrc)
      end if

      ! Clean up regrid route-handle cache
      call catchem_regrid_cleanup(emis_regrid_cache, rc=localrc)

      call ESMF_LogWrite(trim(pName)//': Emission data finalized', &
         ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_finalize

   !> \brief Parse emission category properties from configuration
   !!
   !! Reads additional emission category properties from the ConfigManager
   !! and applies them to the ExtEmisCategoryType.
   !!
   !! \param[inout] category ExtEmisCategoryType to populate with properties
   !! \param[in] config_manager Already loaded CATChem configuration manager
   !! \param[in] category_name Name of the category
   !! \param[out] rc Return code
   subroutine parse_emission_category(category, config_manager, category_name, rc, diag_species)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ConfigManagerType), intent(inout) :: config_manager
      character(len=*), intent(in) :: category_name
      integer, intent(out) :: rc
      character(len=64), optional, allocatable, intent(out) :: diag_species(:)  ! Array for diagnostic species

      ! Local variables
      integer :: localrc
      character(len=EMIS_MAXSTR) :: config_path
      !character(len=*), parameter :: pName = 'parse_emission_category'

      rc = CC_SUCCESS

      ! Build configuration path for this category
      write(config_path, '(A,A)') 'processes/extemis/', trim(category_name)

      ! Read all properties directly into category fields
      call config_manager%get_string(trim(config_path)//'/source_file', category%source_file, localrc, '')
      call config_manager%get_string(trim(config_path)//'/format', category%format, localrc, '')
      call config_manager%get_string(trim(config_path)//'/frequency', category%frequency, localrc, '')
      call config_manager%get_logical(trim(config_path)//'/gridded', category%gridded, localrc, .true.)
      call config_manager%get_logical(trim(config_path)//'/is_2d', category%is_2d, localrc, .true.)
      call config_manager%get_logical(trim(config_path)//'/diagnostics', category%diagnostic, localrc, .false.)
      call config_manager%get_real(trim(config_path)//'/scale_factor', category%global_scale, localrc, 1.0_fp)

      ! Read coordinate names
      call config_manager%get_string(trim(config_path)//'/lat_name', category%latname, localrc, '')
      call config_manager%get_string(trim(config_path)//'/lon_name', category%lonname, localrc, '')
      call config_manager%get_string(trim(config_path)//'/regrid_method', category%regrid_method, localrc, 'none')
      call config_manager%get_string(trim(config_path)//'/time_interpolation', category%time_interpolation, localrc, 'none')
      call config_manager%get_string(trim(config_path)//'/vertical_dist', category%vertical_dist, localrc, 'none')
      call config_manager%get_logical(trim(config_path)//'/reverse_vertical', category%reverse_vertical, localrc, .false.)

      ! Read stack parameter names (for point sources)
      call config_manager%get_string(trim(config_path)//'/stack_diameter', category%stkdmname, localrc, '')
      call config_manager%get_string(trim(config_path)//'/stack_height', category%stkhtname, localrc, '')
      call config_manager%get_string(trim(config_path)//'/stack_temperature', category%stktkname, localrc, '')
      call config_manager%get_string(trim(config_path)//'/stack_velocity', category%stkvename, localrc, '')

      ! Read topfraction and plume rise (for fire/point sources)
      call config_manager%get_real(trim(config_path)//'/topfraction', category%topfraction, localrc, -1.0_fp)
      call config_manager%get_string(trim(config_path)//'/plume_rise', category%plumerise, localrc, '')

      ! Read diagnostic species list using get_array
      call config_manager%get_array(trim(config_path)//'/diag_list', diag_species, localrc, default_values=["All"])

      ! Carbon emission factor (Mie-based BB AOT limiter)
      call config_manager%get_logical(trim(config_path)//'/use_oc_fbb', &
         category%use_oc_fbb, localrc, .false.)

      ! Diurnal biomass burning cycle (following GOCART2G Chem_BiomassDiurnal)
      call config_manager%get_logical(trim(config_path)//'/diurnal_bb', &
         category%diurnal_bb, localrc, .false.)

      ! Apply method: 'add' (default, accumulate onto concentration) or 'replace' (overwrite)
      call config_manager%get_string(trim(config_path)//'/apply_method', &
         category%apply_method, localrc, 'add')

   end subroutine parse_emission_category

   !> \brief Populate emission category in ExtEmisDataType
   !!
   !! Creates ExtEmisCategoryType and ExtEmisFieldType objects
   !! and adds them to the ExtEmisDataType structure.
   !!
   !! \param[inout] ext_emis_data External emission data container
   !! \param[in] category_mapping Emission category mapping from ConfigDataType
   !! \param[in] config_manager Already loaded CATChem configuration manager for reading additional properties
   !! \param[in] grid ESMF grid for field creation
   !! \param[out] rc Return code
   subroutine catchem_emis_populate_category(ext_emis_data, category_mapping, config_manager, nx, ny, nlev, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      type(EmissionCategoryMapping), intent(in) :: category_mapping
      type(ConfigManagerType), intent(inout) :: config_manager
      integer, intent(in) :: nx, ny, nlev
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, ispec, i_diag
      character(len=EMIS_MAXSTR) :: msg, field_name
      type(ExtEmisCategoryType) :: new_category
      type(ExtEmisFieldType) :: new_field
      character(len=64), allocatable :: diag_species_list(:)  ! Array for diagnostic species
      character(len=*), parameter :: pName = 'catchem_emis_populate_category'

      rc = CC_SUCCESS

      ! Initialize new category
      call new_category%init(category_mapping%category_name, 0, & !category_mapping%n_emission_species, &
         'Emission category: '//trim(category_mapping%category_name), localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to initialize category'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ! Set category properties from mapping
      new_category%is_active = category_mapping%is_active

      ! Parse additional properties from configuration using ConfigManager functions
      call parse_emission_category(new_category, config_manager, category_mapping%category_name, localrc, diag_species_list)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A,A)') trim(pName), ': Failed to parse category properties: ', &
            trim(category_mapping%category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
         ! Continue anyway with default properties
      end if

      ! Create emission fields from species mappings
      do ispec = 1, category_mapping%n_emission_species
         field_name = trim(category_mapping%species_mappings(ispec)%emission_field)

         ! Initialize field with default dimensions (would get from file metadata in practice)
         call new_field%init(field_name, nx, ny, nlev, 1, &  ! assuming 1 time step
            trim(category_mapping%species_mappings(ispec)%units), localrc)
         if (localrc == CC_SUCCESS) then
            new_field%long_name = trim(category_mapping%species_mappings(ispec)%long_name)

            ! Check if diagnostics should be enabled for this field
            ! Must meet all conditions: global diagnostics, category diagnostics, and field in diag_list
            if (ext_emis_data%diagnostic .and. new_category%diagnostic) then
               ! Check if field_name is in the diagnostic list array
               if ( allocated(diag_species_list)) then
                  ! if save out all species in this category
                  if (size(diag_species_list) == 1 .and. trim(diag_species_list(1)) == 'All') then
                     new_field%diagnostic = .true.
                  else
                     do i_diag = 1, size(diag_species_list)
                        if (trim(field_name) == trim(diag_species_list(i_diag))) then
                           new_field%diagnostic = .true.
                           exit
                        end if
                     end do
                  end if
               end if
            end if

            call new_category%add_field(new_field, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': Failed to add field: ', trim(field_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
            end if
         end if
      end do

      call ext_emis_data%add_category(new_category, localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to add category to ExtEmisDataType'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      write(msg, '(A,A,A)') trim(pName), ': Successfully populated category ', &
         trim(category_mapping%category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_populate_category

   !> \brief Initialize emission timing for one category
   !!
   !! Pre-loads time coordinates from the NetCDF file (when available) and
   !! sets the initial irec so the first catchem_emis_update reads the correct
   !! time slice.  Reads are driven by period_key comparison in
   !! catchem_emis_update — no ESMF alarms are used.
   !!
   !! \param[inout] category  Emission category to initialise
   !! \param[in]   clock      Model clock (provides startTime / currTime)
   !! \param[out]  rc         Return code
   subroutine catchem_emis_setup_timing(category, clock, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_Clock),          intent(in)    :: clock
      integer,                   intent(out)   :: rc

      integer            :: localrc
      integer            :: curr_month, curr_year, start_month, start_year
      type(ESMF_Time)         :: startTime, currTime
      type(ESMF_TimeInterval) :: timeInterval
      character(len=*), parameter :: pName = 'catchem_emis_setup_timing'

      rc = CC_SUCCESS

      call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! For non-template static-path files, pre-load time coordinates so that
      ! find_time_index can determine the correct initial irec even when the file
      ! contains more time slices than the arithmetic assumption (e.g. 14-month files).
      if (trim(category%frequency) /= 'static' .and. &
         index(trim(category%source_file), '%') == 0 .and. &
         len_trim(category%source_file) > 0 .and. &
         category%n_times == 0) then
         call catchem_emis_read_time_coord(trim(category%source_file), category, localrc)
         ! Non-fatal: if time coord read fails, fall through to arithmetic below
         category%last_resolved_file = trim(category%source_file)
      end if

      if (category%n_times > 0) then
         call catchem_emis_find_time_index(category, currTime, category%frequency, category%irec, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      else if (trim(category%frequency) == "monthly") then
         call ESMF_TimeGet(currTime, mm=curr_month, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
         category%irec = max(0, curr_month - 1)
      else if (trim(category%frequency) == "yearmonth") then
         call ESMF_TimeGet(currTime, yy=curr_year, mm=curr_month, rc=localrc)
         call ESMF_TimeGet(startTime, yy=start_year, mm=start_month, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
         category%irec = max(0, curr_month - 1) + (curr_year - start_year) * 12
      else if (trim(category%frequency) /= 'static') then
         ! Arithmetic fallback for non-template files whose time variable was unreadable.
         ! Computes elapsed periods since simulation start so restarts resume correctly.
         select case (trim(category%frequency))
          case ('hourly');  call ESMF_TimeIntervalSet(timeInterval, h=1,   rc=localrc)
          case ('weekly');  call ESMF_TimeIntervalSet(timeInterval, d=7,   rc=localrc)
          case default;     call ESMF_TimeIntervalSet(timeInterval, d=1,   rc=localrc)
         end select
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
         category%irec = int((currTime - startTime) / timeInterval)
      end if

   end subroutine catchem_emis_setup_timing

   !> \brief Map emission species to CATChem chemical species
   !!
   !! Maps emission field names to CATChem species names using
   !! the emission mapping configuration from ConfigDataType.
   !!
   !! \param[in] emission_mapping Emission mapping configuration from ConfigDataType
   !! \param[in] category_name Name of emission category
   !! \param[in] emis_field_name Emission field name from file
   !! \param[out] catchem_species Array of CATChem species names
   !! \param[out] scale_factors Array of scaling factors for each species
   !! \param[out] species_indices Array of chemical species indices in chem_state
   !! \param[out] n_species Number of mapped species
   !! \param[out] rc Return code
   subroutine catchem_emis_map_species(emission_mapping, category_name, emis_field_name, &
      catchem_species, scale_factors, species_indices, n_species, rc)
      implicit none

      type(EmissionMappingConfig), intent(in) :: emission_mapping
      character(len=*), intent(in) :: category_name
      character(len=*), intent(in) :: emis_field_name
      character(len=64), intent(out) :: catchem_species(:)
      real(fp), intent(out) :: scale_factors(:)
      integer, intent(out) :: species_indices(:)
      integer, intent(out) :: n_species
      integer, intent(out) :: rc

      ! Local variables
      integer ::  j, icat, ispec, localrc
      character(len=EMIS_MAXSTR) :: msg
      character(len=*), parameter :: pName = 'catchem_emis_map_species'

      rc = CC_SUCCESS
      n_species = 0
      catchem_species = ''
      scale_factors = 0.0_fp
      species_indices = 0

      ! Find the category in emission mapping
      do icat = 1, emission_mapping%n_categories
         if (trim(emission_mapping%categories(icat)%category_name) == trim(category_name)) then
            ! Find the species mapping in this category
            do ispec = 1, emission_mapping%categories(icat)%n_emission_species
               if (trim(emission_mapping%categories(icat)%species_mappings(ispec)%emission_field) == trim(emis_field_name)) then
                  ! Found the mapping - copy data
                  n_species = emission_mapping%categories(icat)%species_mappings(ispec)%n_mappings
                  do j = 1, min(n_species, size(catchem_species))
                     catchem_species(j) = emission_mapping%categories(icat)%species_mappings(ispec)%map(j)
                     scale_factors(j) = emission_mapping%categories(icat)%species_mappings(ispec)%scale(j)
                     species_indices(j) = emission_mapping%categories(icat)%species_mappings(ispec)%index(j)
                  end do
                  return
               end if
            end do
            exit  ! Found category but no matching field
         end if
      end do

      ! If we get here, no mapping was found - use fallback
      ! Note: For fallback cases, species indices will be 0 and need to be resolved later
      select case (trim(emis_field_name))
       case ('EMIS_NO', 'NO')
         n_species = 1
         catchem_species(1) = 'NO'
         scale_factors(1) = 1.0_fp
         species_indices(1) = 0  ! Will need lookup
       case ('EMIS_NO2', 'NO2')
         n_species = 1
         catchem_species(1) = 'NO2'
         scale_factors(1) = 1.0_fp
         species_indices(1) = 0  ! Will need lookup
       case ('EMIS_SO2', 'SO2')
         n_species = 1
         catchem_species(1) = 'SO2'
         scale_factors(1) = 1.0_fp
         species_indices(1) = 0  ! Will need lookup
       case ('EMIS_CO', 'CO')
         n_species = 1
         catchem_species(1) = 'CO'
         scale_factors(1) = 1.0_fp
         species_indices(1) = 0  ! Will need lookup
       case default
         ! Unknown mapping
         write(msg, '(A,A,A,A,A)') trim(pName), ': No mapping found for field: ', &
            trim(emis_field_name), ' in category: ', trim(category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
      end select

   end subroutine catchem_emis_map_species

   !> \brief Compute a calendar period key that changes each time data should be re-read
   !!
   !! Returns an integer whose value changes when the model time crosses a period
   !! boundary for the given emission frequency.  The initial value -1 (stored in
   !! last_period_key) always differs from a real key, forcing the first read.
   !!
   !! Key encoding:
   !!   hourly   -> yyyymmddhh
   !!   daily    -> yyyymmdd
   !!   weekly   -> yyyy * 1000 + week_of_year  (1-based, ISO-like)
   !!   monthly / yearmonth -> yyyymm
   !!   static   -> 0  (constant; first read triggered by last_period_key=-1)
   subroutine catchem_emis_period_key(frequency, curr_time, key, rc)
      character(len=*), intent(in)  :: frequency
      type(ESMF_Time),  intent(in)  :: curr_time
      integer,          intent(out) :: key
      integer,          intent(out) :: rc

      integer :: localrc, yy, mm, dd, hh, doy
      character(len=*), parameter :: pName = 'catchem_emis_period_key'

      rc = CC_SUCCESS
      key = 0

      call ESMF_TimeGet(curr_time, yy=yy, mm=mm, dd=dd, h=hh, &
         dayOfYear=doy, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      select case (trim(frequency))
       case ('hourly')
         key = yy*1000000 + mm*10000 + dd*100 + hh
       case ('daily')
         key = yy*10000 + mm*100 + dd
       case ('weekly')
         key = yy*1000 + (doy - 1) / 7   ! integer week-of-year
       case ('monthly', 'yearmonth')
         key = yy*100 + mm
       case ('static')
         key = 0   ! never changes; first read triggered by last_period_key = -1
       case default
         key = yy*10000 + mm*100 + dd   ! treat unknown as daily
      end select

   end subroutine catchem_emis_period_key

   !> \brief Replace all occurrences of old_str with new_str in str (in-place)
   subroutine str_replace_all(str, old_str, new_str)
      character(len=*), intent(inout) :: str
      character(len=*), intent(in)    :: old_str, new_str

      integer :: pos, olen, slen
      character(len=EMIS_MAXSTR) :: tmp

      olen = len_trim(old_str)
      if (olen == 0) return
      do
         pos = index(trim(str), trim(old_str))
         if (pos == 0) exit
         slen = len_trim(str)
         tmp = str(1:pos-1) // trim(new_str) // str(pos+olen:slen)
         str = tmp
      end do
   end subroutine str_replace_all

   !> \brief Substitute date tokens in a filename template using the current model time
   !!
   !! Supported tokens (GEOS ExtData convention):
   !! %y4 = 4-digit year, %m2 = 2-digit month, %d2 = 2-digit day,
   !! %h2 = 2-digit hour, %j3 = 3-digit Julian day-of-year.
   subroutine resolve_filename_template(template, curr_time, filename, rc)
      character(len=*), intent(in)  :: template
      type(ESMF_Time),  intent(in)  :: curr_time
      character(len=*), intent(out) :: filename
      integer,          intent(out) :: rc

      integer :: localrc
      integer :: year, month, day, hour, dayOfYear
      character(len=4) :: y4
      character(len=2) :: m2, d2, h2
      character(len=3) :: j3
      character(len=*), parameter :: pName = 'resolve_filename_template'

      rc = CC_SUCCESS

      call ESMF_TimeGet(curr_time, yy=year, mm=month, dd=day, h=hour, &
         dayOfYear=dayOfYear, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      write(y4, '(I4.4)') year
      write(m2, '(I2.2)') month
      write(d2, '(I2.2)') day
      write(h2, '(I2.2)') hour
      write(j3, '(I3.3)') dayOfYear

      filename = trim(template)
      call str_replace_all(filename, '%y4', y4)
      call str_replace_all(filename, '%m2', m2)
      call str_replace_all(filename, '%d2', d2)
      call str_replace_all(filename, '%h2', h2)
      call str_replace_all(filename, '%j3', j3)

   end subroutine resolve_filename_template

   !> \brief Read and cache the time coordinate from a NetCDF emission file
   !!
   !! Reads the 'time' variable and its CF-standard 'units' attribute, converts each
   !! time value to (yyyymmdd, seconds-of-day), and stores the result in
   !! category%tc_dates / category%tc_secs / category%n_times.
   !! Supports "days since", "hours since", "minutes since", and "seconds since" units.
   !! On any error (no time variable, unrecognised units, etc.) the routine returns
   !! silently with n_times=0, causing the caller to fall back to arithmetic irec.
   subroutine catchem_emis_read_time_coord(filename, category, rc)
      character(len=*),          intent(in)    :: filename
      type(ExtEmisCategoryType), intent(inout) :: category
      integer,                   intent(out)   :: rc

      integer :: localrc, nt
      integer, allocatable :: dates(:), secs(:)
      character(len=EMIS_MAXSTR) :: msg
      character(len=*), parameter :: pName = 'catchem_emis_read_time_coord'

      rc = CC_SUCCESS
      category%n_times = 0

      call AQMIO_ReadTimeCoord(trim(filename), nt, dates, secs, rc=localrc)
      if (localrc /= ESMF_SUCCESS) then
         write(msg, '(A,A,A)') trim(pName), ': Failed reading time coord from: ', trim(filename)
         call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_WARNING, rc=localrc)
         return
      end if

      if (nt < 1) return

      ! Populate category cache arrays
      if (allocated(category%tc_dates)) deallocate(category%tc_dates)
      if (allocated(category%tc_secs))  deallocate(category%tc_secs)
      allocate(category%tc_dates(nt), category%tc_secs(nt))
      category%tc_dates(:) = dates(:)
      category%tc_secs(:)  = secs(:)
      category%n_times = nt

      deallocate(dates, secs)

      write(msg, '(A,A,I0,A,A)') trim(pName), ': Cached ', nt, ' time slices from ', trim(filename)
      call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_read_time_coord

   !> \brief Find the NetCDF time-slice index matching the current model time
   !!
   !! For "monthly" frequency: searches by calendar month, ignoring year
   !! (correct for climatological files reused across simulation years).
   !! For all other frequencies: returns the 1-based lower-bound index
   !! (largest i where tc_dates(i) <= curr_date, or same date with tc_secs(i) <= curr_secs).
   subroutine catchem_emis_find_time_index(category, curr_time, frequency, irec, rc)
      type(ExtEmisCategoryType), intent(in)  :: category
      type(ESMF_Time),           intent(in)  :: curr_time
      character(len=*),          intent(in)  :: frequency
      integer,                   intent(out) :: irec
      integer,                   intent(out) :: rc

      integer :: localrc, i, best
      integer :: curr_yy, curr_mm, curr_dd, curr_hh, curr_mn, curr_ss
      integer :: curr_date, curr_secs, tc_date_i, tc_secs_i, slice_month
      integer :: target_year, target_month
      real(fp) :: frac_dummy
      character(len=*), parameter :: pName = 'catchem_emis_find_time_index'

      rc = CC_SUCCESS
      irec = 1
      if (category%n_times == 0) return

      call ESMF_TimeGet(curr_time, yy=curr_yy, mm=curr_mm, dd=curr_dd, &
         h=curr_hh, m=curr_mn, s=curr_ss, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (trim(frequency) == 'monthly') then
         ! For linear mid-month interpolation on a single multi-record climatology,
         ! the "current" slice (irec) is the LOWER bracketing month — read_regrid
         ! uses irec+1 (with Dec->Jan wrap) as the upper. For non-interpolated or
         ! template data, select the slice for the current calendar month.
         if (trim(category%time_interpolation) == 'linear' .and. category%n_times >= 2) then
            call catchem_emis_month_bracket(curr_time, target_year, target_month, frac_dummy, localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) return
         else
            target_month = curr_mm
         end if
         ! Match by calendar month only — year-independent climatological use
         do i = 1, category%n_times
            slice_month = mod(category%tc_dates(i) / 100, 100)
            if (slice_month == target_month) then
               irec = i
               return
            end if
         end do
         ! Month not found in file (unusual) — fall back to 1-based month index
         irec = max(1, min(target_month, category%n_times))
      else
         ! Lower-bound search: largest i whose time <= curr_time
         curr_date = curr_yy*10000 + curr_mm*100 + curr_dd
         curr_secs = curr_hh*3600  + curr_mn*60  + curr_ss
         best = 1
         do i = 1, category%n_times
            tc_date_i = category%tc_dates(i)
            tc_secs_i = category%tc_secs(i)
            if (tc_date_i < curr_date .or. &
               (tc_date_i == curr_date .and. tc_secs_i <= curr_secs)) then
               best = i
            else
               exit  ! time array is monotonically increasing
            end if
         end do
         irec = best
      end if

   end subroutine catchem_emis_find_time_index

   !> \brief Recompute temporal interpolation weights and blend cached time slices
   !!
   !! Called every timestep for categories with needs_time_blend=.true.
   !! Recomputes weights from the current clock time and blends interp_data_t1/t2
   !! into emission_data.
   subroutine catchem_emis_blend_time(category, curr_time, rc)
      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_Time),           intent(in)    :: curr_time
      integer,                   intent(out)   :: rc

      integer :: localrc, ifield
      integer :: curr_yy, curr_mm, curr_dd, curr_hh, curr_mn, curr_ss
      integer :: dim_days, nk_blend
      integer :: blo_year, blo_month
      real(fp) :: w_next, w_curr

      rc = CC_SUCCESS

      ! Get current time components
      call ESMF_TimeGet(curr_time, yy=curr_yy, mm=curr_mm, dd=curr_dd, &
         h=curr_hh, m=curr_mn, s=curr_ss, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! Compute weight based on frequency
      select case (trim(category%frequency))
       case ('monthly')
         if (category%n_times >= 2) then
            ! Single multi-record climatology: mid-month interpolation. The value
            ! passes through each monthly mean exactly at the middle of its month,
            ! so the monthly mean is preserved (GEOS/GOCART ExtData convention).
            call catchem_emis_month_bracket(curr_time, blo_year, blo_month, w_next, localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__, rcToReturn=rc)) return
         else
            ! Multi-file template (one record per file): legacy start-of-month ramp
            dim_days = days_in_month_func(curr_yy, curr_mm)
            w_next = (real(curr_dd - 1, fp) + real(curr_hh, fp)/24.0_fp + &
               real(curr_mn, fp)/1440.0_fp + real(curr_ss, fp)/86400.0_fp) / &
               real(dim_days, fp)
         end if
       case ('daily')
         w_next = (real(curr_hh, fp) + real(curr_mn, fp)/60.0_fp + &
            real(curr_ss, fp)/3600.0_fp) / 24.0_fp
       case ('hourly')
         w_next = (real(curr_mn, fp) + real(curr_ss, fp)/60.0_fp) / 60.0_fp
       case default
         w_next = 0.0_fp
      end select
      w_curr = 1.0_fp - w_next

      ! Blend stored time slices for each field
      do ifield = 1, category%n_fields
         if (.not. allocated(category%fields(ifield)%interp_data_t1) .or. &
            .not. allocated(category%fields(ifield)%interp_data_t2)) cycle

         nk_blend = size(category%fields(ifield)%interp_data_t1, 3)
         category%fields(ifield)%emission_data(:,:,1:nk_blend,1) = &
            w_curr * category%fields(ifield)%interp_data_t1(:,:,:,1) + &
            w_next * category%fields(ifield)%interp_data_t2(:,:,:,1)
      end do

   end subroutine catchem_emis_blend_time

   !> \brief Compute the mid-month interpolation bracket for monthly-mean data
   !!
   !! Monthly-mean values are treated as valid at the MIDDLE of their month
   !! (GEOS/GOCART ExtData convention). For a given model time this returns the
   !! lower bracketing month (lo_year/lo_month) and the linear weight `frac`
   !! (0..1) of the upper month, where the bracket endpoints are the midpoints of
   !! consecutive months:
   !!   - second half of current month -> bracket [current, next]
   !!   - first half  of current month -> bracket [previous, current]
   !! so emission = (1-frac)*M_lo + frac*M_up, which equals M_m exactly at the
   !! middle of month m and therefore preserves the monthly mean.
   subroutine catchem_emis_month_bracket(curr_time, lo_year, lo_month, frac, rc)
      type(ESMF_Time), intent(in)  :: curr_time
      integer,         intent(out) :: lo_year, lo_month
      real(fp),        intent(out) :: frac
      integer,         intent(out) :: rc

      integer  :: localrc, yy, mm, dd, hh, mn, ss
      integer  :: up_year, up_month, dim_curr, dim_lo, dim_up
      real(fp) :: pos, mid_curr, span
      character(len=*), parameter :: pName = 'catchem_emis_month_bracket'

      rc = CC_SUCCESS
      frac = 0.0_fp

      call ESMF_TimeGet(curr_time, yy=yy, mm=mm, dd=dd, h=hh, m=mn, s=ss, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      dim_curr = days_in_month_func(yy, mm)
      ! Elapsed days since 00:00 on the 1st (0-based, fractional)
      pos = real(dd - 1, fp) + real(hh, fp)/24.0_fp + &
         real(mn, fp)/1440.0_fp + real(ss, fp)/86400.0_fp
      mid_curr = real(dim_curr, fp) / 2.0_fp

      if (pos >= mid_curr) then
         ! Second half of the month: interpolate current -> next
         lo_year = yy;  lo_month = mm
         up_year = yy;  up_month = mm + 1
         if (up_month > 12) then
            up_month = 1;  up_year = yy + 1
         end if
         dim_lo = dim_curr
         dim_up = days_in_month_func(up_year, up_month)
         span   = real(dim_lo, fp)/2.0_fp + real(dim_up, fp)/2.0_fp
         frac   = (pos - mid_curr) / span
      else
         ! First half of the month: interpolate previous -> current
         lo_year = yy;  lo_month = mm - 1
         if (lo_month < 1) then
            lo_month = 12;  lo_year = yy - 1
         end if
         dim_lo = days_in_month_func(lo_year, lo_month)
         dim_up = dim_curr
         span   = real(dim_lo, fp)/2.0_fp + real(dim_up, fp)/2.0_fp
         ! Elapsed from mid(previous): (second half of previous) + pos in current
         frac   = (real(dim_lo, fp)/2.0_fp + pos) / span
      end if

      ! Numerical safety: clamp to [0,1]
      if (frac < 0.0_fp) frac = 0.0_fp
      if (frac > 1.0_fp) frac = 1.0_fp
   end subroutine catchem_emis_month_bracket

   !> \brief Return the number of days in a given month/year
   pure function days_in_month_func(year, month) result(ndays)
      integer, intent(in) :: year, month
      integer :: ndays
      integer, parameter :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      logical :: is_leap

      if (month < 1 .or. month > 12) then
         ndays = 30
         return
      end if
      ndays = mdays(month)
      if (month == 2) then
         is_leap = (mod(year,4)==0 .and. mod(year,100)/=0) .or. (mod(year,400)==0)
         if (is_leap) ndays = 29
      end if
   end function days_in_month_func

end module catchem_emis_mod
