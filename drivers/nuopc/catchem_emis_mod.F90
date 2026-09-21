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

module catchem_nuopc_emis_mod

   use iso_c_binding, only: c_ptr, c_char, c_int, c_double, c_null_char, c_associated, c_f_pointer
   use ESMF
   use NUOPC
   use aqmio
   use netcdf
   use catchem_regrid_mod, only: RegridCache, catchem_regrid_field, catchem_regrid_cleanup
   use catchem_bridge_precision, only: fp, is_exact_zero
   use catchem_bridge_error, only: CC_SUCCESS, CC_FAILURE
   use catchem_nuopc_emis_data_mod, only: ExtEmisDataType, ExtEmisCategoryType, ExtEmisFieldType
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none
   private

   ! Public interfaces
   public :: catchem_emis_init
   public :: catchem_emis_update
   public :: catchem_emis_finalize
   public :: catchem_emis_write_diagnostics
   public :: catchem_map_points_to_grid

   real(c_double), parameter :: g0 = 9.80665_c_double
   real(c_double), parameter :: AIRMW = 28.9644_c_double
   real(c_double), parameter :: AVO = 6.02214076e23_c_double

   !> \brief Parameters for emission handling
   integer, parameter :: EMIS_MAXSTR = 256
   integer, parameter :: EMIS_MAXFIELDS = 100
   real(fp), parameter :: EMIS_MISSING = -999.0_fp

   real(fp), parameter :: EMIS_ACCEPT = 1.e+15_fp ! Same as MAPL library "undefval"
   ! Observed ocean surface DMS rarely exceeds ~100 nmol/L; this rejects
   ! unmasked land/fill-value contamination smeared into coastal cells by
   ! conservative regridding, which otherwise feeds DMSemission unguarded.
   real(fp), parameter :: DMS_OCEAN_NMOLL_MAX = 1.e+3_fp

   interface
      integer(c_int) function catchem_config_has_emission_mapping(core_ptr) bind(C, name="catchem_config_has_emission_mapping")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      integer(c_int) function catchem_config_get_emission_category_count(core_ptr) bind(C, name="catchem_config_get_emission_category_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: core_ptr
      end function

      subroutine catchem_config_get_emission_category_name_at(core_ptr, index, name_out, max_len) &
         bind(C, name="catchem_config_get_emission_category_name_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: name_out(*)
         integer(c_int), value :: max_len
      end subroutine

      integer(c_int) function catchem_config_is_emission_category_active(core_ptr, category_name) &
         bind(C, name="catchem_config_is_emission_category_active")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
      end function

      integer(c_int) function catchem_config_get_emission_field_count(core_ptr, category_name) &
         bind(C, name="catchem_config_get_emission_field_count")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
      end function

      subroutine catchem_config_get_emission_field_name_at(core_ptr, category_name, field_idx, name_out, max_len) &
         bind(C, name="catchem_config_get_emission_field_name_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
         integer(c_int), value :: field_idx
         character(kind=c_char), intent(out) :: name_out(*)
         integer(c_int), value :: max_len
      end subroutine

      subroutine catchem_config_get_emission_field_units(core_ptr, category_name, field_name, units_out, max_len) &
         bind(C, name="catchem_config_get_emission_field_units")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
         character(kind=c_char), intent(in) :: field_name(*)
         character(kind=c_char), intent(out) :: units_out(*)
         integer(c_int), value :: max_len
      end subroutine

      integer(c_int) function catchem_config_get_emission_species_map_count(core_ptr, category_name, field_name) &
         bind(C, name="catchem_config_get_emission_species_map_count")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: category_name(*)
         character(kind=c_char), intent(in) :: field_name(*)
      end function

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
      end subroutine

      integer(c_int) function catchem_config_get_yaml_bool(core_ptr, yaml_path, default_val) &
         bind(C, name="catchem_config_get_yaml_bool")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: yaml_path(*)
         integer(c_int), value :: default_val
      end function

      real(c_double) function catchem_config_get_yaml_double(core_ptr, yaml_path, default_val) &
         bind(C, name="catchem_config_get_yaml_double")
         import :: c_ptr, c_char, c_double
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: yaml_path(*)
         real(c_double), value :: default_val
      end function

      subroutine catchem_config_get_yaml_string(core_ptr, yaml_path, val_out, max_len, default_val) &
         bind(C, name="catchem_config_get_yaml_string")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: yaml_path(*)
         character(kind=c_char), intent(out) :: val_out(*)
         integer(c_int), value :: max_len
         character(kind=c_char), intent(in) :: default_val(*)
      end subroutine

      subroutine catchem_config_find_fengsha_static_file(core_ptr, val_out, max_len) &
         bind(C, name="catchem_config_find_fengsha_static_file")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(out) :: val_out(*)
         integer(c_int), value :: max_len
      end subroutine

      integer(c_int) function catchem_config_get_yaml_list_count(core_ptr, yaml_path) &
         bind(C, name="catchem_config_get_yaml_list_count")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: yaml_path(*)
      end function

      subroutine catchem_config_get_yaml_list_at(core_ptr, yaml_path, index, val_out, max_len) &
         bind(C, name="catchem_config_get_yaml_list_at")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: core_ptr
         character(kind=c_char), intent(in) :: yaml_path(*)
         integer(c_int), value :: index
         character(kind=c_char), intent(out) :: val_out(*)
         integer(c_int), value :: max_len
      end subroutine

      type(c_ptr) function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function

      type(c_ptr) function catchem_state_get_pointer_3d(state_ptr, name) bind(C, name="catchem_state_get_pointer_3d")
         import :: c_ptr, c_char
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function

      type(c_ptr) function catchem_state_get_species_conc_pointer(state_ptr, index) &
         bind(C, name="catchem_state_get_species_conc_pointer")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      integer(c_int) function catchem_state_get_species_count(state_ptr) bind(C, name="catchem_state_get_species_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
      end function

      integer(c_int) function catchem_state_get_species_index(state_ptr, name) bind(C, name="catchem_state_get_species_index")
         import :: c_ptr, c_char, c_int
         type(c_ptr), value :: state_ptr
         character(kind=c_char), intent(in) :: name(*)
      end function

      integer(c_int) function catchem_state_is_species_gas(state_ptr, index) bind(C, name="catchem_state_is_species_gas")
         import :: c_ptr, c_int
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function

      real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: state_ptr
         integer(c_int), value :: index
      end function
   end interface

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
   subroutine catchem_emis_init(ext_emis_data, core_ptr, nx, ny, nlev, clock, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      type(c_ptr), intent(in) :: core_ptr
      integer, intent(in) :: nx, ny, nlev
      type(ESMF_Clock), intent(in) :: clock
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, icat, n_categories, active_category_index
      logical :: extemis_activate, category_active
      character(len=EMIS_MAXSTR) :: msg, category_name
      character(len=*), parameter :: pName = 'catchem_emis_init'

      ! Initialize
      rc = CC_SUCCESS

      if (.not. c_associated(core_ptr)) return

      ! Check top-level processes/extemis/activate switch
      extemis_activate = (catchem_config_get_yaml_bool(core_ptr, 'processes/extemis/activate' // c_null_char, 1_c_int) /= 0)
      if (.not. extemis_activate) then
#ifdef CATCHEM_TRACE_NUOPC
         write(*,'(A)') '[CATCHEM DEBUG] catchem_emis_init: extemis disabled by processes/extemis/activate=false'
         call flush(6)
#endif
         call ESMF_LogWrite(trim(pName)//': External emissions disabled (processes/extemis/activate=false)', &
            ESMF_LOGMSG_INFO, rc=localrc)
      end if

      ! Initialize ExtEmisDataType with 0 to allow push-back population
      call ext_emis_data%init(0, 'CATChem NUOPC Emission Data', localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to initialize ExtEmisDataType'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ext_emis_data%diagnostic = (catchem_config_get_yaml_bool(core_ptr, 'processes/extemis/global_diagnostics' // c_null_char, 1_c_int) /= 0)

      if (extemis_activate) then
         ! Check if emission mapping is loaded in C++ ConfigManager
         if (catchem_config_has_emission_mapping(core_ptr) == 0) then
            write(msg, '(A,A)') trim(pName), ': Emission mapping not loaded in C++ ConfigManager'
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         else
            n_categories = catchem_config_get_emission_category_count(core_ptr)
            do icat = 0, n_categories - 1
               call catchem_config_get_emission_category_name_at(core_ptr, icat, category_name, 64_c_int)
               call clean_c_string(category_name)
               category_active = (catchem_config_is_emission_category_active(core_ptr, trim(category_name) // c_null_char) /= 0)

               ! Purely config-driven: a category (including dust/fengsha) is
               ! populated only when the emission config marks it active and
               ! declares its source_file and fields.  There is intentionally
               ! NO fallback that force-activates dust or injects default field
               ! names -- missing configuration must fail loudly, not be
               ! silently synthesized (which risks diverging from the intended
               ! inputs).
               if (category_active) then
                  call catchem_emis_populate_category(ext_emis_data, core_ptr, category_name, nx, ny, nlev, localrc)
                  if (localrc /= CC_SUCCESS) then
                     write(msg, '(A,A,A)') trim(pName), ': Failed to populate category ', trim(category_name)
                     call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
                     rc = CC_FAILURE
                     return
                  end if

                  active_category_index = ext_emis_data%n_categories
                  call catchem_emis_setup_timing(ext_emis_data%categories(active_category_index), clock, localrc)
                  if (localrc /= CC_SUCCESS) then
                     write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: Failed timing setup for category: ', trim(category_name)
                     call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR, rc=localrc)
                     rc = CC_FAILURE
                     return
                  end if
               end if
            end do
         end if
      end if

#ifdef CATCHEM_TRACE_NUOPC
      write(*,'(A,I0)') '[CATCHEM DEBUG] catchem_emis_init: n_categories=', ext_emis_data%n_categories
      call flush(6)
#endif

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
   subroutine catchem_emis_update(ext_emis_data, core_ptr, current_time, nlev, IO, grid, dt, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      type(c_ptr), intent(in) :: core_ptr
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(in) :: nlev
      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, i, period_key
      integer :: blo_year, blo_month
      real(fp) :: bfrac
      character(len=EMIS_MAXSTR) :: msg, timeString
      character(len=*), parameter :: pName = 'catchem_emis_update'

      rc = CC_SUCCESS

      ! Skip if no emission categories were initialized (e.g. extemis disabled)
      if (ext_emis_data%n_categories == 0) return

      ! Loop through all emission categories and check if updates are needed
      do i = 1, ext_emis_data%n_categories
         if (.not. ext_emis_data%categories(i)%is_active) cycle

         ! Determine the current calendar period key for this category's frequency.
         call catchem_emis_period_key(ext_emis_data%categories(i)%frequency, &
            current_time, period_key, localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return

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
               " @ "//trim(timeString), ESMF_LOGMSG_DEBUG, rc=localrc)

            if (ext_emis_data%categories(i)%n_times == 0 .and. &
               index(trim(ext_emis_data%categories(i)%source_file), '%') == 0) then
               ext_emis_data%categories(i)%irec = ext_emis_data%categories(i)%irec + 1
            end if

            call catchem_emis_read(ext_emis_data%categories(i), ext_emis_data%regrid_cache, IO, grid, &
               nlev, current_time, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: Failed to read data for category: ', &
                  trim(ext_emis_data%categories(i)%category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
               rc = CC_FAILURE
               return
            end if

            ext_emis_data%categories(i)%last_period_key = period_key
         end if

         if (ext_emis_data%categories(i)%needs_time_blend) then
            call catchem_emis_blend_time(ext_emis_data%categories(i), current_time, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: Failed to blend time for category: ', &
                  trim(ext_emis_data%categories(i)%category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
               rc = CC_FAILURE
               return
            end if
         end if

         call catchem_emis_apply(ext_emis_data%categories(i), ext_emis_data%global_scale, core_ptr, dt, current_time, localrc)
         if (localrc /= CC_SUCCESS) then
            write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: Failed to apply emissions for category: ', &
               trim(ext_emis_data%categories(i)%category_name)
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
            rc = CC_FAILURE
            return
         end if
      end do

      call ESMF_LogWrite(trim(pName)//': Emission data updated', &
         ESMF_LOGMSG_DEBUG, rc=localrc)

   end subroutine catchem_emis_update

   !> Inspect the source variables rather than assuming that every field in a
   !! category has the category-level rank.  AQMIO files commonly combine
   !! surface fields with fields carrying a record dimension, and some files
   !! also contain native vertical fields in the same category.
   subroutine catchem_emis_detect_field_ranks(category, filename, rc)
      type(ExtEmisCategoryType), intent(inout) :: category
      character(len=*),          intent(in)    :: filename
      integer,                   intent(out)   :: rc

      integer :: nc_status, ncid, ifield, varid, ndims, uid, spatial_ndims, nlev_var
      integer, allocatable :: dimids(:)
      character(len=NF90_MAX_NAME) :: dim_name
      logical :: is_time

      rc = CC_SUCCESS
      do ifield = 1, category%n_fields
         category%fields(ifield)%is_2d = category%is_2d
         category%fields(ifield)%nlev_file = 1
      end do

      nc_status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      if (nc_status /= NF90_NOERR) return

      uid = -1
      nc_status = nf90_inquire(ncid, unlimitedDimId=uid)
      do ifield = 1, category%n_fields
         nc_status = nf90_inq_varid(ncid, trim(category%fields(ifield)%field_name), varid)
         if (nc_status /= NF90_NOERR) cycle
         nc_status = nf90_inquire_variable(ncid, varid, ndims=ndims)
         if (nc_status /= NF90_NOERR .or. ndims < 1) cycle

         allocate(dimids(ndims))
         nc_status = nf90_inquire_variable(ncid, varid, dimIds=dimids)
         if (nc_status /= NF90_NOERR) then
            deallocate(dimids)
            cycle
         end if

         spatial_ndims = ndims
         is_time = .false.
         if (uid /= -1 .and. dimids(ndims) == uid) then
            is_time = .true.
         else
            dim_name = ''
            nc_status = nf90_inquire_dimension(ncid, dimids(ndims), name=dim_name)
            if (nc_status == NF90_NOERR) then
               dim_name = emis_lower(dim_name)
               is_time = index(dim_name, 'time') > 0 .or. index(dim_name, 'month') > 0 .or. &
                  index(dim_name, 'record') > 0
            end if
         end if
         if (is_time) spatial_ndims = ndims - 1

         category%fields(ifield)%is_2d = (spatial_ndims <= 2)
         if (spatial_ndims >= 3) then
            nc_status = nf90_inquire_dimension(ncid, dimids(spatial_ndims), len=nlev_var)
            if (nc_status == NF90_NOERR) category%fields(ifield)%nlev_file = nlev_var
         end if
         deallocate(dimids)
      end do
      nc_status = nf90_close(ncid)
   end subroutine catchem_emis_detect_field_ranks

   !> Size storage to the native vertical extent of a detected 3D field.
   subroutine catchem_emis_size_field_vertical(field, nlev_default, nlev_file)
      type(ExtEmisFieldType), intent(inout) :: field
      integer,                intent(in)    :: nlev_default
      integer,                intent(out)   :: nlev_file
      integer :: nx, ny, nt

      nlev_file = nlev_default
      if (field%is_2d) return
      if (field%nlev_file > 0) nlev_file = field%nlev_file
      if (.not. allocated(field%emission_data)) return
      if (size(field%emission_data, 3) == nlev_file) return

      nx = size(field%emission_data, 1)
      ny = size(field%emission_data, 2)
      nt = size(field%emission_data, 4)
      deallocate(field%emission_data)
      allocate(field%emission_data(nx, ny, nlev_file, nt))
      field%emission_data = 0.0_fp
      field%nz = nlev_file
   end subroutine catchem_emis_size_field_vertical

   !> \brief Read emission data from files
   !!
   !! Reads emission data from NetCDF files using AQMIO module
   !! and populates the ExtEmisFieldType structures.
   !!
   !! \param[inout] ext_emis_data External emission data container
   !! \param[in] category_name Name of emission category to read
   !! \param[out] rc Return code
   subroutine catchem_emis_read(category, regrid_cache, IO, grid, nlev, curr_time, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(RegridCache), intent(inout) :: regrid_cache
      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      integer, intent(in) :: nlev
      type(ESMF_Time), intent(in) :: curr_time
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc,   ifield, nlev_f
      character(len=EMIS_MAXSTR) :: msg, filename
      character(len=64) :: category_name
      type(ESMF_Field) :: esmf_field
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      real(ESMF_KIND_R4), pointer :: field_data_3d(:,:,:) => null()
      character(len=*), parameter :: pName = 'catchem_emis_read'
      logical :: use_regrid
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

      if (is_null_filename(filename)) then
         write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: No valid source file specified for category: ', trim(category_name)
         call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      ! Check that the file exists before attempting I/O.
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: Source file not found: ', trim(filename)
         call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
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

      ! Keep the source-file rank authoritative.  The category-level is only a
      ! fallback for variables that cannot be inspected.
      call catchem_emis_detect_field_ranks(category, filename, localrc)
#ifdef CATCHEM_TRACE_NUOPC
      do ifield = 1, category%n_fields
         write(*,'(A,A,A,A,A,A,A,L1,A,I0)') '[CATCHEM DEBUG] AQMIO rank category=', trim(category_name), &
            ' file=', trim(filename), &
            ' field=', trim(category%fields(ifield)%field_name), &
            ' is_2d=', category%fields(ifield)%is_2d, &
            ' nlev_file=', category%fields(ifield)%nlev_file
      end do
      call flush(6)
#endif

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
         call catchem_emis_read_regrid(category, regrid_cache, grid, nlev, filename, curr_time, rc)
         return
      end if

      !open file (Note: although AQMIO_Read can open file in its source code, it gives zeros for some reason.
      !           So we have to open it here first.)
      call AQMIO_Open(IO, filename, iomode="read", iofmt=AQMIO_FMT_NETCDF, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

      do ifield = 1, category%n_fields
         !create field to receive data
         if (category%fields(ifield)%is_2d) then
            esmf_field = ESMF_FieldCreate(grid, name=trim(category%fields(ifield)%field_name), &
               typekind=ESMF_TYPEKIND_R4, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
         else !3D field
            call catchem_emis_size_field_vertical(category%fields(ifield), nlev, nlev_f)
            esmf_field = ESMF_FieldCreate(grid, name=trim(category%fields(ifield)%field_name), &
               typekind=ESMF_TYPEKIND_R4, ungriddedLBound=(/1/), ungriddedUBound=(/nlev_f/), rc=localrc)
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

         if (category%fields(ifield)%is_2d) then
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
               category%fields(ifield)%emission_data(:,:,:,1) = real(field_data_3d(:,:,nlev_f:1:-1), fp)  !reverse vertical level
            else
               category%fields(ifield)%emission_data(:,:,:,1) = real(field_data_3d(:,:,:), fp)
            end if
         end if

         category%fields(ifield)%is_loaded = .true.   !set to true; otherwise diagnostics will not be saved.
#ifdef CATCHEM_TRACE_NUOPC
         write(*,'(A,A,A,A,A,L1,A,L1,A,ES12.4,A,ES12.4)') '[CATCHEM DEBUG] AQMIO read category=', trim(category_name), &
            ' field=', trim(category%fields(ifield)%field_name), &
            ' emission_data=', allocated(category%fields(ifield)%emission_data), &
            ' interp_t1=', allocated(category%fields(ifield)%interp_data_t1), &
            ' min=', minval(category%fields(ifield)%emission_data), &
            ' max=', maxval(category%fields(ifield)%emission_data)
         call flush(6)
#endif

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
      call ESMF_LogWrite(msg, ESMF_LOGMSG_DEBUG, rc=localrc)

   end subroutine catchem_emis_read

   !> \brief Read emission data with runtime regridding from lat-lon to model grid
   !!
   !! Reads global lat-lon emission data and regrids it onto the model
   !! grid using ESMF bilinear regridding.  Route handles are cached in
   !! the component-owned regrid cache so weights are computed once per instance.
   subroutine catchem_emis_read_regrid(category, regrid_cache, grid, nlev, filename, curr_time, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(RegridCache),         intent(inout) :: regrid_cache
      type(ESMF_Grid),          intent(in)    :: grid
      integer,                  intent(in)    :: nlev
      character(len=*),         intent(in)    :: filename
      type(ESMF_Time),          intent(in)    :: curr_time
      integer,                  intent(out)   :: rc

      ! Local variables
      integer :: localrc, ifield, klev, nlev_f
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
               call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_DEBUG, rc=localrc)
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
            call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_DEBUG, rc=localrc)
         end if
      end if

      do ifield = 1, category%n_fields
         ! Create 2D destination field on the model grid
         esmf_field = ESMF_FieldCreate(grid, &
            name=trim(category%fields(ifield)%field_name), &
            typekind=ESMF_TYPEKIND_R4, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         if (category%fields(ifield)%is_2d) then
            ! --- 2D field ---
            call catchem_regrid_field( &
               cache     = regrid_cache, &
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
                     cache     = regrid_cache, &
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
                     cache     = regrid_cache, &
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
            call catchem_emis_size_field_vertical(category%fields(ifield), nlev, nlev_f)
            do klev = 1, nlev_f
               call catchem_regrid_field( &
                  cache     = regrid_cache, &
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
                     allocate(category%fields(ifield)%interp_data_t1(nx, ny, nlev_f, 1))
                  end if
                  category%fields(ifield)%interp_data_t1(:,:,klev,1) = real(field_data_2d(:,:), fp)

                  if (multi_file_interp) then
                     call catchem_regrid_field( &
                        cache     = regrid_cache, &
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
                        cache     = regrid_cache, &
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
                     allocate(category%fields(ifield)%interp_data_t2(nx, ny, nlev_f, 1))
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
                  category%fields(ifield)%emission_data(:,:,nlev_f:1:-1,1)
               if (do_time_interp) then
                  category%fields(ifield)%interp_data_t1(:,:,:,1) = &
                     category%fields(ifield)%interp_data_t1(:,:,nlev_f:1:-1,1)
                  category%fields(ifield)%interp_data_t2(:,:,:,1) = &
                     category%fields(ifield)%interp_data_t2(:,:,nlev_f:1:-1,1)
               end if
            end if
         end if

         category%fields(ifield)%is_loaded = .true.
#ifdef CATCHEM_TRACE_NUOPC
         write(*,'(A,A,A,A,A,L1,A,L1)') '[CATCHEM DEBUG] AQMIO regrid read category=', trim(category_name), &
            ' field=', trim(category%fields(ifield)%field_name), &
            ' emission_data=', allocated(category%fields(ifield)%emission_data), &
            ' interp_t1=', allocated(category%fields(ifield)%interp_data_t1)
         call flush(6)
#endif

         call ESMF_FieldDestroy(esmf_field, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         field_data_2d => null()
      end do

      write(msg, '(A,A,A)') trim(pName), &
         ': Successfully read & regridded emission data for category ', trim(category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_DEBUG, rc=localrc)

   end subroutine catchem_emis_read_regrid

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
   subroutine distribute_emissions_vertical(emission_flux, delp, airden, pblh, vertical_dist, nx, ny, nz)
      implicit none

      real(fp), intent(inout) :: emission_flux(:,:,:)
      real(fp), intent(in) :: delp(:,:,:)
      real(fp), intent(in) :: airden(:,:,:)
      real(fp), intent(in) :: pblh(:,:)
      character(len=*), intent(in) :: vertical_dist
      integer, intent(in) :: nx, ny, nz

      ! Local variables
      integer :: i, j, k
      real(fp) :: ps, p0, p1, z0_col, z1_col, dz, deltaz, deltap
      real(fp) :: p100, p500, pPBL, p9000, p10000, zpbl
      logical :: found100, found500, foundPBL, found9000, found10000
      real(fp) :: f_dist, emis_sfc
      real(fp) :: p_top, p_bot  ! pressure range for distribution

      ! Hardcoded aviation emission layers [m] following GOCART2G convention:
      !   aviation_layers = [LTO_bot, LTO_top/CDS_bot, CDS_top/CRS_bot, CRS_top]
      !   LTO (Landing/Take-Off):     0 -   100 m
      !   CDS (Climb/Descent):      100 -  9000 m
      !   CRS (Cruise):            9000 - 10000 m
      ! Only the CDS/CRS tops are needed numerically; the LTO bounds (0-100 m)
      ! are the p100 level computed below.
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
            if (is_exact_zero(emis_sfc)) cycle

            ! Compute surface pressure by summing all layer thicknesses
            ps = 0.0_fp
            do k = 1, nz
               ps = ps + delp(i, j, k)
            end do

            ! Find pressure at target altitudes by walking from surface (k=1) upward (k=nz)
            p0 = ps
            z0_col = 0.0_fp
            p100   = 0.0_fp;  found100   = .false.
            p500   = 0.0_fp;  found500   = .false.
            pPBL   = 0.0_fp;  foundPBL   = .false.
            p9000  = 0.0_fp;  found9000  = .false.
            p10000 = 0.0_fp;  found10000 = .false.

            do k = 1, nz
               p1 = p0 - delp(i, j, k)
               dz = delp(i, j, k) / (airden(i, j, k) * g0)
               z1_col = z0_col + dz

               if (.not. found100 .and. z0_col < 100.0_fp .and. z1_col >= 100.0_fp) then
                  deltaz = z1_col - 100.0_fp
                  deltap = deltaz * airden(i, j, k) * g0
                  p100 = p1 + deltap
                  found100 = .true.
               end if

               if (.not. found500 .and. z0_col < 500.0_fp .and. z1_col >= 500.0_fp) then
                  deltaz = z1_col - 500.0_fp
                  deltap = deltaz * airden(i, j, k) * g0
                  p500 = p1 + deltap
                  found500 = .true.
               end if

               zpbl = max(pblh(i, j), 100.0_fp)
               if (.not. foundPBL .and. z0_col < zpbl .and. z1_col >= zpbl) then
                  deltaz = z1_col - zpbl
                  deltap = deltaz * airden(i, j, k) * g0
                  pPBL = p1 + deltap
                  foundPBL = .true.
               end if

               if (.not. found9000 .and. z0_col < AVN_CDS_TOP .and. z1_col >= AVN_CDS_TOP) then
                  deltaz = z1_col - AVN_CDS_TOP
                  deltap = deltaz * airden(i, j, k) * g0
                  p9000 = p1 + deltap
                  found9000 = .true.
               end if

               if (.not. found10000 .and. z0_col < AVN_CRS_TOP .and. z1_col >= AVN_CRS_TOP) then
                  deltaz = z1_col - AVN_CRS_TOP
                  deltap = deltaz * airden(i, j, k) * g0
                  p10000 = p1 + deltap
                  found10000 = .true.
               end if

               p0 = p1
               z0_col = z1_col
            end do

            ! Fallback: if target height was never reached, use top-of-atmosphere pressure
            if (.not. found100)   p100   = p0
            if (.not. found500)   p500   = p0
            if (.not. foundPBL)   pPBL   = p0
            if (.not. found9000)  p9000  = p0
            if (.not. found10000) p10000 = p0

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
               p1 = p0 - delp(i, j, k)

               ! Compute fractional overlap of this model layer with the target pressure range
               ! p0 = pressure at layer bottom (higher pressure, lower altitude)
               ! p1 = pressure at layer top (lower pressure, higher altitude)
               f_dist = 0.0_fp

               if (p0 <= p_bot .and. p1 >= p_top) then
                  ! Layer fully within target range
                  f_dist = delp(i, j, k) / (p_bot - p_top)
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
   !! Placeholder: the Mie-based AOT check is not implemented yet, so the
   !! factor is 1 everywhere. The inputs are carried so the intended
   !! interface is visible at the call site.
   !!
   !! \param[in]  emission_flux  3D emission flux after vertical distribution [kg/m2/s]
   !! \param[in]  scale_factor   Species-specific scale factor from mapping
   !! \param[in]  dt             Time step [s]
   !! \param[out] f_bb           2D scaling factor [0..1] per column
   !! \param[out] rc             Return code
   subroutine compute_bb_emission_factor(emission_flux, scale_factor, dt, &
      f_bb, rc)
      implicit none

      real(fp), intent(in)    :: emission_flux(:,:,:)
      real(fp), intent(in)    :: scale_factor
      real(fp), intent(in)    :: dt
      real(fp), intent(out)   :: f_bb(:,:)
      integer, intent(out)    :: rc

      ! Inputs are not used by the placeholder; reference them so the
      ! signature stays stable without unused-argument warnings (the
      ! assumed-shape array via size(), scalars via associate).
      associate(unused_scale => scale_factor, unused_dt => dt); end associate
      if (size(emission_flux) >= 0) continue

      rc = CC_SUCCESS
      f_bb = 1.0_fp
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
            if (is_exact_zero(emission_2d(i,j))) cycle

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
   subroutine catchem_emis_apply(category, global_scale, core_ptr, dt, current_time, rc)
      use catchem_bridge_constants, only: g0, AIRMW
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      real(fp), intent(in) :: global_scale
      type(c_ptr), intent(in), optional :: core_ptr
      real(fp), intent(in) :: dt
      type(ESMF_Time), intent(in) :: current_time
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, ifield, ispec, n_mapped_species
      integer :: nx, ny, nz, i, j, k
      character(len=EMIS_MAXSTR) :: msg, field_name, category_name
      character(len=64) :: mapped_species_name
      real(c_double) :: scale_factor
      integer(c_int) :: species_index
      character(len=*), parameter :: pName = 'catchem_emis_apply'

      real(fp), allocatable :: emission_flux(:,:,:)
      real(c_double) :: converter, dqa
      logical :: is_gas

      type(c_ptr) :: state_ptr, c_conc, c_airden, c_pedge, c_pblh, c_lon, c_lat
      real(c_double), pointer :: f_conc(:,:,:)
      real(c_double), pointer :: f_airden(:,:,:)
      real(c_double), pointer :: f_pedge(:,:,:)
      real(c_double), pointer :: f_pblh(:,:)
      real(c_double), pointer :: f_lon(:,:)
      real(c_double), pointer :: f_lat(:,:)
      real(fp), allocatable :: f_delp(:,:,:)
      real(fp), allocatable :: f_bb(:,:)

      interface
         type(c_ptr) function catchem_core_get_state_manager(core_ptr) bind(C, name="catchem_core_get_state_manager")
            import :: c_ptr
            type(c_ptr), value :: core_ptr
         end function
         integer(c_int) function catchem_state_is_species_gas(state_ptr, index) bind(C, name="catchem_state_is_species_gas")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
            integer(c_int), value :: index
         end function
         real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: state_ptr
            integer(c_int), value :: index
         end function
         type(c_ptr) function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d")
            import :: c_ptr, c_char
            type(c_ptr), value :: state_ptr
            character(kind=c_char), intent(in) :: name(*)
         end function
         type(c_ptr) function catchem_state_get_pointer_3d(state_ptr, name) bind(C, name="catchem_state_get_pointer_3d")
            import :: c_ptr, c_char
            type(c_ptr), value :: state_ptr
            character(kind=c_char), intent(in) :: name(*)
         end function
         integer(c_int) function catchem_state_get_species_conc_pointer_checked( &
            state_ptr, index, dim1, dim2, ptr_out) bind(C, name="catchem_state_get_species_conc_pointer_checked")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
            integer(c_int), value :: index, dim1, dim2
            type(c_ptr), intent(out) :: ptr_out
         end function
         integer(c_int) function catchem_state_mark_chem_host_modified(state_ptr) &
            bind(C, name="catchem_state_mark_chem_host_modified")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
         end function
      end interface

      rc = CC_SUCCESS
      if (.not. present(core_ptr)) return
      if (.not. c_associated(core_ptr)) return

      if (is_point_category(category)) then
         call catchem_emis_apply_points(category, global_scale, core_ptr, dt, rc)
         return
      end if

      state_ptr = catchem_core_get_state_manager(core_ptr)
      if (.not. c_associated(state_ptr)) return

      ! Get dimensions from the first loaded field
      nx = 0; ny = 0; nz = 0
      do ifield = 1, category%n_fields
         if (category%fields(ifield)%is_loaded .and. allocated(category%fields(ifield)%emission_data)) then
            nx = size(category%fields(ifield)%emission_data, 1)
            ny = size(category%fields(ifield)%emission_data, 2)
            nz = size(category%fields(ifield)%emission_data, 3)
            exit
         end if
      end do

      if (nx == 0 .or. ny == 0 .or. nz == 0) return

      ! Get meteorological state pointers
      c_pedge = catchem_state_get_pointer_3d(state_ptr, 'PEDGE' // c_null_char)
      c_airden = catchem_state_get_pointer_3d(state_ptr, 'AIRDEN_DRY' // c_null_char)
      c_pblh = catchem_state_get_pointer_2d(state_ptr, 'PBLH' // c_null_char)
      c_lon = catchem_state_get_pointer_2d(state_ptr, 'LON' // c_null_char)
      c_lat = catchem_state_get_pointer_2d(state_ptr, 'LAT' // c_null_char)

      if (.not. c_associated(c_pedge) .or. .not. c_associated(c_airden)) then
         write(msg, '(A,A)') trim(pName), ': Failed to get PEDGE or AIRDEN_DRY pointers'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      call c_f_pointer(c_pedge, f_pedge, [nx, ny, nz+1])
      call c_f_pointer(c_airden, f_airden, [nx, ny, nz])
      if (c_associated(c_pblh)) call c_f_pointer(c_pblh, f_pblh, [nx, ny])
      if (c_associated(c_lon)) call c_f_pointer(c_lon, f_lon, [nx, ny])
      if (c_associated(c_lat)) call c_f_pointer(c_lat, f_lat, [nx, ny])

      allocate(emission_flux(nx, ny, nz))
      allocate(f_delp(nx, ny, nz))

      ! Compute f_delp from f_pedge (PEDGE goes from top to surface)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               f_delp(i, j, k) = abs(f_pedge(i, j, k) - f_pedge(i, j, k+1))
            end do
         end do
      end do

      category_name = trim(category%category_name)

      ! Loop through all fields in this category
      do ifield = 1, category%n_fields
         if (.not. category%fields(ifield)%is_loaded .or. .not. allocated(category%fields(ifield)%emission_data)) cycle

         field_name = trim(category%fields(ifield)%field_name)
         emission_flux(:,:,:) = category%fields(ifield)%emission_data(:,:,:,1)
         emission_flux = emission_flux * category%global_scale * global_scale

         if (category%diurnal_bb .and. c_associated(c_lon) .and. c_associated(c_lat)) then
            call apply_biomass_diurnal(emission_flux(:,:,1), real(f_lon, fp), real(f_lat, fp), current_time, nx, ny, localrc)
         end if

         if (trim(category%vertical_dist) /= 'none' .and. trim(category%vertical_dist) /= '') then
            if (c_associated(c_pblh)) then
               call distribute_emissions_vertical(emission_flux, f_delp, real(f_airden, fp), real(f_pblh, fp), &
                  category%vertical_dist, nx, ny, nz)
            else
               write(msg, '(A,A)') trim(pName), ': Missing PBLH for vertical distribution'
               call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
            end if
         end if

         n_mapped_species = catchem_config_get_emission_species_map_count(core_ptr, trim(category_name) // c_null_char, trim(field_name) // c_null_char)

         do ispec = 1, n_mapped_species
            call catchem_config_get_emission_species_map_at(core_ptr, trim(category_name) // c_null_char, trim(field_name) // c_null_char, ispec - 1, &
               mapped_species_name, 64_c_int, scale_factor, species_index)
            call clean_c_string(mapped_species_name)

            if (species_index <= 0) then
               ! Try MET alias if it maps to MET_
               if (len_trim(mapped_species_name) > 4 .and. (trim(mapped_species_name(1:4)) == 'MET_' .or. trim(mapped_species_name(1:4)) == 'met_')) then
                  ! Meteorological fields are already bound in StateManager during NUOPC import
                  cycle
               end if
               cycle
            end if

            if (catchem_state_get_species_conc_pointer_checked(state_ptr, species_index, &
               int(nx * ny, c_int), int(nz, c_int), c_conc) /= 0_c_int) then
               rc = CC_FAILURE
               return
            end if

            call c_f_pointer(c_conc, f_conc, [nx, ny, nz])

            is_gas = (catchem_state_is_species_gas(state_ptr, species_index) /= 0)

            if (is_gas) then
               converter = AIRMW / catchem_state_get_species_mw(state_ptr, species_index) * 1.0e6_c_double
            else
               converter = 1.0e9_c_double
            end if

            ! Apply BB emission factor if needed (only OC/BC aerosols)
            if (category%use_oc_fbb .and. .not. is_gas .and. &
               (mapped_species_name(1:2) == 'oc' .or. mapped_species_name(1:2) == 'OC' .or. &
               mapped_species_name(1:2) == 'br' .or. mapped_species_name(1:2) == 'BR')) then
               if (.not. allocated(f_bb)) allocate(f_bb(nx, ny))
               call compute_bb_emission_factor(emission_flux, real(scale_factor, fp), dt, f_bb, localrc)
               if (localrc == CC_SUCCESS) then
                  do k = 1, nz
                     emission_flux(:,:,k) = emission_flux(:,:,k) * f_bb(:,:)
                  end do
               end if
            end if

            ! Preserve the legacy CATChem replacement contract.  The legacy
            ! adapter constructed a zeroed full-domain tendency and assigned
            ! it wholesale for a "replace" category.  Clearing first is
            ! therefore required even when a source cell is zero, missing, or
            ! rejected (for example an invalid ocean-DMS climatology value).
            ! Without this, non-advected inputs such as dms_in retain stale
            ! values in the persistent unified chemistry buffer and continue
            ! to drive subsequent process timesteps.
            if (trim(category%apply_method) == 'replace') f_conc = 0.0_c_double

            do k = 1, nz
               do j = 1, ny
                  do i = 1, nx
                     if (emission_flux(i,j,k) > 0.0_fp) then
                        select case (trim(category%fields(ifield)%units))
                         case('nmol/l', 'nmol/L', 'NMOL/L')
                           ! Reject unmasked fill-value/land contamination before it reaches DMSemission.
                           if (.not. ieee_is_finite(emission_flux(i,j,k)) .or. &
                              emission_flux(i,j,k) > DMS_OCEAN_NMOLL_MAX) cycle
                           dqa = emission_flux(i,j,k) * scale_factor
                         case ('1/cm3', '1/cm^3', '#/cm3', 'molec/cm3')
                           dqa = emission_flux(i,j,k) * scale_factor / AVO * AIRMW / f_airden(i,j,k) * 1.e3_c_double
                         case ('mol/mol', 'MOL/MOL')
                           dqa = emission_flux(i,j,k) * scale_factor * 1.e6_c_double
                         case ('kg/m2/s', 'KG/M2/S')
                           if (1.01_fp * emission_flux(i,j,k) / category%global_scale / global_scale > EMIS_ACCEPT) cycle
                           dqa = emission_flux(i,j,k) * scale_factor * dt * g0 / f_delp(i,j,k) * converter
                         case default
                           write(msg, '(A,A,A)') trim(pName), ': Unrecognized emission field units: ', &
                              trim(category%fields(ifield)%units)
                           call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
                           dqa = 0.0_c_double
                        end select

                        if (trim(category%apply_method) == 'replace') then
                           f_conc(i,j,k) = dqa
                        else
                           f_conc(i,j,k) = f_conc(i,j,k) + dqa
                        end if
                     end if
                  end do
               end do
            end do

         end do
      end do

      if (catchem_state_mark_chem_host_modified(state_ptr) /= 0_c_int) then
         rc = CC_FAILURE
         return
      end if

      deallocate(f_delp, emission_flux)
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
      call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_DEBUG, rc=localrc)

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
   subroutine catchem_emis_apply_points(category, global_scale, core_ptr, dt, rc)
      use catchem_bridge_constants, only: g0, AIRMW
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      real(fp), intent(in) :: global_scale
      type(c_ptr), intent(in) :: core_ptr
      real(fp), intent(in) :: dt
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, nx, ny, nz, ifield, ispec, it, k, i, j
      integer :: npts, ksel, n_mapped_species
      real(c_double) :: area, fluxcol, hlow, hup, dzv, zb, zt, ovlp, frac
      real(c_double) :: converter, scale_factor, dmr
      character(len=64) :: mapped_species_name
      character(len=EMIS_MAXSTR) :: category_name, field_name
      character(len=*), parameter :: pName = 'catchem_emis_apply_points'

      integer(c_int) :: species_index
      logical :: is_gas

      type(c_ptr) :: state_ptr, c_conc, c_area, c_pedge, c_z, c_lat, c_lon
      real(c_double), pointer :: f_conc(:,:,:)
      real(c_double), pointer :: f_area(:,:)
      real(c_double), pointer :: f_pedge(:,:,:)
      real(c_double), pointer :: f_z(:,:,:)
      real(c_double), pointer :: f_lat(:,:)
      real(c_double), pointer :: f_lon(:,:)
      real(fp), allocatable :: f_delp(:,:,:)

      interface
         type(c_ptr) function catchem_core_get_state_manager(core_ptr) bind(C, name="catchem_core_get_state_manager")
            import :: c_ptr
            type(c_ptr), value :: core_ptr
         end function
         integer(c_int) function catchem_state_is_species_gas(state_ptr, index) bind(C, name="catchem_state_is_species_gas")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
            integer(c_int), value :: index
         end function
         real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: state_ptr
            integer(c_int), value :: index
         end function
         type(c_ptr) function catchem_state_get_pointer_2d(state_ptr, name) bind(C, name="catchem_state_get_pointer_2d")
            import :: c_ptr, c_char
            type(c_ptr), value :: state_ptr
            character(kind=c_char), intent(in) :: name(*)
         end function
         type(c_ptr) function catchem_state_get_pointer_3d(state_ptr, name) bind(C, name="catchem_state_get_pointer_3d")
            import :: c_ptr, c_char
            type(c_ptr), value :: state_ptr
            character(kind=c_char), intent(in) :: name(*)
         end function
         integer(c_int) function catchem_state_get_species_conc_pointer_checked( &
            state_ptr, index, dim1, dim2, ptr_out) bind(C, name="catchem_state_get_species_conc_pointer_checked")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
            integer(c_int), value :: index, dim1, dim2
            type(c_ptr), intent(out) :: ptr_out
         end function
         integer(c_int) function catchem_state_mark_chem_host_modified(state_ptr) &
            bind(C, name="catchem_state_mark_chem_host_modified")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
         end function
         integer(c_int) function catchem_state_get_nx(state_ptr) bind(C, name="catchem_state_get_nx")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
         end function
         integer(c_int) function catchem_state_get_ny(state_ptr) bind(C, name="catchem_state_get_ny")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
         end function
         integer(c_int) function catchem_state_get_nz(state_ptr) bind(C, name="catchem_state_get_nz")
            import :: c_ptr, c_int
            type(c_ptr), value :: state_ptr
         end function
      end interface

      rc = CC_SUCCESS
      if (.not. c_associated(core_ptr)) return

      state_ptr = catchem_core_get_state_manager(core_ptr)
      if (.not. c_associated(state_ptr)) return

      nx = catchem_state_get_nx(state_ptr)
      ny = catchem_state_get_ny(state_ptr)
      nz = catchem_state_get_nz(state_ptr)

      c_lat = catchem_state_get_pointer_2d(state_ptr, 'LAT' // c_null_char)
      c_lon = catchem_state_get_pointer_2d(state_ptr, 'LON' // c_null_char)
      c_area = catchem_state_get_pointer_2d(state_ptr, 'AREA_M2' // c_null_char)
      c_pedge = catchem_state_get_pointer_3d(state_ptr, 'PEDGE' // c_null_char)
      c_z = catchem_state_get_pointer_3d(state_ptr, 'Z' // c_null_char)

      if (.not. c_associated(c_lat) .or. .not. c_associated(c_lon) .or. .not. c_associated(c_area) .or. &
         .not. c_associated(c_pedge) .or. .not. c_associated(c_z)) then
         call ESMF_LogWrite(trim(pName)//': Missing pointers for points emission.', ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      call c_f_pointer(c_lat, f_lat, [nx, ny])
      call c_f_pointer(c_lon, f_lon, [nx, ny])
      call c_f_pointer(c_area, f_area, [nx, ny])
      call c_f_pointer(c_pedge, f_pedge, [nx, ny, nz+1])
      call c_f_pointer(c_z, f_z, [nx, ny, nz+1])

      allocate(f_delp(nx, ny, nz))
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               f_delp(i, j, k) = abs(f_pedge(i, j, k) - f_pedge(i, j, k+1))
            end do
         end do
      end do

      ! Map points to the local grid once per read
      do ifield = 1, category%n_fields
         if (category%fields(ifield)%npts <= 0) cycle
         if (.not. allocated(category%fields(ifield)%ip)) then
            call catchem_map_points_to_grid(category%fields(ifield)%lat, &
               category%fields(ifield)%lon, category%fields(ifield)%npts, &
               real(f_lat, fp), real(f_lon, fp), category%fields(ifield)%ip, &
               category%fields(ifield)%jp, localrc)
            if (localrc /= CC_SUCCESS) then
               call ESMF_LogWrite(trim(pName)//': point-to-grid mapping failed', ESMF_LOGMSG_ERROR, rc=localrc)
               rc = CC_FAILURE
               return
            end if
         end if
      end do

      category_name = trim(category%category_name)

      do ifield = 1, category%n_fields
         if (.not. category%fields(ifield)%is_loaded) cycle
         npts = category%fields(ifield)%npts
         if (npts <= 0) cycle

         field_name = trim(category%fields(ifield)%field_name)
         n_mapped_species = catchem_config_get_emission_species_map_count(core_ptr, trim(category_name) // c_null_char, trim(field_name) // c_null_char)

         do ispec = 1, n_mapped_species
            call catchem_config_get_emission_species_map_at(core_ptr, trim(category_name) // c_null_char, trim(field_name) // c_null_char, ispec - 1, &
               mapped_species_name, 64_c_int, scale_factor, species_index)
            call clean_c_string(mapped_species_name)

            if (species_index <= 0) cycle

            if (catchem_state_get_species_conc_pointer_checked(state_ptr, species_index, &
               int(nx * ny, c_int), int(nz, c_int), c_conc) /= 0_c_int) then
               rc = CC_FAILURE
               return
            end if

            call c_f_pointer(c_conc, f_conc, [nx, ny, nz])

            is_gas = (catchem_state_is_species_gas(state_ptr, species_index) /= 0)

            if (is_gas) then
               converter = AIRMW / catchem_state_get_species_mw(state_ptr, species_index) * 1.0e6_c_double
            else
               converter = 1.0e9_c_double
            end if

            do it = 1, npts
               i = category%fields(ifield)%ip(it)
               j = category%fields(ifield)%jp(it)
               ! The distributed locator assigns each global point to exactly
               ! one PET and deliberately stores (-1,-1) on every non-owner.
               ! Those points are not errors on this local tile.
               if (i == -1 .and. j == -1) cycle
               if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) then
                  call ESMF_LogWrite(trim(pName)//': invalid local point-source mapping', &
                     ESMF_LOGMSG_ERROR, rc=localrc)
                  rc = CC_FAILURE
                  return
               end if

               area = f_area(i,j)
               if (area <= 1.0_c_double) cycle

               fluxcol = category%fields(ifield)%pemis(it) / area * &
                  scale_factor * category%global_scale * global_scale
               if (fluxcol <= 0.0_c_double) cycle

               hlow = category%fields(ifield)%pbot(it)
               hup  = category%fields(ifield)%ptop(it)

               if (hup > hlow) then
                  ! Explosive plume
                  hlow = hup - (hup - hlow) / 3.0_c_double
                  dzv  = max(hup - hlow, tiny(1.0_c_double))
                  do k = 1, nz
                     zb = min(f_z(i,j,k), f_z(i,j,k+1))
                     zt = max(f_z(i,j,k), f_z(i,j,k+1))
                     ovlp = min(zt, hup) - max(zb, hlow)
                     if (ovlp <= 0.0_c_double) cycle
                     frac = ovlp / dzv
                     dmr = fluxcol * frac * dt * g0 / f_delp(i,j,k)
                     f_conc(i,j,k) = f_conc(i,j,k) + dmr * converter
                  end do
               else
                  ! Degassing
                  ksel = find_point_layer(real(f_z(i,j,:), fp), real(hlow, fp), nz)
                  dmr = fluxcol * dt * g0 / f_delp(i,j,ksel)
                  f_conc(i,j,ksel) = f_conc(i,j,ksel) + dmr * converter
               end if
            end do
         end do
      end do

      deallocate(f_delp)

      if (catchem_state_mark_chem_host_modified(state_ptr) /= 0_c_int) rc = CC_FAILURE

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
            ESMF_LOGMSG_DEBUG, rc=localrc)
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
            if (ext_emis_data%categories(icat)%gridded .and. &
               ext_emis_data%categories(icat)%fields(ifield)%is_2d) then
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
               call ESMF_LogWrite(msg, ESMF_LOGMSG_DEBUG, rc=localrc)
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
      call catchem_regrid_cleanup(ext_emis_data%regrid_cache, rc=localrc)

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
   subroutine parse_emission_category(category, core_ptr, category_name, rc, diag_species)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(c_ptr), intent(in) :: core_ptr
      character(len=*), intent(in) :: category_name
      integer, intent(out) :: rc
      character(len=64), optional, allocatable, intent(out) :: diag_species(:)

      integer :: i, n_diag
      character(len=EMIS_MAXSTR) :: config_path, item_path, c_buf, clean_cat_name

      rc = CC_SUCCESS

      clean_cat_name = category_name
      call clean_c_string(clean_cat_name)

      config_path = 'processes/extemis/' // trim(clean_cat_name)

      ! source_file
      item_path = trim(config_path) // '/source_file'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%source_file = trim(c_buf)

      ! format
      item_path = trim(config_path) // '/format'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%format = trim(c_buf)

      ! frequency
      item_path = trim(config_path) // '/frequency'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%frequency = trim(c_buf)

      ! gridded
      item_path = trim(config_path) // '/gridded'
      category%gridded = (catchem_config_get_yaml_bool(core_ptr, trim(item_path) // c_null_char, 1_c_int) /= 0)

      ! is_2d
      item_path = trim(config_path) // '/is_2d'
      category%is_2d = (catchem_config_get_yaml_bool(core_ptr, trim(item_path) // c_null_char, 1_c_int) /= 0)

      ! diagnostics
      item_path = trim(config_path) // '/diagnostics'
      category%diagnostic = (catchem_config_get_yaml_bool(core_ptr, trim(item_path) // c_null_char, 0_c_int) /= 0)

      ! scale_factor
      item_path = trim(config_path) // '/scale_factor'
      category%global_scale = catchem_config_get_yaml_double(core_ptr, trim(item_path) // c_null_char, 1.0_c_double)

      ! lat_name, lon_name, regrid_method, time_interpolation, vertical_dist
      item_path = trim(config_path) // '/lat_name'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%latname = trim(c_buf)

      item_path = trim(config_path) // '/lon_name'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%lonname = trim(c_buf)

      item_path = trim(config_path) // '/regrid_method'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, 'none' // c_null_char)
      call clean_c_string(c_buf)
      category%regrid_method = trim(c_buf)

      item_path = trim(config_path) // '/time_interpolation'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, 'none' // c_null_char)
      call clean_c_string(c_buf)
      category%time_interpolation = trim(c_buf)

      item_path = trim(config_path) // '/vertical_dist'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, 'none' // c_null_char)
      call clean_c_string(c_buf)
      category%vertical_dist = trim(c_buf)

      item_path = trim(config_path) // '/reverse_vertical'
      category%reverse_vertical = (catchem_config_get_yaml_bool(core_ptr, trim(item_path) // c_null_char, 0_c_int) /= 0)

      ! Stack parameters
      item_path = trim(config_path) // '/stack_diameter'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%stkdmname = trim(c_buf)

      item_path = trim(config_path) // '/stack_height'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%stkhtname = trim(c_buf)

      item_path = trim(config_path) // '/stack_temperature'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%stktkname = trim(c_buf)

      item_path = trim(config_path) // '/stack_velocity'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%stkvename = trim(c_buf)

      ! Topfraction and plume rise
      item_path = trim(config_path) // '/topfraction'
      category%topfraction = catchem_config_get_yaml_double(core_ptr, trim(item_path) // c_null_char, -1.0_c_double)

      item_path = trim(config_path) // '/plume_rise'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, '' // c_null_char)
      call clean_c_string(c_buf)
      category%plumerise = trim(c_buf)

      ! Diurnal / OC / apply_method
      item_path = trim(config_path) // '/use_oc_fbb'
      category%use_oc_fbb = (catchem_config_get_yaml_bool(core_ptr, trim(item_path) // c_null_char, 0_c_int) /= 0)

      item_path = trim(config_path) // '/diurnal_bb'
      category%diurnal_bb = (catchem_config_get_yaml_bool(core_ptr, trim(item_path) // c_null_char, 0_c_int) /= 0)

      item_path = trim(config_path) // '/apply_method'
      call catchem_config_get_yaml_string(core_ptr, trim(item_path) // c_null_char, c_buf, 256_c_int, 'add' // c_null_char)
      call clean_c_string(c_buf)
      category%apply_method = trim(c_buf)

      ! Diagnostic species list
      item_path = trim(config_path) // '/diag_list'
      n_diag = catchem_config_get_yaml_list_count(core_ptr, trim(item_path) // c_null_char)
      if (present(diag_species)) then
         if (n_diag > 0) then
            allocate(diag_species(n_diag))
            do i = 1, n_diag
               call catchem_config_get_yaml_list_at(core_ptr, trim(item_path) // c_null_char, int(i - 1, c_int), c_buf, 64_c_int)
               call clean_c_string(c_buf)
               diag_species(i) = trim(c_buf)
            end do
         else
            allocate(diag_species(1))
            diag_species(1) = "All"
         end if
      end if
   end subroutine parse_emission_category

   subroutine catchem_emis_populate_category(ext_emis_data, core_ptr, category_name, nx, ny, nlev, rc)
      implicit none

      type(ExtEmisDataType), intent(inout) :: ext_emis_data
      type(c_ptr), intent(in) :: core_ptr
      character(len=*), intent(in) :: category_name
      integer, intent(in) :: nx, ny, nlev
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, ispec, i_diag, n_fields, n_species_fields
      character(len=EMIS_MAXSTR) :: field_name, field_units, species_list_path
      type(ExtEmisCategoryType) :: new_category
      type(ExtEmisFieldType) :: new_field
      character(len=64), allocatable :: diag_species_list(:)

      rc = CC_SUCCESS

      call new_category%init(category_name, 0, 'Emission category: '//trim(category_name), localrc)
      if (localrc /= CC_SUCCESS) return

      new_category%is_active = (catchem_config_is_emission_category_active(core_ptr, trim(category_name) // c_null_char) /= 0)

      call parse_emission_category(new_category, core_ptr, category_name, localrc, diag_species_list)

      n_fields = catchem_config_get_emission_field_count(core_ptr, trim(category_name) // c_null_char)

      do ispec = 0, n_fields - 1
         call catchem_config_get_emission_field_name_at(core_ptr, trim(category_name) // c_null_char, ispec, field_name, 64_c_int)
         call clean_c_string(field_name)
         ! Field units are part of the emission mapping contract.  Most sources
         ! are mass fluxes, but some inputs (e.g., an ocean concentration) are
         ! state fields and must not inherit the mass-flux default.
         call catchem_config_get_emission_field_units(core_ptr, trim(category_name) // c_null_char, &
            trim(field_name) // c_null_char, field_units, int(EMIS_MAXSTR, c_int))
         call clean_c_string(field_units)
         if (len_trim(field_units) == 0) field_units = 'kg/m2/s'
         call new_field%init(field_name, nx, ny, nlev, 1, trim(field_units), localrc)
         if (localrc == CC_SUCCESS) then
            new_field%long_name = trim(field_name)
#ifdef CATCHEM_TRACE_NUOPC
            write(*,'(A,A,A,A)') '[CATCHEM DEBUG] AQMIO populate category=', trim(category_name), &
               ' field=', trim(field_name)
            call flush(6)
#endif
            if (ext_emis_data%diagnostic .and. new_category%diagnostic) then
               if (allocated(diag_species_list)) then
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
         end if
      end do

      species_list_path = 'processes/extemis/' // trim(category_name) // '/species'
      n_species_fields = catchem_config_get_yaml_list_count(core_ptr, trim(species_list_path) // c_null_char)
      if (n_species_fields == 0) then
         species_list_path = 'process/extemis/' // trim(category_name) // '/species'
         n_species_fields = catchem_config_get_yaml_list_count(core_ptr, trim(species_list_path) // c_null_char)
      end if

      do ispec = 0, n_species_fields - 1
         call catchem_config_get_yaml_list_at(core_ptr, trim(species_list_path) // c_null_char, &
            int(ispec, c_int), field_name, 64_c_int)
         call clean_c_string(field_name)
         if (len_trim(field_name) == 0) cycle
         if (new_category%find_field(trim(field_name)) > 0) cycle

         call catchem_config_get_emission_field_units(core_ptr, trim(category_name) // c_null_char, &
            trim(field_name) // c_null_char, field_units, int(EMIS_MAXSTR, c_int))
         call clean_c_string(field_units)
         if (len_trim(field_units) == 0) field_units = 'kg/m2/s'
         call new_field%init(field_name, nx, ny, nlev, 1, trim(field_units), localrc)
         if (localrc == CC_SUCCESS) then
            new_field%long_name = trim(field_name)
#ifdef CATCHEM_TRACE_NUOPC
            write(*,'(A,A,A,A)') '[CATCHEM DEBUG] AQMIO populate species-list category=', trim(category_name), &
               ' field=', trim(field_name)
            call flush(6)
#endif
            call new_category%add_field(new_field, localrc)
         end if
      end do

      call ext_emis_data%add_category(new_category, localrc)
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

      rc = CC_SUCCESS

      call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ! For non-template static-path files, pre-load time coordinates so that
      ! find_time_index can determine the correct initial irec even when the file
      ! contains more time slices than the arithmetic assumption (e.g. 14-month files).
      if (trim(category%frequency) /= 'static' .and. &
         index(trim(category%source_file), '%') == 0 .and. &
         .not. is_null_filename(category%source_file) .and. &
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
      logical :: file_exists
      character(len=*), parameter :: pName = 'catchem_emis_read_time_coord'

      rc = CC_SUCCESS
      category%n_times = 0

      if (is_null_filename(filename)) then
         write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: No valid filename specified for category: ', trim(category%category_name)
         call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
         write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: File not found for time coord: ', trim(filename)
         call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
         return
      end if

      call AQMIO_ReadTimeCoord(trim(filename), nt, dates, secs, rc=localrc)
      if (localrc /= ESMF_SUCCESS) then
         write(msg, '(A,A,A)') trim(pName), ': FATAL ERROR: Failed reading time coord from: ', trim(filename)
         call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR, rc=localrc)
         rc = CC_FAILURE
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
      call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_DEBUG, rc=localrc)

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


   !> Helper: Clean C-string null terminator and pad with spaces
   subroutine clean_c_string(str)
      character(len=*), intent(inout) :: str
      integer :: idx
      idx = index(str, c_null_char)
      if (idx > 0) then
         str(idx:) = ' '
      end if
   end subroutine clean_c_string

   !> Helper: Check if filename is empty, null, or none
   elemental logical function is_null_filename(fn)
      character(len=*), intent(in) :: fn
      character(len=EMIS_MAXSTR) :: s
      s = adjustl(fn)
      is_null_filename = (len_trim(s) == 0 .or. &
         trim(s) == 'null' .or. trim(s) == 'NULL' .or. trim(s) == 'Null' .or. &
         trim(s) == 'none' .or. trim(s) == 'NONE' .or. trim(s) == 'None')
   end function is_null_filename
end module catchem_nuopc_emis_mod
