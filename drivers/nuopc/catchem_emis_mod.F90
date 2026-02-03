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
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use ConfigManager_Mod, only: ConfigManagerType, ConfigDataType, EmissionCategoryMapping, &
      EmisSpeciesMappingEntry, EmissionMappingConfig
   use StateManager_Mod, only: StateManagerType
   use ChemState_Mod, only: ChemStateType
   use MetState_Mod, only: MetStateType
   use ExtEmisData_Mod, only: ExtEmisDataType, ExtEmisCategoryType, ExtEmisFieldType

   implicit none
   private

   ! Public interfaces
   public :: catchem_emis_init
   public :: catchem_emis_update
   public :: catchem_emis_finalize
   public :: catchem_emis_write_diagnostics


   !> \brief Parameters for emission handling
   integer, parameter :: EMIS_MAXSTR = 256
   integer, parameter :: EMIS_MAXFIELDS = 100
   real(fp), parameter :: EMIS_MISSING = -999.0_fp
   real(fp), parameter :: EMIS_ACCEPT = 1.e+15_fp

   !> \brief Emission timing and alarm information
   type :: EmissionTimingType
      type(ESMF_Alarm) :: alarm
      type(ESMF_TimeInterval) :: time_interval
      integer :: current_record = 1
      logical :: needs_update = .true.
      character(len=64) :: frequency = 'hourly'
      character(len=64) :: category_name = ''  ! Category name for identification
   end type EmissionTimingType

   !> \brief Module-level storage for emission timing information
   !! This array parallels the categories in ExtEmisDataType but stays in this module
   !! to avoid ESMF dependencies in the core data structures
   type(EmissionTimingType), allocatable, save :: category_timings(:)
   integer, save :: n_category_timings = 0

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
      character(len=EMIS_MAXSTR) :: msg
      character(len=*), parameter :: pName = 'catchem_emis_init'

      ! Initialize
      rc = CC_SUCCESS

      ! Check if emission mapping is loaded
      if (.not. config_manager%config_data%emission_mapping%is_loaded) then
         write(msg, '(A,A)') trim(pName), ': Emission mapping not loaded in ConfigManager'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         rc = CC_FAILURE
         return
      end if

      ! Initialize ExtEmisDataType with number of categories from config
      call ext_emis_data%init(config_manager%config_data%emission_mapping%n_categories, &
         'CATChem NUOPC Emission Data', localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to initialize ExtEmisDataType'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         rc = CC_FAILURE
         return
      end if

      ! Initialize parallel timing storage for categories
      n_category_timings = config_manager%config_data%emission_mapping%n_categories
      if (allocated(category_timings)) deallocate(category_timings)
      allocate(category_timings(n_category_timings))

      ! Populate emission categories from already-loaded configuration
      do icat = 1, config_manager%config_data%emission_mapping%n_categories

         !debug
         write(msg, '(A,A,A)') trim(pName), ': Populate category ', &
            trim(config_manager%config_data%emission_mapping%categories(icat)%category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         !end debug

         if (config_manager%config_data%emission_mapping%categories(icat)%is_active) then
            call catchem_emis_populate_category(ext_emis_data, &
               config_manager%config_data%emission_mapping%categories(icat), &
               config_manager, nx, ny, nlev, localrc)
            if (localrc /= CC_SUCCESS) then
               write(msg, '(A,A,A)') trim(pName), ': Failed to populate category ', &
                  trim(config_manager%config_data%emission_mapping%categories(icat)%category_name)
               call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
               rc = CC_FAILURE
               return
            end if

            ! Initialize timing information for this category
            category_timings(icat)%category_name = config_manager%config_data%emission_mapping%categories(icat)%category_name

            !debug: Check what we're accessing from ext_emis_data
            write(msg, '(A,A,I0,A,A)') trim(pName), ': Accessing ext_emis_data%categories(', icat, &
               ') with name: "', trim(ext_emis_data%categories(icat)%category_name)//'"'
            call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
            !end debug

            category_timings(icat)%frequency = trim(ext_emis_data%categories(icat)%frequency)  ! Get frequency from parsed category
            category_timings(icat)%current_record = 0
            category_timings(icat)%needs_update = .false.
            !set up alarm for this category
            call catchem_emis_setup_timing(ext_emis_data%categories(icat), clock, &
               category_timings(icat)%frequency, localrc)

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
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      integer :: localrc, i
      logical :: alarm_ringing
      character(len=EMIS_MAXSTR) :: msg, timeString
      character(len=*), parameter :: pName = 'catchem_emis_update'

      rc = CC_SUCCESS

      ! Get managers from state manager
      config_manager => state_manager%get_config_ptr()
      met_state => state_manager%get_met_state_ptr()
      chem_state => state_manager%get_chem_state_ptr()

      ! Loop through all emission categories and check if updates are needed
      do i = 1, ext_emis_data%n_categories
         if (.not. ext_emis_data%categories(i)%is_active) cycle

         ! Check timing information from parallel storage
         if (i <= n_category_timings) then
            ! Check if emission timing alarm is ringing
            if (allocated(category_timings)) then
               !check alarm ringing
               alarm_ringing = ESMF_AlarmIsRinging(category_timings(i)%alarm, rc=localrc)
               if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

               category_timings(i)%needs_update = alarm_ringing ! Set update flag

               if (alarm_ringing) then
                  !write infor to log
                  call ESMF_TimeGet(current_time, timeString=timeString, rc=localrc)
                  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                     line=__LINE__,  file=__FILE__,  rcToReturn=rc))  return  ! bail out

                  call ESMF_LogWrite(trim(pName)//': reading emission for '//trim(ext_emis_data%categories(i)%category_name)//&
                     " @ "//trim(timeString), ESMF_LOGMSG_INFO, rc=localrc)

                  ! Read new emission data
                  ext_emis_data%categories(i) % irec = ext_emis_data%categories(i) % irec + 1 !time slice one timestep forward
                  call catchem_emis_read(ext_emis_data%categories(i), IO, grid, localrc)
                  if (localrc /= CC_SUCCESS) then
                     write(msg, '(A,A,A)') trim(pName), ': Failed to read data for category: ', &
                        trim(ext_emis_data%categories(i)%category_name)
                     call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=rc)
                  end if

                  !turn off alrm
                  call ESMF_AlarmRingerOff(category_timings(i)%alarm, rc=localrc)
                  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                     line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
                  ! Update timing state
                  category_timings(i)%current_record = category_timings(i)%current_record + 1
                  category_timings(i)%needs_update = .false.  ! Reset until next alarm

                  ! Apply emissions to chemical state
                  call catchem_emis_apply(ext_emis_data%categories(i), i, ext_emis_data%global_scale, config_manager, chem_state, met_state, dt, localrc)
                  if (localrc /= CC_SUCCESS) then
                     write(msg, '(A,A,A)') trim(pName), ': Failed to apply emissions for category: ', &
                        trim(ext_emis_data%categories(i)%category_name)
                     call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=rc)
                  end if
               end if
            end if
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
   subroutine catchem_emis_read(category, IO, grid, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_GridComp), intent(inout) :: IO
      type(ESMF_Grid), intent(in) :: grid
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc,   ifield
      character(len=EMIS_MAXSTR) :: msg, filename
      character(len=64) :: category_name
      type(ESMF_Field) :: esmf_field
      real(ESMF_KIND_R4), pointer :: field_data_2d(:,:) => null()
      character(len=*), parameter :: pName = 'catchem_emis_read'

      rc = CC_SUCCESS

      category_name = trim(category%category_name)

      ! Get filename from category configuration
      filename = trim(category%source_file)
      if (len_trim(filename) == 0) then
         write(msg, '(A,A,A)') trim(pName), ': No source file specified for category: ', trim(category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         rc = CC_FAILURE
         return
      end if

      do ifield = 1, category%n_fields
         !create field to receive data
         esmf_field = ESMF_FieldCreate(grid, name=trim(category%fields(ifield)%field_name), &
            typekind=ESMF_TYPEKIND_R4, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

         !read data into field
         call AQMIO_Read(IO, (/ esmf_field /), fileName=filename, timeSlice=category % irec, &
            iofmt=AQMIO_FMT_NETCDF, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) then
            ! Clean up field before returning
            call ESMF_FieldDestroy(esmf_field, rc=localrc)
            return  ! bail out
         end if

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

         ! Clean up ESMF field after data transfer
         call ESMF_FieldDestroy(esmf_field, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

         ! Nullify pointer for safety
         field_data_2d => null()

      end do

      write(msg, '(A,A,A)') trim(pName), ': Successfully read emission data for category ', &
         trim(category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_read

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
   subroutine catchem_emis_apply(category, icat, global_scale, config_manager, chem_state, met_state, dt, rc)
      use Constants, only: g0, AIRMW  ! Gravitational acceleration and air molecular weight
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      integer, intent(in) :: icat !category index in the ext_emis_data
      real(fp), intent(in) :: global_scale
      type(ConfigManagerType), intent(in) :: config_manager
      type(ChemStateType), intent(inout) :: chem_state
      type(MetStateType), intent(in) :: met_state
      real(fp), intent(in) :: dt
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
      real(fp) :: converter

      rc = CC_SUCCESS

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
               ! Fallback case - need to lookup species index
               species_idx = chem_state%find_species(trim(mapped_species_name))
               if (species_idx <= 0) then
                  write(msg, '(A,A,A)') trim(pName), ': Species not found in chem_state: ', &
                     trim(mapped_species_name)
                  call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=localrc)
                  cycle
               end if
            end if

            ! Convert to final concentration units
            ! For gas species: convert kg/kg to ppmv
            ! For aerosol species: convert kg/kg to ug/kg
            if (chem_state%ChemSpecies(species_idx)%is_gas) then
               converter = AIRMW / chem_state%ChemSpecies(species_idx)%mw_g * 1.0e6_fp
            else
               converter = 1.0e9_fp
            end if
            species_tendency = 0.0_fp

            do j = 1, ny
               do i = 1, nx
                  if (emission_flux(i,j,k) > 0.0_fp) then
                     do k = 1, nz
                        ! Step 1: Convert to mass mixing ratio change (kg/kg) from emission (kg/m2/s)
                        ! Step 2: Convert to kg/kg or ppmv using converter calculated above
                        species_tendency(i,j,k) = emission_flux(i,j,k) * scale_factor *dt * g0 / met_state%DELP(i,j,k) * converter
                     end do
                  end if
               end do
            end do

            ! Add tendency to concentrations
            ! Only apply to first vertical level (k=1) since emissions are surface-based
            concentrations(:,:,:,species_idx) = concentrations(:,:,:,species_idx) + species_tendency(:,:,:)

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
   end subroutine catchem_emis_apply

   !> \brief Write emission diagnostics to NetCDF file
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
            if (.not. ext_emis_data%categories(icat)%fields(ifield)%diagnostic) cycle
            if (.not. ext_emis_data%categories(icat)%fields(ifield)%is_loaded) cycle
            if (.not. allocated(ext_emis_data%categories(icat)%fields(ifield)%emission_data)) cycle

            field_name = trim(ext_emis_data%categories(icat)%fields(ifield)%field_name)
            field_name = "Emis_" // trim(category_name) // trim(field_name)  ! Prefix for diagnostics
            description = trim(ext_emis_data%categories(icat)%fields(ifield)%long_name)
            units = trim(ext_emis_data%categories(icat)%fields(ifield)%units)

            ! Write field based on whether it's gridded (2D) or not (3D)
            if (ext_emis_data%categories(icat)%gridded) then
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
      integer :: i, localrc
      character(len=*), parameter :: pName = 'catchem_emis_finalize'

      rc = CC_SUCCESS

      ! Clean up ExtEmisDataType
      call ext_emis_data%cleanup(localrc)
      if (localrc /= CC_SUCCESS) then
         call ESMF_LogWrite(trim(pName)//': Warning - ExtEmisDataType cleanup failed', &
            ESMF_LOGMSG_WARNING, rc=localrc)
      end if

      ! Clean up module-level timing storage
      if (allocated(category_timings)) then
         ! Destroy ESMF alarms before deallocating
         do i = 1, n_category_timings
            call ESMF_AlarmDestroy(category_timings(i)%alarm, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
         end do
         deallocate(category_timings)
      end if
      n_category_timings = 0

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
      call config_manager%get_logical(trim(config_path)//'/diagnostics', category%diagnostic, localrc, .false.)
      call config_manager%get_real(trim(config_path)//'/scale_factor', category%global_scale, localrc, 1.0_fp)

      ! Read coordinate names
      call config_manager%get_string(trim(config_path)//'/lat_name', category%latname, localrc, '')
      call config_manager%get_string(trim(config_path)//'/lon_name', category%lonname, localrc, '')

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
      call new_category%init(category_mapping%category_name, category_mapping%n_emission_species, &
         'Emission category: '//trim(category_mapping%category_name), localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to initialize category'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         rc = CC_FAILURE
         return
      end if

      !debug
      write(msg, '(A,A)') trim(pName), ': Populating category '//trim(new_category%category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
      !end debug

      ! Set category properties from mapping
      new_category%is_active = category_mapping%is_active

      ! Parse additional properties from configuration using ConfigManager functions
      call parse_emission_category(new_category, config_manager, category_mapping%category_name, localrc, diag_species_list)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A,A)') trim(pName), ': Failed to parse category properties: ', &
            trim(category_mapping%category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=rc)
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
               call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=rc)
            end if
         end if
      end do

      ! Add category to ExtEmisDataType
      !debug: Check category name before adding
      write(msg, '(A,A,A)') trim(pName), ': Adding category with name: "', &
         trim(new_category%category_name)//'"'
      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
      !end debug

      call ext_emis_data%add_category(new_category, localrc)
      if (localrc /= CC_SUCCESS) then
         write(msg, '(A,A)') trim(pName), ': Failed to add category to ExtEmisDataType'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         rc = CC_FAILURE
         return
      end if

      !debug: Check category name after adding
      if (ext_emis_data%n_categories > 0) then
         write(msg, '(A,A,A)') trim(pName), ': After adding, category name is: "', &
            trim(ext_emis_data%categories(ext_emis_data%n_categories)%category_name)//'"'
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
      end if
      !end debug

      write(msg, '(A,A,A)') trim(pName), ': Successfully populated category ', &
         trim(category_mapping%category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)

   end subroutine catchem_emis_populate_category

   !> \brief Setup emission timing alarms
   !!
   !! Creates ESMF alarms for each emission category based on update frequency.
   !! Stores alarms in module-level category_timings array.
   !!
   !! \param[in] category_name Name of emission category
   !! \param[in] clock Model clock for alarm creation
   !! \param[in] frequency Update frequency string (e.g., 'hourly', 'daily')
   !! \param[out] rc Return code
   subroutine catchem_emis_setup_timing(category, clock, frequency, rc)
      implicit none

      type(ExtEmisCategoryType), intent(inout) :: category
      type(ESMF_Clock), intent(in) :: clock
      character(len=*), intent(in) :: frequency
      integer, intent(out) :: rc

      ! Local variables
      integer :: localrc, cat_idx
      type(ESMF_Time) :: startTime, currTime
      type(ESMF_TimeInterval) :: timeInterval
      character(len=EMIS_MAXSTR) ::  msg
      character(len=*), parameter :: pName = 'catchem_emis_setup_timing'

      rc = CC_SUCCESS

      !debug: Check what category name we received
      write(msg, '(A,A,A)') trim(pName), ': Received category with name: "', &
         trim(category%category_name)//'"'
      call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=localrc)
      !end debug

      ! Find the category index in timing storage
      cat_idx = 0
      if (allocated(category_timings)) then
         do cat_idx = 1, n_category_timings
            if (trim(category_timings(cat_idx)%category_name) == trim(category%category_name)) then
               exit
            end if
         end do
         if (cat_idx > n_category_timings) cat_idx = 0
      end if

      if (cat_idx == 0) then
         write(msg, '(A,A,A)') trim(pName), ': Category not found in timing storage: ', trim(category%category_name)
         call ESMF_LogWrite(msg, ESMF_LOGMSG_ERROR, rc=rc)
         rc = CC_FAILURE
         return
      end if

      ! Set ring interval based on frequency
      select case (trim(frequency))
       case ("hourly")
         call ESMF_TimeIntervalSet(timeInterval, h=1, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
       case ("daily")
         call ESMF_TimeIntervalSet(timeInterval, d=1, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
       case ("weekly")
         call ESMF_TimeIntervalSet(timeInterval, d=7, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
       case ("monthly")
         call ESMF_TimeIntervalSet(timeInterval, mm=1, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
       case ("static")
         call ESMF_TimeIntervalSet(timeInterval, rc=localrc)
         if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out
       case default
         call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="- unknown emission frequency: "//trim(frequency), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      end select

      call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out

      ! -- set input time record according to start type (startup/continue)
      category % irec = int( (currTime - startTime) / timeInterval )

      category_timings(cat_idx)%alarm = ESMF_AlarmCreate(clock, ringTime=startTime, &
         ringInterval=timeInterval, name=trim(category%category_name)//"_alarm", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__,  file=__FILE__,  rcToReturn=rc)) return  ! bail out


      ! Store timing interval
      category_timings(cat_idx)%time_interval = timeInterval

      write(msg, '(A,A,A,A,A)') trim(pName), ': Created alarm ', trim(category%category_name)//"_alarm", &
         ' for category ', trim(category%category_name)
      call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=localrc)

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
      integer ::  j, icat, ispec
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
         call ESMF_LogWrite(msg, ESMF_LOGMSG_WARNING, rc=rc)
      end select

   end subroutine catchem_emis_map_species

end module catchem_emis_mod
