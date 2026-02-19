!> \file SO4chemScheme_GOCART_Mod.F90
!! \brief GOCART SO2 to SO4 production scheme
!!
!! Pure science kernel for gocart scheme in so4chem process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_gocart (search for "TODO")
!! 2. Add scheme-specific helper subroutines as needed
!! 3. Update physical constants for your scheme
!! 4. Customize the environmental response functions
!!
!! INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
!! - Parameter initialization and validation
!! - Input array validation and error handling
!! - Memory management and array allocation
!! - Integration with host model time stepping
!!
!! Generated on: 2026-02-11T13:30:17.194715
!! Author: Wei Li
!! Reference: GOCART2G process library SulfateChemDriver function
module SO4chemScheme_GOCART_Mod

   use precision_mod, only: fp, rae
   use SO4chemCommon_Mod, only: SO4chemSchemeGOCARTConfig
   use error_mod, only: CC_SUCCESS, CC_Error
   use GOCART2G_Process, only: SulfateUpdateOxidants, SulfateChemDriver

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: plid = 0.01_fp    ! Pressure lid [hPa]
   real(fp), parameter :: undefval= 1.0e+15_fp    ! Same as MAPL library

contains

   !> Pure science computation for gocart scheme
   !!
   !! This is a pure computational kernel implementing GOCART SO2 to SO4 production scheme.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  g0    Required constant from Constants module
   !! @param[in]  Cpd    Required constant from Constants module
   !! @param[in]  AVO    Required constant from Constants module
   !! @param[in]  VON_KARMAN    Required constant from Constants module
   !! @param[in]  AIRMW    Required constant from Constants module
   !! @param[in]  PI    Required constant from Constants module
   !! @param[in]  year    Time parameter from TimeState (year)
   !! @param[in]  month    Time parameter from TimeState (month)
   !! @param[in]  day    Time parameter from TimeState (day)
   !! @param[in]  hour    Time parameter from TimeState (hour)
   !! @param[in]  minute    Time parameter from TimeState (minute)
   !! @param[in]  second    Time parameter from TimeState (second)
   !! @param[in]  airden    AIRDEN field [appropriate units]
   !! @param[in]  cldf    CLDF field [appropriate units]
   !! @param[in]  delp    DELP field [appropriate units]
   !! @param[in]  hflux    HFLUX field [appropriate units]
   !! @param[in]  lat    LAT field [appropriate units]
   !! @param[in]  lon    LON field [appropriate units]
   !! @param[in]  lwi    LWI field [appropriate units]
   !! @param[in]  pblh    PBLH field [appropriate units]
   !! @param[in]  pmid    PMID field [appropriate units]
   !! @param[in]  t    T field [appropriate units]
   !! @param[in]  tstep    Time step [s] - retrieved from process interface
   !! @param[in]  ustar    USTAR field [appropriate units]
   !! @param[in]  z    Z field [appropriate units]
   !! @param[in]  z0h    Z0H field [appropriate units]
   !! @param[in]  species_mw_g    Species mw_g property
   !! @param[in]  species_short_name    Species short_name property
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] PSO4_from_SO2_per_level    total sulfate production rate from SO2 per level [kg/kg/s] (num_layers)
   !! @param[inout] PSO4_from_gaseous_SO2_per_level    sulfate production rate from gaseous SO2 per level [kg/kg/s] (num_layers)
   !! @param[inout] PSO4_from_aqueous_SO2_per_level    sulfate production rate from aqueous SO2 per level [kg/kg/s] (num_layers)
   !! @param[inout] PSO2_from_DMS_per_level    SO2 production rate from DMS per level [kg/kg/s] (num_layers)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   subroutine compute_gocart( &
      num_layers, &
      num_species, &
      params, &
      g0, &
      Cpd, &
      AVO, &
      VON_KARMAN, &
      AIRMW, &
      PI, &
      year, &
      month, &
      day, &
      hour, &
      minute, &
      second, &
      airden, &
      cldf, &
      delp, &
      hflux, &
      lat, &
      lon, &
      lwi, &
      pblh, &
      pmid, &
      t, &
      tstep, &
      ustar, &
      z, &
      z0h, &
      species_mw_g, &
      species_short_name, &
      species_conc, &
      species_tendencies, &
      Production_rate_per_species_per_level, &
      PSO4_from_gaseous_SO2_per_level, &
      PSO4_from_aqueous_SO2_per_level, &
      diagnostic_species_id &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(SO4chemSchemeGOCARTConfig), intent(in) :: params
      real(fp), intent(in) :: g0  ! Required constant from Constants module
      real(fp), intent(in) :: Cpd  ! Required constant from Constants module
      real(fp), intent(in) :: AVO  ! Required constant from Constants module
      real(fp), intent(in) :: VON_KARMAN  ! Required constant from Constants module
      real(fp), intent(in) :: AIRMW  ! Required constant from Constants module
      real(fp), intent(in) :: PI  ! Required constant from Constants module
      integer, intent(in) :: year  ! Time parameter from TimeState
      integer, intent(in) :: month  ! Time parameter from TimeState
      integer, intent(in) :: day  ! Time parameter from TimeState
      integer, intent(in) :: hour  ! Time parameter from TimeState
      integer, intent(in) :: minute  ! Time parameter from TimeState
      integer, intent(in) :: second  ! Time parameter from TimeState
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: cldf(num_layers)  ! Surface field - scalar
      real(fp), intent(in) :: delp(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: hflux  ! Surface field - scalar
      real(fp), intent(in) :: lat  ! Surface field - scalar
      real(fp), intent(in) :: lon  ! Surface field - scalar
      integer, intent(in) :: lwi  ! Surface field - scalar
      real(fp), intent(in) :: pblh  ! Surface field - scalar
      real(fp), intent(in) :: pmid(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: t(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      real(fp), intent(in) :: ustar  ! Surface field - scalar
      real(fp), intent(in) :: z(num_layers+1)  ! Edge field - requires nz+1 dimensions
      real(fp), intent(in) :: z0h  ! Surface field - scalar
      real(fp), intent(in) :: species_mw_g(:)  ! Species mw_g property
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: Production_rate_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: PSO4_from_gaseous_SO2_per_level(:)
      real(fp), intent(inout), optional :: PSO4_from_aqueous_SO2_per_level(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: klid = 1 !since the layer is reversed, we give 1 here, which is the top layer
      integer :: diag_idx  ! For diagnostic species indexing
      integer :: species_idx
      integer :: nDMS= -1, nSO2= -1, nSO4= -1, nMSA= -1 ! index position of sulfates
      integer :: nOH= -1, nNO3= -1, nH2O2= -1 ! index position of oxidants
      integer :: nymd, nhms   !YYYYMMDD, HHMMSS time formats
      real(fp), allocatable :: latRad(:,:), lonRad(:,:)
      real(fp) :: fMassMSA, fMassDMS, fMassSO2, fMassSO4 ! gram molecular weights of species
      real(fp) :: rad2deg,  deg2rad  ! PI cannot be used here
      ! Local Variables
      real(fp), pointer :: GOCART_tmpu(:,:,:)
      real(fp), pointer :: GOCART_rhoa(:,:,:)
      real(fp), pointer :: GOCART_HGHTE(:,:,:)
      real(fp), pointer :: GOCART_DELP(:,:,:)
      real(fp), pointer :: GOCART_PRESS(:,:,:)
      real(fp), pointer :: GOCART_cloud(:,:,:)
      real(fp), pointer :: GOCART_LWI(:,:)
      real(fp), pointer :: GOCART_USTAR(:,:)
      real(fp), pointer :: GOCART_PBLH(:,:)
      real(fp), pointer :: GOCART_HFLUX(:,:)
      real(fp), pointer :: GOCART_Z0H(:,:)
      !some chem variables to be populated
      !Monthly climatology of these three oxidenats from GMI is read in and we store them in chem_state arrays.
      real(fp), pointer, dimension(:,:,:) :: oh_clim    !volume mixing ratio
      real(fp), pointer, dimension(:,:,:) :: h2o2_clim  ![# cm-3]
      real(fp), pointer, dimension(:,:,:) :: no3_clim   ![# cm-3]
      !OH and NO3 will go through diurnal variation scaling based on solar zenith angle, while H2O2 is reset to
      !climatology every three hours
      real(fp), dimension(:,:,:), allocatable :: xoh, xno3, xh2o2   !kg/kg
      real(fp), dimension(:,:,:), allocatable :: dms, so2, so4 !kg/kg
      real(fp), pointer, dimension(:,:,:) :: msa  !kg/kg
      real(fp), pointer, dimension(:,:,:) :: SU_dep  ! Sulfate Dry Deposition All Bins [kg/m2/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO2 ! vertical sum of SO2 Prod from DMS oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:) :: SU_PMSA ! vertical sum of MSA Prod from DMS oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO4 ! vertical sum of SO4 Prod from all SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO4g ! vertical sum of SO4 Prod from gaseous SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO4aq ! vertical sum of SO4 Prod from aqueous SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso2  ! SO2 Prod from DMS oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pmsa  ! MSA Prod from DMS oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso4  ! SO4 Prod from all SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso4g  ! SO4 Prod from gaseous SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso4aq  ! SO4 Prod from aqueous SO2 oxidation [kg/kg/s]
      real(fp), dimension(:,:), allocatable :: drydepfrequency
      !some clock variables for oxidants updates
      logical, save :: firsttime = .true.
      integer, save :: nymd_last = -1  ! NYMD of last H2O2 update
      real(fp), save :: tstep_counter = 0.0_fp ! counter (in seconds) for tracking when to update H2O2 based on 3-hourly climatology
      real(fp), allocatable, save :: xh2o2_init(:,:,:) ! initial H2O2 from last time step
      logical :: recycle_h2o2  = .false.
      !error information
      integer :: RC
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at compute_gocart (in SO4chemScheme_GOCART_Mod.F90)'
      !RC = CC_SUCCESS
      RC = 0 !try not to rely on CC_SUCCESS
      !drydepf = 0.0_fp

      rad2deg = 180.0_fp/PI
      deg2rad = PI/180.0_fp
      !construct time in yyyymmdd and hhmmss formats for use in gocart
      nymd = year*10000 + month*100 + day
      nhms = hour*10000 + minute*100 + second

      !get species indices for use in gocart
      do species_idx = 1, num_species
         if (species_short_name(species_idx) == 'SO2' .or. species_short_name(species_idx) == 'so2') then
            nSO2 = species_idx
         else if (species_short_name(species_idx) == 'SO4' .or. species_short_name(species_idx) == 'so4') then
            nSO4 = species_idx
         else if (species_short_name(species_idx) == 'DMS' .or. species_short_name(species_idx) == 'dms') then
            nDMS = species_idx
         else if (species_short_name(species_idx) == 'MSA' .or. species_short_name(species_idx) == 'msa') then
            nMSA = species_idx
         else if (species_short_name(species_idx) == 'OH' .or. species_short_name(species_idx) == 'oh') then
            nOH = species_idx
         else if (species_short_name(species_idx) == 'NO3' .or. species_short_name(species_idx) == 'no3') then
            nNO3 = species_idx
         else if (species_short_name(species_idx) == 'H2O2' .or. species_short_name(species_idx) == 'h2o2') then
            nH2O2 = species_idx
         end if
      end do

      if (nSO2 == -1 .or. nSO4 == -1 .or. nDMS == -1 .or. nMSA == -1 .or. nOH == -1 .or. nNO3 == -1 .or. nH2O2 == -1) then
         errMsg = 'Error in compute_gocart: SO2, SO4, DMS, MSA, OH, NO3, and H2O2 must be present in species list.'
         !call CC_Error(trim(errMsg), RC, thisLoc)
         write(*,'(A)') trim(errMsg)
         return
      end if

      !allocate arrays
      allocate(oh_clim(1,1,num_layers), h2o2_clim(1,1,num_layers), no3_clim(1,1,num_layers), &
         xoh(1,1,num_layers), xno3(1,1,num_layers), xh2o2(1,1,num_layers), &
         dms(1,1,num_layers), so2(1,1,num_layers), so4(1,1,num_layers), msa(1,1,num_layers), &
         SU_dep(1, 1, num_species), SU_PSO2(1, 1), SU_PMSA(1, 1), SU_PSO4(1, 1), SU_PSO4g(1, 1), SU_PSO4aq(1, 1), &
         pso2(1, 1, num_layers), pmsa(1, 1, num_layers), pso4(1, 1, num_layers), pso4g(1, 1, num_layers), &
         pso4aq(1, 1, num_layers), drydepfrequency(1, 1), latRad(1,1), lonRad(1,1))


      !retrieve climatology fields; remember to reverse the vertical layer (TODO: double check the input files for this)
      oh_clim(1,1,:) = species_conc(num_layers:1:-1, nOH)  !TODO: the unit is not in kg/kg and it should not be exchanged with NUOPC although it is in the chem_state.
      no3_clim(1,1,:) = species_conc(num_layers:1:-1, nNO3)
      h2o2_clim(1,1,:) = species_conc(num_layers:1:-1, nH2O2)
      ! Initialize some variables for the first time
      if (firsttime) then
         nymd_last = nymd
         if (.not. allocated(xh2o2_init)) then
            allocate(xh2o2_init(1,1,num_layers))
         end if
         xh2o2_init = h2o2_clim  ! initialize H2O2 to climatology at first time step
         firsttime = .false.
      end if

      !update tstep_counter
      tstep_counter = tstep_counter + tstep
      if (tstep_counter >= 10800) then
         recycle_h2o2 = .true.
         tstep_counter = 0
      end if

      ! transform data for GOCART DryDeposition call
      call PrepMetVarsForGOCART(num_layers,     &
         t,               &
         airden,          &
         z,               &
         cldf,            &
         delp,            &
         lwi,             &
         ustar,           &
         pblh,            &
         pmid,            &
         hflux,           &
         z0h,             &
         GOCART_tmpu,     &
         GOCART_RHOA,     &
         GOCART_HGHTE,    &
         GOCART_cloud,    &
         GOCART_DELP,     &
         GOCART_LWI,      &
         GOCART_USTAR,    &
         GOCART_PBLH,     &
         GOCART_PRESS,    &
         GOCART_HFLUX,    &
         GOCART_Z0H)

      !update oxidants based on climatology and diurnal cycle
      xoh = 0.0_fp; xno3 = 0.0_fp; xh2o2 = xh2o2_init
      latRad(1,1) = lat * deg2rad; lonRad(1,1) = lon * deg2rad
      call SulfateUpdateOxidants(nymd, nhms, lonRad, latRad, GOCART_rhoa, num_layers, tstep, nymd_last, &
         undefval, rad2deg, AVO, PI, AIRMW, oh_clim, no3_clim, h2o2_clim, xoh, xno3, xh2o2, recycle_h2o2, RC)

      if (RC /= 0) then
         ErrMsg = 'Error in compute_gocart: Failed in updating oxidants in GOCART So4chem process.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(ErrMsg)
         return
      end if

      !get pressure lid index
      call findKlid(klid, plid, GOCART_PRESS(:,:,:), RC)
      !if (RC /= CC_SUCCESS) then
      if (RC /= 0) then
         ErrMsg = 'Error in compute_gocart: Failed in finding pressure lid index in GOCART So4chem process.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(ErrMsg)
         return
      end if

      !retrieve sulfate species concentrations
      dms(1,1,:) = species_conc(num_layers:1:-1, nDMS)
      so2(1,1,:) = species_conc(num_layers:1:-1, nSO2)
      so4(1,1,:) = species_conc(num_layers:1:-1, nSO4)
      msa(1,1,:) = species_conc(num_layers:1:-1, nMSA)
      fMassMSA = species_mw_g(nMSA)
      fMassDMS = species_mw_g(nDMS)
      fMassSO2 = species_mw_g(nSO2)
      fMassSO4 = species_mw_g(nSO4)

      !debug
      write(*,*) 'GOCART chemistry driver input:', 'nymd=', nymd, 'nhms=', nhms, 'lat=', lat, 'lon=', lon, &
         'dms=', dms(1,1,num_layers), 'so2=', so2(1,1,num_layers), 'so4=', so4(1,1,num_layers), 'msa=', msa(1,1,num_layers), &
         'oh=', xoh(1,1,num_layers), 'no3=', xno3(1,1,num_layers), 'h2o2=', xh2o2(1,1,num_layers)

      !call GOCART sulfate chemistry driver
      call SulfateChemDriver(num_layers, klid, tstep, PI, rad2deg, VON_KARMAN, AIRMW, AVO, Cpd, g0, fMassMSA,fMassDMS,fMassSO2,fMassSO4,&
         nymd, nhms, lonRad, latRad, dms, so2, so4, msa, nDMS, nSO2, nSO4, nMSA, xoh, xno3, xh2o2, xh2o2_init, GOCART_DELP, GOCART_tmpu, GOCART_cloud, &
         GOCART_rhoa, GOCART_HGHTE, GOCART_USTAR, GOCART_HFLUX, GOCART_LWI, GOCART_PBLH, GOCART_Z0H, SU_dep, SU_PSO2, SU_PMSA, SU_PSO4, SU_PSO4g, &
         SU_PSO4aq, pso2, pmsa, pso4, pso4g, pso4aq, drydepfrequency, RC)

      if (RC /= 0) then
         ErrMsg = 'Error in compute_gocart: Failed in GOCART sulfate chemistry driver.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(ErrMsg)
         return
      end if

      !debug
      write(*,*) 'GOCART chemistry driver output:', 'dryfre=', drydepfrequency, &
         'dms=', dms(1,1,num_layers), 'so2=', so2(1,1,num_layers), 'so4=', so4(1,1,num_layers), 'msa=', msa(1,1,num_layers), &
         'pmsa=', pmsa(1,1,num_layers), 'pso2=', pso2(1,1,num_layers), 'pso4=', pso4(1,1,num_layers), 'pso4g=', pso4g(1,1,num_layers), 'pso4aq=', pso4aq(1,1,num_layers)

      !assign to output tendencies; remember to reverse the vertical layer back to original order
      if (params%update_so2) then !since the chem driver has drydep in it, not sure if we should update so2 chem array here.
         species_tendencies(:, nSO2) = so2(1,1,num_layers:1:-1)
      end if
      species_tendencies(:, nSO4) = so4(1,1,num_layers:1:-1)
      species_tendencies(:, nMSA) = msa(1,1,num_layers:1:-1)
      !species_tendencies(:, nDMS) = dms(1,1,num_layers:1:-1) !!!TODO: DMS is not updated for now since it is read in through monthly files. We should uncommont this once the emission is in.

      ! Per-species-per-level diagnostic: 2D array (levels, species)
      if (present(Production_rate_per_species_per_level) .and. present(diagnostic_species_id)) then
         ! Find position of this species in diagnostic_species_id array
         do diag_idx = 1, size(diagnostic_species_id)
            if (diagnostic_species_id(diag_idx) == nMSA) then
               ! Add your custom production rate (dms to so2, dms to msa, so2 to so4) per species per level calculation
               Production_rate_per_species_per_level(:, diag_idx) = pmsa(1,1,num_layers:1:-1)
            end if
            if (diagnostic_species_id(diag_idx) == nSO2) then
               ! Add your custom production rate (dms to so2, dms to msa, so2 to so4) per species per level calculation
               Production_rate_per_species_per_level(:, diag_idx) = pso2(1,1,num_layers:1:-1)
            end if
            if (diagnostic_species_id(diag_idx) == nSO4) then
               ! Add your custom production rate (dms to so2, dms to msa, so2 to so4) per species per level calculation
               Production_rate_per_species_per_level(:, diag_idx) = pso4(1,1,num_layers:1:-1)
            end if
         end do
      end if

      if (present(PSO4_from_gaseous_SO2_per_level)) then
         PSO4_from_gaseous_SO2_per_level = PSO4g(1,1,num_layers:1:-1)
      end if

      if (present(PSO4_from_aqueous_SO2_per_level)) then
         PSO4_from_aqueous_SO2_per_level = PSO4aq(1,1,num_layers:1:-1)
      end if


      !cleanup pointers
      if (associated(GOCART_TMPU)) nullify(GOCART_TMPU)
      if (associated(GOCART_RHOA)) nullify(GOCART_RHOA)
      if (associated(GOCART_HGHTE)) nullify(GOCART_HGHTE)
      if (associated(GOCART_DELP)) nullify(GOCART_DELP)
      if (associated(GOCART_cloud)) nullify(GOCART_cloud)
      if (associated(GOCART_PRESS)) nullify(GOCART_PRESS)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)
      if (associated(GOCART_USTAR)) nullify(GOCART_USTAR)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)
      if (associated(GOCART_HFLUX)) nullify(GOCART_HFLUX)
      if (associated(GOCART_Z0H)) nullify(GOCART_Z0H)
      if (associated(SU_dep)) nullify(SU_dep)
      if (associated(SU_PSO2)) nullify(SU_PSO2)
      if (associated(SU_PMSA)) nullify(SU_PMSA)
      if (associated(SU_PSO4)) nullify(SU_PSO4)
      if (associated(SU_PSO4g)) nullify(SU_PSO4g)
      if (associated(SU_PSO4aq)) nullify(SU_PSO4aq)
      if (associated(pso2)) nullify(pso2)
      if (associated(pmsa)) nullify(pmsa)
      if (associated(pso4)) nullify(pso4)
      if (associated(pso4g)) nullify(pso4g)
      if (associated(pso4aq)) nullify(pso4aq)
      if (associated(msa)) nullify(msa)
      if (associated(oh_clim)) nullify(oh_clim)
      if (associated(no3_clim)) nullify(no3_clim)
      if (associated(h2o2_clim)) nullify(h2o2_clim)
      !cleanup array allocations
      deallocate( xoh, xno3, xh2o2, dms, so2, so4, drydepfrequency, latRad, lonRad)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines

   !>
   !! \brief PrepMetVarsForGOCART - Prep the meteorological variables for GOCART DryDeposition scheme
   !!
   !! \param [INOUT] metstate
   !! \param [INOUT] tmpu
   !! \param [INOUT] rhoa
   !! \param [INOUT] hghte
   !! \param [INOUT] oro
   !! \param [INOUT] ustar
   !! \param [INOUT] pblh
   !! \param [INOUT] shflux
   !! \param [INOUT] z0h
   !! \param [INOUT] u10m
   !! \param [INOUT] v10m
   !! \param [INOUT] fraclake
   !! \param [INOUT] gwettop
   !! \param [OUT] rc
   !!
   !! \ingroup core_modules
   !!!>
   pure subroutine PrepMetVarsForGOCART(km,              &
      tmpu,            &
      rhoa,            &
      hghte,           &
      cldfrc,          &
      delp,            &
      lwi,             &
      ustar,           &
      pblh,            &
      pmid,            &
      hflux,           &
      z0h,             &
      GOCART_tmpu,     &
      GOCART_RHOA,     &
      GOCART_HGHTE,    &
      GOCART_cloud,    &
      GOCART_DELP,     &
      GOCART_LWI,      &
      GOCART_USTAR,    &
      GOCART_PBLH,     &
      GOCART_PRESS,    &
      GOCART_HFLUX,    &
      GOCART_Z0H)



      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      INTEGER,  intent(in)                    :: lwi                                    ! orography flag; Land, ocean, ice mask
      REAL(fp),  intent(in), DIMENSION(:), target :: tmpu   ! Temperature [K]
      REAL(fp),  intent(in), DIMENSION(:), target :: rhoa   ! Air density [kg/m^3]
      REAL(fp),  intent(in), DIMENSION(:), target :: hghte  ! Height [m]
      REAL(fp),  intent(in), DIMENSION(:), target :: cldfrc  ! Cloud fraction [1]
      REAL(fp),  intent(in), DIMENSION(:), target :: delp    ! Pressure thickness [Pa]
      REAL(fp),  intent(in), DIMENSION(:), target :: pmid    ! Pressure at mid-layer [Pa]
      REAL(fp),  intent(in), target               :: ustar                                 ! friction speed [m/sec]
      REAL(fp),  intent(in), target              :: pblh                                  ! PBL height [m]
      REAL(fp),  intent(in), target              :: hflux                                 ! sfc. sens. heat flux [W m-2]
      REAL(fp),  intent(in), target              :: z0h                                   ! rough height, sens. heat [m]

      ! INPUT/OUTPUTS
      REAL(fp), intent(inout), pointer :: GOCART_TMPU(:,:,:)   !< temperature [K]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA   !< air density [kg/m^3]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE  !< geometric height [m]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_cloud  !< cloud fraction [1]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_DELP    !< pressure thickness [Pa]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_PRESS    !< pressure at mid-layer [Pa]
      real(fp), intent(inout), pointer :: GOCART_LWI(:,:)                 !< orography flag; Land, ocean, ice mask
      REAL(fp), intent(inout), pointer :: GOCART_USTAR(:,:)               !< friction speed [m/sec]
      REAL(fp), intent(inout), pointer :: GOCART_PBLH(:,:)                !< PBL height [m]
      REAL(fp), intent(inout), pointer :: GOCART_HFLUX(:,:)               !< sfc. sens. heat flux [W m-2]
      REAL(fp), intent(inout), pointer :: GOCART_Z0H(:,:)                 !< rough height, sens. heat [m]

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(GOCART_TMPU(1, 1, km))
      allocate(GOCART_RHOA(1, 1, km))
      allocate(GOCART_HGHTE(1, 1, 0:km))
      allocate(GOCART_cloud(1, 1, km))
      allocate(GOCART_DELP(1, 1, km))
      allocate(GOCART_PRESS(1, 1, km))
      allocate(GOCART_LWI(1, 1))
      allocate(GOCART_USTAR(1, 1))
      allocate(GOCART_PBLH(1, 1))
      allocate(GOCART_HFLUX(1, 1))
      allocate(GOCART_Z0H(1, 1))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)

      GOCART_TMPU(1,1,:) = tmpu(size(tmpu):1:-1) ! temperature [K]
      GOCART_RHOA(1,1,:) = rhoa(size(rhoa):1:-1) ! air density [kg/m^3]
      GOCART_cloud(1,1,:) = cldfrc(size(cldfrc):1:-1) ! cloud fraction [1]
      GOCART_DELP(1,1,:) = delp(size(delp):1:-1) ! pressure thickness [Pa]
      GOCART_HGHTE(1,1,:) = hghte(size(hghte):1:-1)    ! top of layer geopotential height [m]
      GOCART_PRESS(1,1,:) = pmid(size(pmid):1:-1)    ! pressure at mid-layer [Pa]
      GOCART_LWI = real(LWI, fp)     ! orography flag; Land, ocean, ice mask
      GOCART_USTAR  = ustar

      ! friction speed [m/sec]
      GOCART_PBLH   = pblh      ! PBL height [m]
      GOCART_HFLUX = hflux     ! sfc. sens. heat flux [W m-2]
      GOCART_Z0H    = z0h       ! rough height, sens. heat [m]


   end subroutine PrepMetVarsForGOCART

   !>
   !! \brief findKlid - Finds corresponding vertical index for defined pressure lid
   !!
   !! \param [INOUT] klid
   !! \param [IN] plid
   !! \param [IN] ple
   !! \param [OUT] rc
   !!!>
   subroutine findKlid (klid, plid, ple, rc)

      implicit NONE
      ! !INPUT PARAMETERS:
      integer, intent(inout) :: klid ! index for pressure lid
      real(fp), intent(in)       :: plid ! pressure lid [hPa]; default is 0.01 hPa
      real(fp), dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]
      ! !OUTPUT PARAMETERS:
      integer, intent(out) :: rc ! return code; 0 - all is good; 1 - bad
      ! !Reference to gocart: https://github.com/GEOS-ESM/GOCART/blob/9ff3df9545dd582f415f682d3297e8c6c841e5cb/ESMF/Shared/Chem_AeroGeneric.F90#L316
      ! !Local Variables
      integer :: k, j, i
      real(fp) :: plid_, diff, refDiff
      real(fp), allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]
      !EOP
      !----------------------------------------------------------------------------------
      !  Begin...
      klid = 1
      rc = 0

      !  convert from hPa to Pa
      plid_ = plid*100.0_fp

      allocate(pres(ubound(ple,3)))

      !  find pressure at each model level
      do k = 1, ubound(ple,3)
         pres(k) = ple(1,1,k)
      end do

      !  find smallest absolute difference between plid and average pressure at each model level
      refDiff = 150000.0_fp
      do k = 1, ubound(ple,3)
         diff = abs(pres(k) - plid_)
         if (diff < refDiff) then
            klid = k
            refDiff = diff
         end if
      end do

      !  Check to make sure that all pressures at (i,j) were the same
      do j = 1, ubound(ple,2)
         do i = 1, ubound(ple,1)
            !if (pres(klid) /= ple(i,j,klid)) then !This gives a warning for floating point comparison. Use rae instead
            if (.not. rae(pres(klid), ple(i,j,klid))) then
               rc = 1
               return
            end if
         end do
      end do

   end subroutine findKlid

end module SO4chemScheme_GOCART_Mod
