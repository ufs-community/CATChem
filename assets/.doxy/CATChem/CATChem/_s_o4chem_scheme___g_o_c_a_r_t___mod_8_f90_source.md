

# File SO4chemScheme\_GOCART\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**so4chem**](dir_fb8fc0df5ebe1b02f5e46b98d91cbc63.md) **>** [**schemes**](dir_429bccfa51a729cf5e11bef5cf73cbc7.md) **>** [**SO4chemScheme\_GOCART\_Mod.F90**](_s_o4chem_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the documentation of this file](_s_o4chem_scheme___g_o_c_a_r_t___mod_8_f90.md)


```Fortran

module so4chemscheme_gocart_mod

   use precision_mod, only: fp, rae
   use so4chemcommon_mod, only: so4chemschemegocartconfig
   use error_mod, only: cc_success, cc_error
   use gocart2g_process, only: sulfateupdateoxidants, sulfatechemdriver, dmsemission

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: plid = 0.01_fp    ! Pressure lid [hPa]
   real(fp), parameter :: undefval= 1.0e+15_fp    ! Same as MAPL library

contains

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
      u10m, &
      ustar, &
      v10m, &
      z, &
      z0h, &
      species_mw_g, &
      species_short_name, &
      species_conc, &
      species_tendencies, &
      firsttime, &
      nymd_last, &
      nhms_last_recycle, &
      xh2o2_init, &
      Production_rate_per_species_per_level, &
      PSO4_from_gaseous_SO2_per_level, &
      PSO4_from_aqueous_SO2_per_level, &
      DMS_emission_flux, &
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
      real(fp), intent(in) :: u10m  ! Surface field - scalar
      real(fp), intent(in) :: ustar  ! Surface field - scalar
      real(fp), intent(in) :: v10m  ! Surface field - scalar
      real(fp), intent(in) :: z(num_layers+1)  ! Edge field - requires nz+1 dimensions
      real(fp), intent(in) :: z0h  ! Surface field - scalar
      real(fp), intent(in) :: species_mw_g(:)  ! Species mw_g property
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      ! Per-column persistent state variables
      logical, intent(inout) :: firsttime  ! flag for first time step
      integer, intent(inout) :: nymd_last  ! last day of H2O2 update
      integer, intent(inout) :: nhms_last_recycle  ! last time step of H2O2 recycle
      real(fp), intent(inout), allocatable :: xh2o2_init(:)  ! H2O2 column initialization
      real(fp), intent(inout), optional :: Production_rate_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: PSO4_from_gaseous_SO2_per_level(:)
      real(fp), intent(inout), optional :: PSO4_from_aqueous_SO2_per_level(:)
      real(fp), intent(inout), optional :: DMS_emission_flux
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: klid = 1 !since the layer is reversed, we give 1 here, which is the top layer
      integer :: diag_idx  ! For diagnostic species indexing
      integer :: species_idx
      integer :: nDMS= -1, nso2= -1, nso4= -1, nmsa= -1, ndms_in= -1 ! index position of sulfates
      integer :: nOH= -1, nno3= -1, nh2o2= -1 ! index position of oxidants
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
      real(fp), pointer :: GOCART_U10M(:,:)
      real(fp), pointer :: GOCART_V10M(:,:)
      !some chem variables to be populated
      !Monthly climatology of these three oxidenats from GMI is read in and we store them in chem_state arrays.
      real(fp), pointer, dimension(:,:,:) :: oh_clim    !volume mixing ratio [mol/mol]
      real(fp), pointer, dimension(:,:,:) :: h2o2_clim  !volume mixing ratio [mol/mol]
      real(fp), pointer, dimension(:,:,:) :: no3_clim   !volume mixing ratio [mol/mol]
      !OH and NO3 will go through diurnal variation scaling based on solar zenith angle, while H2O2 is reset to
      !climatology every three hours and every new day
      real(fp), dimension(:,:,:), allocatable :: xoh, xno3, xh2o2   !kg/kg
      real(fp), dimension(:,:,:), allocatable :: dms, so2, so4 !kg/kg
      real(fp), pointer, dimension(:,:,:) :: msa  !kg/kg
      real(fp), pointer, dimension(:,:,:) :: SU_dep  ! Sulfate Dry Deposition All Bins [kg/m2/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO2 ! vertical sum of SO2 Prod from DMS oxidation [kg/m2/s]
      real(fp), pointer, dimension(:,:) :: SU_PMSA ! vertical sum of MSA Prod from DMS oxidation [kg/m2/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO4 ! vertical sum of SO4 Prod from all SO2 oxidation [kg/m2/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO4g ! vertical sum of SO4 Prod from gaseous SO2 oxidation [kg/m2/s]
      real(fp), pointer, dimension(:,:) :: SU_PSO4aq ! vertical sum of SO4 Prod from aqueous SO2 oxidation [kg/m2/s]
      real(fp), pointer, dimension(:,:,:) :: SU_emis   ! DMS emissions in kg/m2/s
      real(fp), pointer, dimension(:,:,:) :: pso2  ! SO2 Prod from DMS oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pmsa  ! MSA Prod from DMS oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso4  ! SO4 Prod from all SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso4g  ! SO4 Prod from gaseous SO2 oxidation [kg/kg/s]
      real(fp), pointer, dimension(:,:,:) :: pso4aq  ! SO4 Prod from aqueous SO2 oxidation [kg/kg/s]
      real(fp), dimension(:,:), allocatable :: drydepfrequency
      real(fp), dimension(:,:), allocatable :: dmso_conc !DMS source concentration in ocean water [nmol/L]
      ! h2o2_init is reused from last time step
      real(fp), allocatable :: xh2o2_init_gocart(:,:,:) ! initial H2O2 from last time step
      logical :: recycle_h2o2
      !error information
      integer :: RC
      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errmsg = ''
      thisloc = ' -> at compute_gocart (in SO4chemScheme_GOCART_Mod.F90)'
      !RC = CC_SUCCESS
      rc = 0 !try not to rely on CC_SUCCESS
      !drydepf = 0.0_fp

      rad2deg = 180.0_fp/pi
      deg2rad = pi/180.0_fp
      !construct time in yyyymmdd and hhmmss formats for use in gocart
      nymd = year*10000 + month*100 + day
      nhms = hour*10000 + minute*100 + second

      !get species indices for use in gocart
      do species_idx = 1, num_species
         if (species_short_name(species_idx) == 'SO2' .or. species_short_name(species_idx) == 'so2') then
            nso2 = species_idx
         else if (species_short_name(species_idx) == 'SO4' .or. species_short_name(species_idx) == 'so4') then
            nso4 = species_idx
         else if (species_short_name(species_idx) == 'DMS' .or. species_short_name(species_idx) == 'dms') then
            ndms = species_idx
         else if (species_short_name(species_idx) == 'DMS_IN' .or. species_short_name(species_idx) == 'dms_in') then
            ndms_in = species_idx
         else if (species_short_name(species_idx) == 'MSA' .or. species_short_name(species_idx) == 'msa') then
            nmsa = species_idx
         else if (species_short_name(species_idx) == 'OH' .or. species_short_name(species_idx) == 'oh') then
            noh = species_idx
         else if (species_short_name(species_idx) == 'NO3' .or. species_short_name(species_idx) == 'no3') then
            nno3 = species_idx
         else if (species_short_name(species_idx) == 'H2O2' .or. species_short_name(species_idx) == 'h2o2') then
            nh2o2 = species_idx
         end if
      end do

      if (nso2 == -1 .or. nso4 == -1 .or. ndms == -1 .or. ndms_in == -1 .or. nmsa == -1 .or. noh == -1 .or. nno3 == -1 .or. nh2o2 == -1) then
         errmsg = 'Error in compute_gocart: SO2, SO4, DMS, DMS_IN, MSA, OH, NO3, and H2O2 must be present in species list.'
         !call CC_Error(trim(errMsg), RC, thisLoc)
         write(*,'(A)') trim(errmsg)
         return
      end if

      !allocate arrays
      allocate(oh_clim(1,1,num_layers), h2o2_clim(1,1,num_layers), no3_clim(1,1,num_layers), &
         xoh(1,1,num_layers), xno3(1,1,num_layers), xh2o2(1,1,num_layers), xh2o2_init_gocart(1,1,num_layers), &
         dms(1,1,num_layers), so2(1,1,num_layers), so4(1,1,num_layers), msa(1,1,num_layers), &
         su_dep(1, 1, num_species), su_emis(1, 1, num_species), su_pso2(1, 1), su_pmsa(1, 1), su_pso4(1, 1), su_pso4g(1, 1), su_pso4aq(1, 1), &
         pso2(1, 1, num_layers), pmsa(1, 1, num_layers), pso4(1, 1, num_layers), pso4g(1, 1, num_layers), &
         pso4aq(1, 1, num_layers), drydepfrequency(1, 1), latrad(1,1), lonrad(1,1), dmso_conc(1,1))


      !retrieve climatology fields; remember to reverse the vertical layer (TODO: double check the input files for this)
      oh_clim(1,1,:) = species_conc(num_layers:1:-1, noh) * 1.0e-6_fp !change from ppm to mol/mol.
      no3_clim(1,1,:) = species_conc(num_layers:1:-1, nno3) * 1.0e-6_fp !change from ppm to mol/mol.
      h2o2_clim(1,1,:) = species_conc(num_layers:1:-1, nh2o2) * 1.0e-6_fp !change from ppm to mol/mol.
      ! Initialize some variables for the first time
      if (firsttime) then
         ! IMPORTANT: nymd_last must NOT equal nymd_current so that the
         ! "if (nymd_last == nymd_current)" block inside SulfateUpdateOxidants
         ! does NOT fire every timestep. In the original GOCART, nymd_oxidants
         ! is initialized to -1 and never updated (the condition is never true).
         ! H2O2 is only recycled via recycle_h2o2 every 3 hours.
         nymd_last = -1
         ! First time, set initial recycle time
         nhms_last_recycle = nhms
         !allocate and initialize xh2o2_init to climatology for the first time step
         if (.not. allocated(xh2o2_init)) then
            allocate(xh2o2_init(num_layers))
         end if
         xh2o2_init = h2o2_clim(1,1,:)  ! initialize H2O2 to climatology at first time step
         firsttime = .false.
      end if

      ! Recycle H2O2 every 3 hours (matching GOCART's daily_alarm(clock,30000) behavior).
      ! Do NOT update nymd_last - it must stay at -1 to prevent SulfateUpdateOxidants
      ! from resetting xh2o2 to climatology every timestep.
      recycle_h2o2 = .false.
      if ((nhms - nhms_last_recycle >= 30000) .or. &
         (nhms < nhms_last_recycle)) then  ! handles day rollover (e.g., 230000 -> 010000)
         nhms_last_recycle = nhms
         recycle_h2o2 = .true.
      end if

      ! transform data for GOCART DryDeposition call
      call prepmetvarsforgocart(num_layers,     &
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
         u10m,            &
         v10m,            &
         z0h,             &
         gocart_tmpu,     &
         gocart_rhoa,     &
         gocart_hghte,    &
         gocart_cloud,    &
         gocart_delp,     &
         gocart_lwi,      &
         gocart_ustar,    &
         gocart_pblh,     &
         gocart_press,    &
         gocart_hflux,    &
         gocart_u10m,    &
         gocart_v10m,    &
         gocart_z0h)

      !update oxidants based on climatology and diurnal cycle
      xoh = 0.0_fp; xno3 = 0.0_fp; xh2o2_init_gocart(1,1,:)= xh2o2_init; xh2o2 = xh2o2_init_gocart
      latrad(1,1) = lat * deg2rad; lonrad(1,1) = lon * deg2rad
      call sulfateupdateoxidants(nymd, nhms, lonrad, latrad, gocart_rhoa, num_layers, tstep, nymd_last, &
         undefval, rad2deg, avo, pi, airmw, oh_clim, no3_clim, h2o2_clim, xoh, xno3, xh2o2, recycle_h2o2, rc)

      if (rc /= 0) then
         errmsg = 'Error in compute_gocart: Failed in updating oxidants in GOCART So4chem process.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(errmsg)
         return
      end if

      !get pressure lid index
      call findklid(klid, plid, gocart_press(:,:,:), rc)
      !if (RC /= CC_SUCCESS) then
      if (rc /= 0) then
         errmsg = 'Error in compute_gocart: Failed in finding pressure lid index in GOCART So4chem process.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(errmsg)
         return
      end if

      !retrieve sulfate species concentrations
      fmassmsa = species_mw_g(nmsa)
      fmassdms = species_mw_g(ndms)
      fmassso2 = species_mw_g(nso2)
      fmassso4 = species_mw_g(nso4)
      !dms(1,1,:) = species_conc(num_layers:1:-1, nDMS) * 1.0e-9_fp  !ug/kg ==> kg/kg
      dms(1,1,:) = species_conc(num_layers:1:-1, ndms) * 1.0e-6_fp * fmassdms / airmw  !ppm ==> kg/kg
      so2(1,1,:) = species_conc(num_layers:1:-1, nso2) * 1.0e-6_fp * fmassso2 / airmw  ! ppm ==> kg/kg
      so4(1,1,:) = species_conc(num_layers:1:-1, nso4) * 1.0e-9_fp  !ug/kg ==> kg/kg
      !msa(1,1,:) = species_conc(num_layers:1:-1, nMSA) * 1.0e-6_fp * fMassMSA / AIRMW  ! ppm ==> kg/kg
      msa(1,1,:) = species_conc(num_layers:1:-1, nmsa) * 1.0e-9_fp  ! ug/kg ==> kg/kg

      !run DMS emission scheme
      dmso_conc = species_conc(1, ndms_in) !in [nmol/L]. Note this is a special unit case since it is not atmospheric composition.
      su_emis = 0.0_fp
      call dmsemission (num_layers, tstep, g0, gocart_tmpu, gocart_u10m, gocart_v10m, gocart_lwi, &
         gocart_delp, fmassdms, dmso_conc, dms, su_emis, ndms, rc)
      if (rc /= 0) then
         errmsg = 'Error in compute_gocart: Failed in GOCART DMSemission.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(errmsg)
         return
      end if

      !call GOCART sulfate chemistry driver
      !force dz to be a big value (negative will not work depending on compiler) at the surface to make drydep frequency equal zero.
      !https://github.com/GEOS-ESM/GOCART/blob/9ff3df9545dd582f415f682d3297e8c6c841e5cb/Process_Library/GOCART2G_Process.F90#L3124
      !Five functions need to be customized here if we want to turn it off compleltely.
      !This is to ensure dry deposition does not run twice for SO2 and SO4
      gocart_hghte(:,:,num_layers - 1) = gocart_hghte(:,:,num_layers) + 1.0e38_fp
      call sulfatechemdriver(num_layers, klid, tstep, pi, rad2deg, von_karman, airmw, avo, cpd, g0, fmassmsa,fmassdms,fmassso2,fmassso4,&
         nymd, nhms, lonrad, latrad, dms, so2, so4, msa, ndms, nso2, nso4, nmsa, xoh, xno3, xh2o2, xh2o2_init_gocart, gocart_delp, gocart_tmpu, gocart_cloud, &
         gocart_rhoa, gocart_hghte, gocart_ustar, gocart_hflux, gocart_lwi, gocart_pblh, gocart_z0h, su_dep, su_pso2, su_pmsa, su_pso4, su_pso4g, &
         su_pso4aq, pso2, pmsa, pso4, pso4g, pso4aq, drydepfrequency, rc)

      if (rc /= 0) then
         errmsg = 'Error in compute_gocart: Failed in GOCART sulfate chemistry driver.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(errmsg)
         return
      end if

      !save H2O2 initialization for next time step
      xh2o2_init = xh2o2_init_gocart(1,1,:)

      !assign to output tendencies; remember to reverse the vertical layer back to original order
      if (params%update_so2) then !since the chem driver has drydep in it, not sure if we should update so2 chem array here.
         species_tendencies(:, nso2) = so2(1,1,num_layers:1:-1) * 1.0e6_fp * airmw / fmassso2  ! kg/kg ==> ppm
      else
         species_tendencies(:, nso2) = species_conc(:, nso2)  !keep SO2 unchanged.
      end if
      species_tendencies(:, nso4) = so4(1,1,num_layers:1:-1) * 1.0e9_fp  !kg/kg ==> ug/kg
      species_tendencies(:, nmsa) = msa(1,1,num_layers:1:-1) * 1.0e9_fp  ! kg/kg ==> ug/kg
      species_tendencies(:, ndms) = dms(1,1,num_layers:1:-1) * 1.0e6_fp * airmw / fmassdms  ! kg/kg ==> ppm
      species_tendencies(:, ndms_in) = species_conc(:, ndms_in)  !Note: DMS in ocean is unchanged since it is read in through monthly files.
      species_tendencies(:, noh) = species_conc(:, noh) !keep OH and NO3 oxidants unchanged (no cross-process consumption modeled)
      species_tendencies(:, nno3) = species_conc(:, nno3)
      !H2O2: write the post-chem (afterchem) H2O2 depleted by aqueous SO2 oxidation back into the shared
      !array so that the downstream wet-deposition process sees the same H2O2 already consumed here, as in
      !GEOS-Chem/GOCART. The host (catchem_emis_mod) re-imports the time-interpolated GMI climatology into
      !species_conc(nH2O2) at the start of the next timestep, so this per-step overwrite does NOT corrupt the
      !climatology baseline; the cross-step/3-hourly H2O2 depletion memory is carried by xh2o2_init above.
      species_tendencies(:, nh2o2) = xh2o2_init_gocart(1,1,num_layers:1:-1) * 1.0e6_fp  ! mol/mol ==> ppm

      ! Per-species-per-level diagnostic: 2D array (levels, species)
      if (present(production_rate_per_species_per_level) .and. present(diagnostic_species_id)) then
         ! Find position of this species in diagnostic_species_id array
         do diag_idx = 1, size(diagnostic_species_id)
            if (diagnostic_species_id(diag_idx) == nmsa) then
               ! Add your custom production rate (dms to so2, dms to msa, so2 to so4) per species per level calculation
               production_rate_per_species_per_level(:, diag_idx) = pmsa(1,1,num_layers:1:-1)
            end if
            if (diagnostic_species_id(diag_idx) == nso2) then
               ! Add your custom production rate (dms to so2, dms to msa, so2 to so4) per species per level calculation
               production_rate_per_species_per_level(:, diag_idx) = pso2(1,1,num_layers:1:-1)
            end if
            if (diagnostic_species_id(diag_idx) == nso4) then
               ! Add your custom production rate (dms to so2, dms to msa, so2 to so4) per species per level calculation
               production_rate_per_species_per_level(:, diag_idx) = pso4(1,1,num_layers:1:-1)
            end if
         end do
      end if

      if (present(pso4_from_gaseous_so2_per_level)) then
         pso4_from_gaseous_so2_per_level = pso4g(1,1,num_layers:1:-1)
      end if

      if (present(pso4_from_aqueous_so2_per_level)) then
         pso4_from_aqueous_so2_per_level = pso4aq(1,1,num_layers:1:-1)
      end if

      if (present(dms_emission_flux)) then
         dms_emission_flux = su_emis(1,1,ndms)
      end if


      !cleanup pointers
      if (associated(gocart_tmpu)) deallocate(gocart_tmpu); nullify(gocart_tmpu)
      if (associated(gocart_rhoa)) deallocate(gocart_rhoa); nullify(gocart_rhoa)
      if (associated(gocart_hghte)) deallocate(gocart_hghte); nullify(gocart_hghte)
      if (associated(gocart_delp)) deallocate(gocart_delp); nullify(gocart_delp)
      if (associated(gocart_cloud)) deallocate(gocart_cloud); nullify(gocart_cloud)
      if (associated(gocart_press)) deallocate(gocart_press); nullify(gocart_press)
      if (associated(gocart_lwi)) deallocate(gocart_lwi); nullify(gocart_lwi)
      if (associated(gocart_ustar)) deallocate(gocart_ustar); nullify(gocart_ustar)
      if (associated(gocart_hflux)) deallocate(gocart_hflux); nullify(gocart_hflux)
      if (associated(gocart_u10m)) deallocate(gocart_u10m); nullify(gocart_u10m)
      if (associated(gocart_v10m)) deallocate(gocart_v10m); nullify(gocart_v10m)
      if (associated(gocart_z0h)) deallocate(gocart_z0h); nullify(gocart_z0h)
      if (associated(su_dep)) deallocate(su_dep); nullify(su_dep)
      if (associated(su_pso2)) deallocate(su_pso2); nullify(su_pso2)
      if (associated(su_pmsa)) deallocate(su_pmsa); nullify(su_pmsa)
      if (associated(su_pso4)) deallocate(su_pso4); nullify(su_pso4)
      if (associated(su_pso4g)) deallocate(su_pso4g); nullify(su_pso4g)
      if (associated(su_pso4aq)) deallocate(su_pso4aq); nullify(su_pso4aq)
      if (associated(su_emis)) deallocate(su_emis); nullify(su_emis)
      if (associated(pso2)) deallocate(pso2); nullify(pso2)
      if (associated(pmsa)) deallocate(pmsa); nullify(pmsa)
      if (associated(pso4)) deallocate(pso4); nullify(pso4)
      if (associated(pso4g)) deallocate(pso4g); nullify(pso4g)
      if (associated(pso4aq)) deallocate(pso4aq); nullify(pso4aq)
      if (associated(msa)) deallocate(msa); nullify(msa)
      if (associated(oh_clim)) deallocate(oh_clim); nullify(oh_clim)
      if (associated(no3_clim)) deallocate(no3_clim); nullify(no3_clim)
      if (associated(h2o2_clim)) deallocate(h2o2_clim); nullify(h2o2_clim)
      !cleanup array allocations
      deallocate( xoh, xno3, xh2o2, xh2o2_init_gocart, dms, so2, so4, drydepfrequency, latrad, lonrad, dmso_conc)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines

   subroutine prepmetvarsforgocart(km,              &
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
      u10m,            &
      v10m,            &
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
      GOCART_U10M,    &
      GOCART_V10M,    &
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
      REAL(fp),  intent(in), target              :: u10m                                  ! 10m wind speed [m/sec]
      REAL(fp),  intent(in), target              :: v10m                                  ! 10m wind speed [m/sec]
      REAL(fp),  intent(in), target              :: z0h                                   ! rough height, sens. heat [m]

      ! INPUT/OUTPUTS
      REAL(fp), intent(inout), pointer :: GOCART_TMPU(:,:,:)
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_cloud
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_DELP
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_PRESS
      real(fp), intent(inout), pointer :: GOCART_LWI(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_USTAR(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_PBLH(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_HFLUX(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_U10M(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_V10M(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_Z0H(:,:)

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(gocart_tmpu(1, 1, km))
      allocate(gocart_rhoa(1, 1, km))
      allocate(gocart_hghte(1, 1, 0:km))
      allocate(gocart_cloud(1, 1, km))
      allocate(gocart_delp(1, 1, km))
      allocate(gocart_press(1, 1, km))
      allocate(gocart_lwi(1, 1))
      allocate(gocart_ustar(1, 1))
      allocate(gocart_pblh(1, 1))
      allocate(gocart_hflux(1, 1))
      allocate(gocart_u10m(1, 1))
      allocate(gocart_v10m(1, 1))
      allocate(gocart_z0h(1, 1))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)

      gocart_tmpu(1,1,:) = tmpu(size(tmpu):1:-1) ! temperature [K]
      gocart_rhoa(1,1,:) = rhoa(size(rhoa):1:-1) ! air density [kg/m^3]
      gocart_cloud(1,1,:) = cldfrc(size(cldfrc):1:-1) ! cloud fraction [1]
      gocart_delp(1,1,:) = delp(size(delp):1:-1) ! pressure thickness [Pa]
      gocart_hghte(1,1,:) = hghte(size(hghte):1:-1)    ! top of layer geopotential height [m]
      gocart_press(1,1,:) = pmid(size(pmid):1:-1)    ! pressure at mid-layer [Pa]
      gocart_lwi = real(lwi, fp)     ! orography flag; Land, ocean, ice mask
      gocart_ustar  = ustar

      ! friction speed [m/sec]
      gocart_pblh   = pblh      ! PBL height [m]
      gocart_hflux = hflux     ! sfc. sens. heat flux [W m-2]
      gocart_u10m = u10m       ! 10m wind speed [m/sec]
      gocart_v10m = v10m       ! 10m wind speed [m/sec]
      gocart_z0h    = z0h       ! rough height, sens. heat [m]


   end subroutine prepmetvarsforgocart

   subroutine findklid (klid, plid, ple, rc)

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
      refdiff = 150000.0_fp
      do k = 1, ubound(ple,3)
         diff = abs(pres(k) - plid_)
         if (diff < refdiff) then
            klid = k
            refdiff = diff
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

   end subroutine findklid

end module so4chemscheme_gocart_mod
```


