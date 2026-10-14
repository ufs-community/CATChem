!> \file WetDepScheme_GOCART_Mod.F90
!! \brief GOCART2G wet removal scheme: SU_Wet_Removal for sulfate species (DMS/SO2/SO4/MSA) and WetRemovalUFS for all other species
!!
!! Pure science kernel for gocart scheme in wetdep process.
!! This module wraps the GOCART2G process-library wet-removal routines. Rather than
!! re-implementing the science, it prepares the column meteorology in the layout the
!! GOCART routines expect (top-to-bottom ordering, (1,1,k) shaped pointers) and calls
!! them directly:
!!   - SU_Wet_Removal : the sulfate group (DMS/SO2/SO4/MSA) in a single call.
!!   - WetRemovalUFS  : every other aerosol species, one call per species.
!!
!! Author: Wei Li
!! Reference: GOCART2G Process Library (GEOS-ESM/GOCART): SU_Wet_Removal and WetRemovalUFS (Jacob et al. [2000] Harvard scheme)
module WetDepScheme_GOCART_Mod

   use precision_mod, only: fp, rae
   use WetDepCommon_Mod, only: WetDepSchemeGOCARTConfig
   use Constants, only: g0, AIRMW  ! gravity and dry-air molecular weight for GOCART calls / unit conversion
   use GOCART2G_Process, only: SU_Wet_Removal, WetRemovalUFS

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: plid = 0.01_fp    ! Pressure lid [hPa]

contains

   !> Science computation for gocart wet-removal scheme.
   !!
   !! Routes sulfate species (SO4/SO2/DMS/MSA) through GOCART's SU_Wet_Removal and all other
   !! aerosol species through GOCART's WetRemovalUFS. H2O2 is used only as an auxiliary input
   !! to SU_Wet_Removal (to cap soluble SO2) and is not itself wet-deposited.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  mairden    MAIRDEN field (moist air density) [kg/m3]
   !! @param[in]  pedge    PEDGE field (level-edge air pressure) [Pa] (num_layers+1)
   !! @param[in]  pfilsan    PFILSAN field (3D flux of ice nonconvective precip) [kg/m2/s] (num_layers+1)
   !! @param[in]  pfllsan    PFLLSAN field (3D flux of liquid nonconvective precip) [kg/m2/s] (num_layers+1)
   !! @param[in]  preccon    PRECCON field (convective precip at ground) [m] - divided by tstep -> [m/s]
   !! @param[in]  preclsc    PRECLSC field (large-scale precip at ground) [m] - divided by tstep -> [m/s]
   !! @param[in]  t    T field (temperature) [K]
   !! @param[in]  tstep    Time step [s] - retrieved from process interface
   !! @param[in]  species_is_aerosol    Species is_aerosol property
   !! @param[in]  species_short_name    Species short_name property
   !! @param[in]  species_henry_cr    Species henry_cr property (unused; GOCART routines carry their own)
   !! @param[in]  species_henry_k0    Species henry_k0 property (unused)
   !! @param[in]  species_henry_pKa    Species henry_pKa property (unused)
   !! @param[in]  species_wd_retfactor    Species wd_retfactor property (unused)
   !! @param[in]  species_wd_LiqAndGas    Species wd_LiqAndGas property (unused)
   !! @param[in]  species_wd_convfacI2G    Species wd_convfacI2G property (unused)
   !! @param[in]  species_wd_rainouteff    Species wd_rainouteff property (3 temperature-dependent efficiencies) (num_species,3)
   !! @param[in]  species_wd_reevap_frac    Species wd_reevap_frac property (unused)
   !! @param[in]  species_radius    Species radius property [um]
   !! @param[in]  species_mw_g    Species mw_g property
   !! @param[in]  species_conc   Species concentrations [ppm for gases, ug/kg for aerosols] (num_layers, num_species)
   !! @param[inout] species_tendencies  Updated species concentrations (replacement mode) (num_layers, num_species)
   !! @param[inout] wetdep_mass_per_species_per_level    Wet deposition mass loss per species per level [kg/m2] (num_layers, num_species)
   !! @param[inout] wetdep_flux_per_species_per_level    Wet deposition flux per species per level [kg/m2/s] (num_layers, num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   subroutine compute_gocart( &
      num_layers, &
      num_species, &
      params, &
      mairden, &
      pedge, &
      pfilsan, &
      pfllsan, &
      preccon, &
      preclsc, &
      t, &
      tstep, &
      species_is_aerosol, &
      species_short_name, &
      species_henry_cr, &
      species_henry_k0, &
      species_henry_pKa, &
      species_wd_retfactor, &
      species_wd_LiqAndGas, &
      species_wd_convfacI2G, &
      species_wd_rainouteff, &
      species_wd_reevap_frac, &
      species_radius, &
      species_mw_g, &
      species_conc, &
      species_tendencies, &
      wetdep_mass_per_species_per_level, &
      wetdep_flux_per_species_per_level, &
      diagnostic_species_id &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(WetDepSchemeGOCARTConfig), intent(in) :: params
      real(fp), intent(in) :: mairden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: pedge(num_layers+1)  ! Edge field - requires nz+1 dimensions
      real(fp), intent(in) :: pfilsan(num_layers+1)  ! Edge field - requires nz+1 dimensions
      real(fp), intent(in) :: pfllsan(num_layers+1)  ! Edge field - requires nz+1 dimensions
      real(fp), intent(in) :: preccon  ! Surface field - scalar
      real(fp), intent(in) :: preclsc  ! Surface field - scalar
      real(fp), intent(in) :: t(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      logical, intent(in) :: species_is_aerosol(:)  ! Species is_aerosol property
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      real(fp), intent(in) :: species_henry_cr(:)  ! Species henry_cr property
      real(fp), intent(in) :: species_henry_k0(:)  ! Species henry_k0 property
      real(fp), intent(in) :: species_henry_pKa(:)  ! Species henry_pKa property
      real(fp), intent(in) :: species_wd_retfactor(:)  ! Species wd_retfactor property
      logical, intent(in) :: species_wd_LiqAndGas(:)  ! Species wd_LiqAndGas property
      real(fp), intent(in) :: species_wd_convfacI2G(:)  ! Species wd_convfacI2G property
      real(fp), intent(in) :: species_wd_rainouteff(:,:)  ! Species wd_rainouteff property (num_species,3)
      real(fp), intent(in) :: species_wd_reevap_frac(:)  ! Species wd_reevap_frac property
      real(fp), intent(in) :: species_radius(:)  ! Species radius property
      real(fp), intent(in) :: species_mw_g(:)  ! Species mw_g property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: wetdep_mass_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: wetdep_flux_per_species_per_level(:,:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! GOCART expects top-to-bottom ordering with (1,1,k) shaped POINTER arrays (i=j=1 column).
      real(fp), pointer :: GOCART_tmpu(:,:,:)     ! temperature [K]
      real(fp), pointer :: GOCART_rhoa(:,:,:)     ! moist air density [kg/m3]
      real(fp), pointer :: GOCART_ple(:,:,:)      ! level-edge pressure [Pa] (0:km)
      real(fp), pointer :: GOCART_pfllsan(:,:,:)  ! liquid nonconvective precip flux [kg/m2/s] (0:km)
      real(fp), pointer :: GOCART_pfilsan(:,:,:)  ! ice nonconvective precip flux [kg/m2/s] (0:km)
      real(fp), pointer :: GOCART_precc(:,:)      ! convective precip rate [m/s]
      real(fp), pointer :: GOCART_precl(:,:)      ! large-scale precip rate [m/s]

      ! Sulfate working arrays (kg/kg, top-first) and derived quantities
      real(fp) :: GOCART_delp(1,1,num_layers)   ! pressure thickness [Pa], derived from PEDGE
      real(fp) :: GOCART_press(1,1,num_layers)  ! mid-layer pressure [Pa], for the pressure-lid search
      real(fp) :: dms(1,1,num_layers), so2(1,1,num_layers), so4(1,1,num_layers)
      real(fp) :: dms0(1,1,num_layers), so20(1,1,num_layers), so40(1,1,num_layers)
      real(fp) :: h2o2_int(1,1,num_layers)
      real(fp) :: aerosol(1,1,num_layers), aerosol0(1,1,num_layers)
      real(fp), pointer :: msa(:,:,:), msa0(:,:,:)

      ! GOCART diagnostic pointers (unused here -> left disassociated; GOCART guards with associated())
      real(fp), pointer :: su_flux(:,:,:), su_pso4(:,:,:), su_pso4wet(:,:,:)
      real(fp), pointer :: su_pso4col(:,:), su_pso4wetcol(:,:)
      real(fp), pointer :: ufs_flux(:,:,:)

      ! Per-species removed mass mixing ratio [kg/kg] in GOCART (top-first) order, for diagnostics
      real(fp) :: removed_g(num_layers, num_species)
      real(fp) :: mass_col(num_layers)

      integer :: species_idx, k, diag_idx
      integer :: nDMS, nSO2, nSO4, nMSA, nH2O2  ! species-array indices (by name)
      integer :: klid, rc
      real(fp) :: fMassSO4, fMassSO2, fMassDMS

      ! Local bin indices used INSIDE SU_Wet_Removal to index its DC/fluxout arrays (not species-array indices)
      integer, parameter :: iDMS = 1, iSO2 = 2, iSO4 = 3, iMSA = 4, nbins_su = 4

      rc = 0
      removed_g(:,:) = 0.0_fp

      ! Replacement mode -> every species keeps its current value unless overwritten below.
      species_tendencies(:,:) = species_conc(:,:)

      ! Identify sulfate group + H2O2 by short name.
      nDMS = -1; nSO2 = -1; nSO4 = -1; nMSA = -1; nH2O2 = -1
      do species_idx = 1, num_species
         select case (trim(species_short_name(species_idx)))
          case ('SO2', 'so2')
            nSO2 = species_idx
          case ('SO4', 'so4')
            nSO4 = species_idx
          case ('DMS', 'dms')
            nDMS = species_idx
          case ('MSA', 'msa')
            nMSA = species_idx
          case ('H2O2', 'h2o2')
            nH2O2 = species_idx
         end select
      end do

      ! Prepare column meteorology in GOCART top-to-bottom order.
      allocate(GOCART_tmpu(1,1,num_layers), GOCART_rhoa(1,1,num_layers))
      allocate(GOCART_ple(1,1,0:num_layers), GOCART_pfllsan(1,1,0:num_layers), GOCART_pfilsan(1,1,0:num_layers))
      allocate(GOCART_precc(1,1), GOCART_precl(1,1))

      GOCART_tmpu(1,1,:)    = t(num_layers:1:-1)
      GOCART_rhoa(1,1,:)    = mairden(num_layers:1:-1)
      GOCART_ple(1,1,:)     = pedge(num_layers+1:1:-1)
      GOCART_pfllsan(1,1,:) = pfllsan(num_layers+1:1:-1)
      GOCART_pfilsan(1,1,:) = pfilsan(num_layers+1:1:-1)

      ! GOCART only uses precc/precl for a "> 0" test; divide by tstep to reproduce the UFS host operation
      ! (rain amount [m] -> rain rate [m/s]). CATChem already separates large-scale vs convective components.
      GOCART_precc(1,1) = preccon / tstep
      GOCART_precl(1,1) = preclsc / tstep

      ! Layer pressure thickness and mid-layer pressure from edge pressure (top-first: pressure increases downward)
      do k = 1, num_layers
         GOCART_delp(1,1,k)  = GOCART_ple(1,1,k) - GOCART_ple(1,1,k-1)
         GOCART_press(1,1,k) = 0.5_fp * (GOCART_ple(1,1,k-1) + GOCART_ple(1,1,k))
      end do

      ! Pressure-lid index (level nearest plid), matching the other GOCART schemes.
      call findKlid(klid, plid, GOCART_press, rc)

      nullify(msa, msa0, su_flux, su_pso4, su_pso4wet, su_pso4col, su_pso4wetcol, ufs_flux)

      ! --------------------------------------------------------------------------------------------------
      ! (1) Sulfate group -> SU_Wet_Removal (single call for DMS/SO2/SO4/MSA); requires H2O2.
      ! DMS/MSA are optional: SU_Wet_Removal never removes DMS (DC(nDMS)=0), so DMS need not be carried in
      ! the wetdep list. Requiring it would force is_wetdep on DMS and alter the Jacob scheme's behavior.
      ! --------------------------------------------------------------------------------------------------
      if (nSO2 > 0 .and. nSO4 > 0 .and. nH2O2 > 0) then
         fMassSO2 = species_mw_g(nSO2)
         fMassSO4 = species_mw_g(nSO4)
         if (nDMS > 0) then
            fMassDMS = species_mw_g(nDMS)
         else
            fMassDMS = 62.0_fp  ! GOCART DMS molar mass [g/mol]; DMS array stays zero (never removed)
         end if

         ! Convert to GOCART units (kg/kg mass mixing ratio), top-first.
         if (nDMS > 0) then
            dms(1,1,:) = species_conc(num_layers:1:-1, nDMS) * 1.0e-6_fp * fMassDMS / AIRMW  ! ppm -> kg/kg
         else
            dms(1,1,:) = 0.0_fp
         end if
         so2(1,1,:) = species_conc(num_layers:1:-1, nSO2) * 1.0e-6_fp * fMassSO2 / AIRMW  ! ppm -> kg/kg
         so4(1,1,:) = species_conc(num_layers:1:-1, nSO4) * 1.0e-9_fp                     ! ug/kg -> kg/kg
         ! H2O2 read fresh from the shared array each step (so4chem runs before wetdep and writes the
         ! post-chemistry H2O2 here); used only as vmr auxiliary to cap soluble SO2, not written back.
         h2o2_int(1,1,:) = species_conc(num_layers:1:-1, nH2O2) * 1.0e-6_fp               ! ppm -> mol/mol (vmr)

         dms0 = dms; so20 = so2; so40 = so4  ! keep pre-removal copies for diagnostics

         if (nMSA > 0) then
            allocate(msa(1,1,num_layers), msa0(1,1,num_layers))
            msa(1,1,:) = species_conc(num_layers:1:-1, nMSA) * 1.0e-9_fp  ! ug/kg -> kg/kg
            msa0 = msa
         end if

         call SU_Wet_Removal(num_layers, nbins_su, klid, tstep, .true., g0, AIRMW, &
            GOCART_delp, fMassSO4, fMassSO2, &
            h2o2_int, GOCART_ple, GOCART_rhoa, GOCART_precc, GOCART_precl, &
            GOCART_pfllsan, GOCART_pfilsan, GOCART_tmpu, &
            iDMS, iSO2, iSO4, iMSA, dms, so2, so4, msa, &
            su_flux, su_pso4col, su_pso4wetcol, su_pso4, su_pso4wet, rc)

         ! Convert back to CATChem units (replacement mode); record kg/kg removed for diagnostics.
         species_tendencies(:, nSO2) = so2(1,1,num_layers:1:-1) * 1.0e6_fp * AIRMW / fMassSO2  ! kg/kg -> ppm
         species_tendencies(:, nSO4) = so4(1,1,num_layers:1:-1) * 1.0e9_fp                     ! kg/kg -> ug/kg
         removed_g(:, nSO2) = so20(1,1,:) - so2(1,1,:)
         removed_g(:, nSO4) = so40(1,1,:) - so4(1,1,:)
         if (nDMS > 0) then
            species_tendencies(:, nDMS) = dms(1,1,num_layers:1:-1) * 1.0e6_fp * AIRMW / fMassDMS  ! kg/kg -> ppm
            removed_g(:, nDMS) = dms0(1,1,:) - dms(1,1,:)
         end if
         if (nMSA > 0) then
            species_tendencies(:, nMSA) = msa(1,1,num_layers:1:-1) * 1.0e9_fp  ! kg/kg -> ug/kg
            removed_g(:, nMSA) = msa0(1,1,:) - msa(1,1,:)
         end if
         ! H2O2 stays unchanged (auxiliary input only); species_tendencies(:,nH2O2) already = species_conc.
      end if

      ! --------------------------------------------------------------------------------------------------
      ! (2) All other aerosol species -> WetRemovalUFS (one call each). Gases (incl. H2O2) are skipped.
      ! --------------------------------------------------------------------------------------------------
      do species_idx = 1, num_species
         if (species_idx == nSO2 .or. species_idx == nSO4 .or. species_idx == nDMS .or. &
            species_idx == nMSA .or. species_idx == nH2O2) cycle
         if (.not. species_is_aerosol(species_idx)) cycle  ! WetRemovalUFS supports aerosols (and NH3) only

         aerosol(1,1,:) = species_conc(num_layers:1:-1, species_idx) * 1.0e-9_fp  ! ug/kg -> kg/kg
         aerosol0 = aerosol

         call WetRemovalUFS(num_layers, klid, 1, tstep, trim(species_short_name(species_idx)), .true., g0, &
            species_radius(species_idx), species_wd_rainouteff(species_idx, :), &
            params%washout_tuning, params%radius_threshold, aerosol, &
            GOCART_ple, GOCART_tmpu, GOCART_rhoa, GOCART_pfllsan, GOCART_pfilsan, ufs_flux, rc)

         species_tendencies(:, species_idx) = aerosol(1,1,num_layers:1:-1) * 1.0e9_fp  ! kg/kg -> ug/kg
         removed_g(:, species_idx) = aerosol0(1,1,:) - aerosol(1,1,:)
      end do

      ! --------------------------------------------------------------------------------------------------
      ! Per-level wet-deposition diagnostics: mass = removed[kg/kg] * delp/g0 [kg/m2]; flux = mass/tstep.
      ! Results are written back in CATChem (surface-first) order.
      ! --------------------------------------------------------------------------------------------------
      if (present(diagnostic_species_id) .and. &
         (present(wetdep_mass_per_species_per_level) .or. present(wetdep_flux_per_species_per_level))) then
         do diag_idx = 1, size(diagnostic_species_id)
            species_idx = diagnostic_species_id(diag_idx)
            if (species_idx < 1 .or. species_idx > num_species) cycle
            do k = 1, num_layers
               mass_col(k) = removed_g(k, species_idx) * GOCART_delp(1,1,k) / g0
            end do
            if (present(wetdep_mass_per_species_per_level)) &
               wetdep_mass_per_species_per_level(:, diag_idx) = mass_col(num_layers:1:-1)
            if (present(wetdep_flux_per_species_per_level)) &
               wetdep_flux_per_species_per_level(:, diag_idx) = mass_col(num_layers:1:-1) / tstep
         end do
      end if

      ! Cleanup
      if (associated(GOCART_tmpu))    deallocate(GOCART_tmpu);    nullify(GOCART_tmpu)
      if (associated(GOCART_rhoa))    deallocate(GOCART_rhoa);    nullify(GOCART_rhoa)
      if (associated(GOCART_ple))     deallocate(GOCART_ple);     nullify(GOCART_ple)
      if (associated(GOCART_pfllsan)) deallocate(GOCART_pfllsan); nullify(GOCART_pfllsan)
      if (associated(GOCART_pfilsan)) deallocate(GOCART_pfilsan); nullify(GOCART_pfilsan)
      if (associated(GOCART_precc))   deallocate(GOCART_precc);   nullify(GOCART_precc)
      if (associated(GOCART_precl))   deallocate(GOCART_precl);   nullify(GOCART_precl)
      if (associated(msa))            deallocate(msa);            nullify(msa)
      if (associated(msa0))           deallocate(msa0);           nullify(msa0)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================

   !>
   !! \brief findKlid - Finds corresponding vertical index for defined pressure lid
   !!
   !! \param [INOUT] klid
   !! \param [IN] plid
   !! \param [IN] ple
   !! \param [OUT] rc
   !!!>
   subroutine findKlid(klid, plid, ple, rc)

      implicit none
      ! !INPUT PARAMETERS:
      integer, intent(inout) :: klid ! index for pressure lid
      real(fp), intent(in) :: plid ! pressure lid [hPa]; default is 0.01 hPa
      real(fp), dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]
      ! !OUTPUT PARAMETERS:
      integer, intent(out) :: rc ! return code; 0 - all is good; 1 - bad
      ! !Reference to gocart: https://github.com/GEOS-ESM/GOCART/blob/9ff3df9545dd582f415f682d3297e8c6c841e5cb/ESMF/Shared/Chem_AeroGeneric.F90#L316
      ! !Local Variables
      integer :: k, j, i
      real(fp) :: plid_, diff, refDiff
      real(fp), allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]

      klid = 1
      rc = 0

      !  convert from hPa to Pa
      plid_ = plid*100.0_fp

      allocate(pres(ubound(ple,3)))

      !  find pressure at each model level
      do k = 1, ubound(ple,3)
         pres(k) = ple(1,1,k)
      end do

      !  find smallest absolute difference between plid and pressure at each model level
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
            if (.not. rae(pres(klid), ple(i,j,klid))) then
               rc = 1
               return
            end if
         end do
      end do

   end subroutine findKlid

end module WetDepScheme_GOCART_Mod
