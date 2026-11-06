

# File WetDepScheme\_GOCART\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**schemes**](dir_8ca87c5e2f5cf830ab1a41055168a46b.md) **>** [**WetDepScheme\_GOCART\_Mod.F90**](_wet_dep_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the documentation of this file](_wet_dep_scheme___g_o_c_a_r_t___mod_8_f90.md)


```Fortran

module wetdepscheme_gocart_mod

   use catchem_bridge_precision, only: fp, rae
   use catchem_bridge_error, only: cc_warning, cc_success
   use wetdepcommon_mod, only: wetdepschemegocartconfig
   use catchem_bridge_constants, only: g0, airmw  ! gravity and dry-air molecular weight for GOCART calls / unit conversion
   use gocart2g_process, only: su_wet_removal, wetremovalufs

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: plid = 0.01_fp    ! Pressure lid [hPa]

contains

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
      species_is_wetdep, &
      species_short_name, &
      species_wd_rainouteff, &
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
      logical, intent(in) :: species_is_wetdep(:)   ! Species is_wetdep participation flag
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      real(fp), intent(in) :: species_wd_rainouteff(:,:)  ! Species wd_rainouteff property (num_species,3)
      real(fp), intent(in) :: species_radius(:)  ! Species radius property [um]
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
      real(fp) :: rainout_eff(3)
      character(len=512) :: thisLoc, errMsg

      ! Local bin indices used INSIDE SU_Wet_Removal to index its DC/fluxout arrays (not species-array indices)
      integer, parameter :: iDMS = 1, iso2 = 2, iso4 = 3, imsa = 4, nbins_su = 4

      rc = 0
      removed_g(:,:) = 0.0_fp

      ! Replacement mode -> every species keeps its current value unless overwritten below.
      ! This echo is also the participation gate: species with is_wetdep = .false. never
      ! reach an overwrite branch (see route gates), so they return bit-exactly unchanged
      ! (same contract as the Jacob kernel).
      species_tendencies(:,:) = species_conc(:,:)

      ! Identify sulfate group + H2O2 by short name.
      ndms = -1; nso2 = -1; nso4 = -1; nmsa = -1; nh2o2 = -1
      do species_idx = 1, num_species
         select case (trim(species_short_name(species_idx)))
          case ('SO2', 'so2')
            nso2 = species_idx
          case ('SO4', 'so4')
            nso4 = species_idx
          case ('DMS', 'dms')
            ndms = species_idx
          case ('MSA', 'msa')
            nmsa = species_idx
          case ('H2O2', 'h2o2')
            nh2o2 = species_idx
         end select
      end do

      ! Prepare column meteorology in GOCART top-to-bottom order.
      allocate(gocart_tmpu(1,1,num_layers), gocart_rhoa(1,1,num_layers))
      allocate(gocart_ple(1,1,0:num_layers), gocart_pfllsan(1,1,0:num_layers), gocart_pfilsan(1,1,0:num_layers))
      allocate(gocart_precc(1,1), gocart_precl(1,1))

      gocart_tmpu(1,1,:)    = t(num_layers:1:-1)
      gocart_rhoa(1,1,:)    = mairden(num_layers:1:-1)
      gocart_ple(1,1,:)     = pedge(num_layers+1:1:-1)
      gocart_pfllsan(1,1,:) = pfllsan(num_layers+1:1:-1)
      gocart_pfilsan(1,1,:) = pfilsan(num_layers+1:1:-1)

      ! GOCART only uses precc/precl for a "> 0" test; divide by tstep to reproduce the UFS host operation
      ! (rain amount [m] -> rain rate [m/s]). CATChem already separates large-scale vs convective components.
      gocart_precc(1,1) = preccon / tstep
      gocart_precl(1,1) = preclsc / tstep

      ! Layer pressure thickness and mid-layer pressure from edge pressure (top-first: pressure increases downward)
      do k = 1, num_layers
         gocart_delp(1,1,k)  = gocart_ple(1,1,k) - gocart_ple(1,1,k-1)
         gocart_press(1,1,k) = 0.5_fp * (gocart_ple(1,1,k-1) + gocart_ple(1,1,k))
      end do

      ! Pressure-lid index (level nearest plid), matching the other GOCART schemes.
      call findklid(klid, plid, gocart_press, rc)

      nullify(msa, msa0, su_flux, su_pso4, su_pso4wet, su_pso4col, su_pso4wetcol, ufs_flux)

      ! --------------------------------------------------------------------------------------------------
      ! (1) Sulfate group -> SU_Wet_Removal (single call for DMS/SO2/SO4/MSA); requires H2O2.
      ! DMS/MSA are optional: SU_Wet_Removal never removes DMS (DC(nDMS)=0), so DMS need not be carried in
      ! the wetdep list. Requiring it would force is_wetdep on DMS and alter the Jacob scheme's behavior.
      ! Participation gate: the coupled route runs only when the core sulfur species are enrolled in
      ! wet deposition; optional partners (DMS/MSA) participate only when enrolled.
      ! H2O2 is a hard input for the coupled solver (it caps soluble SO2).  When the sulfur pair is
      ! enrolled but H2O2 is absent, the route is skipped with a loud CC_Warning rather than silently
      ! producing zero, so a mis-configured mechanism stays observable in the log (this mirrors the
      ! Jacob kernel's missing-species warnings).  SO2/SO4 echo unchanged in that case.
      ! --------------------------------------------------------------------------------------------------
      if (nso2 > 0 .and. nso4 > 0 .and. species_is_wetdep(nso2) .and. species_is_wetdep(nso4) .and. &
         nh2o2 <= 0) then
         thisloc = 'WetDepScheme_GOCART_Mod.compute_gocart>su_wet_removal'
         errmsg = 'H2O2 is not a chemical species in the model. GOCART sulfate wet removal is skipped ' // &
            'for SO2/SO4 this step; species echoed unchanged (no silent zero).'
         rc = cc_success
         CALL cc_warning(errmsg, rc, thisloc)
      end if
      if (nso2 > 0 .and. nso4 > 0 .and. nh2o2 > 0 .and. &
         species_is_wetdep(nso2) .and. species_is_wetdep(nso4)) then
         fmassso2 = species_mw_g(nso2)
         fmassso4 = species_mw_g(nso4)
         if (ndms > 0) then
            fmassdms = species_mw_g(ndms)
         else
            fmassdms = 62.0_fp  ! GOCART DMS molar mass [g/mol]; DMS array stays zero (never removed)
         end if

         ! Convert to GOCART units (kg/kg mass mixing ratio), top-first.
         if (ndms > 0 .and. species_is_wetdep(ndms)) then
            dms(1,1,:) = species_conc(num_layers:1:-1, ndms) * 1.0e-6_fp * fmassdms / airmw  ! ppm -> kg/kg
         else
            dms(1,1,:) = 0.0_fp
         end if
         so2(1,1,:) = species_conc(num_layers:1:-1, nso2) * 1.0e-6_fp * fmassso2 / airmw  ! ppm -> kg/kg
         so4(1,1,:) = species_conc(num_layers:1:-1, nso4) * 1.0e-9_fp                     ! ug/kg -> kg/kg
         ! H2O2 read fresh from the shared array each step (so4chem runs before wetdep and writes the
         ! post-chemistry H2O2 here); used only as vmr auxiliary to cap soluble SO2, not written back.
         h2o2_int(1,1,:) = species_conc(num_layers:1:-1, nh2o2) * 1.0e-6_fp               ! ppm -> mol/mol (vmr)

         dms0 = dms; so20 = so2; so40 = so4  ! keep pre-removal copies for diagnostics

         if (nmsa > 0 .and. species_is_wetdep(nmsa)) then
            allocate(msa(1,1,num_layers), msa0(1,1,num_layers))
            msa(1,1,:) = species_conc(num_layers:1:-1, nmsa) * 1.0e-9_fp  ! ug/kg -> kg/kg
            msa0 = msa
         end if

         call su_wet_removal(num_layers, nbins_su, klid, tstep, .true., g0, airmw, &
            gocart_delp, fmassso4, fmassso2, &
            h2o2_int, gocart_ple, gocart_rhoa, gocart_precc, gocart_precl, &
            gocart_pfllsan, gocart_pfilsan, gocart_tmpu, &
            idms, iso2, iso4, imsa, dms, so2, so4, msa, &
            su_flux, su_pso4col, su_pso4wetcol, su_pso4, su_pso4wet, rc)

         ! Convert back to CATChem units (replacement mode); record kg/kg removed for diagnostics.
         species_tendencies(:, nso2) = so2(1,1,num_layers:1:-1) * 1.0e6_fp * airmw / fmassso2  ! kg/kg -> ppm
         species_tendencies(:, nso4) = so4(1,1,num_layers:1:-1) * 1.0e9_fp                     ! kg/kg -> ug/kg
         removed_g(:, nso2) = so20(1,1,:) - so2(1,1,:)
         removed_g(:, nso4) = so40(1,1,:) - so4(1,1,:)
         if (ndms > 0 .and. species_is_wetdep(ndms)) then
            species_tendencies(:, ndms) = dms(1,1,num_layers:1:-1) * 1.0e6_fp * airmw / fmassdms  ! kg/kg -> ppm
            removed_g(:, ndms) = dms0(1,1,:) - dms(1,1,:)
         end if
         if (nmsa > 0 .and. species_is_wetdep(nmsa)) then
            species_tendencies(:, nmsa) = msa(1,1,num_layers:1:-1) * 1.0e9_fp  ! kg/kg -> ug/kg
            removed_g(:, nmsa) = msa0(1,1,:) - msa(1,1,:)
         end if
         ! H2O2 stays unchanged (auxiliary input only); species_tendencies(:,nH2O2) already = species_conc.
      end if

      ! --------------------------------------------------------------------------------------------------
      ! (2) All other aerosol species -> WetRemovalUFS (one call each). Gases (incl. H2O2) are skipped.
      ! Participation gate: species not enrolled in wet deposition keep their echoed input values.
      ! --------------------------------------------------------------------------------------------------
      do species_idx = 1, num_species
         if (species_idx == nso2 .or. species_idx == nso4 .or. species_idx == ndms .or. &
            species_idx == nmsa .or. species_idx == nh2o2) cycle
         if (.not. species_is_aerosol(species_idx)) cycle  ! WetRemovalUFS supports aerosols (and NH3) only
         if (.not. species_is_wetdep(species_idx)) cycle   ! participation gate (is_wetdep flag)

         aerosol(1,1,:) = species_conc(num_layers:1:-1, species_idx) * 1.0e-9_fp  ! ug/kg -> kg/kg
         aerosol0 = aerosol
         rainout_eff = species_wd_rainouteff(species_idx, :)

         call wetremovalufs(num_layers, klid, 1, tstep, trim(species_short_name(species_idx)), .true., g0, &
            species_radius(species_idx), rainout_eff, &
            params%washout_tuning, params%radius_threshold, aerosol, &
            gocart_ple, gocart_tmpu, gocart_rhoa, gocart_pfllsan, gocart_pfilsan, ufs_flux, rc)

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
               mass_col(k) = removed_g(k, species_idx) * gocart_delp(1,1,k) / g0
            end do
            if (present(wetdep_mass_per_species_per_level)) &
               wetdep_mass_per_species_per_level(:, diag_idx) = mass_col(num_layers:1:-1)
            if (present(wetdep_flux_per_species_per_level)) &
               wetdep_flux_per_species_per_level(:, diag_idx) = mass_col(num_layers:1:-1) / tstep
         end do
      end if

      ! Cleanup
      if (associated(gocart_tmpu))    deallocate(gocart_tmpu);    nullify(gocart_tmpu)
      if (associated(gocart_rhoa))    deallocate(gocart_rhoa);    nullify(gocart_rhoa)
      if (associated(gocart_ple))     deallocate(gocart_ple);     nullify(gocart_ple)
      if (associated(gocart_pfllsan)) deallocate(gocart_pfllsan); nullify(gocart_pfllsan)
      if (associated(gocart_pfilsan)) deallocate(gocart_pfilsan); nullify(gocart_pfilsan)
      if (associated(gocart_precc))   deallocate(gocart_precc);   nullify(gocart_precc)
      if (associated(gocart_precl))   deallocate(gocart_precl);   nullify(gocart_precl)
      if (associated(msa))            deallocate(msa);            nullify(msa)
      if (associated(msa0))           deallocate(msa0);           nullify(msa0)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================

   subroutine findklid(klid, plid, ple, rc)

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
            if (.not. rae(pres(klid), ple(i,j,klid))) then
               rc = 1
               return
            end if
         end do
      end do

   end subroutine findklid

end module wetdepscheme_gocart_mod
```


