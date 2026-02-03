

# File WetDepScheme\_JACOB\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md) **>** [**schemes**](dir_8ca87c5e2f5cf830ab1a41055168a46b.md) **>** [**WetDepScheme\_JACOB\_Mod.F90**](_wet_dep_scheme___j_a_c_o_b___mod_8_f90.md)

[Go to the documentation of this file](_wet_dep_scheme___j_a_c_o_b___mod_8_f90.md)


```Fortran

module wetdepscheme_jacob_mod

   use precision_mod, only: fp, zero, one, rae, tiny_
   use error_mod, only: cc_warning, cc_success !CC_Error
   use wetdepcommon_mod, only: wetdepschemejacobconfig
   use constants, only: g0, airmw  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_jacob

   ! -- local parameters
   real(fp), parameter :: density_ice = 917.0_fp                 ! density of ice in kg m-3
   real(fp), parameter :: density_liq = 1.e+03_fp                ! density of liquid water in kg m-3
   real(fp), parameter :: m_to_cm  = 100.0_fp                    ! conversion factor from m to cm
   real(fp), parameter :: kg_to_cm3_liq = m_to_cm / density_liq  ! conversion factor from kg to cm3 for liquid water
   real(fp), parameter :: kg_to_cm3_ice = m_to_cm / density_ice  ! conversion factor from kg to cm3 for ice
   real(fp), parameter :: qq_thr   = 0.0_fp                      ! cm3 (h2o) / cm3 (air) / s
   real(fp), parameter :: pdwn_thr = 0.0_fp                      ! cm3 (h2o) / cm2 (air) / s
   real(fp), parameter :: k_min = 1.e-04_fp ! s-1
   real(fp), parameter :: cwc   = 1.e-06_fp ! s-1 (recommended by Qiaoqiao Wang et al., 2014. Originally 1.5e-6, see Jacob et al., 2000)


contains

   subroutine compute_jacob( &
      num_layers, &
      num_species, &
      params, &
      airden_dry, &
      mairden, &
      pedge, &
      pfilsan, &
      pfllsan, &
      reevapls, &
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
      type(WetDepSchemeJACOBConfig), intent(in) :: params
      real(fp), intent(in) :: airden_dry(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: mairden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: pedge(num_layers+1)  ! Edge field - requires nz+1 dimensions
      real(fp), intent(in) :: pfilsan(num_layers+1)    ! 3D atmospheric field
      real(fp), intent(in) :: pfllsan(num_layers+1)    ! 3D atmospheric field
      real(fp), intent(in) :: reevapls(num_layers)    ! 3D atmospheric field
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
      real(fp), intent(in) :: species_wd_rainouteff(:,:)  ! Species wd_rainouteff property
      real(fp), intent(in) :: species_radius(:)  ! Species radius property
      real(fp), intent(in) :: species_mw_g(:)  ! Species mw_g property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: wetdep_mass_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: wetdep_flux_per_species_per_level(:,:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: k, species_idx, km1, ktop, kbot, so2_id, so4_id, h2o2_id, rc
      integer :: diag_idx  ! For diagnostic species indexing
      ! local physical variables
      real(fp)     :: delp       ! pressure thickness [Pa]
      real(fp)     :: dqls       ! liquid water flux gradient [kg/(m^2 s)]
      real(fp)     :: dqis       ! ice water flux gradient [kg/(m^2 s)]
      real(fp)     :: dqls_kgm3s ! liquid water flux gradient [kg/(m^3 s)]
      real(fp)     :: dqis_kgm3s ! ice water flux gradient [kg/(m^3 s)]
      real(fp)     :: f          ! total precipitation fraction (f_rainout + f_washout) [1]
      real(fp)     :: ftop       ! top of grid box rainout fraction [1]
      real(fp)     :: f_prime    ! rainout fraction in middle layers [1]
      real(fp)     :: f_rainout  ! rainout fraction [1]
      real(fp)     :: f_washout  ! washout fraction [1]
      real(fp)     :: k_rain     ! rainout rate [m^3/s]
      logical      :: kin        ! kinetic process flag [kinetic or equilibrium]
      real(fp)     :: dt         ! chemistry model time-step [sec]
      real(fp)     :: lossfrac   ! loss fraction
      real(fp)     :: qdwn       ! cm3 (h2o) / cm2 (air) / s
      real(fp)     :: press      ! pressure [Pa]
      real(fp)     :: delz       ! thickness of layer [m]
      real(fp)     :: efficiency(3)  ! efficiency factors for rainout
      real(fp), dimension(:), allocatable :: qq      ! precipatitng water rate [cm3 (h2o) / cm2 (air) / s]
      real(fp), dimension(:), allocatable :: pdwn    ! preciptation rate at top of grid cells [cm3 (h2o) / cm2 (air) / s]
      real(fp), dimension(:), allocatable :: dpog    ! pressure thickness of grid cells divided by gravity [Pa / (m/s^2)]
      real(fp), dimension(:), allocatable :: conc    ! concentration [kg/m2] converted from conc_in [kg/kg]
      real(fp), dimension(:), allocatable :: SO2     ! concentration of SO2 [kg/kg]; converted in rainout and washout, not here
      real(fp), dimension(:), allocatable :: SO4     ! concentration of SO4 [kg/m2]; converted from input SO4_in [kg/kg]
      real(fp), dimension(:), allocatable :: H2O2    ! concentration of H2O2 [kg/kg]; converted in rainout and washout, not here
      real(fp), dimension(:), allocatable :: dconc   ! concentration loss kg/m2
      real(fp), dimension(:), allocatable :: c_h2o   ! concentration of h2o
      real(fp), dimension(:), allocatable :: cldice  ! ice concentration
      real(fp), dimension(:), allocatable :: cldliq  ! liquid water concentration
      real(fp), dimension(:), allocatable :: reevap  ! evaporation rate [cm3 (h2o) / cm2 (air) / s]
      real(fp), dimension(:), allocatable :: delz_cm ! thickness of layer [cm]
      !
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg
      errmsg  = ''
      thisloc = ' -> at compute_jacob (in process/wetdep/scheme/WetDepScheme_JACOB_Mod.F90)'

      ! -- begin
      rc = cc_success

      !initialize variables
      ktop = num_layers
      kbot = 1
      km1 = 2  !This is to depress the warning 'km1 may be used uninitialized'
      dt = tstep

      allocate(qq(kbot:ktop), pdwn(kbot:ktop), conc(kbot:ktop), dconc(kbot:ktop), dpog(kbot:ktop), c_h2o(kbot:ktop), &
         cldice(kbot:ktop), cldliq(kbot:ktop), delz_cm(kbot:ktop), so2(kbot:ktop), so4(kbot:ktop), h2o2(kbot:ktop), reevap(kbot:ktop))

      ! find species indices for SO2, SO4, and H2O2
      so2_id = -1
      so4_id = -1
      h2o2_id = -1
      so2_id = max(find_species_ind(species_short_name, 'SO2'), find_species_ind(species_short_name, 'so2'))
      so4_id = max(find_species_ind(species_short_name, 'SO4'), find_species_ind(species_short_name, 'so4'), &
         find_species_ind(species_short_name, 'aso4j'), find_species_ind(species_short_name, 'ASO4J'))
      h2o2_id = max(find_species_ind(species_short_name, 'H2O2'), find_species_ind(species_short_name, 'h2o2'))
      if (so2_id < 1 ) then
         errmsg = 'SO2 is not a chemical species in the model. Jacob wet deposition scheme will assign zero to it.'
         CALL cc_warning( errmsg, rc, thisloc )
         so2 = zero
      else
         so2 = species_conc(:, so2_id) * species_mw_g(so2_id) * 1.0e-6_fp / airmw !convert from ppmv to kg/kg
      endif
      if (so4_id < 1 ) then
         errmsg = 'SO4 is not a chemical species in the model. Jacob wet deposition scheme will assign zero to it.'
         CALL cc_warning( errmsg, rc, thisloc )
         so4 = zero
      else
         so4 = species_conc(:, so4_id) * 1.e-09_fp  !convert from ug/kg to kg/kg
      endif
      if (h2o2_id < 1 ) then
         errmsg = 'H2O2 is not a chemical species in the model. Jacob wet deposition scheme will assign zero to it.'
         CALL cc_warning( errmsg, rc, thisloc )
         h2o2 = zero
      else
         h2o2 = species_conc(:, h2o2_id) * species_mw_g(h2o2_id) * 1.0e-6_fp / airmw !convert from ppmv to kg/kg
      endif

      ! calculate vertical met first
      do k = kbot, ktop
         km1 = k + 1

         ! -- initialize auxiliary arrays
         !if (k == ktop) then
         !   !TODO: GOCART has an additional index on the model top edge;
         !   dqls = pfllsan(k)
         !   dqis = pfilsan(k)
         !   pdwn(k) = kg_to_cm3_liq * pfllsan(k) + kg_to_cm3_ice * pfilsan(k)
         !else

         !Here we follow GOCART with an additional index; otherwise, uncomment the if else statement above
         ! -- liquid/ice precipitation formation in grid cell (kg/m2/s)
         dqls = pfllsan(k) - pfllsan(km1)
         dqis = pfilsan(k) - pfilsan(km1)
         ! -- precipitation flux from upper level (convert from kg/m2/s to cm3/cm2/s)
         pdwn(k) = kg_to_cm3_liq * pfllsan(km1) + kg_to_cm3_ice * pfilsan(km1)

         !end if ! if (k == ktop)

         delp = pedge(k) - pedge(km1)
         dpog(k) = delp / g0
         delz = dpog(k) / mairden(k) ! thickness of layer [m]
         delz_cm(k) = delz * m_to_cm  ! thickness of layer [cm]

         ! -- liquid/ice precipitation formation in grid cell (kg/m2/s)
         !dqls = pfllsan(k) - pfllsan(km1)
         !dqis = pfilsan(k) - pfilsan(km1)

         ! -- convert from kg/m2/s to kg (H2O) / m3(air) / s
         dqls_kgm3s = dqls / delz
         dqis_kgm3s = dqis / delz

         ! -- total precipitation formation (convert from kg (H2O) / m3(air) / s to cm3 (H2O) / cm3 (air) /s)
         ! -- To convert from kg (H2O) / m3(air) / s to cm3 (H2O) / cm3 (air) / s, divide by the density of
         ! -- the precipitation (ice or liquid)
         qq(k) =  dqls_kgm3s / density_liq +  dqis_kgm3s / density_ice
         reevap(k) = reevapls(k) * (airden_dry(k) / 1000.0_fp) ! convert from kg/kg/s to cm3/cm2/s

         ! -- precipitation flux from upper level (convert from kg/m2/s to cm3/cm2/s)
         !pdwn(k) = kg_to_cm3_liq * pfllsan(km1) + kg_to_cm3_ice * pfilsan(km1)

         ! -- initialize concentrations array, converting from kg/kg to kg/m2
         !this seems for both gas and aerosol
         !SO2(k)  = conc_in(k) !already assigned before the loop
         so4(k)  = so4(k) * dpog(k)

         ! -- compute mixing ratio of saturated water vapour over ice (from SETUP_WETSCAV)
         press     = 0.5_fp * ( pedge(km1) + pedge(k) ) !pressure in grid box
         c_h2o(k) = 10._fp ** (-2663.5_fp / t(k) + 12.537_fp ) / press

         ! -- estimate cloud ice and liquid water content (from SETUP_WETSCAV)
         if ( t(k) >= 268.0_fp ) then
            cldliq(k) = cwc
         else if ( t(k) > 248.0_fp ) then
            cldliq(k) = cwc * ( t(k) - 248.0_fp ) / 20.0_fp
         else
            cldliq(k) = zero
         end if
         cldice(k) = max(cwc - cldliq(k), zero) ! ensure cldice >= 0
      end do

      ! loop each species for wet deposition calculation
      do species_idx = 1, num_species
         !get input concentration for this species
         ! -- initialize concentrations array, converting from ug/kg or ppmv to kg/m2
         if (species_is_aerosol(species_idx)) then
            conc(:) = species_conc(:, species_idx)  * 1.e-09_fp * dpog(:) !convert from ug/kg to kg/kg and then to kg/m2
         else
            conc(:) = species_conc(:, species_idx) * species_mw_g(species_idx) * 1.0e-6_fp / airmw * dpog(:) !convert from ppmv to kg/kg and then to kg/m2
         end if
         ! -- initialize loss array
         dconc(:) = zero
         efficiency(:) = species_wd_rainouteff(species_idx, :)

         ! -- starts at the top
         k = ktop
         f = zero
         if (qq(k) > qq_thr) then
            ! -- compute rainout rate
            k_rain = k_min + qq(k) / cwc
            f = qq(k) / ( k_rain * cwc )

            call rainout(species_is_aerosol(species_idx), efficiency, species_wd_liqandgas(species_idx), &
               species_henry_k0(species_idx), species_henry_cr(species_idx), species_henry_pka(species_idx),              &
               species_wd_convfaci2g(species_idx), species_wd_retfactor(species_idx), f, k_rain, dt, t(k),  c_h2o(k),     &
               cldice(k), cldliq(k), species_short_name(species_idx), lossfrac, so2(k), h2o2(k))

            ! -- compute and apply effective loss fraction
            call rainout_loss( k, lossfrac, conc, dconc )

         end if

         ! -- middle layers
         ftop = f
         do k = ktop-1 , kbot+1, -1
            km1 = k + 1

            f_prime = zero
            ! -- if precipitation is forming in the grid cell
            if (qq(k) > qq_thr) then
               k_rain = k_min + qq(k) / cwc
               f_prime = qq(k) / ( k_rain * cwc )
            end if

            ! -- account for precipitation flux
            f_rainout = zero
            f_washout = zero

            if (pdwn(k) > pdwn_thr) then
               f_rainout = f_prime
               f_washout = max( zero, ftop - f_rainout )
            end if

            f = f_rainout + f_washout

            if ( f > zero ) then
               if ( f_rainout > zero ) then

                  call rainout(species_is_aerosol(species_idx), efficiency, species_wd_liqandgas(species_idx), &
                     species_henry_k0(species_idx), species_henry_cr(species_idx), species_henry_pka(species_idx),              &
                     species_wd_convfaci2g(species_idx), species_wd_retfactor(species_idx), f, k_rain, dt, t(k),  c_h2o(k),     &
                     cldice(k), cldliq(k), species_short_name(species_idx), lossfrac, so2(k), h2o2(k))

                  ! -- compute and apply effective loss fraction
                  call rainout_loss( k, lossfrac, conc, dconc )

               end if
               if ( f_washout > zero ) then
                  if ( f_rainout > zero ) then
                     ! -- washout from precipitation entering from the top
                     qdwn = pdwn(km1)
                     !TODO: is reevap available in GFS? Not used in GOCART version?
                     reevap(k) = max(reevap(k), 0e+0_fp)
                  else
                     ! -- washout from precipitation leaving through the bottom
                     qdwn = pdwn(k)
                  end if

                  call washout(species_radius(species_idx), f, t(k), qdwn, delz_cm(k), dt, species_short_name(species_idx), &
                     species_is_aerosol(species_idx), species_henry_k0(species_idx), species_henry_cr(species_idx), &
                     species_henry_pka(species_idx), params%scale_factor, params%radius_threshold, lossfrac, kin, so2(k), h2o2(k))

                  ! -- compute and apply effective loss fraction
                  call washout_loss( k, lossfrac, kin, f_washout, f_rainout, pdwn, reevap(k), &
                     delz_cm, conc, dconc, species_short_name(species_idx), so4 )

               end if
            else
               ! -- complete resuspension of rainout + washout from level above
               call complete_reevap( k, conc, dconc, species_short_name(species_idx), so4 )

            end if

            ftop = f

         end do

         ! -- surface level
         k = kbot
         if (pdwn(km1) > pdwn_thr) then
            f = ftop
            if ( f > zero ) then
               qdwn = pdwn(km1)

               call washout(species_radius(species_idx), f, t(k), qdwn, delz_cm(k), dt, species_short_name(species_idx), &
                  species_is_aerosol(species_idx), species_henry_k0(species_idx), species_henry_cr(species_idx), &
                  species_henry_pka(species_idx), params%scale_factor, params%radius_threshold, lossfrac, kin, so2(k), h2o2(k))

               ! -- compute and apply effective loss fraction
               call washout_loss( k, lossfrac, kin, f_washout, f_rainout, pdwn, reevap(k), &
                  delz_cm, conc, dconc, species_short_name(species_idx), so4 )

            end if
         end if

         ! calculate vertical met first
         do k = kbot, ktop

            ! -- convert back to ug/kg or ppmv
            if (species_is_aerosol(species_idx)) then
               species_tendencies(k, species_idx) = max(0.0_fp, conc(k)) / dpog(k) * 1.0e9_fp
            else
               species_tendencies(k, species_idx) = max(0.0_fp, conc(k)) / dpog(k) * airmw / species_mw_g(species_idx) * 1.0e6_fp
            end if

            ! Update diagnostic fields here based on your scheme's requirements
            ! Per-species-per-level diagnostic: 2D array (levels, species)
            if (present(wetdep_mass_per_species_per_level) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom wet deposition mass loss per species per level calculation
                     wetdep_mass_per_species_per_level(k, diag_idx) = dconc(k)
                     exit
                  end if
               end do
            end if
            ! Per-species-per-level diagnostic: 2D array (levels, species)
            if (present(wetdep_flux_per_species_per_level) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom wet deposition flux per species per level calculation
                     wetdep_flux_per_species_per_level(k, diag_idx) = dconc(k) / dt
                     exit
                  end if
               end do
            end if
         end do ! End layer loop

      end do ! End species loop

      deallocate(qq, pdwn, conc, dconc, dpog, delz_cm, c_h2o, cldice, cldliq, so2, so4, h2o2, reevap)

   end subroutine compute_jacob

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   function find_species_ind(SpeciesNames, species_name) result(species_index)
      implicit none
      character(len=*), intent(in) :: SpeciesNames(:)
      character(len=*), intent(in) :: species_name
      integer :: species_index

      integer :: i

      species_index = 0  ! Not found

      if (size(speciesnames) > 0 ) then
         do i = 1, size(speciesnames)
            if (trim(speciesnames(i)) == trim(species_name)) then
               species_index = i
               exit
            endif
         enddo
      endif

   end function find_species_ind

   subroutine rainout( is_aero, efficiency, wd_LidAndGas, k0, cr, pKa, cnvI2G, retfac, f, k, dt, tk, c_h2o, cldice, cldliq, spc, lossfrac, SO2, H2O2 )
      IMPLICIT NONE
      ! Parameters
      !-----------
      logical,    intent(in)  :: is_aero
      real(fp),   intent(in)  :: efficiency(3)
      logical,    intent(in)  :: wd_LidAndGas
      real(fp),   intent(in)  :: k0
      real(fp),   intent(in)  :: cr
      real(fp),   intent(in)  :: pKa
      real(fp),   intent(in)  :: cnvI2G
      real(fp),   intent(in)  :: retfac
      real(fp),   intent(in)  :: f
      real(fp),   intent(in)  :: k
      real(fp),   intent(in)  :: dt
      real(fp),   intent(in)  :: tk
      real(fp),   intent(in)  :: c_h2o
      real(fp),   intent(in)  :: cldice
      real(fp),   intent(in)  :: cldliq
      character(len = 20),  intent(in)  :: spc
      real(fp),   intent(out) :: lossfrac
      real(fp),   intent(in) :: SO2
      real(fp),   intent(in) :: H2O2

      ! Local Variables
      !----------------
      real(fp) :: i2g, l2g, c_tot
      real(fp) :: f_i, f_l, ki
      real(fp) :: SO2s, H2O2s   !after uint conversion to mol/mol
      real(fp) :: SO2LOSS
      real(fp), parameter :: kc = 5.e-3_fp    ! conversion rate from cloud condensate to precip (s-1)

      !--------------------------------------------
      ! main function
      !--------------------------------------------
      !=================================================================
      ! %%% SPECIAL CASE %%%
      ! SO2, HNO3 and H2SO4 scavenges like an aerosol although they are
      ! considered to be a gas-phase species elsewhere
      !=================================================================
      if (is_aero .or. spc == 'SO2' .or. spc == 'HNO3' .or. spc == 'H2SO4'  &
         .or. spc == 'so2' .or. spc == 'hno3' .or. spc == 'h2so4') then
         lossfrac = rainfrac( f, k, dt )

         ! -- apply rainout efficiency (simplify from APPLY_RAINOUT_EFF)
         if (tk < 237.0_fp) then
            ! ice
            lossfrac = efficiency(1) * lossfrac
         else if (tk < 258.0_fp) then
            ! snow
            lossfrac = efficiency(2) * lossfrac
         else
            ! liquid rain
            lossfrac = efficiency(3) * lossfrac
         end if

         if (spc == 'SO2' .or. spc == 'so2') then
            !unit conversion to mol/mol using airMW and SO2MW or H2OMW
            so2s = so2 * (airmw / 64.04_fp)
            h2o2s = h2o2 * (airmw / 34.02_fp)
            ! Update SO2 and H2O2
            if ( so2s > tiny_ ) then
               ! Limit lossfrac
               so2loss      = min( so2s, h2o2s )
               lossfrac     = so2loss * lossfrac / so2s
               lossfrac     = max( lossfrac, 0e+0_fp )

               ! Update saved H2O2 concentration
               !TODO:: comment out the depletion since we may not have afterchem H2O2 like in GOES-Chem
               !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
               !H2O2s = H2O2s - ( SO2s * lossfrac )
               !H2O2s = MAX( H2O2s, TINY_ )

            else

               ! If SO2 is not defined then set lossfrac to 0
               lossfrac     = 0.0_fp

            endif

            ! Update saved SO2 concentration
            !TODO:: comment out the depletion since we may not have afterchem SO2 like in GOES-Chem
            !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
            !SO2s     = SO2s * ( 1.0_fp - lossfrac )
            !SO2s     = MAX( SO2s, TINY_ )
         end if

      else !for gases exceptfor SO2, HNO3 and H2SO4

         ! -- compute ice to gas ratio assuming scavenging by co-condensation
         !    simplified from COMPUTE_Ki
         ! wd_LidAndGas is true only for SO2, NH3 and H2O2 for now
         i2g = zero
         if (wd_lidandgas) then
            if ( c_h2o > zero ) i2g = cnvi2g * cldice / c_h2o
         end if
         ! -- compute l2g (adopted from COMPUTE_L2G)
         !TODO: seems an error in GOCART version here
         l2g = liq_to_gas_ratio( k0, cr, pka, tk, cldliq)

         ! -- fraction of species in liquid and ice phases
         c_tot = one + l2g + i2g
         f_l   = l2g / c_tot
         f_i   = i2g / c_tot

         ! -- compute Ki for loss due to scavenging from convective updraft
         if ( tk >= 268.0_fp ) then
            ki = kc * ( f_l + f_i )
         else if ( tk > 248.0_fp ) then
            ki = kc * ( retfac * f_l + f_i )
         else
            ki = kc * f_i
         end if

         ! -- compute rained-out fraction
         lossfrac = rainfrac( f, ki, dt )

      end if

   end subroutine rainout


   subroutine washout( radius, f, tk, qdwn, dz, dt, spc, is_aero, k0, cr, pKa, wtune, radius_thr, washfrac, kin, SO2, H2O2)

      implicit none

      real(fp),  intent(in)  :: radius
      real(fp),  intent(in)  :: f
      real(fp),  intent(in)  :: tk
      real(fp),  intent(in)  :: qdwn
      real(fp),  intent(in)  :: dz
      real(fp),  intent(in)  :: dt
      character(len = 20),  intent(in) :: spc
      logical,   intent(in)  :: is_aero
      real(fp),  intent(in) :: k0
      real(fp),  intent(in) :: cr
      real(fp),  intent(in) :: pKa
      real(fp),  intent(in)  :: wtune
      real(fp),  intent(in)  :: radius_thr
      real(fp),  intent(out) :: washfrac
      logical,   intent(out)  :: kin
      real(fp),  intent(in) :: H2O2
      real(fp),  intent(in) :: SO2
      ! -- local variables
      real(fp)               :: SO2LOSS
      real(fp)               :: SO2s, H2O2s !unit conversion to mol/mol

      ! -- begin
      washfrac = zero

      !=================================================================
      ! HNO3 scavenges like an aerosol although it is considered
      ! to be a gas-phase species elsewhere (e.g. dry deposition)
      !=================================================================

      IF ( spc == 'HNO3' .or. spc == 'hno3' ) THEN  !TODO: better way to check species name?

         ! Washout is a kinetic process
         kin      = .true.

         ! Get washout fraction
         washfrac = washfrac_hno3(  f, tk, qdwn,  dt )

         !=================================================================
         ! SO2 scavenges like an aerosol although it is considered
         ! to be a gas-phase species elsewhere (e.g. dry deposition)
         !=================================================================
      ELSE IF ( spc == 'SO2' .or. spc == 'H2SO4' .or. &
         spc == 'so2' .or. spc == 'h2so4') THEN  !TODO: better way to check species name?

         ! NOTE: Even though SO2 is not an aerosol we treat it as SO4 in
         ! wet scavenging.  When evaporation occurs, it returns to SO4.
         kin      = .true.
         !NOTE: Here we use a dummy radius = 0.5 um for SO2/H2SO4 since
         !      it is considered as fine aerosol here
         washfrac = washfrac_aerosol( 0.5_fp, f, tk, qdwn, dt, wtune, 1.0_fp )

         IF (spc == 'SO2' .or. spc == 'so2') THEN
            ! Use the wet-scavenging following [Chin et al, 1996] such
            ! that a soluble fraction of SO2 is limited by the availability
            ! of H2O2 in the precipitating grid box.  Then scavenge the
            ! soluble SO2 at the same rate as sulfate.

            !unit conversion to mol/mol using airMW and SO2MW or H2OMW
            so2s = so2 * (airmw / 64.04_fp)
            h2o2s = h2o2 * (airmw / 34.02_fp)
            IF ( tk >= 268e+0_fp .AND. so2s > tiny_ ) THEN

               ! Adjust WASHFRAC
               so2loss  = min( so2s, h2o2s )
               washfrac = so2loss * washfrac / so2s
               washfrac = max( washfrac, 0e+0_fp )

               ! Deplete H2O2s the same as SO2s
               !TODO:: comment out the depletion since we may not have afterchem H2O2 like in GOES-Chem
               !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
               !H2O2s = H2O2s - ( SO2s * WASHFRAC )
               !H2O2s = MAX( H2O2s, TINY_ )

            ELSE
               washfrac = 0e+0_fp

            ENDIF

            ! Update saved SO2 concentration
            !TODO:: comment out the depletion since we may not have afterchem SO2 like in GOES-Chem
            !       Otherwise, the depletion will do twice with the later one in washout_loss subroutine
            !SO2s = SO2s * ( 1e+0_fp - WASHFRAC )
            !SO2s = MAX( SO2s, TINY_ )
         END IF

         !-----------------------------------------------------------------
         ! Washout for aerosol species
         !-----------------------------------------------------------------
      ELSE IF ( is_aero ) THEN

         ! Washout is a kinetic process
         kin      = .true.
         washfrac = washfrac_aerosol( radius, f, tk, qdwn, dt, wtune, radius_thr )

         !-----------------------------------------------------------------
         ! Washout for gas-phase species
         ! (except H2SO4, NO3, and SO2, which scavenge like aerosols;
         !-----------------------------------------------------------------
      ELSE
         ! Washout is an equilibrium process
         kin = .false.
         CALL washfrac_liq_gas( f, tk, qdwn, dz, dt, k0, cr, pka, washfrac, kin )

      END IF

   end subroutine washout

   real(fp) function rainfrac( f, k, dt )

      implicit none
      ! INPUT Parameters
      real(fp), intent(in) :: f
      real(fp), intent(in) :: k
      real(fp), intent(in) :: dt

      rainfrac = f * ( one - exp( -k * dt ) )

   end function rainfrac

   real(fp) function liq_to_gas_ratio( k0, cr, pKa, tk, qliq )

      real(fp), intent(in) :: k0
      real(fp), intent(in) :: cr
      real(fp), intent(in) :: pKa
      real(fp), intent(in) :: tk
      real(fp), intent(in) :: qliq

      ! -- local variables
      real(fp) :: h !, t
      ! -- local parameters
      real(fp), parameter :: cloud_pH = 4.5_fp
      real(fp), parameter  :: Tref = 298.15_fp     ! K
      real(fp), parameter  :: R    = 8.3144598_fp  ! J K-1 mol-1
      real(fp), parameter  :: Pref = 101.325_fp    ! mPa

      ! -- compute Henry's law constant for a given temperature
      !t = real(tk, kind=f8) TODO: not sure why this conversion is needed; do not use it for now
      h = k0 * exp( cr * (1._fp/tk - 1._fp/tref) ) * r * tk / pref

      ! -- adjust Henry's law constant for chemical equilibriums in liquid phase
      if ( pka > -100._fp ) h = h * ( one + 10._fp ** ( cloud_ph - pka ) )

      liq_to_gas_ratio = h * qliq

   end function liq_to_gas_ratio

   real(fp) function washfrac_aerosol( radius, f, tk, pdwn, dt, tuning, radius_fine )

      implicit none

      real(fp), intent(in) :: radius
      real(fp), intent(in) :: f
      real(fp), intent(in) :: tk
      real(fp), intent(in) :: pdwn
      real(fp), intent(in) :: dt
      real(fp), intent(in) :: tuning
      real(fp), intent(in) :: radius_fine

      ! -- local variables
      real(fp)          :: dth, pph
      ! -- local parameters
      real(fp), parameter :: k_wash = 1.06e-03_fp
      real(fp), parameter :: h2s = 3600.0_fp ! s-1

      ! -- begin
      washfrac_aerosol = zero

      if ( f > zero ) then
         ! -- convert instant rates (s-1) to hourly rates
         pph = 10.0_fp * pdwn * h2s
         dth = dt / h2s

         if ( radius < radius_fine ) then  !for fine aerosol (simplified from WASHFRAC_FINE_AEROSOL)
            if ( tk >= 268e+0_fp  ) then
               washfrac_aerosol = f * ( one  - exp(-k_wash * tuning * (pph / f ) ** 0.61e+0_fp * dth))
            else
               washfrac_aerosol = f * ( one  - exp(-2.6e+1_fp * k_wash * tuning  * (pph / f ) ** 0.96e+0_fp * dth))
            endif
         else  !for coarse aerosol (simplified from WASHFRAC_COARSE_AEROSOL)
            if ( tk >= 268e+0_fp  ) then
               washfrac_aerosol = f * ( one  - exp(-0.92e+0_fp * tuning * (pph / f ) ** 0.79e+0_fp * dth))
            else
               !TODO: GOCART applied a factor of 0.5 to the tuning factor for coarse aerosol????
               !washfrac_aerosol = F * ( one  - EXP(-1.57e+0_fp / 0.5e+0_fp * tuning * (pph / f ) ** 0.96e+0_fp * dth))
               washfrac_aerosol = f * ( one  - exp(-1.57e+0_fp * tuning * (pph / f ) ** 0.96e+0_fp * dth))
            endif
         endif
      endif

   end function washfrac_aerosol

   real(fp) function washfrac_hno3( f, tk, pdwn, dt )

      implicit none

      real(fp), intent(in) :: f
      real(fp), intent(in) :: tk
      real(fp), intent(in) :: pdwn
      real(fp), intent(in) :: dt

      ! -- local parameters
      real(fp), parameter :: k_wash = 1.0_fp    ! First order washout rate (cm-1)

      ! -- begin
      washfrac_hno3 = zero
      ! -- compute washout fraction only if T >= 268K
      if ( tk >= 268.0_fp .and. f > zero ) then
         washfrac_hno3 = f * ( one - exp( -k_wash * pdwn * dt / f ) )
      end if

   end function washfrac_hno3

   subroutine washfrac_liq_gas( f, tk, pdwn, dz, dt, k0, cr, pKa, washfrac, kin)

      implicit none

      real(fp),  intent(in) :: f
      real(fp),  intent(in) :: tk
      real(fp),  intent(in) :: pdwn
      real(fp),  intent(in) :: dz
      real(fp),  intent(in) :: dt
      real(fp),  intent(in) :: k0
      real(fp),  intent(in) :: cr
      real(fp),  intent(in) :: pKa
      real(fp),  intent(out) :: washfrac
      logical,   intent(out) :: kin

      ! -- local variables
      real(fp) :: qliq, l2g, washfrac_kin

      ! -- begin

      ! Start with the assumption that washout will be an
      ! equilibrium process
      kin = .false.

      if ( tk < 268.0_fp ) then
         ! -- no washout
         washfrac = zero

      else

         ! -- compute L2G
         qliq = pdwn * dt / ( f * dz )

         ! Compute liquid to gas ratio
         l2g = liq_to_gas_ratio( k0, cr, pka, tk, qliq )

         ! -- washout fraction from Henry's Law
         washfrac = l2g / ( one + l2g )

         ! -- washout fraction from kinetic processes (HNO3)
         ! set f = one and call washfrac_hno3 function above
         washfrac_kin = washfrac_hno3( one, tk, pdwn, dt )

         ! -- equilibrium washout must not exceed kinetic washout
         !TODO: GOCART is missing  'washfrac_kin * f' here ??????
         if ( washfrac > washfrac_kin ) then
            washfrac = washfrac_kin * f
            kin = .true. ! washout is a kinetic process
         end if

      end if

   end subroutine washfrac_liq_gas

   subroutine rainout_loss( k, lossfrac, conc, dconc )

      implicit none

      ! -- input/output parameters
      integer,   intent(in)    :: k
      real(fp),  intent(in)    :: lossfrac
      real(fp),  dimension(:), intent(inout) :: conc
      real(fp),  dimension(:), intent(inout) :: dconc

      ! -- local variables
      real(fp)   :: wetloss

      ! -- apply loss (simplified from the DO_WASHOUT_ONLY subroutine)
      wetloss = lossfrac * conc(k)
      conc(k) = conc(k) - wetloss
      if ( k == size(conc) ) then !If it is the top layer; we assume the layer index is not reversed
         ! Dconc is an accumulator array for rained-out species.
         ! The species in Dconc are in the liquid phase and will
         ! precipitate to the levels below until a washout occurs.
         dconc(k) = wetloss
      else
         ! Add to Dconc the species lost to rainout in grid box
         ! (I,J,L) plus the species lost to rainout from grid box
         ! (I,J,L+1), which has by now precipitated down into
         ! grid box (I,J,L).  Dconc will continue to accumulate
         ! rained out species in this manner until a washout
         ! event occurs.
         dconc(k) = dconc(k+1) + wetloss
      end if

   end subroutine rainout_loss

   subroutine washout_loss( k, lossfrac, kin, f_washout, f_rainout, pdwn, reevap, delz_cm, conc, dconc, spc, SO4 )

      implicit none

      ! -- input/output parameters
      integer,   intent(in)    :: k
      real(fp),  intent(inout) :: lossfrac
      logical,   intent(in)    :: kin
      real(fp),  intent(in)    :: f_washout
      real(fp),  intent(in)    :: f_rainout
      real(fp),  dimension(:), intent(in)    :: pdwn
      !real(fp), dimension(:), intent(in)    :: qq      !< precipatitng water rate [cm3 (h2o) / cm2 (air) / s]
      real(fp)                 :: reevap
      real(fp),  dimension(:), intent(in)    :: delz_cm
      real(fp),  dimension(:), intent(inout) :: conc
      real(fp),  dimension(:), intent(inout) :: dconc
      real(fp),  dimension(:), intent(inout) :: SO4
      character(len = 20),  intent(in) :: spc

      ! -- local variables
      integer    :: km1
      real(fp)   :: f
      real(fp)   :: wetloss
      real(fp)   :: alpha
      real(fp)   :: gain
      real(fp)   :: washed

      !begins here (simplified from the DO_WASHOUT_ONLY and DO_WASHOUT_AT_SFC subroutines)
      f = f_washout + f_rainout
      km1 = k + 1
      if ( k == 1 ) then !for the surface layer; note this assumes the layer index starts from surface not the top

         ! -- f is included in lossfrac for aerosols and HNO3
         if ( kin ) then
            wetloss = lossfrac * conc(k)
         else
            wetloss = f * lossfrac * conc(k)
         end if
         conc(k) = conc(k) - wetloss
         dconc(k) = dconc(km1) + wetloss

      else  !for middle layers
         if ( kin ) then
            ! -- adjust loss fraction for aerosols
            lossfrac = lossfrac * f_washout / f

            ! Define ALPHA, the fraction of the raindrops that
            ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
            if ( pdwn(km1) - zero > zero  ) then !avoid divide by zero
               !TODO: is qq(k) right in here?
               !alpha = abs( qq(k) ) * delz_cm(k) / pdwn(km1)
               alpha = abs( reevap ) * delz_cm(k) / pdwn(km1)
            else
               alpha = one
            end if
            ! Restrict ALPHA to be less than 1
            alpha = min( one, alpha )
            ! Assume 50% of the re-evaporated water rains out to aerosols
            ! GAINED is the rained out aerosol coming down from
            ! grid box (I,J,L+1) that will evaporate and re-enter
            ! the atmosphere in the gas phase in grid box (I,J,L).
            gain  = 0.5_fp * alpha * dconc(km1)
            wetloss  = conc(k) * lossfrac - gain
            ! SO2 in sulfate chemistry is wet-scavenged on the
            ! raindrop and converted to SO4 by aqeuous chem.
            ! If evaporation occurs then SO2 comes back as SO4
            if (spc == 'SO2' .or. spc == 'so2') then
               so4(k) = so4(k) + gain * 96e+0_fp / 64e+0_fp
               conc(k) = conc(k) * (1e+0_fp - lossfrac)
            else
               conc(k) = conc(k) - wetloss
            end if

         else !not kinetic process (mainly for gases)

            !Not sure why dconc is not counted here in wetloss calculation
            washed  = f_washout * conc(k) + dconc(km1)
            wetloss = lossfrac * ( washed - dconc(km1) )
            conc(k) = conc(k) - wetloss

         end if ! kinetic or equilibrium process

         ! Add washout losses in grid box (I,J,L) to dconc
         if ( f_rainout > zero ) then
            dconc(k) = dconc(k) + wetloss
         else
            dconc(k) = dconc(km1) + wetloss
         end if

      end if !surface or middle layer

   end subroutine washout_loss

   subroutine complete_reevap( k, conc, dconc, spc, SO4 )

      implicit none

      ! -- input/output parameters
      integer,   intent(in)    :: k
      real(fp),  dimension(:), intent(inout) :: conc
      real(fp),  dimension(:), intent(inout) :: dconc
      character(len = 20),  intent(in) :: spc
      real(fp),  dimension(:), intent(inout)    :: SO4

      ! -- local variables
      integer    :: km1
      real(fp)   :: wetloss

      ! -- begin
      km1 = k + 1
      wetloss = -dconc(km1)
      ! All of the rained-out species coming from grid box
      ! (I,J,L+1) goes back into the gas phase at (I,J,L)
      ! In evap, SO2 comes back as SO4
      if (spc == 'SO2' .or. spc == 'so2') then
         so4(k) = so4(k) - (wetloss * 96e+0_fp / 64e+0_fp )
      else
         conc(k) = conc(k) - wetloss
      end if

      ! There is nothing rained out/washed out in grid box
      ! (I,J,L), so set dconc at grid box (I,J,L) to zero.
      dconc(k) = 0e+0_fp

   end subroutine complete_reevap


end module wetdepscheme_jacob_mod
```


