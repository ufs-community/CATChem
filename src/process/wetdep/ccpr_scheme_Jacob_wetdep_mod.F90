!>
!! \file
!! \brief CCPr Scheme for wet deposition
!!
!!
!! Reference:
!!
!! \author Wei Li
!! \date 04/2025
!!!>
module CCPr_Scheme_Jacob_WetDep_Mod

   implicit none

   private

   public :: CCPr_Scheme_Jacob_WetDep

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param
   !!
   !! \ingroup catchem_wetdep_process
   !!!>
   subroutine CCPr_Scheme_Jacob_WetDep( km, cdt, spc, is_aero, wd_LiqAndGas, k0, cr, pKa, retfac, cvtI2G, grav, radius, rainout_eff, &
      wtune, radius_thr, ple, tmpu, rhoa, pfllsan, pfilsan, qreevap, airden, conc_in, H2O2, SO4_in, fluxout,  rc )

      ! Uses
      use precision_mod, only : fp, zero     !pull in a precision from the PRECISION MODULE
      Use Error_Mod,     Only : CC_SUCCESS   ! Error Check Success
      USE CCPr_WetDep_Common_Mod
      implicit none

      ! !INPUT PARAMETERS:
      integer,                 intent(in)    :: km           !< total model levels
      real(fp),                intent(in)    :: cdt          !< chemistry model time-step [sec]
      character(len=20),       intent(in)    :: spc          !< species name
      logical,                 intent(in)    :: is_aero      !< true for aerosol
      logical,                 intent(in)    :: wd_LiqAndGas !< ice-to-gas ratio is computed by co-condensation?
      real(fp),                intent(in)    :: k0           !< Henry's solubility constant [M/atm]
      real(fp),                intent(in)    :: cr           !< Henry's volatility constant [K]
      real(fp),                intent(in)    :: pKa          !< Henry's pH correction factor [1]
      real(fp),                intent(in)    :: retfac       !< Retention factor [-]
      real(fp),                intent(in)    :: cvtI2G       !< Conversion factor from ice to gas ratio if wd_LidAndGas is true
      real(fp),                intent(in)    :: grav         !< gravity [m/sec^2]
      real(fp),                intent(in)    :: radius       !< Particle radius [um]
      real(fp), dimension(3),  intent(in)    :: rainout_eff  !< temperature-dependent rainout efficiencies TODO: can we read in as a list from species yaml file?
      real(fp),                intent(in)    :: wtune        !< Washout Tuning factor [-]; Note we add this, not from GC
      real(fp),                intent(in)    :: radius_thr   !< Threshold particle radius for washout[um]
      real(fp), dimension(:),  intent(in)    :: ple          !< pressure level thickness [Pa]
      real(fp), dimension(:),  intent(in)    :: tmpu         !< temperature [K]
      real(fp), dimension(:),  intent(in)    :: rhoa         !< moist air density [kg/m^3]
      real(fp), dimension(:),  intent(in)    :: pfllsan      !< 3D flux of liquid nonconvective precipitation [kg/(m^2 sec)]
      real(fp), dimension(:),  intent(in)    :: pfilsan      !< 3D flux of ice nonconvective precipitation [kg/(m^2 sec)]
      real(fp), dimension(:),  intent(in)    :: qreevap      !< Evaporation of precip LS+anvil [kg/kg/s]
      real(fp), dimension(:),  intent(in)    :: airden       !< dry air density [kg/m^3]
      real(fp), dimension(:),  intent(inout) :: conc_in      !< concentrations [kg/kg]
      real(fp), dimension(:),  intent(in)    :: H2O2         !< H2O2 concentration [kg/kg] used for SO2 washout
      real(fp), dimension(:),  intent(inout) :: SO4_in       !< SO4 concentration [kg/kg] used for SO2 washout
      real(fp), dimension(:),  intent(inout) :: fluxout      !< tracer loss flux [kg m-2 s-1]
      ! !OUTPUT PARAMETERS:
      integer,                intent(out)   :: rc           ! Error return code

      ! looping indexes
      integer  :: k, km1, ktop, kbot
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
      logical  :: kin        ! kinetic process flag [kinetic or equilibrium]
      real(fp)     :: dt         ! chemistry model time-step [sec]
      !real(fp)     :: totloss    ! total loss fraction
      real(fp)     :: lossfrac   ! loss fraction
      !real(fp)     :: wetloss    ! wet loss fraction before evaporation
      real(fp)     :: qdwn       ! cm3 (h2o) / cm2 (air) / s
      real(fp)     :: press       ! pressure [Pa]
      !real(fp)     :: alpha      ! ratio of evap. to sublimation
      !real(fp)     :: gain       ! gain fraction
      !real(fp)     :: washed     ! concentration of washed out tracer
      real(fp)     :: delz       ! thickness of layer [m]
      real(fp), dimension(:), allocatable :: qq      ! precipatitng water rate [cm3 (h2o) / cm2 (air) / s]
      real(fp), dimension(:), allocatable :: pdwn    ! preciptation rate at top of grid cells [cm3 (h2o) / cm2 (air) / s]
      real(fp), dimension(:), allocatable :: dpog    ! pressure thickness of grid cells divided by gravity [Pa / (m/s^2)]
      real(fp), dimension(:), allocatable :: conc    ! concentration [kg/m2] converted from conc_in [kg/kg]
      real(fp), dimension(:), allocatable :: SO2     ! concentration of SO2 [kg/kg]; converted in rainout and washout, not here
      real(fp), dimension(:), allocatable :: SO4     ! concentration of SO4 [kg/m2]; converted from input SO4_in [kg/kg]
      real(fp), dimension(:), allocatable :: dconc   ! concentration loss kg/m2
      real(fp), dimension(:), allocatable :: c_h2o   ! concentration of h2o
      real(fp), dimension(:), allocatable :: cldice  ! ice concentration
      real(fp), dimension(:), allocatable :: cldliq  ! liquid water concentration
      real(fp), dimension(:), allocatable :: reevap  ! evaporation rate [cm3 (h2o) / cm2 (air) / s]
      real(fp), dimension(:), allocatable :: delz_cm ! thickness of layer [cm]

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

      ! -- begin
      rc = cc_success

      !initialize variables
      ktop = km
      kbot = 1
      km1 = 2  !This is to depress the warning 'km1 may be used uninitialized'
      dt = cdt

      allocate(qq(kbot:ktop), pdwn(kbot:ktop), conc(kbot:ktop), dconc(kbot:ktop), dpog(kbot:ktop), &
         c_h2o(kbot:ktop), cldice(kbot:ktop), cldliq(kbot:ktop), delz_cm(kbot:ktop), SO2(kbot:ktop), SO4(kbot:ktop), reevap(kbot:ktop))

      ! -- compute column quantities
      do k = kbot, ktop
         km1 = k + 1 !TODO: check if it works when k = ktop

         ! -- initialize auxiliary arrays
         if (k == ktop) then
            !TODO: why GOCART does not have errors here?
            delp = ple(k)
            dqls = pfllsan(k)
            dqis = pfilsan(k)
            pdwn(k) = kg_to_cm3_liq * pfllsan(k) + kg_to_cm3_ice * pfilsan(k)
            press     = 0.5 * ( ZERO + ple(k) )
         else
            delp = ple(k) - ple(km1)
            dpog(k) = delp / grav
            delz = dpog(k) / rhoa(k) ! thickness of layer [m]
            delz_cm(k) = delz * m_to_cm  ! thickness of layer [cm]

            ! -- liquid/ice precipitation formation in grid cell (kg/m2/s)
            dqls = pfllsan(k) - pfllsan(km1)
            dqis = pfilsan(k) - pfilsan(km1)

            ! -- convert from kg/m2/s to kg (H2O) / m3(air) / s
            dqls_kgm3s = dqls / delz
            dqis_kgm3s = dqis / delz

            ! -- total precipitation formation (convert from kg (H2O) / m3(air) / s to cm3 (H2O) / cm3 (air) /s)
            ! -- To convert from kg (H2O) / m3(air) / s to cm3 (H2O) / cm3 (air) / s, divide by the density of
            ! -- the precipitation (ice or liquid)
            qq(k) =  dqls_kgm3s / density_liq +  dqis_kgm3s / density_ice
            reevap(k) = qreevap(k) * (airden(k) / 1000.0_fp) ! convert from kg/kg/s to cm3/cm2/s

            ! -- precipitation flux from upper level (convert from kg/m2/s to cm3/cm2/s)
            pdwn(k) = kg_to_cm3_liq * pfllsan(km1) + kg_to_cm3_ice * pfilsan(km1)

            ! -- initialize concentrations array, converting from kg/kg to kg/m2
            !this seems for both gas and aerosol
            SO2(k)  = conc_in(k) !SO2 is still in kg/kg; only used when spc == 'SO2' so using conc_in is fine
            conc(k) = conc_in(k) * dpog(k)
            SO4(k)  = SO4_in(k) * dpog(k)

            ! -- initialize loss array
            dconc(k) = zero

            ! -- compute mixing ratio of saturated water vapour over ice (from SETUP_WETSCAV)
            press     = 0.5 * ( ple(km1) + ple(k) )
            c_h2o(k) = 10._fp ** (-2663.5_fp / tmpu(k) + 12.537_fp ) / press

            ! -- estimate cloud ice and liquid water content (from SETUP_WETSCAV)
            if ( tmpu(k) >= 268.0_fp ) then
               cldliq(k) = cwc
            else if ( tmpu(k) > 248.0_fp ) then
               cldliq(k) = cwc * ( tmpu(k) - 248.0_fp ) / 20.0_fp
            else
               cldliq(k) = zero
            end if
            cldice(k) = MAX(cwc - cldliq(k), zero) ! ensure cldice >= 0
         end if ! if (k == ktop)
      end do

      ! -- starts at the top
      k = ktop
      f = zero
      if (qq(k) > qq_thr) then
         ! -- compute rainout rate
         k_rain = k_min + qq(k) / cwc
         f = qq(k) / ( k_rain * cwc )

         call rainout(is_aero, rainout_eff, wd_LiqAndGas, k0, cr, pKa, cvtI2G, retfac, f, k_rain, dt, tmpu(k), &
            c_h2o(k), cldice(k), cldliq(k), spc, lossfrac, SO2(k), H2O2(k))

         ! -- compute and apply effective loss fraction
         call rainout_loss( k, lossfrac, conc, dconc )

      end if

      ! -- middle layers
      ftop = f
      do k = ktop-1 , kbot+1, -1
         km1 = k - 1

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

               call rainout(is_aero, rainout_eff, wd_LiqAndGas, k0, cr, pKa, cvtI2G, retfac, f_rainout, k_rain, dt, tmpu(k), &
                  c_h2o(k), cldice(k), cldliq(k), spc, lossfrac, SO2(k), H2O2(k))

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

               call washout(radius, f, tmpu(k), qdwn, delz_cm(k), dt, spc, is_aero, &
                  k0, cr, pKa, wtune, radius_thr, lossfrac, kin, SO2(k) ,H2O2(k))

               ! -- compute and apply effective loss fraction
               call washout_loss( k, lossfrac, kin, f_washout, f_rainout, pdwn, reevap(k), &
                  delz_cm, conc, dconc, spc, SO4 )

            end if
         else
            ! -- complete resuspension of rainout + washout from level above
            call complete_reevap( k, conc, dconc, spc, SO4 )

         end if

         ftop = f

      end do

      ! -- surface level
      k = kbot
      if (pdwn(km1) > pdwn_thr) then
         f = ftop
         if ( f > zero ) then
            qdwn = pdwn(km1)

            call washout(radius, f, tmpu(k), qdwn, delz_cm(k), dt, spc, is_aero, &
               k0, cr, pKa, wtune, radius_thr, lossfrac, kin, SO2(k), H2O2(k))

            ! -- compute and apply effective loss fraction
            call washout_loss( k, lossfrac, kin, f_washout, f_rainout, pdwn, reevap(k), &
               delz_cm, conc, dconc, spc, SO4 )

         end if
      end if

      do k = ktop, kbot
         ! -- convert back to kg/kg
         conc_in(k) = conc(k) / dpog(k)
         SO4_in(k) = SO4(k) / dpog(k)
      end do

      !calculate fluxout
      fluxout = dconc / dt

      deallocate(qq, pdwn, conc, dconc, dpog, delz_cm, c_h2o, cldice, cldliq, SO2, SO4, reevap)

   end subroutine CCPr_Scheme_Jacob_WetDep

end module CCPr_Scheme_Jacob_WetDep_Mod
