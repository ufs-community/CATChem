!>
!! \file SettlingPhysics_Mod.F90
!! \brief Internalized settling physics routines from GOCART2G.
!!
!! All routines operate on 1D column arrays in CATChem's native
!! bottom-to-top vertical ordering. No vertical reversal is needed.
!!
!! Ported from src/external/GOCART/Process_Library/GOCART2G_Process.F90
!!
!! \ingroup settling_process
!<
module SettlingPhysics_Mod

   use precision_mod, only: fp
   implicit none
   private

   ! Public settling routines
   public :: settling_compute
   public :: settling_calc_vsettle
   public :: settling_solver_default
   public :: settling_solver
   public :: settling_particle_swelling
   public :: settling_swelling_fitzgerald
   public :: settling_swelling_gerber
   public :: settling_swelling_gerber_nh4so4
   public :: settling_swelling_pk2007

   ! Physical constants
   real(fp), parameter :: RHOW = 1000.0_fp              ! Density of water [kg/m3]
   real(fp), parameter :: V_UPWARD_MARING = 0.33e-2_fp  ! Maring upward velocity correction [m/s]

   ! Boltzmann constant and air molecule mass for Chem_CalcVsettle
   real(fp), parameter :: KB = 1.3807e-23_fp       ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real(fp), parameter :: M_AIR = 4.8096e-26_fp    ! Mass of avg air molecule [kg]
   real(fp), parameter :: PI = 3.141529265_fp

   ! Derived constants
   real(fp), parameter :: F_VT = 8.0_fp * KB / PI / M_AIR
   real(fp), parameter :: TWO_OVER_NINE = 2.0_fp / 9.0_fp

   ! Drag correction coefficients (Pruppacher and Klett)
   real(fp), parameter :: A0 = -3.18657_fp
   real(fp), parameter :: A1 =  0.992696_fp
   real(fp), parameter :: A2 = -1.53193e-3_fp
   real(fp), parameter :: A3 = -9.870593e-4_fp
   real(fp), parameter :: A4 = -5.78878e-4_fp
   real(fp), parameter :: A5 =  8.55176e-5_fp
   real(fp), parameter :: A6 = -3.27815e-6_fp

contains

   !---------------------------------------------------------------------------
   !> Main settling computation for a single column.
   !! Replaces GOCART2G Chem_Settling, operates on 1D column arrays
   !! in CATChem's native bottom-to-top vertical ordering.
   !!
   !! @param[in]    nlev             Number of vertical levels
   !! @param[in]    klid             Index for pressure lid
   !! @param[in]    cdt              Time step [s]
   !! @param[in]    grav             Gravitational acceleration [m/s2]
   !! @param[in]    radius_dry       Dry particle radius [m]
   !! @param[in]    rhop_dry         Dry particle density [kg/m3]
   !! @param[in]    swelling_flag    Particle swelling method (0=none, 1=Fitzgerald, 2=Gerber, 3=Gerber NH4SO4, 4=PK2007)
   !! @param[inout] conc             Species concentration [kg/kg]
   !! @param[in]    t                Temperature [K]
   !! @param[in]    airden           Air density [kg/m3]
   !! @param[in]    rh               Relative humidity [0-1]
   !! @param[in]    z_edge           Geopotential height at layer edges [m], size nlev+1
   !! @param[in]    delp             Pressure thickness [Pa]
   !! @param[out]   vsettle_out      Settling velocity [m/s] (optional)
   !! @param[out]   fluxout          Surface mass flux lost by settling [kg/m2/s] (optional)
   !! @param[in]    correction_maring Apply Maring et al. 2003 correction (optional)
   !! @param[in]    solver_type      1=default solver, 2=UFS solver (optional, default=1)
   !! @param[out]   rc               Return code (0=success)
   !---------------------------------------------------------------------------
   subroutine settling_compute(nlev, klid, cdt, grav, &
      radius_dry, rhop_dry, swelling_flag, &
      conc, t, airden, rh, z_edge, delp, &
      vsettle_out, fluxout, correction_maring, solver_type, rc)

      integer, intent(in) :: nlev
      integer, intent(in) :: klid
      real(fp), intent(in) :: cdt
      real(fp), intent(in) :: grav
      real(fp), intent(in) :: radius_dry
      real(fp), intent(in) :: rhop_dry
      integer, intent(in) :: swelling_flag
      real(fp), intent(inout) :: conc(nlev)
      real(fp), intent(in) :: t(nlev)
      real(fp), intent(in) :: airden(nlev)
      real(fp), intent(in) :: rh(nlev)
      real(fp), intent(in) :: z_edge(nlev+1)
      real(fp), intent(in) :: delp(nlev)
      real(fp), intent(out), optional :: vsettle_out(nlev)
      real(fp), intent(out), optional :: fluxout
      logical, intent(in), optional :: correction_maring
      integer, intent(in), optional :: solver_type
      integer, intent(out) :: rc

      ! Local variables
      integer :: k, solver
      real(fp) :: ONE_OVER_G
      real(fp) :: dz(nlev)
      real(fp) :: radius(nlev), rhop(nlev)
      real(fp) :: vsettle(nlev)
      real(fp) :: qa(nlev)
      real(fp) :: cmass_before, cmass_after

      rc = 0

      ! Default solver type
      solver = 1
      if (present(solver_type)) solver = solver_type

      ONE_OVER_G = 1.0_fp / grav

      ! If radius <= 0 then nothing to do
      if (radius_dry <= 0.0_fp) then
         rc = 100
         return
      end if

      ! Compute layer thickness from edge heights
      ! z_edge is at layer edges: z_edge(k) is bottom edge of layer k,
      ! z_edge(k+1) is top edge of layer k (bottom-to-top ordering)
      do k = 1, nlev
         dz(k) = z_edge(k+1) - z_edge(k)
      end do

      ! Copy concentration to working array
      qa = conc

      ! Compute column dry mass before sedimentation
      cmass_before = 0.0_fp
      do k = klid, nlev
         cmass_before = cmass_before + qa(k) * delp(k) * ONE_OVER_G
      end do

      ! Particle swelling
      call settling_particle_swelling(nlev, rh, radius_dry, rhop_dry, &
         radius, rhop, swelling_flag)

      ! Settling velocity of the wet particle
      vsettle = 0.0_fp
      do k = klid, nlev
         call settling_calc_vsettle(radius(k), rhop(k), airden(k), &
            t(k), grav, vsettle(k))
      end do

      ! Apply Maring correction if requested
      if (present(correction_maring)) then
         if (correction_maring) then
            do k = 1, nlev
               vsettle(k) = max(1.0e-9_fp, vsettle(k) - V_UPWARD_MARING)
            end do
         end if
      end if

      ! Output settling velocity if requested
      if (present(vsettle_out)) then
         vsettle_out = vsettle
      end if

      ! Time integration
      select case (solver)
       case (1)
         call settling_solver_default(nlev, cdt, delp, dz, vsettle, qa)
       case (2)
         call settling_solver(nlev, cdt, delp, dz, vsettle, qa)
      end select

      ! Compute column dry mass after sedimentation for flux calculation
      cmass_after = 0.0_fp
      do k = klid, nlev
         cmass_after = cmass_after + qa(k) * delp(k) * ONE_OVER_G
      end do

      ! Surface flux
      if (present(fluxout)) then
         fluxout = (cmass_before - cmass_after) / cdt
      end if

      ! Copy back
      conc = qa

   end subroutine settling_compute

   !---------------------------------------------------------------------------
   !> Stokes settling velocity with Cunningham slip correction.
   !! Replaces GOCART2G Chem_CalcVsettle.
   !!
   !! Calculates the aerosol settling velocity following Seinfeld and Pandis,
   !! Pruppacher and Klett, and CARMA (Toon et al., 1988).
   !! Includes drag correction for Reynolds number > 0.01.
   !!
   !! @param[in]  radius   Particle radius [m]
   !! @param[in]  rhop     Particle density [kg/m3]
   !! @param[in]  rhoa     Air density [kg/m3]
   !! @param[in]  tmpu     Temperature [K]
   !! @param[in]  grav     Gravitational acceleration [m/s2]
   !! @param[out] vsettle  Settling velocity [m/s]
   !---------------------------------------------------------------------------
   subroutine settling_calc_vsettle(radius, rhop, rhoa, tmpu, grav, vsettle)

      real(fp), intent(in) :: radius
      real(fp), intent(in) :: rhop
      real(fp), intent(in) :: rhoa
      real(fp), intent(in) :: tmpu
      real(fp), intent(in) :: grav
      real(fp), intent(out) :: vsettle

      ! Local variables
      real(fp) :: rmu    ! Dynamic viscosity [kg m-1 s-1]
      real(fp) :: vt     ! Thermal velocity of air molecule [m/s]
      real(fp) :: rmfp   ! Air molecule mean free path [m]
      real(fp) :: bpm    ! Cunningham slip correction factor
      real(fp) :: rkn    ! Knudsen number
      real(fp) :: re     ! Reynolds number
      real(fp) :: x, y   ! Parameters for drag correction

      ! Dynamic viscosity from corrected Sutherland's Equation
      rmu = 1.8325e-5_fp * (416.16_fp / (tmpu + 120.0_fp)) * (tmpu / 296.16_fp)**1.5_fp

      ! Thermal velocity of air molecule
      vt = sqrt(tmpu * F_VT)

      ! Air molecule mean free path
      rmfp = 2.0_fp * rmu / (rhoa * vt)

      ! Knudsen number
      rkn = rmfp / radius

      ! Cunningham slip correction factor (linearized form, Binkowski and Shankar 1995)
      bpm = 1.0_fp + 1.246_fp * rkn

      ! Fall speed (assumes Reynolds # < 0.01)
      vsettle = TWO_OVER_NINE * rhop * radius * radius * grav * bpm / rmu

      ! Check Reynolds number for drag correction
      re = 2.0_fp * rhoa * radius * vsettle / rmu

      ! If Re > 0.01 apply drag correction (Pruppacher and Klett regime 2, eq. 10-142)
      if (re > 0.01_fp) then
         x = log(24.0_fp * re / bpm)
         y = A0 + x * (A1 + x * (A2 + x * (A3 + x * (A4 + x * (A5 + A6 * x)))))
         re = exp(y) * bpm
         vsettle = 0.5_fp * rmu * re / (rhoa * radius)
      end if

   end subroutine settling_calc_vsettle

   !---------------------------------------------------------------------------
   !> Default settling solver with CFL-based sub-stepping.
   !! Replaces GOCART2G SettlingSolver, operates on 1D column arrays.
   !!
   !! @param[in]    nlev  Number of vertical levels
   !! @param[in]    cdt   Time step [s]
   !! @param[in]    delp  Pressure thickness [Pa]
   !! @param[in]    dz    Layer thickness [m]
   !! @param[in]    vs    Settling velocity [m/s]
   !! @param[inout] qa    Species concentration [kg/kg]
   !---------------------------------------------------------------------------
   subroutine settling_solver_default(nlev, cdt, delp, dz, vs, qa)

      integer, intent(in) :: nlev
      real(fp), intent(in) :: cdt
      real(fp), intent(in) :: delp(nlev)
      real(fp), intent(in) :: dz(nlev)
      real(fp), intent(in) :: vs(nlev)
      real(fp), intent(inout) :: qa(nlev)

      ! Local variables
      integer :: k, iit, nSubSteps
      real(fp) :: tau(nlev)
      real(fp) :: dt, dt_cfl

      ! Compute tau = vs/dz for each level (guard against zero dz)
      do k = 1, nlev
         if (dz(k) > 0.0_fp) then
            tau(k) = vs(k) / dz(k)
         else
            tau(k) = 0.0_fp
         end if
      end do

      ! CFL-based sub-stepping (guard against zero maxval)
      if (maxval(tau) <= 0.0_fp) then
         nSubSteps = 0
         dt = cdt
      else
         dt_cfl = 1.0_fp / maxval(tau)

         if (dt_cfl > cdt) then
            ! No need for time sub-splitting
            nSubSteps = 1
            dt = cdt
         else
            nSubSteps = ceiling(cdt / dt_cfl)
            dt = cdt / real(nSubSteps, fp)
         end if
      end if

      do iit = 1, nSubSteps
         ! Update bottom layer (index 1 in bottom-to-top ordering = surface)
         ! In the original GOCART2G code with top-to-bottom ordering, index 1 is the top.
         ! Here we keep the same mathematical structure but on 1D arrays.
         ! The settling moves mass downward (from higher index to lower index in bottom-to-top).
         ! Level 1 = bottom (surface), level nlev = top.
         ! Mass settles from level k to level k-1 (top to bottom = high index to low index).

         ! Update top layer (only loss, mass settles downward)
         qa(nlev) = qa(nlev) * (1.0_fp - dt * tau(nlev))

         ! Update interior and bottom layers (gain from above, loss downward)
         do k = nlev - 1, 1, -1
            qa(k) = qa(k) + (delp(k+1) / delp(k)) * (dt * tau(k+1)) * qa(k+1) &
               - dt * tau(k) * qa(k)
         end do
      end do

   end subroutine settling_solver_default

   !---------------------------------------------------------------------------
   !> UFS settling solver with enhanced numerical safeguards.
   !! Replaces GOCART2G SettlingSolverUFS, operates on 1D column arrays.
   !!
   !! @param[in]    nlev  Number of vertical levels
   !! @param[in]    cdt   Time step [s]
   !! @param[in]    delp  Pressure thickness [Pa]
   !! @param[in]    dz    Layer thickness [m]
   !! @param[in]    vs    Settling velocity [m/s]
   !! @param[inout] qa    Species concentration [kg/kg]
   !---------------------------------------------------------------------------
   subroutine settling_solver(nlev, cdt, delp, dz, vs, qa)

      integer, intent(in) :: nlev
      real(fp), intent(in) :: cdt
      real(fp), intent(in) :: delp(nlev)
      real(fp), intent(in) :: dz(nlev)
      real(fp), intent(in) :: vs(nlev)
      real(fp), intent(inout) :: qa(nlev)

      ! Local variables
      integer :: k, iit, nSubSteps
      real(fp) :: tau(nlev)
      real(fp) :: qa_old(nlev)
      real(fp) :: dt, dt_cfl, max_tau
      real(fp) :: transfer_factor, loss_factor

      real(fp), parameter :: eps = 1.0e-30_fp     ! Small number to prevent division by zero
      real(fp), parameter :: cfl_factor = 0.1_fp   ! CFL stability factor

      ! Compute tau = vs/dz for each level (guard against zero dz)
      do k = 1, nlev
         if (dz(k) > 0.0_fp) then
            tau(k) = vs(k) / dz(k)
         else
            tau(k) = 0.0_fp
         end if
      end do

      ! CFL-based sub-stepping with stability factor (guard against zero max_tau)
      max_tau = maxval(tau)
      if (max_tau <= 0.0_fp) then
         nSubSteps = 0
         dt = cdt
      else
         dt_cfl = cfl_factor / max_tau

         if (dt_cfl >= cdt) then
            ! No need for time sub-splitting
            nSubSteps = 0
            dt = cdt
         else
            nSubSteps = max(1, ceiling(cdt / dt_cfl))
            dt = cdt / real(nSubSteps, fp)
         end if
      end if

      ! Time integration with numerical safeguards
      ! In bottom-to-top ordering: level nlev = top, level 1 = bottom
      ! Mass settles from higher index to lower index
      do iit = 1, nSubSteps
         ! Store old values for mass transfer
         qa_old = qa

         ! Update top layer (only loss)
         loss_factor = max(0.0_fp, min(1.0_fp, dt * tau(nlev)))
         qa(nlev) = max(0.0_fp, qa(nlev) * (1.0_fp - min(loss_factor, 1.0_fp)))

         ! Update interior and bottom layers
         do k = nlev - 1, 1, -1
            loss_factor = max(0.0_fp, min(1.0_fp, dt * tau(k)))

            ! Check if pressure layers are valid
            if (delp(k+1) > eps .and. delp(k) > eps) then
               transfer_factor = (delp(k+1) / delp(k)) * dt * tau(k+1)
               qa(k) = max(0.0_fp, qa(k) * (1.0_fp - min(loss_factor, 1.0_fp))) + &
                  transfer_factor * qa_old(k+1)
            else
               qa(k) = max(0.0_fp, qa(k) * (1.0_fp - min(loss_factor, 1.0_fp)))
            end if
         end do
      end do

   end subroutine settling_solver

   !---------------------------------------------------------------------------
   !> Hygroscopic growth dispatcher.
   !! Replaces GOCART2G ParticleSwelling, operates on 1D column arrays.
   !!
   !! @param[in]  nlev        Number of vertical levels
   !! @param[in]  rh          Relative humidity [0-1]
   !! @param[in]  radius_dry  Dry particle radius [m]
   !! @param[in]  rhop_dry    Dry particle density [kg/m3]
   !! @param[out] radius      Wet particle radius [m]
   !! @param[out] rhop        Wet particle density [kg/m3]
   !! @param[in]  flag        Swelling method (0=none, 1=Fitzgerald, 2=Gerber, 3=Gerber NH4SO4, 4=PK2007)
   !---------------------------------------------------------------------------
   subroutine settling_particle_swelling(nlev, rh, radius_dry, rhop_dry, &
      radius, rhop, flag)

      integer, intent(in) :: nlev
      integer, intent(in) :: flag
      real(fp), intent(in) :: rh(nlev)
      real(fp), intent(in) :: radius_dry
      real(fp), intent(in) :: rhop_dry
      real(fp), intent(out) :: radius(nlev)
      real(fp), intent(out) :: rhop(nlev)

      select case (flag)
       case (0)
         radius = radius_dry
         rhop = rhop_dry

       case (1)
         call settling_swelling_fitzgerald(nlev, rh, radius_dry, rhop_dry, radius, rhop)

       case (2)
         call settling_swelling_gerber(nlev, rh, radius_dry, rhop_dry, radius, rhop)

       case (3)
         call settling_swelling_gerber_nh4so4(nlev, rh, radius_dry, rhop_dry, radius, rhop)

       case (4)
         call settling_swelling_pk2007(nlev, rh, radius_dry, rhop_dry, radius, rhop)

       case default
         radius = radius_dry
         rhop = rhop_dry
      end select

   end subroutine settling_particle_swelling

   !---------------------------------------------------------------------------
   !> Fitzgerald 1975 hygroscopic growth parameterization.
   !! Replaces GOCART2G ParticleSwelling_Fitzgerald, operates on 1D column arrays.
   !!
   !! Adjusts particle size for relative humidity effects based on
   !! Fitzgerald, Journal of Applied Meteorology, 1975.
   !!
   !! @param[in]  nlev        Number of vertical levels
   !! @param[in]  rh          Relative humidity [0-1]
   !! @param[in]  radius_dry  Dry particle radius [m]
   !! @param[in]  rhop_dry    Dry particle density [kg/m3]
   !! @param[out] radius      Wet particle radius [m]
   !! @param[out] rhop        Wet particle density [kg/m3]
   !---------------------------------------------------------------------------
   subroutine settling_swelling_fitzgerald(nlev, rh, radius_dry, rhop_dry, radius, rhop)

      integer, intent(in) :: nlev
      real(fp), intent(in) :: rh(nlev)
      real(fp), intent(in) :: radius_dry
      real(fp), intent(in) :: rhop_dry
      real(fp), intent(out) :: radius(nlev)
      real(fp), intent(out) :: rhop(nlev)

      ! Local variables
      ! Parameters from Fitzgerald 1975 for seasalt-like particles
      ! real(fp), parameter :: epsilon = 1.0_fp       ! Soluble fraction
      real(fp), parameter :: alphaNaCl = 1.35_fp

      real(fp) :: alpha, alpha1, beta, theta, sat, rrat
      integer :: k

      do k = 1, nlev
         radius(k) = radius_dry
         rhop(k) = rhop_dry

         sat = rh(k)

         if (sat > 0.80_fp) then
            ! Parameterization blows up for RH > 0.995
            sat = min(0.995_fp, sat)

            ! Beta parameter
            beta = exp((0.00077_fp * sat) / (1.009_fp - sat))

            ! Theta parameter
            if (sat <= 0.97_fp) then
               theta = 1.058_fp
            else
               theta = 1.058_fp - (0.0155_fp * (sat - 0.97_fp)) / (1.02_fp - sat**1.4_fp)
            end if

            ! Alpha parameter
            alpha1 = 1.2_fp * exp((0.066_fp * sat) / (theta - sat))
            ! Since epsilon == 1, simplified form:
            alpha = alphaNaCl * alpha1

            radius(k) = alpha * radius_dry**beta

            rrat = radius_dry / radius(k)
            rrat = rrat * rrat * rrat

            rhop(k) = rrat * rhop_dry + (1.0_fp - rrat) * RHOW
         end if
      end do

   end subroutine settling_swelling_fitzgerald

   !---------------------------------------------------------------------------
   !> Gerber 1985 hygroscopic growth parameterization.
   !! Replaces GOCART2G ParticleSwelling_Gerber, operates on 1D column arrays.
   !!
   !! @param[in]  nlev        Number of vertical levels
   !! @param[in]  rh          Relative humidity [0-1]
   !! @param[in]  radius_dry  Dry particle radius [m]
   !! @param[in]  rhop_dry    Dry particle density [kg/m3]
   !! @param[out] radius      Wet particle radius [m]
   !! @param[out] rhop        Wet particle density [kg/m3]
   !---------------------------------------------------------------------------
   subroutine settling_swelling_gerber(nlev, rh, radius_dry, rhop_dry, radius, rhop)

      integer, intent(in) :: nlev
      real(fp), intent(in) :: rh(nlev)
      real(fp), intent(in) :: radius_dry
      real(fp), intent(in) :: rhop_dry
      real(fp), intent(out) :: radius(nlev)
      real(fp), intent(out) :: rhop(nlev)

      ! Local variables
      ! Parameters from Gerber 1985 (units require radius in cm)
      real(fp), parameter :: c1 = 0.7674_fp
      real(fp), parameter :: c2 = 3.079_fp
      real(fp), parameter :: c3 = 2.573e-11_fp
      real(fp), parameter :: c4 = -1.424_fp

      real(fp) :: sat, rrat, rcm
      integer :: k

      do k = 1, nlev
         sat = max(rh(k), tiny(1.0_fp))  ! Avoid zero FPE
         sat = min(0.995_fp, sat)

         rcm = radius_dry * 100.0_fp  ! Radius in cm

         radius(k) = 0.01_fp * (c1 * rcm**c2 / (c3 * rcm**c4 - log10(sat)) &
            + rcm * rcm * rcm)**(1.0_fp / 3.0_fp)

         rrat = radius_dry / radius(k)
         rrat = rrat * rrat * rrat

         rhop(k) = rrat * rhop_dry + (1.0_fp - rrat) * RHOW
      end do

   end subroutine settling_swelling_gerber

   !---------------------------------------------------------------------------
   !> Gerber 1985 parameterization for Ammonium Sulfate.
   !! Replaces GOCART2G ParticleSwelling_Gerber_AmmoniumSulfate, operates on 1D column arrays.
   !!
   !! @param[in]  nlev        Number of vertical levels
   !! @param[in]  rh          Relative humidity [0-1]
   !! @param[in]  radius_dry  Dry particle radius [m]
   !! @param[in]  rhop_dry    Dry particle density [kg/m3]
   !! @param[out] radius      Wet particle radius [m]
   !! @param[out] rhop        Wet particle density [kg/m3]
   !---------------------------------------------------------------------------
   subroutine settling_swelling_gerber_nh4so4(nlev, rh, radius_dry, rhop_dry, radius, rhop)

      integer, intent(in) :: nlev
      real(fp), intent(in) :: rh(nlev)
      real(fp), intent(in) :: radius_dry
      real(fp), intent(in) :: rhop_dry
      real(fp), intent(out) :: radius(nlev)
      real(fp), intent(out) :: rhop(nlev)

      ! Local variables
      ! Parameters for ammonium sulfate from Gerber 1985 (units require radius in cm)
      real(fp), parameter :: SU_c1 = 0.4809_fp
      real(fp), parameter :: SU_c2 = 3.082_fp
      real(fp), parameter :: SU_c3 = 3.110e-11_fp
      real(fp), parameter :: SU_c4 = -1.428_fp

      real(fp) :: sat, rrat, rcm
      integer :: k

      do k = 1, nlev
         sat = max(rh(k), tiny(1.0_fp))  ! Avoid zero FPE
         sat = min(0.995_fp, sat)

         rcm = radius_dry * 100.0_fp  ! Radius in cm

         radius(k) = 0.01_fp * (SU_c1 * rcm**SU_c2 / (SU_c3 * rcm**SU_c4 - log10(sat)) &
            + rcm * rcm * rcm)**(1.0_fp / 3.0_fp)

         rrat = radius_dry / radius(k)
         rrat = rrat * rrat * rrat

         rhop(k) = rrat * rhop_dry + (1.0_fp - rrat) * RHOW
      end do

   end subroutine settling_swelling_gerber_nh4so4

   !---------------------------------------------------------------------------
   !> Petters and Kreidenweis (ACP 2007) hygroscopic growth parameterization.
   !! Replaces GOCART2G ParticleSwelling_PK2007, operates on 1D column arrays.
   !!
   !! @param[in]  nlev        Number of vertical levels
   !! @param[in]  rh          Relative humidity [0-1]
   !! @param[in]  radius_dry  Dry particle radius [m]
   !! @param[in]  rhop_dry    Dry particle density [kg/m3]
   !! @param[out] radius      Wet particle radius [m]
   !! @param[out] rhop        Wet particle density [kg/m3]
   !---------------------------------------------------------------------------
   subroutine settling_swelling_pk2007(nlev, rh, radius_dry, rhop_dry, radius, rhop)

      integer, intent(in) :: nlev
      real(fp), intent(in) :: rh(nlev)
      real(fp), intent(in) :: radius_dry
      real(fp), intent(in) :: rhop_dry
      real(fp), intent(out) :: radius(nlev)
      real(fp), intent(out) :: rhop(nlev)

      ! Local variables
      real(fp) :: sat, rrat
      integer :: k

      do k = 1, nlev
         sat = rh(k)
         sat = min(0.99_fp, sat)

         radius(k) = radius_dry * (1.0_fp + 1.19_fp * sat / (1.0_fp - sat))**(1.0_fp / 3.0_fp)

         rrat = radius_dry / radius(k)
         rrat = rrat * rrat * rrat

         rhop(k) = rrat * rhop_dry + (1.0_fp - rrat) * RHOW
      end do

   end subroutine settling_swelling_pk2007

end module SettlingPhysics_Mod
