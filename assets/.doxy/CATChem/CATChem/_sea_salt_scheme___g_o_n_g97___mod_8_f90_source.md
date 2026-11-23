

# File SeaSaltScheme\_GONG97\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**seasalt**](dir_3e6fd2e121e43ca7d4114b6c0b4e05b6.md) **>** [**schemes**](dir_ec083b49fedbd640552af85049fd7226.md) **>** [**SeaSaltScheme\_GONG97\_Mod.F90**](_sea_salt_scheme___g_o_n_g97___mod_8_f90.md)

[Go to the documentation of this file](_sea_salt_scheme___g_o_n_g97___mod_8_f90.md)


```Fortran
!> \file SeaSaltScheme_GONG97_Mod.F90
!! \brief Gong 1997 sea salt emission scheme
!!
!! Pure science kernel for gong97 scheme in seasalt process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_gong97 (search for "TODO")
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
!! Generated on: 2025-09-16T00:40:10.206610
!! Author: Barry Baker
!! Reference: Gong et al. [1997]
module seasaltscheme_gong97_mod

   use catchem_bridge_precision, only: fp,zero, rae
   use seasaltcommon_mod, only: seasaltschemegong97config
   use catchem_bridge_constants, only: pi  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gong97

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor

contains

   !> Pure science computation for gong97 scheme
   !!
   !! This is a pure computational kernel implementing Gong 1997 sea salt emission scheme.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  PI    Required constant from Constants module
   !! @param[in]  frocean    FROCEAN field [appropriate units]
   !! @param[in]  frseaice    FRSEAICE field [appropriate units]
   !! @param[in]  lat    LAT field [appropriate units]
   !! @param[in]  lon    LON field [appropriate units]
   !! @param[in]  sst    SST field [appropriate units]
   !! @param[in]  u10m    U10M field [appropriate units]
   !! @param[in]  v10m    V10M field [appropriate units]
   !! @param[in]  species_density    Species density property
   !! @param[in]  species_radius    Species radius property
   !! @param[in]  species_lower_radius    Species lower_radius property
   !! @param[in]  species_upper_radius    Species upper_radius property
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] seasalt_mass_emission_total    Total mass emission diagnostic [ug/m2/s]
   !! @param[inout] seasalt_number_emission_total  Total number emission diagnostic [#/m2/s]
   !! @param[inout] seasalt_mass_emission_per_bin   Mass emission per bin diagnostic [kg/m2/s] (num_species)
   !! @param[inout] seasalt_number_emission_per_bin Number emission per bin diagnostic [#/m2/s] (num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   subroutine compute_gong97( &
      num_layers, &
      num_species, &
      params, &
      PI, &
      frocean, &
      frseaice, &
      lat, &
      lon, &
      sst, &
      u10m, &
      v10m, &
      species_density, &
      species_radius, &
      species_lower_radius, &
      species_upper_radius, &
      species_conc, &
      species_tendencies, &
      seasalt_mass_emission_total, &
      seasalt_number_emission_total, &
      seasalt_mass_emission_per_bin, &
      seasalt_number_emission_per_bin, &
      diagnostic_species_id  &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(SeaSaltSchemeGONG97Config), intent(in) :: params
      real(fp), intent(in) :: PI  ! Required constant from Constants module
      real(fp), intent(in) :: frocean  ! Surface field - scalar
      real(fp), intent(in) :: frseaice  ! Surface field - scalar
      real(fp), intent(in) :: lat  ! Surface field - scalar
      real(fp), intent(in) :: lon  ! Surface field - scalar
      real(fp), intent(in) :: sst  ! Surface field - scalar
      real(fp), intent(in) :: u10m  ! Surface field - scalar
      real(fp), intent(in) :: v10m  ! Surface field - scalar
      real(fp), intent(in) :: species_density(num_species)  ! Species density property
      real(fp), intent(in) :: species_radius(num_species)  ! Species radius property
      real(fp), intent(in) :: species_lower_radius(num_species)  ! Species lower_radius property
      real(fp), intent(in) :: species_upper_radius(num_species)  ! Species upper_radius property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: seasalt_mass_emission_total
      real(fp), intent(inout), optional :: seasalt_number_emission_total
      real(fp), intent(inout), optional :: seasalt_mass_emission_per_bin(:)
      real(fp), intent(inout), optional :: seasalt_number_emission_per_bin(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: k, RC
      integer :: diag_idx  ! For diagnostic species indexing
      logical :: do_seasalt                            !< Enable Dust Calculation Flag
      integer :: n, ir                                 !< Loop counter
      real(fp) :: w10m                                 !< 10m wind speed [m/s]
      integer, parameter :: nr = 10                    !< Number of (linear) sub-size bins
      real(fp), parameter    :: r80fac = 1.65_fp           !< ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
      real(fp) :: DryRadius                            !< sub-bin radius         (dry, um)
      real(fp) :: DeltaDryRadius                       !< sub-bin radius spacing (dry, um)
      real(fp) :: rwet, drwet                          !< sub-bin radius spacing (rh=80%, um)
      real(fp) :: NumberEmissions                      !< sub-bin number emission rate [#/m2/s]
      real(fp) :: MassEmissions                        !< sub-bin number emission rate [kg/m2/s]
      real(fp) :: mass_emission_flux(num_layers, num_species)
      real(fp) :: numb_emission_flux(num_layers, num_species)
      real(fp) :: aFac
      real(fp) :: bFac
      real(fp) :: scalefac
      real(fp) :: rpow
      real(fp) :: exppow
      real(fp) :: wpow
      real(fp) :: MassScaleFac
      real(fp) :: gweibull
      real(fp) :: deep_lakes_mask, dummylon
      real(fp) :: fsstemis
      real(fp) :: fhoppel
      real(fp) :: scale

      ! `species_radius`, `species_conc` are part of the shared scheme calling convention and
      ! intentionally unused by this scheme; reference them so the interface
      ! stays uniform without an unused-dummy-argument warning.
      associate(unused_species_radius => species_radius); end associate
      associate(unused_species_conc => species_conc); end associate

      ! Initialize output (pure subroutines must initialize all outputs)
      rc = 0
      mass_emission_flux = 0.0_fp
      numb_emission_flux = 0.0_fp
      massemissions = 0.0_fp
      numberemissions = 0.0_fp
      gweibull = 1.0_fp
      fsstemis = 1.0_fp
      fhoppel = 1.0_fp
      deep_lakes_mask = 1.0_fp
      dummylon = lon

      !initialize diagnostics if present
      if (present(seasalt_mass_emission_total)) seasalt_mass_emission_total = 0.0_fp
      if (present(seasalt_number_emission_total)) seasalt_number_emission_total = 0.0_fp
      if (present(seasalt_mass_emission_per_bin)) seasalt_mass_emission_per_bin = 0.0_fp
      if (present(seasalt_number_emission_per_bin)) seasalt_number_emission_per_bin = 0.0_fp

      do_seasalt = .true. ! Default value for all cases

      ! Don't do Sea Salt over land
      !----------------------------------------------------------------
      scale = frocean - frseaice
      if (scale <= 0.0_fp) then
         do_seasalt = .false.
      endif

      if (do_seasalt) then
         ! Gong 1997 Params
         !-----------------
         scalefac = 3.0_fp
         rpow     = 1.05_fp
         exppow   = 1.19_fp
         wpow     = 3.41_fp

         ! get 10m mean wind speed
         !------------------------
         w10m = sqrt(u10m ** 2 + v10m ** 2)

         ! Weibull Distribution following Fan and Toon 2011 if WeibullFlag
         !----------------------------------------------------------------------------
         call weibulldistribution(gweibull, params%weibull_flag, w10m, rc)
         if (rc /= 0) then
            rc = -1
            print *, 'Error in weibullDistribution' !TODO: this may not allowed in a pure function
            return
         endif

         ! Get Jeagle SST Correction
         call jeaglesstcorrection(fsstemis, sst,1, rc)
         if (rc /= 0) then
            rc = -1
            print *, 'Error in jeagleSSTcorrection'
            return
         endif

         ! Deep Lakes Mask for Great Lakes lon = [93W,75W], lat = [40.5N, 50N]
         if( dummylon < 0.0 ) dummylon = dummylon + 360.0_fp
         if (lat >= 40.5_fp .and. lat <= 50.0_fp .and. dummylon >= 267.0_fp .and. dummylon <= 285.0_fp) then
            deep_lakes_mask = 0.0_fp
         endif
         ! The Caspian Sea: lon = [45.0, 56], lat = 35, 48]
         if (lat >= 35.0_fp .and. lat <= 48.0_fp .and. dummylon >= 45.0_fp .and. dummylon <= 56.0_fp) then
            deep_lakes_mask = 0.0_fp
         endif

         !update scale factor with all corrections
         scale = min(max(0.0_fp, scale * deep_lakes_mask), 1.0_fp) * gweibull * fsstemis * params%scale_factor

         ! Main computation loop
         do k = 1, num_layers

            ! Apply to each species
            do n = 1, num_species
               ! delta dry radius
               !-----------------
               deltadryradius = (species_upper_radius(n) - species_lower_radius(n) )/ nr

               ! Dry Radius Substep
               !-------------------
               dryradius = species_lower_radius(n) + 0.5_fp * deltadryradius

               do ir = 1, nr ! SubSteps

                  ! Mass scale fcator
                  massscalefac = scalefac * 4._fp/3._fp*pi*species_density(n)*(dryradius**3._fp) * 1.e-18_fp

                  ! Effective Wet Radius in Sub Step
                  rwet  = r80fac * dryradius

                  ! Effective Delta Wet Radius
                  drwet = r80fac * deltadryradius

                  afac = 3.0_fp
                  bfac = (0.380_fp-log10(rwet))/0.650_fp

                  ! Number emissions flux (# m-2 s-1)
                  numberemissions = numberemissions + seasaltemissiongong( rwet, drwet, w10m, scalefac, &
                     afac, bfac, rpow, exppow, wpow )

                  ! Mass emissions flux (kg m-2 s-1)
                  massemissions = massemissions + seasaltemissiongong( rwet, drwet, w10m, massscalefac, &
                     afac, bfac, rpow, exppow, wpow )

                  dryradius = dryradius + deltadryradius

               enddo ! ir loop

               mass_emission_flux(k, n) = massemissions * scale
               numb_emission_flux(k, n) = numberemissions * scale

               ! Reset for next species
               massemissions = 0.0_fp
               numberemissions = 0.0_fp

               ! Ensure non-negative emissions
               species_tendencies(k, n) = max(0.0_fp, mass_emission_flux(k, n))

               ! TODO: Update diagnostic fields here based on your scheme's requirements
               ! Each process should implement custom diagnostic calculations
               ! Example patterns:
               if (present(seasalt_mass_emission_total)) then
                  seasalt_mass_emission_total = seasalt_mass_emission_total + mass_emission_flux(k, n)
               end if
               if (present(seasalt_number_emission_total)) then
                  seasalt_number_emission_total = seasalt_number_emission_total + numb_emission_flux(k, n)
               end if
               if (present(seasalt_mass_emission_per_bin) .and. present(diagnostic_species_id)) then
                  ! Find position of this species in diagnostic_species_id array
                  do diag_idx = 1, size(diagnostic_species_id)
                     if (diagnostic_species_id(diag_idx) == n) then
                        ! Add your custom sea salt mass emission flux per bin calculation
                        seasalt_mass_emission_per_bin(diag_idx) = mass_emission_flux(k, n)
                        exit
                     end if
                  end do
               end if
               if (present(seasalt_number_emission_per_bin) .and. present(diagnostic_species_id)) then
                  ! Find position of this species in diagnostic_species_id array
                  do diag_idx = 1, size(diagnostic_species_id)
                     if (diagnostic_species_id(diag_idx) == n) then
                        ! Add your custom sea salt mass emission flux per bin calculation
                        seasalt_number_emission_per_bin(diag_idx) = numb_emission_flux(k, n)
                        exit
                     end if
                  end do
               end if
            end do

         end do

      end if

   end subroutine compute_gong97

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !>
   !! \brief Jeagle et al. 2012 SST correction
   !!
   !! Jaeglé, L., Quinn, P. K., Bates, T. S., Alexander, B., and Lin, J.-T.:
   !! Global distribution of sea salt aerosols: new constraints from in situ and remote
   !! sensing observations, Atmos. Chem. Phys., 11, 3137–3157,
   !! https://doi.org/10.5194/acp-11-3137-2011, 2011.
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   subroutine jeaglesstcorrection(fsstemis, sst, sstFlag, rc)

      ! !USES:
      implicit NONE

      ! !INPUT/OUTPUT PARAMETERS:
      real(fp), intent(inout) :: fsstemis     !
      real(fp), intent(in)  :: sst  ! surface temperature (K)
      integer, intent(in) :: sstFlag

      ! !OUTPUT PARAMETERS:
      integer, optional, intent(out) :: rc
      !EOP

      ! !Local Variables
      real(fp) :: tskin_c
      !EOP
      !-------------------------------------------------------------------------
      !  Begin...
      rc = -1 ! Error code
      fsstemis = 1.0_fp

      fsstemis = zero
      tskin_c  = sst - 273.15_fp
      if (sstflag .eq. 1) then
         fsstemis = max(0.0_fp,(0.3_fp + 0.1_fp*tskin_c - 0.0076_fp*tskin_c**2 + 0.00021_fp*tskin_c**3))
      else
         ! temperature range (0, 36) C
         tskin_c = max(-0.1_fp, tskin_c)
         tskin_c = min(36.0_fp, tskin_c)

         fsstemis = (-1.107211_fp -0.010681_fp * tskin_c -0.002276_fp * tskin_c**2.0_fp &
            + 60.288927_fp*1.0_fp/(40.0_fp - tskin_c))
         fsstemis = max(0.0_fp, fsstemis)
         fsstemis = min(7.0_fp, fsstemis)
      endif

      rc = 0
   end subroutine jeaglesstcorrection

   !>
   !! \brief Function to compute sea salt emissions following the Gong style parameterization.
   !!
   !! Functional form is from Gong 2003:
   !!   \f$dN/dr = scalefac * 1.373 * (w^wpow) * (r^-aFac) * (1+0.057*r^rpow) * 10^(exppow*exp(-bFac^2))\f$
   !! where r is the particle radius at 80% RH, dr is the size bin width at 80% RH, and w is the wind speed
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   function seasaltemissiongong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

      real(fp), intent(in) :: r         !< Wet particle radius [um]
      real(fp), intent(in) :: dr        !< Wet particle bin width [um]
      real(fp), intent(in) :: w         !< Grid box mean wind speed [m s-1] (10-m or ustar wind)
      real(fp), intent(in) :: scalefac  !< scale factor
      real(fp), intent(in) :: aFac
      real(fp), intent(in) :: bFac
      real(fp), intent(in) :: rpow
      real(fp), intent(in) :: exppow
      real(fp), intent(in) :: wpow
      real(fp)             :: SeasaltEmissionGong

      !  Initialize
      seasaltemissiongong = 0.

      !  Particle size distribution function
      seasaltemissiongong = scalefac * 1.373_fp*r**(-afac)*(1._fp+0.057_fp*r**rpow) &
         *10._fp**(exppow*exp(-bfac**2._fp))*dr
      !  Apply wind speed function
      seasaltemissiongong = w**wpow * seasaltemissiongong

   end function seasaltemissiongong

   !>
   !! \brief Calculate the weibull distribution for 10m wind speed (u10, v10)
   !!
   !! The Weibull distribution correction ends up being a multiplicative constant
   !! (g) times our present source function (see Eq. 12 in Fan & Toon, 2011 and notes for
   !! (9/22/11). This constant is derived from the incomplete and complete forms of the gamma
   !! function, hence the utilities pasted below.  The Weibull function and shape
   !! parameters (k, c) assumed are from Justus 1978.
   !!
   !! \param[inout] gweibull Multiplicative constant
   !! \param[in]    weibullFlag  Flag for weibull correction
   !! \param[in]    wm 10m wind speed
   !! \param[out]   RC Return Code
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   subroutine weibulldistribution(gweibull, weibullFlag, wm, RC)

      implicit none

      ! Input/Output
      !-------------
      real(fp), intent(inout) :: gweibull

      ! Input
      !------
      logical,  intent(in)    :: weibullFlag
      real(fp), intent(in)    :: wm

      ! Output
      !-------
      integer,  intent(out)   :: RC

      ! Local Variables
      real(fp) :: a, c, k, wt, x
      character(len=256) :: errMsg, thisLoc !  needed for error handling thisLoc
      ! Initialize
      errmsg = ''
      thisloc = ' -> at weibullDistribution (in util/metutils_mod.F90)'
      rc = 0
      gweibull = 1.0_fp

      wt = 4.0_fp

      if (weibullflag) then
         gweibull = 0.0_fp

         if (wm > 0.01_fp) then
            k = 0.94_fp * sqrt(wm)
            c = wm / gamma(1.0_fp + 1.0_fp / k)
            x = (wt / c) ** k
            a = 3.41_fp / k + 1.0_fp
            gweibull = (c / wm) ** 3.41_fp * igamma(a, x, rc)
         endif
      endif


   end subroutine weibulldistribution

   !>
   !! \brief Calculate the incomplete Gamma function
   !!
   !! The incomplete Gamma function is defined as
   !! \int_x^\infty t^{A-1}\exp(-t) dt
   !!
   !! \param[in]    A
   !! \param[in]    X
   !! \param[out]   RC
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   real(fp) function igamma(A, X, RC)

      IMPLICIT NONE

      REAL(fp), INTENT(in) :: A
      REAL(fp), INTENT(IN) :: X
      integer, intent(out) :: rc

      ! LOCAL VARIABLE
      REAL(fp) :: XAM, GIN, S, R, T0
      INTEGER K
      rc = 0
      igamma = 0

      xam=-x+a*log(x)
      IF (xam.GT.700.0_fp.OR.a.GT.170.0_fp) THEN
         WRITE(*,*)'IGAMMA: a and/or x too large, X = ', x
         WRITE(*,*) 'A = ', a
         rc = -1
         return
      ENDIF

      IF (rae(x, 0.0_fp)) THEN
         !IF ( X == 0.0_fp) THEN
         igamma=gamma(a)

      ELSE IF (x.LE.1.0_fp+a) THEN
         s=1.0_fp/a
         r=s
         DO  k=1,60
            r=r*x/(a+k)
            s=s+r
            IF (abs(r/s).LT.1.0e-15_fp) EXIT
         END DO
         gin=exp(xam)*s
         igamma=gamma(a)-gin
      ELSE IF (x.GT.1.0_fp+a) THEN
         t0=0.0_fp
         DO k=60,1,-1
            t0=(k-a)/(1.0_fp+k/(x+t0))
         end do

         igamma=exp(xam)/(x+t0)

      ENDIF

   end function igamma

end module seasaltscheme_gong97_mod
```


