

# File SeaSaltScheme\_GONG97\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**seasalt**](dir_3e6fd2e121e43ca7d4114b6c0b4e05b6.md) **>** [**schemes**](dir_ec083b49fedbd640552af85049fd7226.md) **>** [**SeaSaltScheme\_GONG97\_Mod.F90**](_sea_salt_scheme___g_o_n_g97___mod_8_f90.md)

[Go to the documentation of this file](_sea_salt_scheme___g_o_n_g97___mod_8_f90.md)


```Fortran

module seasaltscheme_gong97_mod

   use precision_mod, only: fp,zero, rae
   use seasaltcommon_mod, only: seasaltschemegong97config
   use constants, only: pi  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gong97

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor

contains

   subroutine compute_gong97( &
      num_layers, &
      num_species, &
      params, &
      frocean, &
      frseaice, &
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
      real(fp), intent(in) :: frocean  ! Surface field - scalar
      real(fp), intent(in) :: frseaice  ! Surface field - scalar
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
      logical :: do_seasalt
      integer :: n, ir
      real(fp) :: w10m
      integer, parameter :: nr = 10
      real, parameter    :: r80fac = 1.65_fp           
      real(fp) :: DryRadius
      real(fp) :: DeltaDryRadius
      real(fp) :: rwet, drwet
      real(fp) :: NumberEmissions
      real(fp) :: MassEmissions
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
      real(fp) :: fsstemis
      real(fp) :: fhoppel
      real(fp) :: scale


      ! Initialize output (pure subroutines must initialize all outputs)
      rc = 0
      mass_emission_flux = 0.0_fp
      numb_emission_flux = 0.0_fp
      massemissions = 0.0_fp
      numberemissions = 0.0_fp
      gweibull = 1.0_fp
      fsstemis = 1.0_fp
      fhoppel = 1.0_fp


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

         ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
         do k = 1, num_layers

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

            scale = scale * gweibull * fsstemis * params%scale_factor

            ! Apply to each species
            do n = 1, num_species
               ! delta dry radius
               !-----------------
               deltadryradius = (species_upper_radius(n) - species_lower_radius(n) )/ nr

               ! Dry Radius Substep
               !-------------------
               dryradius = species_lower_radius(n) + 0.5 * deltadryradius

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

   function seasaltemissiongong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

      real(fp), intent(in) :: r
      real(fp), intent(in) :: dr
      real(fp), intent(in) :: w
      real(fp), intent(in) :: scalefac
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

         if (wm > 0.012_fp) then
            k = 0.94_fp * sqrt(wm)
            c = wm / gamma(1.0_fp + 1.0_fp / k)
            x = (wt / c) ** k
            a = 3.41_fp / k + 1.0_fp
            gweibull = (c / wm) ** 3.41_fp * igamma(a, x, rc)
         endif
      endif


   end subroutine weibulldistribution

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


