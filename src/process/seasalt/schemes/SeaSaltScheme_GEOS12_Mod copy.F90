!> \file SeaSaltScheme_GEOS12_Mod.F90
!! \brief GEOS-Chem 2012 sea salt emission scheme with observational constraints
!!
!! Pure science kernel for geos12 scheme in seasalt process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_geos12 (search for "TODO")
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
!! Generated on: 2025-09-15T17:20:44.139071
!! Author: Barry Baker
!! Reference: Jaeglé et al. [2011]
module SeaSaltScheme_GEOS12_Mod

   use precision_mod, only: fp, zero
   use SeaSaltCommon_Mod, only: SeaSaltSchemeGEOS12Config
   use Constants, only: PI  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_geos12

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor

contains

   !> Pure science computation for geos12 scheme
   !!
   !! This is a pure computational kernel implementing GEOS-Chem 2012 sea salt emission scheme with observational constraints.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  frocean    FROCEAN field [appropriate units]
   !! @param[in]  frseaice    FRSEAICE field [appropriate units]
   !! @param[in]  sst    SST field [appropriate units]
   !! @param[in]  ustar    USTAR field [appropriate units]
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
   pure subroutine compute_geos12( &
      num_layers, &
      num_species, &
      params, &
      frocean, &
      frseaice, &
      sst, &
      ustar, &
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
      type(SeaSaltSchemeGEOS12Config), intent(in) :: params
      real(fp), intent(in) :: frocean  ! Surface field - scalar
      real(fp), intent(in) :: frseaice  ! Surface field - scalar
      real(fp), intent(in) :: sst  ! Surface field - scalar
      real(fp), intent(in) :: ustar  ! Surface field - scalar
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
      integer :: k, species_idx, RC
      integer :: diag_idx  ! For diagnostic species indexing
      real(fp) :: base_emission_factor
      real(fp) :: environmental_factor
      real(fp) :: species_factor

      logical :: do_seasalt                            !< Enable Dust Calculation Flag
      integer :: n, ir                                 !< Loop counter
      integer, parameter :: nr = 10                    !< Number of (linear) sub-size bins
      real(fp), parameter    :: r80fac = 1.65              !< ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
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
      real(fp) :: fsstemis
      real(fp) :: fhoppel
      real(fp) :: scale

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.
      RC = 0
      mass_emission_flux = 0.0_fp
      numb_emission_flux = 0.0_fp
      MassEmissions = 0.0_fp
      NumberEmissions = 0.0_fp
      fsstemis = 1.0_fp
      fhoppel = 1.0_fp

      do_seasalt = .true. ! Default value for all cases

      ! Don't do Sea Salt over land
      !----------------------------------------------------------------
      scale = FROCEAN - FRSEAICE
      if (scale <= 0.0_fp) then
         do_seasalt = .False.
      endif

      if (do_seasalt) then
         ! GEOS 12 Params
         !---------------
         scalefac = 33.0e3_fp
         rpow     = 3.45_fp
         exppow   = 1.607_fp
         wpow     = 3.41_fp - 1._fp

         ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
         do k = 1, num_layers

            ! TODO: Replace this generic implementation with your scheme's algorithm
            ! This is a placeholder that demonstrates the expected structure
            ! Get Jeagle SST Correction
            call jeagleSSTcorrection(fsstemis, SST,1, RC)
            if (RC /= 0) then
               RC = -1
               !print *, 'Error in jeagleSSTcorrection'
               return
            endif

            scale = scale * fsstemis * params%scale_factor

            ! Apply to each species
            do n = 1, num_species
               ! delta dry radius
               !-----------------
               DeltaDryRadius = (species_upper_radius(n) - species_lower_radius(n) )/ nr

               ! Dry Radius Substep
               !-------------------
               DryRadius = species_lower_radius(n) + 0.5_fp * DeltaDryRadius

               ! Mass scale fcator
               MassScaleFac = scalefac * 4._fp/3._fp*PI*species_density(n)*(DryRadius**3._fp) * 1.e-18_fp

               do ir = 1, nr ! SubSteps

                  ! Effective Wet Radius in Sub Step
                  rwet  = r80fac * DryRadius

                  ! Effective Delta Wet Radius
                  drwet = r80fac * DeltaDryRadius

                  aFac = 4.7_fp*(1._fp + 30._fp*rwet)**(-0.017_fp*rwet**(-1.44_fp))
                  bFac = (0.380_fp-log10(rwet))/0.65_fp

                  ! Number emissions flux (# m-2 s-1)
                  NumberEmissions = NumberEmissions + SeasaltEmissionGong( rwet, drwet, USTAR, scalefac, &
                     aFac, bFac, rpow, exppow, wpow )

                  ! Mass emissions flux (kg m-2 s-1)
                  MassEmissions = MassEmissions + SeasaltEmissionGong( rwet, drwet, USTAR, MassScaleFac, &
                     aFac, bFac, rpow, exppow, wpow )

                  DryRadius = DryRadius + DeltaDryRadius

               enddo ! ir loop

               mass_emission_flux(k, n) = MassEmissions * scale
               numb_emission_flux(k, n) = NumberEmissions * scale
               ! Reset for next species
               MassEmissions = 0.0_fp
               NumberEmissions = 0.0_fp

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

      end if ! do_seasalt

   end subroutine compute_geos12

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !> Example helper function for environmental response
   pure function compute_environmental_response_geos12(met_value, reference_value) result(factor)
      real(fp), intent(in) :: met_value       ! Meteorological value
      real(fp), intent(in) :: reference_value ! Reference value
      real(fp) :: factor

      ! Simple exponential response - customize for your scheme
      factor = exp((met_value - reference_value) / reference_value)
      factor = max(0.0_fp, min(10.0_fp, factor))  ! Reasonable bounds
   end function compute_environmental_response_geos12

   !> Example helper function for species-specific scaling
   pure function compute_species_scaling_geos12(species_idx, params) result(scaling)
      integer, intent(in) :: species_idx
      type(SeaSaltSchemeGEOS12Config), intent(in) :: params
      real(fp) :: scaling

      ! Species-specific scaling - customize for your scheme
      select case (species_idx)
       case (1)
         scaling = 1.0_fp    ! First species baseline
       case (2:3)
         scaling = 0.5_fp    ! Reduced emission for species 2-3
       case default
         scaling = 0.1_fp    ! Low emission for other species
      end select

   end function compute_species_scaling_geos12

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
   pure subroutine jeagleSSTcorrection(fsstemis, sst, sstFlag, rc)

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
      RC = -1 ! Error code
      fsstemis = 1.0_fp

      fsstemis = ZERO
      tskin_c  = sst - 273.15_fp
      if (sstFlag .eq. 1) then
         fsstemis = max(0.0_fp,(0.3_fp + 0.1_fp*tskin_c - 0.0076_fp*tskin_c**2 + 0.00021_fp*tskin_c**3))
      else
         ! temperature range (0, 36) C
         tskin_c = max(-0.1_fp, Tskin_c)
         tskin_c = min(36.0_fp, tskin_c)

         fsstemis = (-1.107211_fp -0.010681_fp * tskin_c -0.002276_fp * tskin_c**2.0_fp &
            + 60.288927_fp*1.0_fp/(40.0_fp - tskin_c))
         fsstemis = max(0.0_fp, fsstemis)
         fsstemis = min(7.0_fp, fsstemis)
      endif

      RC = 0
   end subroutine jeagleSSTcorrection

   !>
   !! \brief Function to compute sea salt emissions following the Gong style parameterization.
   !!
   !! Functional form is from Gong 2003:
   !!   \f$dN/dr = scalefac * 1.373 * (w^wpow) * (r^-aFac) * (1+0.057*r^rpow) * 10^(exppow*exp(-bFac^2))\f$
   !! where r is the particle radius at 80% RH, dr is the size bin width at 80% RH, and w is the wind speed
   !!
   !! \ingroup catchem_seasalt_process
   !!!>
   pure function SeasaltEmissionGong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

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
      SeasaltEmissionGong = 0.

      !  Particle size distribution function
      SeasaltEmissionGong = scalefac * 1.373_fp*r**(-aFac)*(1._fp+0.057_fp*r**rpow) &
         *10._fp**(exppow*exp(-bFac**2._fp))*dr
      !  Apply wind speed function
      SeasaltEmissionGong = w**wpow * SeasaltEmissionGong

   end function SeasaltEmissionGong

end module SeaSaltScheme_GEOS12_Mod
