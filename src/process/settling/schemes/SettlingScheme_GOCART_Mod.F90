!> \file SettlingScheme_GOCART_Mod.F90
!! \brief GOCART gravitational settling scheme
!!
!! Pure science kernel for gocart scheme in settling process.
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
!! Generated on: 2025-12-17T15:27:52.203209
!! Author: Wei Li
!! Reference: GOCART2G process library Chem_SettlingSimple function
module SettlingScheme_GOCART_Mod

   use precision_mod, only: fp
   use SettlingCommon_Mod, only: SettlingSchemeGOCARTConfig
   use Constants, only: PI  !load the constants needed for this scheme
   use GOCART2G_MieMod, only: GOCART2G_Mie  ! For Mie data in gocart scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor

contains

   !> Pure science computation for gocart scheme
   !!
   !! This is a pure computational kernel implementing GOCART gravitational settling scheme.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  airden    AIRDEN field [appropriate units]
   !! @param[in]  delp    DELP field [appropriate units]
   !! @param[in]  pmid    PMID field [appropriate units]
   !! @param[in]  rh    RH field [appropriate units]
   !! @param[in]  t    T field [appropriate units]
   !! @param[in]  tstep    Time step [s] - retrieved from process interface
   !! @param[in]  zmid    ZMID field [appropriate units]
   !! @param[in]  species_short_name    Species short_name property
   !! @param[in]  mie_data           Complete Mie data array from ChemState
   !! @param[in]  species_mie_map    Mapping from process species to MieData indices
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] settling_velocity_per_species_per_level    settling velocity per species per level [m/s] (num_layers, num_species)
   !! @param[inout] settling_flux_per_species    settling flux per species across column [kg/m2/s] (num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   pure subroutine compute_gocart( &
      num_layers, &
      num_species, &
      params, &
      airden, &
      delp, &
      pmid, &
      rh, &
      t, &
      tstep, &
      zmid, &
      species_short_name, &
      mie_data, &
      species_mie_map, &
      species_conc, &
      species_tendencies, &
      settling_velocity_per_species_per_level, &
      settling_flux_per_species, &
      diagnostic_species_id &
   )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(SettlingSchemeGOCARTConfig), intent(in) :: params
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: delp(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: pmid(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: rh(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: t(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      real(fp), intent(in) :: zmid(num_layers)    ! 3D atmospheric field
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      type(GOCART2G_Mie), intent(in) :: mie_data(:)  ! Complete Mie data array from ChemState
      integer, intent(in) :: species_mie_map(num_species)  ! Mapping from process species to MieData indices
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: settling_velocity_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: settling_flux_per_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: k, species_idx
      integer :: diag_idx  ! For diagnostic species indexing
      real(fp) :: base_emission_factor
      real(fp) :: environmental_factor
      real(fp) :: species_factor

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! TODO: Replace this generic implementation with your scheme's algorithm
         ! This is a placeholder that demonstrates the expected structure

         ! Initialize environmental factors
         environmental_factor = 1.0_fp

         ! Apply scheme-specific environmental responses based on meteorological fields
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how AIRDEN affects your emissions
         ! environmental_factor = environmental_factor * some_function(airden(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how DELP affects your emissions
         ! environmental_factor = environmental_factor * some_function(delp(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how PMID affects your emissions
         ! environmental_factor = environmental_factor * some_function(pmid(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how RH affects your emissions
         ! environmental_factor = environmental_factor * some_function(rh(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how T affects your emissions
         ! environmental_factor = environmental_factor * some_function(t(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how TSTEP affects your emissions
         ! environmental_factor = environmental_factor * some_function(tstep(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how ZMID affects your emissions
         ! environmental_factor = environmental_factor * some_function(zmid(k))

         ! Apply to each species
         do species_idx = 1, num_species
            ! Base emission factor (customize this for species-specific emissions)
            base_emission_factor = DEFAULT_SCALING

            ! Species-specific factor (customize based on species properties)
            species_factor = 1.0_fp  ! TODO: Add species-specific scaling

            ! Compute emission flux using your scheme's formula
            ! This is a simple example - replace with your actual algorithm
            species_tendencies(k, species_idx) = base_emission_factor * &
                                          environmental_factor * &
                                          species_factor * &
                                          (1.0_fp + species_conc(k, species_idx))

            ! Ensure non-negative emissions
            species_tendencies(k, species_idx) = max(0.0_fp, species_tendencies(k, species_idx))

            ! TODO: Update diagnostic fields here based on your scheme's requirements
            ! Each process should implement custom diagnostic calculations
            ! Example patterns:
            ! Per-species-per-level diagnostic: 2D array (levels, species)
            if (present(settling_velocity_per_species_per_level) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom settling velocity per species per level calculation
                     settling_velocity_per_species_per_level(k, diag_idx) = species_tendencies(k, species_idx) * 1.0_fp  ! TODO: Replace with actual calculation
                     exit
                  end if
               end do
            end if
            ! Per-species diagnostic: only update for diagnostic species
            if (present(settling_flux_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom settling flux per species across column calculation
                     settling_flux_per_species(diag_idx) = species_tendencies(k, species_idx) * 1.0_fp  ! TODO: Replace with actual calculation
                     exit
                  end if
               end do
            end if
         end do

      end do

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !> Example helper function for environmental response
   pure function compute_environmental_response_gocart(met_value, reference_value) result(factor)
      real(fp), intent(in) :: met_value       ! Meteorological value
      real(fp), intent(in) :: reference_value ! Reference value
      real(fp) :: factor

      ! Simple exponential response - customize for your scheme
      factor = exp((met_value - reference_value) / reference_value)
      factor = max(0.0_fp, min(10.0_fp, factor))  ! Reasonable bounds
   end function compute_environmental_response_gocart

   !> Example helper function for species-specific scaling
   pure function compute_species_scaling_gocart(species_idx, params) result(scaling)
      integer, intent(in) :: species_idx
      type(SettlingSchemeGOCARTConfig), intent(in) :: params
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

   end function compute_species_scaling_gocart

end module SettlingScheme_GOCART_Mod