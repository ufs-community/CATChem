!> \file DustScheme_GINOUX_Mod.F90
!! \brief Ginoux dust emission scheme
!!
!! Pure science kernel for ginoux scheme in dust process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_ginoux (search for "TODO")
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
!! Generated on: 2026-04-17T13:57:10.254102
!! Author: Barry Baker & Wei Li
!! Reference: Ginoux et al. [2001]
module DustScheme_GINOUX_Mod

   use precision_mod, only: fp
   use DustCommon_Mod, only: DustSchemeGINOUXConfig

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_ginoux

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter    :: SSM_THRESH  = 1.0E-02_fp  ! Minimum erodibility threshold

contains

   !> Pure science computation for ginoux scheme
   !!
   !! This is a pure computational kernel implementing Ginoux dust emission scheme.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  g0    Required constant from Constants module
   !! @param[in]  airden    AIRDEN field [appropriate units]
   !! @param[in]  frlake    FRLAKE field [appropriate units]
   !! @param[in]  frsno    FRSNO field [appropriate units]
   !! @param[in]  gwettop    GWETTOP field [appropriate units]
   !! @param[in]  lwi    LWI field [appropriate units]
   !! @param[in]  ssm    SSM field [appropriate units]
   !! @param[in]  tskin    TSKIN field [appropriate units]
   !! @param[in]  u10m    U10M field [appropriate units]
   !! @param[in]  v10m    V10M field [appropriate units]
   !! @param[in]  species_density    Species density property
   !! @param[in]  species_radius    Species radius property
   !! @param[in]  species_conc   Species concentrations [ppm or ug/kg] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] dust_emission_total    Total dust emissions for all bins [kg/m2/s]
   !! @param[inout] dust_emission_per_bin    Dust emission flux per bin [kg/m2/s] (num_species)
   !! @param[inout] utar_threshold    Threshold friction velocity to have dust emission [m/s]
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
   pure subroutine compute_ginoux( &
      num_layers, &
      num_species, &
      params, &
      g0, &
      airden, &
      frlake, &
      frsno, &
      gwettop, &
      lwi, &
      ssm, &
      tskin, &
      u10m, &
      v10m, &
      species_density, &
      species_radius, &
      species_conc, &
      species_tendencies, &
      dust_emission_total, &
      dust_emission_per_bin, &
      utar_threshold_per_bin, &
      diagnostic_species_id &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(DustSchemeGINOUXConfig), intent(in) :: params
      real(fp), intent(in) :: g0  ! Required constant from Constants module
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: frlake  ! Surface field - scalar
      real(fp), intent(in) :: frsno  ! Surface field - scalar
      real(fp), intent(in) :: gwettop  ! Surface field - scalar
      integer, intent(in) :: lwi  ! Surface field - scalar
      real(fp), intent(in) :: ssm  ! Surface field - scalar
      real(fp), intent(in) :: tskin  ! Surface field - scalar
      real(fp), intent(in) :: u10m  ! Surface field - scalar
      real(fp), intent(in) :: v10m  ! Surface field - scalar
      real(fp), intent(in) :: species_density(:)  ! Species density property
      real(fp), intent(in) :: species_radius(:)  ! Species radius property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: dust_emission_total
      real(fp), intent(inout), optional :: dust_emission_per_bin(:)
      real(fp), intent(inout), optional :: utar_threshold_per_bin(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: k, species_idx
      integer :: diag_idx  ! For diagnostic species indexing
      logical :: skip  ! Flag to determine if we should skip computation
      real(fp) :: ginoux_scaling                       !< Ginoux scaling
      real(fp) :: u_thresh0                            !< Dry threshold wind speed [m/s]
      real(fp) :: u_thresh                             !< Moisture Corrected threshold wind speed [m/s]
      real(fp) :: w10m                                 !< 10m wind speed [m/s]
      real(fp) :: emission_temp                        !< Temporary variable for emission calculation

      !needs to reinitialize otherwise the skip condition below will cause weird maps.
      if (present(utar_threshold_per_bin)) utar_threshold_per_bin = 0.0_fp
      if (present(dust_emission_total)) dust_emission_total = 0.0_fp
      if (present(dust_emission_per_bin)) dust_emission_per_bin = 0.0_fp

      ! Skip criteria evaluation
      skip = (LWI /= 1)  !land = 1, water = 0, ice = 2

      if (.not. skip) then
         skip = (SSM < SSM_THRESH)
      endif

      ! Don't do dust over frozen soil
      !--------------------------------
      if (TSKIN <= 273.15_fp) then
         skip = .true.
      endif

      ! Don't do dust if surface is wet
      !--------------------------------
      if (gwettop >= 0.5_fp) then
         skip = .true.
      endif

      ! Skip computation if criteria not met
      if (skip) then
         return
      end if

      ! get the scaling factor following Ginoux et al. (2001)
      ! Note the GOCART2G version does not have the SSM factor
      !Note not using (1 - frlake) * (1 - frsno) as GOCART
      !ginoux_scaling = (1 - frlake) * (1 - frsno) * SSM
      ginoux_scaling = min(1.0_fp, max(0.0_fp, 1.0_fp - frlake - frsno) ) * SSM

      ! get 10m mean wind speed
      w10m = sqrt(U10M ** 2 + V10M ** 2)

      ! Main computation loop
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species

            !initialize emission_temp to zero for this species and layer
            emission_temp = 0.0_fp

            ! get threshold friction velocity following MB97
            call MB97_threshold_velocity(species_density(species_idx), AIRDEN(1), species_radius(species_idx), g0, u_thresh0)

            ! add the moisture correction following Ginoux et al. (2001)
            u_thresh = max(0.0_fp, u_thresh0 * (1.2_fp + 0.2_fp*log10(max(1.e-3_fp, GWETTOP))) )

            ! Compute emission flux
            emission_temp = 0.0_fp
            if (w10m .gt. u_thresh) then
               emission_temp = ginoux_scaling * w10m ** 2 * max(0.0_fp,(w10m - u_thresh) )  ! kg/m2/s
               !TODO: Note Chu_DU is used in GOCART2G for the conversion from du_src
               !The Chu_DU list in GOCART goes through the Chem_UtilResVal function, after which all bins have the
               !same value before the 1e-9 conversion.
               !we do not have du_src input and use SSM instead in ginoux_scaling calculation above.
               emission_temp = emission_temp * params%Ch_DU(species_idx) * 1.0e-9
            endif

            species_tendencies(k, species_idx) = max(0.0_fp, emission_temp)

            ! Update scheme-specific diagnostic fields
            ! Per-species diagnostic: only update for diagnostic species
            if (present(dust_emission_total)) then
               ! Add your custom total dust emissions for all bins calculation
               dust_emission_total = dust_emission_total + species_tendencies(k, species_idx)
            end if
            ! Per-species diagnostic: only update for diagnostic species
            if (present(dust_emission_per_bin) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom dust emission flux per bin calculation
                     dust_emission_per_bin(diag_idx) = species_tendencies(k, species_idx)
                     exit
                  end if
               end do
            end if

            ! Update scheme-specific diagnostic fields
            if (present(utar_threshold_per_bin) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom friction velocity threshold per bin to initiate dust emission calculation
                     utar_threshold_per_bin(diag_idx) = u_thresh
                     exit
                  end if
               end do
            end if

         end do ! species_idx loop
      end do ! k loop

   end subroutine compute_ginoux

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !>
   !! \brief Computes the Threshold Friction Velocity from MB97
   !!
   !! Marticorena, B. and Bergametti, G.: Modeling the atmospheric dust cycle:
   !! 1. Design of a soil-derived dust emission scheme,
   !! J. Geophys. Res.-Atmos., 100, 16415–16430, https://doi.org/10.1029/95JD00690, 1995 | TODO fix with correct reference
   !!
   !! \param soil_density soil density
   !! \param air_density air density
   !! \param radius particle radius
   !! \param ustar_threshold threshold friction velocity
   !!
   !! \ingroup catchem_dust_process
   !!!>
   pure subroutine MB97_threshold_velocity(soil_density, air_density, radius, g0, ustar_threshold)
      ! USES
      IMPLICIT NONE

      ! Input Parameters
      !-----------------
      real(fp), intent(in) :: radius       !< particle radius
      real(fp), intent(in) :: soil_density !< soil density
      real(fp), intent(in) :: air_density  !< air density
      real(fp), intent(in) :: g0           !< gravitational acceleration

      ! Output Parameters
      !------------------
      real(fp), intent(out) :: ustar_threshold !< threshold friction velocity

      ! Local Variables
      !-----------------
      real(fp) :: diameter !< diameter of particle [m]

      diameter = 2.0_fp * radius * 1.0e-6_fp !< convert radius to meters
      ustar_threshold = 0.13_fp * sqrt(soil_density*g0*diameter/air_density) &
         * sqrt(1.0_fp + 6.e-7_fp/(soil_density*g0*diameter**2.5_fp)) &
         / sqrt(1.928_fp*(1331.0_fp*(100._fp*diameter)**1.56_fp+0.38_fp)**0.092_fp - 1.0_fp)

   end subroutine MB97_threshold_velocity

end module DustScheme_GINOUX_Mod
