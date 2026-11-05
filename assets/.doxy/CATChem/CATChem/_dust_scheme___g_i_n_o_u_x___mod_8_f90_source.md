

# File DustScheme\_GINOUX\_Mod.F90

[**File List**](files.md) **>** [**dust**](dir_1c14dfbaca1e3f4c2e26e74290119ebd.md) **>** [**schemes**](dir_11b2254edcf6ee5df673de29b129f986.md) **>** [**DustScheme\_GINOUX\_Mod.F90**](_dust_scheme___g_i_n_o_u_x___mod_8_f90.md)

[Go to the documentation of this file](_dust_scheme___g_i_n_o_u_x___mod_8_f90.md)


```Fortran

module dustscheme_ginoux_mod

   use precision_mod, only: fp
   use dustcommon_mod, only: dustschemeginouxconfig

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_ginoux

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter    :: SSM_THRESH  = 1.0e-02_fp  ! Minimum erodibility threshold

contains

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
      real(fp) :: ginoux_scaling
      real(fp) :: u_thresh0
      real(fp) :: u_thresh
      real(fp) :: w10m
      real(fp) :: emission_temp

      !needs to reinitialize otherwise the skip condition below will cause weird maps.
      if (present(utar_threshold_per_bin)) utar_threshold_per_bin = 0.0_fp
      if (present(dust_emission_total)) dust_emission_total = 0.0_fp
      if (present(dust_emission_per_bin)) dust_emission_per_bin = 0.0_fp

      ! Skip criteria evaluation
      skip = (lwi /= 1)  !land = 1, water = 0, ice = 2

      if (.not. skip) then
         skip = (ssm < ssm_thresh)
      endif

      ! Don't do dust over frozen soil
      !--------------------------------
      if (tskin <= 273.15_fp) then
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
      ginoux_scaling = min(1.0_fp, max(0.0_fp, 1.0_fp - frlake - frsno) ) * ssm

      ! get 10m mean wind speed
      w10m = sqrt(u10m ** 2 + v10m ** 2)

      ! Main computation loop
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species

            !initialize emission_temp to zero for this species and layer
            emission_temp = 0.0_fp

            ! get threshold friction velocity following MB97
            call mb97_threshold_velocity(species_density(species_idx), airden(1), species_radius(species_idx), g0, u_thresh0)

            ! add the moisture correction following Ginoux et al. (2001)
            u_thresh = max(0.0_fp, u_thresh0 * (1.2_fp + 0.2_fp*log10(max(1.e-3_fp, gwettop))) )

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

   pure subroutine mb97_threshold_velocity(soil_density, air_density, radius, g0, ustar_threshold)
      ! USES
      IMPLICIT NONE

      ! Input Parameters
      !-----------------
      real(fp), intent(in) :: radius
      real(fp), intent(in) :: soil_density
      real(fp), intent(in) :: air_density
      real(fp), intent(in) :: g0

      ! Output Parameters
      !------------------
      real(fp), intent(out) :: ustar_threshold

      ! Local Variables
      !-----------------
      real(fp) :: diameter

      diameter = 2.0_fp * radius * 1.0e-6_fp 
      ustar_threshold = 0.13_fp * sqrt(soil_density*g0*diameter/air_density) &
         * sqrt(1.0_fp + 6.e-7_fp/(soil_density*g0*diameter**2.5_fp)) &
         / sqrt(1.928_fp*(1331.0_fp*(100._fp*diameter)**1.56_fp+0.38_fp)**0.092_fp - 1.0_fp)

   end subroutine mb97_threshold_velocity

end module dustscheme_ginoux_mod
```


