module DustScienceBridge_Mod
   use iso_c_binding
   use precision_mod, only: fp
   use Constants, only: g0
   use DustScheme_FENGSHA_Mod, only: compute_fengsha
   use DustScheme_GINOUX_Mod, only: compute_ginoux
   use DustCommon_Mod, only: DustSchemeFENGSHAConfig, DustSchemeGINOUXConfig

   implicit none
   private

contains

   subroutine run_dust_science_bridge( &
      n_cols, n_levels, n_species, n_soil, dt, &
      active_scheme, diagnostics, &
      airden, clayfrac, frlake, frsno, gvf, lai, lwi, rdrag, sandfrac, &
      soilm, ssm, tskin, u10m, v10m, ustar, ustar_threshold, z0, &
      species_density, species_radius, species_lower_radius, species_upper_radius, &
      conc, tendency, &
      diag_emission_total, diag_emission_bin, diag_horizontal_flux, diag_moisture_correction, diag_effective_threshold, diag_utar_threshold, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_dust_science_bridge")

      ! C-interoperable metadata
      integer(c_int), value :: n_cols, n_levels, n_species, n_soil
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: active_scheme(*)
      integer(c_int), value :: diagnostics

      ! 2D/3D Pointers from C++ Kokkos state
      type(c_ptr), value :: airden
      type(c_ptr), value :: clayfrac
      type(c_ptr), value :: frlake
      type(c_ptr), value :: frsno
      type(c_ptr), value :: gvf
      type(c_ptr), value :: lai
      type(c_ptr), value :: lwi
      type(c_ptr), value :: rdrag
      type(c_ptr), value :: sandfrac
      type(c_ptr), value :: soilm
      type(c_ptr), value :: ssm
      type(c_ptr), value :: tskin
      type(c_ptr), value :: u10m
      type(c_ptr), value :: v10m
      type(c_ptr), value :: ustar
      type(c_ptr), value :: ustar_threshold
      type(c_ptr), value :: z0

      ! Species properties
      type(c_ptr), value :: species_density
      type(c_ptr), value :: species_radius
      type(c_ptr), value :: species_lower_radius
      type(c_ptr), value :: species_upper_radius

      ! Multi-dimensional Views
      type(c_ptr), value :: conc
      type(c_ptr), value :: tendency

      ! Diagnostics
      type(c_ptr), value :: diag_emission_total
      type(c_ptr), value :: diag_emission_bin
      type(c_ptr), value :: diag_horizontal_flux
      type(c_ptr), value :: diag_moisture_correction
      type(c_ptr), value :: diag_effective_threshold
      type(c_ptr), value :: diag_utar_threshold
      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      ! Local Fortran Pointers
      real(c_double), pointer :: f_airden(:,:)
      real(c_double), pointer :: f_clayfrac(:), f_frlake(:), f_frsno(:), f_gvf(:), f_lai(:)
      integer(c_int), pointer :: f_lwi(:)
      real(c_double), pointer :: f_rdrag(:), f_sandfrac(:), f_soilm(:,:), f_ssm(:), f_tskin(:)
      real(c_double), pointer :: f_u10m(:), f_v10m(:), f_ustar(:), f_ustar_threshold(:), f_z0(:)
      real(c_double), pointer :: f_conc(:,:,:), f_tendency(:,:,:)

      real(c_double), pointer :: f_density(:), f_radius(:), f_lower_radius(:), f_upper_radius(:)

      real(c_double), pointer :: f_diag_emission_total(:), f_diag_emission_bin(:,:)
      real(c_double), pointer :: f_diag_horizontal_flux(:), f_diag_moisture_correction(:)
      real(c_double), pointer :: f_diag_effective_threshold(:), f_diag_utar_threshold(:,:)

      ! Local Slices for computation
      real(fp) :: col_airden(n_levels)
      real(fp) :: col_soilm(n_soil)
      real(fp) :: col_conc(n_levels, n_species)
      real(fp) :: col_tendency(n_levels, n_species)

      real(fp) :: f_species_density(n_species)
      real(fp) :: f_species_radius(n_species)
      real(fp) :: f_species_lower_radius(n_species)
      real(fp) :: f_species_upper_radius(n_species)

      real(fp) :: col_emission_total
      real(fp) :: col_emission_bin(n_species)
      real(fp) :: col_horizontal_flux
      real(fp) :: col_moisture_correction
      real(fp) :: col_effective_threshold
      real(fp) :: col_utar_threshold(n_species)

      type(DustSchemeFENGSHAConfig) :: fengsha_config
      type(DustSchemeGINOUXConfig)  :: ginoux_config
      character(len=32) :: local_scheme
      integer :: icol, i

      ! Map Scheme Name
      local_scheme = ""
      do i = 1, 32
         if (active_scheme(i) == c_null_char) exit
         local_scheme(i:i) = active_scheme(i)
      end do
      local_scheme = trim(local_scheme)

      ! Pointer Associations
      call c_f_pointer(airden, f_airden, [n_cols, n_levels])
      call c_f_pointer(clayfrac, f_clayfrac, [n_cols])
      call c_f_pointer(frlake, f_frlake, [n_cols])
      call c_f_pointer(frsno, f_frsno, [n_cols])
      call c_f_pointer(gvf, f_gvf, [n_cols])
      call c_f_pointer(lai, f_lai, [n_cols])
      call c_f_pointer(lwi, f_lwi, [n_cols])
      call c_f_pointer(rdrag, f_rdrag, [n_cols])
      call c_f_pointer(sandfrac, f_sandfrac, [n_cols])
      call c_f_pointer(soilm, f_soilm, [n_cols, n_soil])
      call c_f_pointer(ssm, f_ssm, [n_cols])
      call c_f_pointer(tskin, f_tskin, [n_cols])
      call c_f_pointer(u10m, f_u10m, [n_cols])
      call c_f_pointer(v10m, f_v10m, [n_cols])
      call c_f_pointer(ustar, f_ustar, [n_cols])
      call c_f_pointer(ustar_threshold, f_ustar_threshold, [n_cols])
      call c_f_pointer(z0, f_z0, [n_cols])

      call c_f_pointer(conc, f_conc, [n_cols, n_levels, n_species])
      call c_f_pointer(tendency, f_tendency, [n_cols, n_levels, n_species])

      call c_f_pointer(species_density, f_density, [n_species])
      call c_f_pointer(species_radius, f_radius, [n_species])
      call c_f_pointer(species_lower_radius, f_lower_radius, [n_species])
      call c_f_pointer(species_upper_radius, f_upper_radius, [n_species])

      if (diagnostics /= 0) then
         call c_f_pointer(diag_emission_total, f_diag_emission_total, [n_cols])
         call c_f_pointer(diag_emission_bin, f_diag_emission_bin, [n_cols, n_species])
         call c_f_pointer(diag_horizontal_flux, f_diag_horizontal_flux, [n_cols])
         call c_f_pointer(diag_moisture_correction, f_diag_moisture_correction, [n_cols])
         call c_f_pointer(diag_effective_threshold, f_diag_effective_threshold, [n_cols])
         call c_f_pointer(diag_utar_threshold, f_diag_utar_threshold, [n_cols, n_species])
      end if

      ! Cast species properties to fp precision
      f_species_density = real(f_density, fp)
      f_species_radius  = real(f_radius, fp)
      f_species_lower_radius = real(f_lower_radius, fp)
      f_species_upper_radius = real(f_upper_radius, fp)

      ! Process Columns
      do icol = 1, n_cols

         col_airden(:) = real(f_airden(icol, :), fp)
         col_soilm(:)  = real(f_soilm(icol, :), fp)
         col_conc(:,:) = real(f_conc(icol, :, :), fp)
         col_tendency(:,:) = 0.0_fp

         col_emission_total = 0.0_fp
         col_emission_bin(:) = 0.0_fp
         col_horizontal_flux = 0.0_fp
         col_moisture_correction = 0.0_fp
         col_effective_threshold = 0.0_fp
         col_utar_threshold(:) = 0.0_fp

         if (local_scheme == "fengsha") then
            call compute_fengsha( &
               n_levels, n_species, fengsha_config, g0, &
               col_airden, real(f_clayfrac(icol), fp), real(f_frlake(icol), fp), real(f_frsno(icol), fp), &
               real(f_gvf(icol), fp), real(f_lai(icol), fp), int(f_lwi(icol)), real(f_rdrag(icol), fp), &
               real(f_sandfrac(icol), fp), col_soilm, real(f_ssm(icol), fp), real(f_tskin(icol), fp), &
               real(f_ustar(icol), fp), real(f_ustar_threshold(icol), fp), real(f_z0(icol), fp), &
               f_species_radius, f_species_lower_radius, f_species_upper_radius, &
               col_conc, col_tendency, &
               dust_emission_total=col_emission_total, &
               dust_emission_per_bin=col_emission_bin, &
               dust_horizontal_flux=col_horizontal_flux, &
               dust_moisture_correction=col_moisture_correction, &
               dust_effective_threshold=col_effective_threshold, &
               diagnostic_species_id=diagnostic_species_id(1:n_diag_species))
         else if (local_scheme == "ginoux") then
            call compute_ginoux( &
               n_levels, n_species, ginoux_config, g0, &
               col_airden, real(f_frlake(icol), fp), real(f_frsno(icol), fp), 0.0_fp, &
               int(f_lwi(icol)), real(f_ssm(icol), fp), real(f_tskin(icol), fp), &
               real(f_u10m(icol), fp), real(f_v10m(icol), fp), &
               f_species_density, f_species_radius, &
               col_conc, col_tendency, &
               dust_emission_total=col_emission_total, &
               dust_emission_per_bin=col_emission_bin, &
               utar_threshold_per_bin=col_utar_threshold, &
               diagnostic_species_id=diagnostic_species_id(1:n_diag_species))
         end if

         ! Write tendencies
         f_tendency(icol, :, :) = f_tendency(icol, :, :) + real(col_tendency, c_double)

         ! Write diagnostics
         if (diagnostics /= 0) then
            f_diag_emission_total(icol) = real(col_emission_total, c_double)
            f_diag_emission_bin(icol, :) = real(col_emission_bin, c_double)
            f_diag_horizontal_flux(icol) = real(col_horizontal_flux, c_double)
            f_diag_moisture_correction(icol) = real(col_moisture_correction, c_double)
            f_diag_effective_threshold(icol) = real(col_effective_threshold, c_double)
            f_diag_utar_threshold(icol, :) = real(col_utar_threshold, c_double)
         end if

      end do

   end subroutine run_dust_science_bridge

end module DustScienceBridge_Mod
