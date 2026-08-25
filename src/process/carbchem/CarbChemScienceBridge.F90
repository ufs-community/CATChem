module CarbChemScienceBridge_Mod
   use iso_c_binding
   use precision_mod, only: fp
   use Constants, only: g0
   use CarbChemScheme_GOCART_Mod, only: compute_gocart
   use CarbChemCommon_Mod, only: CarbChemSchemeGOCARTConfig

   implicit none
   private

contains

   subroutine run_carbchem_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      active_scheme, diagnostics, &
      year, month, day, hour, minute, second, &
      airden, delp, pmid, &
      species_t_chem_loss, species_names_char, &
      conc, tendency, &
      diag_prod_mass, diag_loss_flux, diag_phobic_mass, diag_phobic_flux, &
      diagnostic_species_id, n_diag_species &
      ) bind(C, name="run_carbchem_science_bridge")

      ! C-interoperable dimensions and metadata
      integer(c_int), value :: n_cols, n_levels, n_species
      real(c_double), value :: dt
      character(kind=c_char), intent(in) :: active_scheme(*)
      integer(c_int), value :: diagnostics
      integer(c_int), value :: year, month, day, hour, minute, second

      ! C++ Raw Pointers
      type(c_ptr), value :: airden
      type(c_ptr), value :: delp
      type(c_ptr), value :: pmid
      type(c_ptr), value :: species_t_chem_loss
      type(c_ptr), value :: species_names_char
      type(c_ptr), value :: conc
      type(c_ptr), value :: tendency

      ! Diagnostics Pointers
      type(c_ptr), value :: diag_prod_mass
      type(c_ptr), value :: diag_loss_flux
      type(c_ptr), value :: diag_phobic_mass
      type(c_ptr), value :: diag_phobic_flux

      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(n_diag_species)

      ! Local Fortran Pointers for multidimensional mapping
      real(c_double), pointer :: f_airden(:,:), f_delp(:,:), f_pmid(:,:)
      real(c_double), pointer :: f_t_chem_loss(:)
      character(kind=c_char), pointer :: f_names_char(:,:)
      real(c_double), pointer :: f_conc(:,:,:), f_tendency(:,:,:)
      real(c_double), pointer :: f_diag_prod_mass(:,:,:)
      real(c_double), pointer :: f_diag_loss_flux(:,:)
      real(c_double), pointer :: f_diag_phobic_mass(:,:,:)
      real(c_double), pointer :: f_diag_phobic_flux(:,:)

      ! Solver slices
      real(fp) :: col_airden(n_levels)
      real(fp) :: col_delp(n_levels)
      real(fp) :: col_pmid(n_levels)
      real(fp) :: col_t_chem_loss(n_species)
      character(len=32) :: col_names(n_species)
      real(fp) :: col_conc(n_levels, n_species)
      real(fp) :: col_tendency(n_levels, n_species)

      real(fp) :: col_prod_mass(n_levels, n_species)
      real(fp) :: col_loss_flux(n_species)
      real(fp) :: col_phobic_mass(n_levels, n_species)
      real(fp) :: col_phobic_flux(n_species)

      ! Control structures
      type(CarbChemSchemeGOCARTConfig) :: gocart_config
      character(len=32) :: local_scheme
      integer :: icol, i, j, k

      ! Extract scheme string
      local_scheme = ""
      do i = 1, 32
         if (active_scheme(i) == c_null_char) exit
         local_scheme(i:i) = active_scheme(i)
      end do
      local_scheme = trim(local_scheme)

      ! Map Pointers
      call c_f_pointer(airden, f_airden, [n_cols, n_levels])
      call c_f_pointer(delp, f_delp, [n_cols, n_levels])
      call c_f_pointer(pmid, f_pmid, [n_cols, n_levels])

      call c_f_pointer(species_t_chem_loss, f_t_chem_loss, [n_species])
      call c_f_pointer(species_names_char, f_names_char, [32, n_species])

      call c_f_pointer(conc, f_conc, [n_cols, n_levels, n_species])
      call c_f_pointer(tendency, f_tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0) then
         call c_f_pointer(diag_prod_mass, f_diag_prod_mass, [n_cols, n_levels, n_species])
         call c_f_pointer(diag_loss_flux, f_diag_loss_flux, [n_cols, n_species])
         call c_f_pointer(diag_phobic_mass, f_diag_phobic_mass, [n_cols, n_levels, n_species])
         call c_f_pointer(diag_phobic_flux, f_diag_phobic_flux, [n_cols, n_species])
      end if

      ! Map metadata
      do i = 1, n_species
         col_t_chem_loss(i) = real(f_t_chem_loss(i), fp)
         col_names(i) = ""
         do j = 1, 32
            col_names(i)(j:j) = f_names_char(j, i)
         end do
         col_names(i) = trim(adjustl(col_names(i)))
      end do

      ! Iterate Columns
      do icol = 1, n_cols

         col_airden(:) = real(f_airden(icol, :), fp)
         col_delp(:) = real(f_delp(icol, :), fp)
         col_pmid(:) = real(f_pmid(icol, :), fp)

         col_conc(:, :) = real(f_conc(icol, :, :), fp)
         col_tendency(:, :) = 0.0_fp

         col_prod_mass(:, :) = 0.0_fp
         col_loss_flux(:) = 0.0_fp
         col_phobic_mass(:, :) = 0.0_fp
         col_phobic_flux(:) = 0.0_fp

         if (local_scheme == "gocart") then
            call compute_gocart( &
               n_levels, n_species, gocart_config, &
               g0, year, month, day, hour, minute, second, &
               col_airden, col_delp, col_pmid, real(dt, fp), &
               col_t_chem_loss, col_names, &
               col_conc, col_tendency, &
               Production_mass_per_species_per_level=col_prod_mass, &
               loss_flux_per_species=col_loss_flux, &
               PhobicToPhilic_mass_per_species_per_level=col_phobic_mass, &
               PhobicToPhilic_flux_per_species=col_phobic_flux, &
               diagnostic_species_id=diagnostic_species_id)
         end if

         ! Copy back tendency
         f_tendency(icol, :, :) = f_tendency(icol, :, :) + real(col_tendency(:, :), c_double)

         if (diagnostics /= 0) then
            f_diag_prod_mass(icol, :, :) = real(col_prod_mass(:, :), c_double)
            f_diag_loss_flux(icol, :) = real(col_loss_flux(:), c_double)
            f_diag_phobic_mass(icol, :, :) = real(col_phobic_mass(:, :), c_double)
            f_diag_phobic_flux(icol, :) = real(col_phobic_flux(:), c_double)
         end if

      end do

   end subroutine run_carbchem_science_bridge

end module CarbChemScienceBridge_Mod
