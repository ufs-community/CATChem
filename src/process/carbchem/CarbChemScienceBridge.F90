module CarbChemScienceBridge_Mod
   use iso_c_binding
   use catchem_bridge_precision, only: fp
   use catchem_bridge_constants, only: g0
   use CarbChemScheme_GOCART_Mod, only: compute_gocart
   use CarbChemCommon_Mod, only: CarbChemSchemeGOCARTConfig

   implicit none
   private

contains

   subroutine run_carbchem_science_bridge( &
      n_cols, n_levels, n_species, dt, &
      active_scheme, diagnostics, &
      gocart_time_days_hydrophobic_to_hydrophilic, &
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

      ! Scheme tuning options staged by CarbChemProcess::init from the runtime
      ! YAML.  The C++ layer owns parsing and validation; the bridge only
      ! applies them onto the GOCART configuration type.
      real(c_double), value :: gocart_time_days_hydrophobic_to_hydrophilic

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

      ! Per-process diagnostics.  The diag_* pointers reference
      ! DiagnosticManager field storage; they are null (and n_diag_species is
      ! 0) when diagnostics are disabled or no species are selected, so the
      ! pointers are only c_f_pointer'd when n_diag_species > 0.
      ! diagnostic_species_id holds 1-based GLOBAL catalog positions (the
      ! space the GOCART scheme matches on); the size-1 dummy is never
      ! dereferenced when n_diag_species == 0 (mirrors the other bridges).
      integer(c_int), value :: n_diag_species
      integer(c_int), intent(in) :: diagnostic_species_id(max(n_diag_species,1))

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

      real(fp) :: col_prod_mass(n_levels, max(n_diag_species,1))
      real(fp) :: col_loss_flux(max(n_diag_species,1))
      real(fp) :: col_phobic_mass(n_levels, max(n_diag_species,1))
      real(fp) :: col_phobic_flux(max(n_diag_species,1))

      ! Control structures
      type(CarbChemSchemeGOCARTConfig) :: gocart_config
      character(len=32) :: local_scheme
      integer :: icol, i, j

      ! Extract scheme string
      local_scheme = ""
      do i = 1, 32
         if (active_scheme(i) == c_null_char) exit
         local_scheme(i:i) = active_scheme(i)
      end do
      local_scheme = trim(local_scheme)

      ! Apply the YAML tuning option staged by the C++ process layer so the
      ! scheme no longer runs on compiled defaults alone.
      gocart_config%time_days_hydrophobic_to_hydrophilic = real(gocart_time_days_hydrophobic_to_hydrophilic, fp)

      ! Map Pointers
      call c_f_pointer(airden, f_airden, [n_cols, n_levels])
      call c_f_pointer(delp, f_delp, [n_cols, n_levels])
      call c_f_pointer(pmid, f_pmid, [n_cols, n_levels])

      call c_f_pointer(species_t_chem_loss, f_t_chem_loss, [n_species])
      call c_f_pointer(species_names_char, f_names_char, [32, n_species])

      call c_f_pointer(conc, f_conc, [n_cols, n_levels, n_species])
      call c_f_pointer(tendency, f_tendency, [n_cols, n_levels, n_species])

      if (diagnostics /= 0 .and. n_diag_species > 0) then
         call c_f_pointer(diag_prod_mass, f_diag_prod_mass, [n_cols, n_levels, n_diag_species])
         call c_f_pointer(diag_loss_flux, f_diag_loss_flux, [n_cols, n_diag_species])
         call c_f_pointer(diag_phobic_mass, f_diag_phobic_mass, [n_cols, n_levels, n_diag_species])
         call c_f_pointer(diag_phobic_flux, f_diag_phobic_flux, [n_cols, n_diag_species])
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

         ! GOCART carbon returns updated aerosol concentrations in ug/kg for computed species.
         ! Extract the tendency from the updated concentrations and apply it.
         do i = 1, n_species
            if (any(abs(col_tendency(:, i)) > 1.0e-32_fp)) then
               ! col_tendency contains the NEW concentration in ug/kg.
               ! Calculate rate of change and update conc/tendency in-place.
               f_tendency(icol, :, i) = (real(col_tendency(:, i), c_double) - f_conc(icol, :, i)) / dt
               f_conc(icol, :, i)     = real(col_tendency(:, i), c_double)
            else
               f_tendency(icol, :, i) = 0.0_c_double
            end if
         end do

         if (diagnostics /= 0 .and. n_diag_species > 0) then
            ! The scheme already scattered each species into its diag_idx
            ! slot (1..n_diag_species), so the copy is a straight 1:1 write of
            ! the ndiag-wide buffers into the registered field storage.
            f_diag_prod_mass(icol, :, 1:n_diag_species) = real(col_prod_mass(:, 1:n_diag_species), c_double)
            f_diag_loss_flux(icol, 1:n_diag_species) = real(col_loss_flux(1:n_diag_species), c_double)
            f_diag_phobic_mass(icol, :, 1:n_diag_species) = real(col_phobic_mass(:, 1:n_diag_species), c_double)
            f_diag_phobic_flux(icol, 1:n_diag_species) = real(col_phobic_flux(1:n_diag_species), c_double)
         end if

      end do

   end subroutine run_carbchem_science_bridge

end module CarbChemScienceBridge_Mod
