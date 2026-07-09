

# File CarbChemScheme\_GOCART\_Mod.F90

[**File List**](files.md) **>** [**carbchem**](dir_5dbdd03f815becc4c35f94e0692b4e09.md) **>** [**schemes**](dir_2aa0506a6ee0350ff25dc7952f19a30f.md) **>** [**CarbChemScheme\_GOCART\_Mod.F90**](_carb_chem_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the documentation of this file](_carb_chem_scheme___g_o_c_a_r_t___mod_8_f90.md)


```Fortran

module carbchemscheme_gocart_mod

   use precision_mod, only: fp, rae
   use carbchemcommon_mod, only: carbchemschemegocartconfig
   use gocart2g_process, only: carbonchemloss, phobictophilic, chem_utilidow, chem_utilcdow

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: plid = 0.01_fp    ! Pressure lid [hPa]

contains

   subroutine compute_gocart( &
      num_layers, &
      num_species, &
      params, &
      g0, &
      year, &
      month, &
      day, &
      hour, &
      minute, &
      second, &
      airden, &
      delp, &
      pmid, &
      tstep, &
      species_t_chem_loss, &
      species_short_name, &
      species_conc, &
      species_tendencies, &
      Production_mass_per_species_per_level, &
      loss_flux_per_species, &
      PhobicToPhilic_mass_per_species_per_level, &
      PhobicToPhilic_flux_per_species, &
      diagnostic_species_id &
      )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(CarbChemSchemeGOCARTConfig), intent(in) :: params
      real(fp), intent(in) :: g0  ! Required constant from Constants module
      integer, intent(in) :: year  ! Time parameter from TimeState
      integer, intent(in) :: month  ! Time parameter from TimeState
      integer, intent(in) :: day  ! Time parameter from TimeState
      integer, intent(in) :: hour  ! Time parameter from TimeState
      integer, intent(in) :: minute  ! Time parameter from TimeState
      integer, intent(in) :: second  ! Time parameter from TimeState
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: delp(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: pmid(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      real(fp), intent(in) :: species_t_chem_loss(:)  ! Species t_chem_loss property
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: Production_mass_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: loss_flux_per_species(:)
      real(fp), intent(inout), optional :: PhobicToPhilic_mass_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: PhobicToPhilic_flux_per_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: n, species_idx, phobic_species_idx, philic_species_idx
      integer :: klid, diag_idx, curr_idx  ! For diagnostic species indexing
      real(fp), pointer :: GOCART_RHOA(:,:,:)
      real(fp), pointer :: GOCART_DELP(:,:,:)
      real(fp), pointer :: GOCART_PRESS(:,:,:)
      real(fp), pointer :: flux_toPhilic(:,:) ! For phobic to philic conversion flux [kg/m2/s]
      real(fp), pointer, dimension(:,:,:,:)  :: intPtr_phobic_philic !mass loss [kg/kg]
      real(fp), pointer, dimension(:,:,:)  :: fluxout  !Mass lost by chemistry [kg/m^2/s]
      logical, allocatable :: sepcies_computed(:) !track species have been computed or not
      real(fp), pointer :: tChemLoss(:) ! tChemLoss for each bin [days]
      real(fp), allocatable :: qUpdate(:), delq(:)  !intermediate variables for diagnostics
      integer :: nymd, nhms   !YYYYMMDD, HHMMSS time formats
      character(len=3) :: cdow  ! Character day of week
      integer :: idow  ! Integer day of week
      !myDOW seems to be -1 all the time in gocart. We keep the functionmality here, but the reset to zero never runs below.
      integer, parameter :: myDOW = -1  ! my Day of the week: Sun=1, Mon=2,...,Sat=7
      integer, parameter :: nbins = 2  ! Number of bins for phobic and philic species in GOCART
      logical reset_con !whether to reset concentration based on day of week
      !error information
      integer :: RC
      character(len=256) :: errMsg

      !allocate some arrays
      allocate(flux_tophilic(1,1), sepcies_computed(num_species), intptr_phobic_philic(1,1, num_layers,nbins), &
         tchemloss(nbins), fluxout(1,1,nbins), qupdate(num_layers), delq(num_layers))
      flux_tophilic = 0.0_fp ! Initialize to zero
      sepcies_computed = .false. ! Initialize to false
      reset_con = .false.
      intptr_phobic_philic = 0.0_fp ! Initialize to zero
      tchemloss = 0.0_fp ! Initialize to zero
      fluxout = 0.0_fp ! Initialize to zero
      qupdate = 0.0_fp ! Initialize to zero
      delq = 0.0_fp ! Initialize to zero

      !construct time in yyyymmdd and hhmmss formats for use in gocart
      nymd = year*10000 + month*100 + day
      nhms = hour*10000 + minute*100 + second
      !   Reset tracer to zero at 0Z on specific day of week
      !   --------------------------------------------------
      idow = chem_utilidow(nymd)
      if ( (nhms==0) .and. (idow == mydow) ) then
         reset_con = .true.
         cdow = chem_utilcdow(nymd)
         intptr_phobic_philic = tiny(1.0_fp) ! avoid division by zero
         write(*, '(A, I8, I8)') 'Note: Carbon '//cdow//' tracer being reset to zero on ', nymd, nhms
      end if

      ! transform data for GOCART DryDeposition call
      call prepmetvarsforgocart(num_layers, airden, delp, pmid, gocart_rhoa, gocart_press, gocart_delp)

      !get pressure lid index
      call findklid(klid, plid, gocart_press(:,:,:), rc)
      !if (RC /= CC_SUCCESS) then
      if (rc /= 0) then
         errmsg = 'Error in compute_gocart: Failed in finding pressure lid index in GOCART So4chem process.'
         !call CC_Error(trim(ErrMsg), RC, thisLoc)
         write(*,'(A)') trim(errmsg)
         return
      end if

      do species_idx = 1, num_species

         if (sepcies_computed(species_idx)) cycle ! Skip if already computed for this species

         ! Identify corresponding hydrophobic and hydrophilic species indices for this species
         select case (species_short_name(species_idx))
          case('oc1', 'OC1')  ! Example short names for OC species - used in gocart
            phobic_species_idx = species_idx
            ! Find corresponding hydrophilic species index
            philic_species_idx = max(find_species_ind(species_short_name, 'oc2'), find_species_ind(species_short_name, 'OC2'))
            if (philic_species_idx == 0) then
               write (*,'(A)') 'Error: No corresponding hydrophilic species found for OC1 species in GOCART scheme.'
               return
            end if
          case('bc1', 'BC1')  ! Example short names for BC species - used in gocart
            phobic_species_idx = species_idx
            ! Find corresponding hydrophilic species index
            philic_species_idx = max(find_species_ind(species_short_name, 'bc2'), find_species_ind(species_short_name, 'BC2'))
            if (philic_species_idx == 0) then
               write (*,'(A)') 'Error: No corresponding hydrophilic species found for BC1 species in GOCART scheme.'
               return
            end if
          case('br1', 'BR1')  ! Example short names for brown carbon species - used in gocart
            phobic_species_idx = species_idx
            ! Find corresponding hydrophilic species index
            philic_species_idx = max(find_species_ind(species_short_name, 'br2'), find_species_ind(species_short_name, 'Br2'))
            if (philic_species_idx == 0) then
               write (*,'(A)') 'Error: No corresponding hydrophilic species found for BR1 species in GOCART scheme.'
               return
            end if
          case('oc2', 'OC2')  ! Example short names for OC species - used in gocart
            philic_species_idx = species_idx
            ! Find corresponding hydrophobic species index
            phobic_species_idx = max(find_species_ind(species_short_name, 'oc1'), find_species_ind(species_short_name, 'OC1'))
            if (phobic_species_idx == 0) then
               write (*,'(A)') 'Error: No corresponding hydrophobic species found for OC2 species in GOCART scheme.'
               return
            end if
          case('bc2', 'BC2')  ! Example short names for BC species - used in gocart
            philic_species_idx = species_idx
            ! Find corresponding hydrophobic species index
            phobic_species_idx = max(find_species_ind(species_short_name, 'bc1'), find_species_ind(species_short_name, 'BC1'))
            if (phobic_species_idx == 0) then
               write (*,'(A)') 'Error: No corresponding hydrophobic species found for BC2 species in GOCART scheme.'
               return
            end if
          case('br2', 'BR2')  ! Example short names for brown carbon species - used in gocart
            philic_species_idx = species_idx
            ! Find corresponding hydrophobic species index
            phobic_species_idx = max(find_species_ind(species_short_name, 'br1'), find_species_ind(species_short_name, 'Br1'))
            if (phobic_species_idx == 0) then
               write (*,'(A)') 'Error: No corresponding hydrophobic species found for BR2 species in GOCART scheme.'
               return
            end if
          case default
            cycle  ! Skip species that are not BC or OC in this example
         end select

         !get concentration for this species after flipping vertical levels for gocart
         ! Unit conversion: model state [ug/kg] -> GOCART internal [kg/kg] (multiply by 1e-9)
         ! Vertical flip: model convention (surface=1) -> GOCART convention (top=1)
         if (.not. reset_con) then
            intptr_phobic_philic(1,1,:, 1) = species_conc(num_layers:1:-1, phobic_species_idx) * 1.0e-9_fp ! [ug/kg] -> [kg/kg]
            intptr_phobic_philic(1,1,:, 2) = species_conc(num_layers:1:-1, philic_species_idx) * 1.0e-9_fp ! [ug/kg] -> [kg/kg]
         end if

         !for diagnostics only; have to reproduce the calculation here to save out the mass in addition to the flux
         qupdate = intptr_phobic_philic(1,1,:,1)*exp(-tstep/(params%time_days_hydrophobic_to_hydrophilic*86400.0_fp))
         qupdate = max(qupdate,1.0e-32_fp)
         delq = max(0.0_fp, intptr_phobic_philic(1,1,:,1) - qupdate)

         !Ad Hoc transfer of hydrophobic to hydrophilic aerosols
         !Rate controlled in RC file; tConvPhobicToPhilic < 0 means no transfer
         call phobictophilic (intptr_phobic_philic(:,:,:,1), intptr_phobic_philic(:,:,:,2), flux_tophilic, &
            params%time_days_hydrophobic_to_hydrophilic, num_layers, tstep, g0, gocart_delp, rc)
         if (rc /= 0) then
            errmsg = 'Error in compute_gocart: Failed in GOCART phobicToPhilic.'
            !call CC_Error(trim(ErrMsg), RC, thisLoc)
            write(*,'(A)') trim(errmsg)
            return
         end if

         ! Per-species-per-level diagnostic: 2D array (levels, species)
         if (present(phobictophilic_mass_per_species_per_level) .and. present(diagnostic_species_id)) then
            ! Find position of this species in diagnostic_species_id array
            do diag_idx = 1, size(diagnostic_species_id)
               !note we give the same valuers to phobic and philic species in this case of conversion between the two
               if (diagnostic_species_id(diag_idx) == phobic_species_idx .or. diagnostic_species_id(diag_idx) == philic_species_idx) then
                  ! Add your custom conversion mass from hydrophobic to hydrophilic per species per level calculation
                  phobictophilic_mass_per_species_per_level(:, diag_idx) = delq(num_layers:1:-1) !flip the layers [kg/kg]
                  !exit !comment out to give the same value for both phobic and philic species in this case of conversion between the two
               end if
            end do
         end if
         ! Per-species diagnostic: only update for diagnostic species
         if (present(phobictophilic_flux_per_species) .and. present(diagnostic_species_id)) then
            ! Find position of this species in diagnostic_species_id array
            do diag_idx = 1, size(diagnostic_species_id)
               !note we give the same valuers to phobic and philic species in this case of conversion between the two
               if (diagnostic_species_id(diag_idx) == phobic_species_idx .or. diagnostic_species_id(diag_idx) == philic_species_idx) then
                  ! Add your custom conversion flux from hydrophobic to hydrophilic per species calculation
                  phobictophilic_flux_per_species(diag_idx) = flux_tophilic(1,1) !column total flux [kg/m2/s]
                  !exit !comment out to give the same value for both phobic and philic species in this case of conversion between the two
               end if
            end do
         end if

         !retrieve tChemLoss for this species
         tchemloss(1) = species_t_chem_loss(phobic_species_idx)
         tchemloss(2) = species_t_chem_loss(philic_species_idx)

         !Ad Hoc chemical destruction of carbon
         !This applies a simple exponential decay to both hydrophobic and
         ! hydrophilic modes with the time constant tChemLoss (e-folding time in days)
         do n = 1, nbins
            !for diagnostics only; have to reproduce the calculation here to save out the mass in addition to the flux
            qupdate = intptr_phobic_philic(1, 1, :, n)*exp(-tstep/(tchemloss(n)*86400.0_fp))
            qupdate = max(qupdate,1.e-32_fp)
            delq = max(0.0_fp,intptr_phobic_philic(1, 1, :, n)-qupdate)

            call carbonchemloss (num_layers, klid, n, tstep, g0, gocart_delp, &
               tchemloss(n), intptr_phobic_philic(:, :, :, n), fluxout, rc)
            if (rc /= 0) then
               errmsg = 'Error in compute_gocart: Failed in GOCART carbonChemLoss.'
               !call CC_Error(trim(ErrMsg), RC, thisLoc)
               write(*,'(A)') trim(errmsg)
               return
            end if

            !get current index for diagnostics
            if (n == 1) then
               curr_idx = phobic_species_idx
            else
               curr_idx = philic_species_idx
            end if

            ! Per-species-per-level diagnostic: 2D array (levels, species)
            if (present(production_mass_per_species_per_level) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == curr_idx) then
                     ! Add your custom production mass (loss here) per species per level calculation
                     production_mass_per_species_per_level(:, diag_idx) = delq(num_layers:1:-1) !note it is always loss here
                     exit
                  end if
               end do
            end if

            ! Per-species diagnostic: only update for diagnostic species
            if (present(loss_flux_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == curr_idx) then
                     ! Add your custom chemical loss flux per species calculation
                     loss_flux_per_species(diag_idx) = fluxout(1,1,n) !column total flux [kg/m2/s]
                     exit
                  end if
               end do
            end if

         end do

         ! Unit conversion: GOCART internal [kg/kg] -> model state [ug/kg] (multiply by 1e9)
         ! Vertical flip: GOCART convention (top=1) -> model convention (surface=1)
         ! Note: 1e-9 (input) and 1e9 (output) are symmetric, preserving mass consistency
         species_tendencies(:, phobic_species_idx) =  intptr_phobic_philic(1,1,num_layers:1:-1,1) * 1.0e9_fp  ! [kg/kg] -> [ug/kg]
         species_tendencies(:, philic_species_idx) =  intptr_phobic_philic(1,1,num_layers:1:-1,2) * 1.0e9_fp  ! [kg/kg] -> [ug/kg]

         !set computed species to true
         sepcies_computed(phobic_species_idx) = .true.
         sepcies_computed(philic_species_idx) = .true.

      end do ! End of loop over species

      if (associated(gocart_rhoa)) deallocate(gocart_rhoa); nullify(gocart_rhoa)
      if (associated(gocart_delp)) deallocate(gocart_delp); nullify(gocart_delp)
      if (associated(gocart_press)) deallocate(gocart_press); nullify(gocart_press)
      if (associated(flux_tophilic)) deallocate(flux_tophilic); nullify(flux_tophilic)
      if (associated(intptr_phobic_philic)) deallocate(intptr_phobic_philic); nullify(intptr_phobic_philic)
      if (associated(tchemloss)) deallocate(tchemloss); nullify(tchemloss)
      if (associated(fluxout)) deallocate(fluxout); nullify(fluxout)
      !other arrays
      deallocate(sepcies_computed, qupdate, delq)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   subroutine prepmetvarsforgocart(km,              &
      rhoa,            &
      delp,            &
      pmid,            &
      GOCART_RHOA,     &
      GOCART_PRESS,     &
      GOCART_DELP)



      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      REAL(fp),  intent(in), DIMENSION(:), target :: rhoa   ! Air density [kg/m^3]
      REAL(fp),  intent(in), DIMENSION(:), target :: delp    ! Pressure thickness [Pa]
      REAL(fp),  intent(in), DIMENSION(:), target :: pmid    ! Pressure at mid-point of layer [Pa]

      ! INPUT/OUTPUTS
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_PRESS
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_DELP

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(gocart_rhoa(1, 1, km))
      allocate(gocart_press(1, 1, km))
      allocate(gocart_delp(1, 1, km))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)
      gocart_rhoa(1,1,:) = rhoa(size(rhoa):1:-1) ! air density [kg/m^3]
      gocart_delp(1,1,:) = delp(size(delp):1:-1) ! pressure thickness [Pa]
      gocart_press(1,1,:) = pmid(size(pmid):1:-1) ! air pressure [Pa]

   end subroutine prepmetvarsforgocart

   subroutine findklid (klid, plid, ple, rc)

      implicit NONE
      ! !INPUT PARAMETERS:
      integer, intent(inout) :: klid ! index for pressure lid
      real(fp), intent(in)       :: plid ! pressure lid [hPa]; default is 0.01 hPa
      real(fp), dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]
      ! !OUTPUT PARAMETERS:
      integer, intent(out) :: rc ! return code; 0 - all is good; 1 - bad
      ! !Reference to gocart: https://github.com/GEOS-ESM/GOCART/blob/9ff3df9545dd582f415f682d3297e8c6c841e5cb/ESMF/Shared/Chem_AeroGeneric.F90#L316
      ! !Local Variables
      integer :: k, j, i
      real(fp) :: plid_, diff, refDiff
      real(fp), allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]
      !EOP
      !----------------------------------------------------------------------------------
      !  Begin...
      klid = 1
      rc = 0

      !  convert from hPa to Pa
      plid_ = plid*100.0_fp

      allocate(pres(ubound(ple,3)))

      !  find pressure at each model level
      do k = 1, ubound(ple,3)
         pres(k) = ple(1,1,k)
      end do

      !  find smallest absolute difference between plid and average pressure at each model level
      refdiff = 150000.0_fp
      do k = 1, ubound(ple,3)
         diff = abs(pres(k) - plid_)
         if (diff < refdiff) then
            klid = k
            refdiff = diff
         end if
      end do

      !  Check to make sure that all pressures at (i,j) were the same
      do j = 1, ubound(ple,2)
         do i = 1, ubound(ple,1)
            !if (pres(klid) /= ple(i,j,klid)) then !This gives a warning for floating point comparison. Use rae instead
            if (.not. rae(pres(klid), ple(i,j,klid))) then
               rc = 1
               return
            end if
         end do
      end do

   end subroutine findklid

   function find_species_ind(SpeciesNames, species_name) result(species_index)
      implicit none
      character(len=*), intent(in) :: SpeciesNames(:)
      character(len=*), intent(in) :: species_name
      integer :: species_index

      integer :: i

      species_index = 0  ! Not found

      if (size(speciesnames) > 0 ) then
         do i = 1, size(speciesnames)
            if (trim(speciesnames(i)) == trim(species_name)) then
               species_index = i
               exit
            endif
         enddo
      endif

   end function find_species_ind

end module carbchemscheme_gocart_mod
```


