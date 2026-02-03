

# File SettlingScheme\_GOCART\_Mod.F90

[**File List**](files.md) **>** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md) **>** [**settling**](dir_1a0bba2ffdf6e6637fcb76856471cb75.md) **>** [**schemes**](dir_34df91cc26d24067840a7381fe21b817.md) **>** [**SettlingScheme\_GOCART\_Mod.F90**](_settling_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the documentation of this file](_settling_scheme___g_o_c_a_r_t___mod_8_f90.md)


```Fortran

module settlingscheme_gocart_mod

   use precision_mod, only: fp, rae
   use settlingcommon_mod, only: settlingschemegocartconfig
   use error_mod, only: cc_success, cc_error
   use constants, only: g0  !load the constants needed for this scheme
   use gocart2g_miemod, only: gocart2g_mie  ! For Mie data in gocart scheme

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
      airden, &
      delp, &
      pmid, &
      rh, &
      t, &
      tstep, &
      z, &
      species_short_name, &
      mie_data, &
      species_mie_map, &
      species_radius, &
      species_density, &
      species_conc, &
      species_tendencies, &
      settling_velocity_per_species_per_level, &
      settling_flux_per_species, &
      diagnostic_species_id &
      )
      ! Uses
      USE gocart2g_process, only: chem_settlingsimple, chem_settling
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
      real(fp), intent(in) :: z(num_layers+1)    ! 3D atmospheric field
      character(len=32), intent(in) :: species_short_name(:)  ! Species short_name property
      type(GOCART2G_Mie), intent(in) :: mie_data(:)  ! Complete Mie data array from ChemState
      integer, intent(in) :: species_mie_map(num_species)  ! Mapping from process species to MieData indices
      real(fp), intent(in) :: species_radius(:)  ! Species radius property
      real(fp), intent(in) :: species_density(:)  ! Species density property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      real(fp), intent(inout), optional :: settling_velocity_per_species_per_level(:,:)
      real(fp), intent(inout), optional :: settling_flux_per_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: rc, species_idx, p
      integer :: diag_idx  ! For diagnostic species indexing
      integer :: bin  ! For bin index
      integer :: klid  ! For pressure lid index
      ! Local Variables
      real(fp), pointer :: GOCART_tmpu(:,:,:)
      real(fp), pointer :: GOCART_rhoa(:,:,:)
      real(fp), pointer :: GOCART_HGHTE(:,:,:)
      real(fp), pointer :: GOCART_RH(:,:,:)
      real(fp), pointer :: GOCART_PRESS(:,:,:)
      real(fp), pointer :: GOCART_DELP(:,:,:)
      real(fp) :: qa(1,1,num_layers)  ! concentration in [kg/kg]
      real(fp), pointer :: SD(:,:,:)  ! settling velocity [m/s]
      real(fp), pointer :: fluxout(:,:,:)  ! flux out across column [kg/m2/s]
      real(fp), pointer :: fluxout_temp(:,:)  ! flux out across column [kg/m2/s]
      !error information
      character(len=255)       :: thisLoc
      character(len=512)       :: ErrMsg
      errmsg  = ''
      thisloc = ' -> at compute_gocart (in process/settling/schemes/SettlingScheme_GOCART_Mod.F90)'
      ! Initialize
      rc = cc_success
      klid = 1 !since the layer is reversed, we give 1 here, which is the top layer
      qa = 0.0_fp
      allocate(sd(1, 1, num_layers))
      allocate(fluxout_temp(1,1))
      sd = 0.0_fp
      fluxout_temp = 0.0_fp

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.

      !reverse met variables
      call prepmetvarsforgocart(num_layers,              &
         t,               &
         airden,          &
         z,               &
         rh,              &
         pmid,            &
         delp,            &
         gocart_tmpu,     &
         gocart_rhoa,     &
         gocart_hghte,    &
         gocart_rh,       &
         gocart_press,    &
         gocart_delp)

      !get pressure lid index
      call findklid(klid, plid, gocart_press(:,:,:), rc)
      if (rc /= cc_success) then
         errmsg = 'Error in finding pressure lid index in GOCART settling scheme.'
         call cc_error(trim(errmsg), rc, thisloc)
         return
      end if

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      ! Apply to each species
      do species_idx = 1, num_species
         ! get bin based on name ending with 1, 2,3,4,5 or a letter
         if (len_trim(species_short_name(species_idx)) >= 1) then
            p = len_trim(species_short_name(species_idx))
            select case (species_short_name(species_idx)(p:p))
             case ('1')
               bin = 1
             case ('2')
               bin = 2
             case ('3')
               bin = 3
             case ('4')
               bin = 4
             case ('5')
               bin = 5
             case default
               bin = 1  ! default bin if no match
            end select
         else
            bin = 1  ! default bin if name is too short
         end if

         !initialize fluxout
         allocate(fluxout(1, 1, bin))
         fluxout = 0.0_fp

         !reverse vertical layer and convert from ug/kg to kg/kg
         qa(1,1,:) = species_conc(num_layers:1:-1, species_idx) * 1.0e-9_fp  ! from ug/kg to kg/kg

         if (params%simple_scheme) then !call gocart simple settling function with mie data provided
            !check mie data is available
            if (species_mie_map(species_idx) <= 0) then
               errmsg = 'Invalid Mie data mapping found. Check if proper mie files are provided.'
               rc = 1
               call cc_error(trim(errmsg), rc, thisloc)
               return
            end if
            call chem_settlingsimple (num_layers, klid, mie_data(species_mie_map(species_idx)), bin, tstep, g0, &
               qa, gocart_tmpu, gocart_rhoa, gocart_rh, gocart_hghte, gocart_delp, fluxout_temp, &
               vsettleout=sd, correctionmaring=params%correction_maring, settling_scheme=2, rc=rc) !hardcode settling_scheme=2 for GOCART scheme
            if (rc /= cc_success) then
               errmsg = 'Error in running GOCART Chem_SettlingSimple scheme.'
               call cc_error(trim(errmsg), rc, thisloc)
               return
            end if
            fluxout(:,:,bin) = fluxout_temp(:,:)
         else  !call gocart simple settling function with internal mie calculation
            call chem_settling (num_layers, klid, bin, params%swelling_method, tstep, g0, species_radius(species_idx)*1.0e-6_fp,  &  !um to m
               species_density(species_idx), qa, gocart_tmpu, gocart_rhoa, gocart_rh, gocart_hghte, gocart_delp, fluxout, &
               vsettleout=sd, correctionmaring=params%correction_maring, settling_scheme=2, rc=rc) !hardcode settling_scheme=2 for GOCART scheme
            if (rc /= cc_success) then
               errmsg = 'Error in running GOCART Chem_Settling scheme.'
               call cc_error(trim(errmsg), rc, thisloc)
               return
            end if
         end if

         ! convert concentrations back to original order and units
         species_tendencies(:, species_idx) = max(0.0_fp, qa(1, 1, num_layers:1:-1) * 1.0e9_fp)   ! from kg/kg to ug/kg

         ! TODO: Update diagnostic fields here based on your scheme's requirements
         ! Each process should implement custom diagnostic calculations
         ! Example patterns:
         ! Per-species-per-level diagnostic: 2D array (levels, species)
         if (present(settling_velocity_per_species_per_level) .and. present(diagnostic_species_id)) then
            ! Find position of this species in diagnostic_species_id array
            do diag_idx = 1, size(diagnostic_species_id)
               if (diagnostic_species_id(diag_idx) == species_idx) then
                  ! Add your custom settling velocity per species per level calculation
                  settling_velocity_per_species_per_level(:, diag_idx) = sd(1,1,num_layers:1:-1)  !reverse layers
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
                  settling_flux_per_species(diag_idx) = fluxout(1,1, bin)
                  exit
               end if
            end do
         end if
      end do

      !cleanup pointers
      if (associated(gocart_tmpu)) nullify(gocart_tmpu)
      if (associated(gocart_rhoa)) nullify(gocart_rhoa)
      if (associated(gocart_hghte)) nullify(gocart_hghte)
      if (associated(gocart_rh)) nullify(gocart_rh)
      if (associated(gocart_press)) nullify(gocart_press)
      if (associated(gocart_delp)) nullify(gocart_delp)
      if (associated(sd)) nullify(sd)
      if (associated(fluxout)) nullify(fluxout)
      if (associated(fluxout_temp)) nullify(fluxout_temp)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

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

   subroutine prepmetvarsforgocart(km,              &
      tmpu,            &
      rhoa,            &
      hghte,           &
      rh,              &
      press,           &
      delp,            &
      GOCART_tmpu,     &
      GOCART_RHOA,     &
      GOCART_HGHTE,    &
      GOCART_RH,       &
      GOCART_PRESS,    &
      GOCART_DELP)

      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      REAL(fp),  intent(in), DIMENSION(:), target :: tmpu   ! Temperature [K]
      REAL(fp),  intent(in), DIMENSION(:), target :: rhoa   ! Air density [kg/m^3]
      REAL(fp),  intent(in), DIMENSION(:), target :: hghte  ! Geopotential Height [m]
      REAL(fp),  intent(in), DIMENSION(:), target :: rh     ! Relative Humidity [decimal; not %]
      REAL(fp),  intent(in), DIMENSION(:), target :: press  ! Pressure [Pa]
      REAL(fp),  intent(in), DIMENSION(:), target :: delp  ! Pressure thickness [Pa]

      ! INPUT/OUTPUTS
      REAL(fp), intent(inout), pointer :: GOCART_TMPU(:,:,:)
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RH
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_PRESS
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_DELP


      allocate(gocart_tmpu(1, 1, km))
      allocate(gocart_rhoa(1, 1, km))
      allocate(gocart_hghte(1,1, 0:km))
      allocate(gocart_rh(1,1, km))
      allocate(gocart_press(1,1, km))
      allocate(gocart_delp(1,1, km))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)
      gocart_tmpu(1,1,:) = tmpu(size(tmpu):1:-1)       ! temperature [K]
      gocart_rhoa(1,1,:) = rhoa(size(rhoa):1:-1)       ! air density [kg/m^3]
      gocart_hghte(1,1,:) = hghte(size(hghte):1:-1)    ! geopotential height [m]
      gocart_rh(1,1,:) = rh(size(rh):1:-1)             ! relative humidity [decimal; not %]
      gocart_press(1,1,:) = press(size(press):1:-1)    ! pressure [Pa]
      gocart_delp(1,1,:) = delp(size(delp):1:-1)       ! pressure thickness [Pa]

   end subroutine prepmetvarsforgocart

end module settlingscheme_gocart_mod
```


