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

   use catchem_bridge_precision, only: fp, rae
   use SettlingCommon_Mod, only: SettlingSchemeGOCARTConfig
   use catchem_bridge_error, only: CC_SUCCESS, CC_Error
   use catchem_bridge_constants, only: g0  !load the constants needed for this scheme
   use GOCART2G_MieMod, only: GOCART2G_Mie  ! For Mie data in gocart scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: plid = 0.01_fp    ! Pressure lid [hPa]

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
   !! @param[in]  z     Z field [appropriate units]
   !! @param[in]  species_short_name    Species short_name property
   !! @param[in]  mie_data           Complete Mie data array from ChemState
   !! @param[in]  species_mie_map    Mapping from process species to MieData indices
   !! @param[in]  species_radius    Species radius property
   !! @param[in]  species_density    Species density property
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] settling_velocity_per_species_per_level    settling velocity per species per level [m/s] (num_layers, num_species)
   !! @param[inout] settling_flux_per_species    settling flux per species across column [kg/m2/s] (num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)
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
      species_is_dust, &
      species_conc, &
      species_tendencies, &
      settling_velocity_per_species_per_level, &
      settling_flux_per_species, &
      diagnostic_species_id &
      )
      ! Uses
      USE GOCART2G_Process, only: Chem_SettlingSimple, Chem_Settling
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
      logical, intent(in) :: species_is_dust(:)  ! Species is_dust property
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
      logical :: apply_maring  ! per-species Maring correction switch
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
      ErrMsg  = ''
      ThisLoc = ' -> at compute_gocart (in process/settling/schemes/SettlingScheme_GOCART_Mod.F90)'
      ! Initialize
      RC = CC_SUCCESS
      klid = 1 !since the layer is reversed, we give 1 here, which is the top layer
      qa = 0.0_fp
      allocate(SD(1, 1, num_layers))
      allocate(fluxout_temp(1,1))
      SD = 0.0_fp
      fluxout_temp = 0.0_fp
      ! fluxout is (re)allocated each species iteration below; give it a defined
      ! disassociated status here so the first associated()/deallocate() guard is
      ! valid. A pointer without => null()/nullify has UNDEFINED association status
      ! on entry, so associated() is illegal and compiler-dependent (ifort happened
      ! to return .false., gfortran .true. -> deallocate of a bogus pointer -> crash).
      nullify(fluxout)

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.

      !reverse met variables
      call PrepMetVarsForGOCART(num_layers,              &
         t,               &
         airden,          &
         z,               &
         rh,              &
         pmid,            &
         delp,            &
         GOCART_tmpu,     &
         GOCART_RHOA,     &
         GOCART_HGHTE,    &
         GOCART_RH,       &
         GOCART_PRESS,    &
         GOCART_DELP)

      !get pressure lid index
      call findKlid(klid, plid, GOCART_PRESS(:,:,:), RC)
      if (RC /= CC_SUCCESS) then
         ErrMsg = 'Error in finding pressure lid index in GOCART settling scheme.'
         call CC_Error(trim(ErrMsg), RC, thisLoc)
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

         !Check if it's a multi-bin species (like dust4) vs chemical name (like SO4, which has two mie bins with MSA)
         if (size(mie_data) > 0 .and. species_mie_map(species_idx) >= 1 .and. &
            species_mie_map(species_idx) <= size(mie_data)) then
            if (bin > mie_data(species_mie_map(species_idx))%nbin) then
               bin = 1
            end if
         end if

         !initialize fluxout
         !NOTE: fluxout is a pointer allocated inside this per-species loop; free any
         !prior iteration's allocation first to avoid orphaning it (memory leak).
         if (associated(fluxout)) deallocate(fluxout)
         allocate(fluxout(1, 1, bin))
         fluxout = 0.0_fp

         !reverse vertical layer and convert from ug/kg to kg/kg
         qa(1,1,:) = species_conc(num_layers:1:-1, species_idx) * 1.0e-9_fp  ! from ug/kg to kg/kg

         ! Maring (2003) correction was derived for dust; GOCART applies it to dust only, not sea salt
         ! (mirrors drydep dust_resuspension_only logic)
         apply_maring = (params%maring_dust_only .and. species_is_dust(species_idx)) .or. &
            (.not. params%maring_dust_only .and. params%correction_maring)

         if (params%simple_scheme) then !call gocart simple settling function with mie data provided
            !check mie data is available
            if (species_mie_map(species_idx) <= 0) then
               ErrMsg = 'Invalid Mie data mapping found. Check if proper mie files are provided.'
               RC = 1
               call CC_Error(trim(ErrMsg), RC, thisLoc)
               return
            end if
            call Chem_SettlingSimple (num_layers, klid, mie_data(species_mie_map(species_idx)), bin, tstep, g0, &
               qa, GOCART_tmpu, GOCART_rhoa, GOCART_RH, GOCART_HGHTE, GOCART_DELP, fluxout_temp, &
               vsettleOut=SD, correctionMaring=apply_maring, settling_scheme=2, rc=RC) !hardcode settling_scheme=2 for GOCART scheme
            if (RC /= CC_SUCCESS) then
               ErrMsg = 'Error in running GOCART Chem_SettlingSimple scheme.'
               call CC_Error(trim(ErrMsg), RC, thisLoc)
               return
            end if
            fluxout(:,:,bin) = fluxout_temp(:,:)
         else  !call gocart simple settling function with internal mie calculation
            call Chem_Settling (num_layers, klid, bin, params%swelling_method, tstep, g0, species_radius(species_idx)*1.0e-6_fp,  &  !um to m
               species_density(species_idx), qa, GOCART_tmpu, GOCART_RHOA, GOCART_RH, GOCART_HGHTE, GOCART_DELP, fluxout, &
               vsettleOut=SD, correctionMaring=apply_maring, settling_scheme=2, rc=RC) !hardcode settling_scheme=2 for GOCART scheme
            if (RC /= CC_SUCCESS) then
               ErrMsg = 'Error in running GOCART Chem_Settling scheme.'
               call CC_Error(trim(ErrMsg), RC, thisLoc)
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
                  settling_velocity_per_species_per_level(:, diag_idx) = SD(1,1,num_layers:1:-1)  !reverse layers
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
      if (associated(GOCART_TMPU)) deallocate(GOCART_TMPU); nullify(GOCART_TMPU)
      if (associated(GOCART_RHOA)) deallocate(GOCART_RHOA); nullify(GOCART_RHOA)
      if (associated(GOCART_HGHTE)) deallocate(GOCART_HGHTE); nullify(GOCART_HGHTE)
      if (associated(GOCART_RH)) deallocate(GOCART_RH); nullify(GOCART_RH)
      if (associated(GOCART_PRESS)) deallocate(GOCART_PRESS); nullify(GOCART_PRESS)
      if (associated(GOCART_DELP)) deallocate(GOCART_DELP); nullify(GOCART_DELP)
      if (associated(SD)) deallocate(SD); nullify(SD)
      if (associated(fluxout)) deallocate(fluxout); nullify(fluxout)
      if (associated(fluxout_temp)) deallocate(fluxout_temp); nullify(fluxout_temp)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !>
   !! \brief findKlid - Finds corresponding vertical index for defined pressure lid
   !!
   !! \param [INOUT] klid
   !! \param [IN] plid
   !! \param [IN] ple
   !! \param [OUT] rc
   !!!>
   subroutine findKlid (klid, plid, ple, rc)

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
      refDiff = 150000.0_fp
      do k = 1, ubound(ple,3)
         diff = abs(pres(k) - plid_)
         if (diff < refDiff) then
            klid = k
            refDiff = diff
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

   end subroutine findKlid

   !>
   !! \brief PrepMetVarsForGOCART - Prep the meteorological variables for GOCART settling scheme
   !!
   !! \param [IN]    km
   !! \param [IN] tmpu
   !! \param [IN] rhoa
   !! \param [IN] hghte
   !! \param [IN] rh
   !! \param [IN] press
   !! \param [IN] delp
   !! \param [INOUT] GOCART_tmpu
   !! \param [INOUT] GOCART_RHOA
   !! \param [INOUT] GOCART_HGHTE
   !! \param [INOUT] GOCART_RH
   !! \param [INOUT] GOCART_PRESS
   !! \param [INOUT] GOCART_DELP
   !!!>
   subroutine PrepMetVarsForGOCART(km,              &
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
      REAL(fp), intent(inout), pointer :: GOCART_TMPU(:,:,:)   !< temperature [K]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA   !< air density [kg/m^3]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE  !< geometric height [m]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RH     !< relative humidity [decimal; not %]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_PRESS  !< pressure [Pa]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_DELP   !< pressure thickness [Pa]


      allocate(GOCART_TMPU(1, 1, km))
      allocate(GOCART_RHOA(1, 1, km))
      allocate(GOCART_HGHTE(1,1, 0:km))
      allocate(GOCART_RH(1,1, km))
      allocate(GOCART_PRESS(1,1, km))
      allocate(GOCART_DELP(1,1, km))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)
      GOCART_TMPU(1,1,:) = tmpu(size(tmpu):1:-1)       ! temperature [K]
      GOCART_RHOA(1,1,:) = rhoa(size(rhoa):1:-1)       ! air density [kg/m^3]
      GOCART_HGHTE(1,1,:) = hghte(size(hghte):1:-1)    ! geopotential height [m]
      GOCART_RH(1,1,:) = rh(size(rh):1:-1)             ! relative humidity [decimal; not %]
      GOCART_PRESS(1,1,:) = press(size(press):1:-1)    ! pressure [Pa]
      GOCART_DELP(1,1,:) = delp(size(delp):1:-1)       ! pressure thickness [Pa]

   end subroutine PrepMetVarsForGOCART

end module SettlingScheme_GOCART_Mod
