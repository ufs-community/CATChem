!> \file DryDepScheme_GOCART_Mod.F90
!! \brief GOCART-2G aerosol dry deposition scheme
!!
!! Pure science kernel for gocart scheme in drydep process.
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
!! Generated on: 2025-11-14T22:58:26.525574
!! Author: Wei Li & Lacey Holland
!! Reference: Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS)
!! Allison B. Collow, Peter R. Colarco, Arlindo M. da Silva, Virginie Buchard,
!! Huisheng Bian, M Chin, Sampa Das, Ravi Govindaraju, Dongchul Kim, and Valentina Aquila,
!! Geosci. Model Development, 17, 14431468, 2024
!! https://doi.org/10.5194/gmd-17-1443-2024
module DryDepScheme_GOCART_Mod

   use precision_mod, only: fp
   use DryDepCommon_Mod, only: DryDepSchemeGOCARTConfig
   use error_mod, only: CC_SUCCESS, CC_Error
   use Constants, only: Cp, g0, VON_KARMAN  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

contains

   !> Pure science computation for gocart scheme
   !!
   !! This is a pure computational kernel implementing GOCART-2G aerosol dry deposition scheme.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  airden    AIRDEN field [appropriate units]
   !! @param[in]  frlake    FRLAKE field [appropriate units]
   !! @param[in]  gwettop    GWETTOP field [appropriate units]
   !! @param[in]  hflux    HFLUX field [appropriate units]
   !! @param[in]  lwi    LWI field [appropriate units]
   !! @param[in]  pblh    PBLH field [appropriate units]
   !! @param[in]  t    T field [appropriate units]
   !! @param[in]  tstep    Time step [s] - retrieved from process interface
   !! @param[in]  u10m    U10M field [appropriate units]
   !! @param[in]  ustar    USTAR field [appropriate units]
   !! @param[in]  v10m    V10M field [appropriate units]
   !! @param[in]  z0h    Z0H field [appropriate units]
   !! @param[in]  z    Z field [appropriate units]
   !! @param[in]  species_density    Species density property
   !! @param[in]  species_radius    Species radius property
   !! @param[in]  species_is_seasalt    Species is_seasalt property
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[inout] species_tendencies  Species tendency terms [mol/mol/s] (num_layers, num_species)
   !! @param[inout] drydep_con_per_species    Dry deposition concentration per species [ug/kg or ppm] (num_species)
   !! @param[inout] drydep_velocity_per_species    Dry deposition velocity [m/s] (num_species)
   !! @param[in] diagnostic_species_id Indices mapping diagnostic species to species array (optional, for per-species diagnostics)

   subroutine compute_gocart( &
      num_layers, &
      num_species, &
      params, &
      airden, &
      frlake, &
      gwettop, &
      hflux, &
      lwi, &
      pblh, &
      t, &
      tstep, &
      u10m, &
      ustar, &
      v10m, &
      z, &
      z0h, &
      species_density, &
      species_radius, &
      species_is_seasalt, &
      species_conc, &
      species_tendencies, &
      is_gas, &
      drydep_con_per_species, &
      drydep_velocity_per_species, &
      diagnostic_species_id &
      )

      ! Uses
      USE GOCART2G_Process, only: DryDeposition
      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(DryDepSchemeGOCARTConfig), intent(in) :: params
      real(fp), intent(in) :: airden(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: frlake  ! Surface field - scalar
      real(fp), intent(in) :: gwettop  ! Surface field - scalar
      real(fp), intent(in) :: hflux  ! Surface field - scalar
      integer, intent(in) :: lwi  ! Surface field - scalar
      real(fp), intent(in) :: pblh  ! Surface field - scalar
      real(fp), intent(in) :: t(num_layers)    ! 3D atmospheric field
      real(fp), intent(in) :: tstep  ! Time step [s] - from process interface
      real(fp), intent(in) :: u10m  ! Surface field - scalar
      real(fp), intent(in) :: ustar  ! Surface field - scalar
      real(fp), intent(in) :: v10m  ! Surface field - scalar
      real(fp), intent(in) :: z0h  ! Surface field - scalar
      real(fp), intent(in) :: z(num_layers+1)    ! 3D atmospheric field
      real(fp), intent(in) :: species_density(num_species)  ! Species density property
      real(fp), intent(in) :: species_radius(num_species)  ! Species radius property
      logical, intent(in) :: species_is_seasalt(num_species)  ! Species is seasalt property
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(inout) :: species_tendencies(num_layers, num_species)
      logical, intent(in) :: is_gas(num_species)  ! Species type flags (true=gas, false=aerosol)
      real(fp), intent(inout), optional :: drydep_con_per_species(:)
      real(fp), intent(inout), optional :: drydep_velocity_per_species(:)
      integer, intent(in), optional :: diagnostic_species_id(:)  ! Indices mapping diagnostic species to species array

      ! Local variables
      integer :: rc, k, species_idx
      integer :: diag_idx  ! For diagnostic species indexing
      real(fp) :: VD
      real(fp) :: drydepf(1,1)
      ! Local Variables
      real(fp), pointer :: GOCART_tmpu(:,:,:)
      real(fp), pointer :: GOCART_rhoa(:,:,:)
      real(fp), pointer :: GOCART_HGHTE(:,:,:)
      real(fp), pointer :: GOCART_LWI(:,:)
      real(fp), pointer :: GOCART_USTAR(:,:)
      real(fp), pointer :: GOCART_PBLH(:,:)

      real(fp), pointer :: GOCART_HFLUX(:,:)
      real(fp), pointer :: GOCART_Z0H(:,:)
      real(fp), pointer :: GOCART_U10(:,:)
      real(fp), pointer :: GOCART_V10(:,:)
      real(fp), pointer :: GOCART_FRACLAKE(:,:)
      real(fp), pointer :: GOCART_GWETTOP(:,:)

      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at compute_gocart (in DryDepScheme_GOCART_Mod.F90)'
      RC = CC_SUCCESS
      VD = 0.0_fp
      drydepf = 0.0_fp

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.

      ! transform data for GOCART DryDeposition call
      call PrepMetVarsForGOCART(num_layers,     &
         t,            &
         airden,            &
         z,            &
         u10m,             &
         v10m,             &
         frlake,        &
         gwettop,         &
         lwi,             &
         ustar,           &
         pblh,            &
         hflux,           &
         z0h,             &
         GOCART_tmpu,     &
         GOCART_RHOA,     &
         GOCART_HGHTE,    &
         GOCART_U10,      &
         GOCART_V10,      &
         GOCART_FRACLAKE, &
         GOCART_GWETTOP,  &
         GOCART_LWI,      &
         GOCART_USTAR,    &
         GOCART_PBLH,     &
         GOCART_HFLUX,    &
         GOCART_Z0H)

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species
            ! Skip species that don't match scheme type (gas vs aerosol)
            if (is_gas(species_idx)) cycle

            if (params%resuspension) then
               call DryDeposition(num_layers, GOCART_TMPU, GOCART_RHOA, GOCART_HGHTE, GOCART_LWI, GOCART_USTAR, &
                  GOCART_PBLH, GOCART_HFLUX, von_karman, cp, g0, GOCART_Z0H, drydepf, RC, &
                  species_radius(species_idx)*1e-6_fp, species_density(species_idx), GOCART_U10, GOCART_V10, &
                  GOCART_FRACLAKE, GOCART_GWETTOP)
            else
               call DryDeposition(num_layers, GOCART_TMPU, GOCART_RHOA, GOCART_HGHTE, GOCART_LWI, GOCART_USTAR, &
                  GOCART_PBLH, GOCART_HFLUX, von_karman, cp, g0, GOCART_Z0H, drydepf, RC)
            endif

            ! Ensure non-negative values
            species_tendencies(k, species_idx) = max(0.0_fp, drydepf(1,1)) * params%scale_factor
            !increase drydep frequency by factor of 5 for seasalt species to match GOCART2G. See the codes in:
            !https://github.com/GEOS-ESM/GOCART/blob/9ff3df9545dd582f415f682d3297e8c6c841e5cb/ESMF/GOCART2G_GridComp/SS2G_GridComp/SS2G_GridCompMod.F90#L820
            if (species_is_seasalt(species_idx) .and. abs(LWI - LAND) < 0.5) then
               species_tendencies(k, species_idx) = species_tendencies(k, species_idx) * 5.0_fp
            end if

            VD = max(species_tendencies(k, species_idx) * (z(k+1) -z(k) ), 1.e-4_fp)


            ! TODO: Update diagnostic fields here based on your scheme's requirements
            ! Each process should implement custom diagnostic calculations
            ! Example patterns:
            ! Per-species diagnostic: only update for diagnostic species
            if (present(drydep_con_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom dry deposition concentration per species calculation
                     drydep_con_per_species(diag_idx) =  &
                        MAX(0.0_fp, species_conc(k,species_idx) * (1.0_fp - exp(-1.0_fp * species_tendencies(k, species_idx) * tstep)))
                     exit
                  end if
               end do
            end if
            ! Per-species diagnostic: only update for diagnostic species
            if (present(drydep_velocity_per_species) .and. present(diagnostic_species_id)) then
               ! Find position of this species in diagnostic_species_id array
               do diag_idx = 1, size(diagnostic_species_id)
                  if (diagnostic_species_id(diag_idx) == species_idx) then
                     ! Add your custom dry deposition velocity calculation
                     drydep_velocity_per_species(diag_idx) = VD
                     exit
                  end if
               end do
            end if
         end do

      end do

      !cleanup pointers
      if (associated(GOCART_TMPU)) nullify(GOCART_TMPU)
      if (associated(GOCART_RHOA)) nullify(GOCART_RHOA)
      if (associated(GOCART_HGHTE)) nullify(GOCART_HGHTE)
      if (associated(GOCART_U10)) nullify(GOCART_U10)
      if (associated(GOCART_V10)) nullify(GOCART_V10)
      if (associated(GOCART_FRACLAKE)) nullify(GOCART_FRACLAKE)
      if (associated(GOCART_GWETTOP)) nullify(GOCART_GWETTOP)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)
      if (associated(GOCART_USTAR)) nullify(GOCART_USTAR)
      if (associated(GOCART_LWI)) nullify(GOCART_LWI)
      if (associated(GOCART_HFLUX)) nullify(GOCART_HFLUX)
      if (associated(GOCART_Z0H)) nullify(GOCART_Z0H)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines

   !>
   !! \brief PrepMetVarsForGOCART - Prep the meteorological variables for GOCART DryDeposition scheme
   !!
   !! \param [INOUT] metstate
   !! \param [INOUT] tmpu
   !! \param [INOUT] rhoa
   !! \param [INOUT] hghte
   !! \param [INOUT] oro
   !! \param [INOUT] ustar
   !! \param [INOUT] pblh
   !! \param [INOUT] shflux
   !! \param [INOUT] z0h
   !! \param [INOUT] u10m
   !! \param [INOUT] v10m
   !! \param [INOUT] fraclake
   !! \param [INOUT] gwettop
   !! \param [OUT] rc
   !!
   !! \ingroup core_modules
   !!!>
   subroutine PrepMetVarsForGOCART(km,              &
      tmpu,            &
      rhoa,            &
      hghte,           &
      u10m,             &
      v10m,             &
      fraclake,        &
      gwettop,         &
      lwi,             &
      ustar,           &
      pblh,            &
      hflux,           &
      z0h,             &
      GOCART_tmpu,     &
      GOCART_RHOA,     &
      GOCART_HGHTE,    &
      GOCART_U10,      &
      GOCART_V10,      &
      GOCART_FRACLAKE, &
      GOCART_GWETTOP,  &
      GOCART_LWI,      &
      GOCART_USTAR,    &
      GOCART_PBLH,     &
      GOCART_HFLUX,    &
      GOCART_Z0H)



      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      INTEGER,  intent(in)                    :: lwi                                    ! orography flag; Land, ocean, ice mask
      REAL(fp),  intent(in), DIMENSION(:), target :: tmpu   ! Temperature [K]
      REAL(fp),  intent(in), DIMENSION(:), target :: rhoa   ! Air density [kg/m^3]
      REAL(fp),  intent(in), DIMENSION(:), target :: hghte  ! Height [m]
      REAL(fp),  intent(in), target               :: ustar                                 ! friction speed [m/sec]
      REAL(fp),  intent(in), target              :: pblh                                  ! PBL height [m]
      REAL(fp),  intent(in), target              :: hflux                                 ! sfc. sens. heat flux [W m-2]
      REAL(fp),  intent(in), target              :: z0h                                   ! rough height, sens. heat [m]
      REAL(fp),  intent(in), target :: u10m                   ! 10-m u-wind component [m/sec]
      REAL(fp),  intent(in), target :: v10m                   ! 10-m v-wind component [m/sec]
      REAL(fp),  intent(in), target :: fraclake               ! fraction covered by water [1]
      REAL(fp),  intent(in), target :: gwettop                ! fraction soil moisture [1]

      ! INPUT/OUTPUTS
      REAL(fp), intent(inout), pointer :: GOCART_TMPU(:,:,:)   !< temperature [K]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA   !< air density [kg/m^3]
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE  !< geometric height [m]
      REAL(fp), intent(inout), pointer :: GOCART_U10(:,:)                 !< 10-m u-wind component [m/sec]
      REAL(fp), intent(inout), pointer :: GOCART_V10 (:,:)                !< 10-m v-wind component [m/sec]
      REAL(fp), intent(inout), pointer :: GOCART_FRACLAKE(:,:)            !< fraction covered by water [1]
      REAL(fp), intent(inout), pointer :: GOCART_GWETTOP(:,:)             !< fraction soil moisture [1]
      real(fp), intent(inout), pointer :: GOCART_LWI(:,:)                 !< orography flag; Land, ocean, ice mask
      REAL(fp), intent(inout), pointer :: GOCART_USTAR(:,:)               !< friction speed [m/sec]
      REAL(fp), intent(inout), pointer :: GOCART_PBLH(:,:)                !< PBL height [m]
      REAL(fp), intent(inout), pointer :: GOCART_HFLUX(:,:)               !< sfc. sens. heat flux [W m-2]
      REAL(fp), intent(inout), pointer :: GOCART_Z0H(:,:)                 !< rough height, sens. heat [m]

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(GOCART_TMPU(1, 1, km))
      allocate(GOCART_RHOA(1, 1, km))
      allocate(GOCART_HGHTE(1, 1, 0:km))
      allocate(GOCART_U10(1, 1))
      allocate(GOCART_V10(1, 1))
      allocate(GOCART_FRACLAKE(1, 1))
      allocate(GOCART_GWETTOP(1, 1))
      allocate(GOCART_LWI(1, 1))
      allocate(GOCART_USTAR(1, 1))
      allocate(GOCART_PBLH(1, 1))
      allocate(GOCART_HFLUX(1, 1))
      allocate(GOCART_Z0H(1, 1))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)

      GOCART_TMPU(1,1,:) = tmpu(size(tmpu):1:-1) ! temperature [K]
      GOCART_RHOA(1,1,:) = rhoa(size(rhoa):1:-1) ! air density [kg/m^3]
      GOCART_HGHTE(1,1,:) = hghte(size(hghte):1:-1)    ! top of layer geopotential height [m]
      GOCART_LWI = real(LWI, fp)     ! orography flag; Land, ocean, ice mask
      GOCART_USTAR  = ustar

      ! friction speed [m/sec]
      GOCART_PBLH   = pblh      ! PBL height [m]
      GOCART_HFLUX = hflux     ! sfc. sens. heat flux [W m-2]
      GOCART_Z0H    = z0h       ! rough height, sens. heat [m]
      GOCART_U10 = u10m         ! zonal wind component (E/W) [m/s]
      GOCART_V10 = v10m         ! meridional wind component (N/S) [m/s]
      GOCART_FRACLAKE = fraclake   ! unitless, lake fraction (0-1)
      GOCART_GWETTOP = gwettop     ! unitless, soil moisture fraction (0-1)


   end subroutine PrepMetVarsForGOCART

end module DryDepScheme_GOCART_Mod
