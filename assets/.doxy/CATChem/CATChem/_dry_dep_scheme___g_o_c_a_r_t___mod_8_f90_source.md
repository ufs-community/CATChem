

# File DryDepScheme\_GOCART\_Mod.F90

[**File List**](files.md) **>** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md) **>** [**schemes**](dir_5a3c86e36f17958630366ebc2b7ca21b.md) **>** [**DryDepScheme\_GOCART\_Mod.F90**](_dry_dep_scheme___g_o_c_a_r_t___mod_8_f90.md)

[Go to the documentation of this file](_dry_dep_scheme___g_o_c_a_r_t___mod_8_f90.md)


```Fortran

module drydepscheme_gocart_mod

   use precision_mod, only: fp
   use drydepcommon_mod, only: drydepschemegocartconfig
   use error_mod, only: cc_success, cc_error
   use constants, only: cp, g0, von_karman  !load the constants needed for this scheme

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gocart

   ! Additional physical constants (modify as needed for your scheme)
   real(fp), parameter :: OCEAN=0.0, land = 1.0, sea_ice = 2.0

contains


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
      USE gocart2g_process, only: drydeposition
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
      errmsg = ''
      thisloc = ' -> at compute_gocart (in DryDepScheme_GOCART_Mod.F90)'
      rc = cc_success
      vd = 0.0_fp
      drydepf = 0.0_fp

      ! Note: species_tendencies and diagnostic arrays are already initialized
      ! by the host ProcessInterface before calling this subroutine.
      ! Do not re-initialize them here.

      ! transform data for GOCART DryDeposition call
      call prepmetvarsforgocart(num_layers,     &
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
         gocart_tmpu,     &
         gocart_rhoa,     &
         gocart_hghte,    &
         gocart_u10,      &
         gocart_v10,      &
         gocart_fraclake, &
         gocart_gwettop,  &
         gocart_lwi,      &
         gocart_ustar,    &
         gocart_pblh,     &
         gocart_hflux,    &
         gocart_z0h)

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! Apply to each species
         do species_idx = 1, num_species
            ! Skip species that don't match scheme type (gas vs aerosol)
            if (is_gas(species_idx)) cycle

            if (params%resuspension) then
               call drydeposition(num_layers, gocart_tmpu, gocart_rhoa, gocart_hghte, gocart_lwi, gocart_ustar, &
                  gocart_pblh, gocart_hflux, von_karman, cp, g0, gocart_z0h, drydepf, rc, &
                  species_radius(species_idx)*1e-6_fp, species_density(species_idx), gocart_u10, gocart_v10, &
                  gocart_fraclake, gocart_gwettop)
            else
               call drydeposition(num_layers, gocart_tmpu, gocart_rhoa, gocart_hghte, gocart_lwi, gocart_ustar, &
                  gocart_pblh, gocart_hflux, von_karman, cp, g0, gocart_z0h, drydepf, rc)
            endif

            ! Ensure non-negative values
            species_tendencies(k, species_idx) = max(0.0_fp, drydepf(1,1)) * params%scale_factor
            !increase drydep frequency by factor of 5 for seasalt species to match GOCART2G. See the codes in:
            !https://github.com/GEOS-ESM/GOCART/blob/9ff3df9545dd582f415f682d3297e8c6c841e5cb/ESMF/GOCART2G_GridComp/SS2G_GridComp/SS2G_GridCompMod.F90#L820
            if (species_is_seasalt(species_idx) .and. abs(lwi - land) < 0.5) then
               species_tendencies(k, species_idx) = species_tendencies(k, species_idx) * 5.0_fp
            end if

            vd = max(species_tendencies(k, species_idx) * (z(k+1) -z(k) ), 1.e-4_fp)


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
                        max(0.0_fp, species_conc(k,species_idx) * (1.0_fp - exp(-1.0_fp * species_tendencies(k, species_idx) * tstep)))
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
                     drydep_velocity_per_species(diag_idx) = vd
                     exit
                  end if
               end do
            end if
         end do

      end do

      !cleanup pointers
      if (associated(gocart_tmpu)) nullify(gocart_tmpu)
      if (associated(gocart_rhoa)) nullify(gocart_rhoa)
      if (associated(gocart_hghte)) nullify(gocart_hghte)
      if (associated(gocart_u10)) nullify(gocart_u10)
      if (associated(gocart_v10)) nullify(gocart_v10)
      if (associated(gocart_fraclake)) nullify(gocart_fraclake)
      if (associated(gocart_gwettop)) nullify(gocart_gwettop)
      if (associated(gocart_lwi)) nullify(gocart_lwi)
      if (associated(gocart_ustar)) nullify(gocart_ustar)
      if (associated(gocart_lwi)) nullify(gocart_lwi)
      if (associated(gocart_hflux)) nullify(gocart_hflux)
      if (associated(gocart_z0h)) nullify(gocart_z0h)

   end subroutine compute_gocart

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines

   subroutine prepmetvarsforgocart(km,              &
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
      REAL(fp), intent(inout), pointer :: GOCART_TMPU(:,:,:)
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_RHOA
      REAL(fp), intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE
      REAL(fp), intent(inout), pointer :: GOCART_U10(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_V10 (:,:)
      REAL(fp), intent(inout), pointer :: GOCART_FRACLAKE(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_GWETTOP(:,:)
      real(fp), intent(inout), pointer :: GOCART_LWI(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_USTAR(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_PBLH(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_HFLUX(:,:)
      REAL(fp), intent(inout), pointer :: GOCART_Z0H(:,:)

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(gocart_tmpu(1, 1, km))
      allocate(gocart_rhoa(1, 1, km))
      allocate(gocart_hghte(1, 1, 0:km))
      allocate(gocart_u10(1, 1))
      allocate(gocart_v10(1, 1))
      allocate(gocart_fraclake(1, 1))
      allocate(gocart_gwettop(1, 1))
      allocate(gocart_lwi(1, 1))
      allocate(gocart_ustar(1, 1))
      allocate(gocart_pblh(1, 1))
      allocate(gocart_hflux(1, 1))
      allocate(gocart_z0h(1, 1))

      !Note: GOCART scheme expects vertical levels in reverse order (top to bottom)

      gocart_tmpu(1,1,:) = tmpu(size(tmpu):1:-1) ! temperature [K]
      gocart_rhoa(1,1,:) = rhoa(size(rhoa):1:-1) ! air density [kg/m^3]
      gocart_hghte(1,1,:) = hghte(size(hghte):1:-1)    ! top of layer geopotential height [m]
      gocart_lwi = real(lwi, fp)     ! orography flag; Land, ocean, ice mask
      gocart_ustar  = ustar

      ! friction speed [m/sec]
      gocart_pblh   = pblh      ! PBL height [m]
      gocart_hflux = hflux     ! sfc. sens. heat flux [W m-2]
      gocart_z0h    = z0h       ! rough height, sens. heat [m]
      gocart_u10 = u10m         ! zonal wind component (E/W) [m/s]
      gocart_v10 = v10m         ! meridional wind component (N/S) [m/s]
      gocart_fraclake = fraclake   ! unitless, lake fraction (0-1)
      gocart_gwettop = gwettop     ! unitless, soil moisture fraction (0-1)


   end subroutine prepmetvarsforgocart

end module drydepscheme_gocart_mod
```


