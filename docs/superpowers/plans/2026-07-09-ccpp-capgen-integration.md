# CCPP Capgen Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize the CATChem CCPP driver to incorporate modern CCPP capgen dynamic constituent registration and mapping capabilities.

**Architecture:** We introduce the standard `_register` phase (`ccpp_catchem_interface_register`) to dynamically declare the CATChem chemical species to CCPP. During the `_init` phase, we resolve each registered species standard name to its global host index via CCPP's standard `ccpp_const_get_idx` utility. During the `_run` phase, we extract non-contiguous global tracer mixing ratios into a contiguous local subset, bind it directly to Kokkos/C++ unmanaged views, execute the solver in-place, and synchronize the results back to the host model.

**Tech Stack:** Fortran 2003/2008, CCPP Framework APIs, Kokkos/C++ Core API bindings.

## Global Constraints
- Target File: `drivers/ccpp/ccpp_catchem_interface.meta`
- Target File: `drivers/ccpp/ccpp_catchem_interface.F90`
- Zero legacy Fortran tracer indexing dependency (remove `ntrac`, `ntchs`, `ntchm`, `chemarr_phys`, `chemarr`).
- Compliant with modern `ccpp-capgen` metadata specification.

---

### Task 1: Modernize CCPP Driver Metadata Specification

We replace the outdated contiguous 3D array interface with standard dynamic constituent metadata in the companion `.meta` file.

**Files:**
- Modify: `drivers/ccpp/ccpp_catchem_interface.meta`

**Interfaces:**
- Consumes: None (Metadata declaration only)
- Produces: Metadata parsed by CCPP `ccpp-capgen` to generate the calling CAP layer for `ccpp_catchem_interface_register`, `ccpp_catchem_interface_init`, and `ccpp_catchem_interface_run`.

- [ ] **Step 1: Write the updated metadata content**

Replace the entire contents of `drivers/ccpp/ccpp_catchem_interface.meta` with the following CCPP-compliant specifications:

```ini
[ccpp-table-properties]
  name = ccpp_catchem_interface
  type = scheme
  dependencies = ./catchem_types.F90,./catchem_wrapper_utils.F90,../../src/api/catchem.F90

########################################################################
[ccpp-arg-table]
  name = ccpp_catchem_interface_register
  type = scheme
[constituent_props]
  standard_name = dynamic_constituents_for_catchem
  units = none
  dimensions = (:)
  allocatable = True
  type = ccpp_constituent_properties_t
  intent = out
[errmsg]
  standard_name = ccpp_error_message
  units = none
  type = character | kind = len=*
  dimensions = ()
  intent = out
[errflg]
  standard_name = ccpp_error_code
  units = 1
  type = integer
  dimensions = ()
  intent = out

########################################################################
[ccpp-arg-table]
  name = ccpp_catchem_interface_init
  type = scheme
[im]
  standard_name = horizontal_dimension
  units = count
  dimensions = ()
  type = integer
  intent = in
[do_catchem]
  standard_name = do_catchem_coupling
  units = flag
  dimensions = ()
  type = logical
  intent = in
[catchem_configfile_in]
  standard_name = catchem_configfile
  units = none
  dimensions = ()
  type = character | kind = len=*
  intent = in
[constituent_props_ptr]
  standard_name = ccpp_constituent_properties
  units = none
  type = ccpp_constituent_prop_ptr_t
  dimensions = (number_of_ccpp_constituents)
  intent = in
[errmsg]
  standard_name = ccpp_error_message
  units = none
  type = character | kind = len=*
  dimensions = ()
  intent = out
[errflg]
  standard_name = ccpp_error_code
  units = 1
  type = integer
  dimensions = ()
  intent = out

########################################################################
[ccpp-arg-table]
  name = ccpp_catchem_interface_finalize
  type = scheme
[do_catchem]
  standard_name = do_catchem_coupling
  units = flag
  dimensions = ()
  type = logical
  intent = in
[errmsg]
  standard_name = ccpp_error_message
  units = none
  type = character | kind = len=*
  dimensions = ()
  intent = out
[errflg]
  standard_name = ccpp_error_code
  units = 1
  type = integer
  dimensions = ()
  intent = out

########################################################################
[ccpp-arg-table]
  name = ccpp_catchem_interface_run
  type = scheme
[im]
  standard_name = horizontal_loop_extent
  units = count
  dimensions = ()
  type = integer
  intent = in
[kte]
  standard_name = vertical_layer_dimension
  units = count
  dimensions = ()
  type = integer
  intent = in
[kme]
  standard_name = vertical_interface_dimension
  units = count
  dimensions = ()
  type = integer
  intent = in
[garea]
  standard_name = cell_area
  units = m2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[nsoil]
  standard_name = vertical_dimension_of_soil
  units = count
  dimensions = ()
  type = integer
  intent = in
[nlndcat]
  standard_name = number_of_vegetation_categories
  units = count
  dimensions = ()
  type = integer
  intent = in
[nsoilcat]
  standard_name = number_of_soil_categories
  units = count
  dimensions = ()
  type = integer
  intent = in
[lat]
  standard_name = latitude_in_degree
  units = degree_north
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[lon]
  standard_name = longitude_in_degree
  units = degree_east
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[do_catchem]
  standard_name = do_catchem_coupling
  units = flag
  dimensions = ()
  type = logical
  intent = in
[dt]
  standard_name = timestep_for_physics
  units = s
  dimensions = ()
  type = real | kind = kind_phys
  intent = in
[jdate]
  standard_name = current_date_and_time
  units = yyyymmdd_and_utc
  dimensions = (8)
  type = integer
  intent = in
[xcosz]
  standard_name = cosine_of_solar_zenith_angle
  units = none
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[lwi]
  standard_name = land_sea_mask
  units = index
  dimensions = (horizontal_loop_extent)
  type = integer
  intent = in
[frlanduse]
  standard_name = fraction_of_vegetation_type
  units = fraction
  dimensions = (horizontal_loop_extent,number_of_vegetation_categories)
  type = real | kind = kind_phys
  intent = in
[gvf]
  standard_name = green_vegetation_fraction
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[seaicefrac]
  standard_name = sea_ice_concentration
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[oceanfrac]
  standard_name = fraction_of_ocean
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[lakefrac]
  standard_name = fraction_of_lake
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[landfrac]
  standard_name = fraction_of_land
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[stype]
  standard_name = soil_type_index
  units = index
  dimensions = (horizontal_loop_extent)
  type = integer
  intent = in
[vtype]
  standard_name = vegetation_type_index
  units = index
  dimensions = (horizontal_loop_extent)
  type = integer
  intent = in
[snowdepth]
  standard_name = physical_snow_depth
  units = m
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[frsnow]
  standard_name = snow_cover_fraction
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[lai]
  standard_name = leaf_area_index
  units = m2 m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[frsoil]
  standard_name = fraction_of_soil_type
  units = fraction
  dimensions = (horizontal_loop_extent,number_of_soil_categories)
  type = real | kind = kind_phys
  intent = in
[pores]
  standard_name = dry_soil_porosity
  units = none
  dimensions = (30)
  type = real | kind = kind_phys
  intent = in
[resid]
  standard_name = residual_soil_moisture
  units = none
  dimensions = (30)
  type = real | kind = kind_phys
  intent = in
[ustar]
  standard_name = friction_velocity
  units = m s-1
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[u10m]
  standard_name = zonal_wind_at_10m
  units = m s-1
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[v10m]
  standard_name = meridional_wind_at_10m
  units = m s-1
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[tskin]
  standard_name = ground_skin_temperature
  units = K
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[ts]
  standard_name = surface_temperature_for_air_sea_fluxes
  units = K
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[hf2d]
  standard_name = sensible_heat_flux
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[lf2d]
  standard_name = latent_heat_flux
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[znt]
  standard_name = roughness_length_for_momentum
  units = m
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[prsfc]
  standard_name = surface_air_pressure
  units = Pa
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[pblh]
  standard_name = planetary_boundary_layer_height
  units = m
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[dswsfc]
  standard_name = downward_shortwave_flux_at_surface
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[nirbmdi]
  standard_name = downward_near_infrared_direct_flux
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[nirdfdi]
  standard_name = downward_near_infrared_diffuse_flux
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[visbmdi]
  standard_name = downward_visible_direct_flux
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[visdfdi]
  standard_name = downward_visible_diffuse_flux
  units = W m-2
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[sfc_alb_nir_dir]
  standard_name = surface_albedo_due_to_NIR_direct
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[sfc_alb_nir_dif]
  standard_name = surface_albedo_due_to_NIR_diffuse
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[sfc_alb_uvvis_dir]
  standard_name = surface_albedo_due_to_UV_and_VIS_direct
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[sfc_alb_uvvis_dif]
  standard_name = surface_albedo_due_to_UV_and_VIS_diffuse
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[soilmoist]
  standard_name = soil_moisture_fraction
  units = m3 m-3
  dimensions = (horizontal_loop_extent,vertical_dimension_of_soil)
  type = real | kind = kind_phys
  intent = in
[pr3d]
  standard_name = air_pressure_at_interface
  units = Pa
  dimensions = (horizontal_loop_extent,vertical_interface_dimension)
  type = real | kind = kind_phys
  intent = in
[phl3d]
  standard_name = geopotential_height_at_interface
  units = m
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[prl3d]
  standard_name = air_pressure
  units = Pa
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[tk3d]
  standard_name = air_temperature
  units = K
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[q3d]
  standard_name = water_vapor_mixing_ratio_wrt_dry_air
  units = kg kg-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[us3d]
  standard_name = eastward_wind
  units = m s-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[vs3d]
  standard_name = northward_wind
  units = m s-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[rh]
  standard_name = relative_humidity
  units = percent
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[delp]
  standard_name = air_pressure_thickness
  units = Pa
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[airden]
  standard_name = dry_air_density
  units = kg m-3
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[pfl_lsan]
  standard_name = cloud_liquid_water_mixing_ratio
  units = kg kg-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[pfl_isan]
  standard_name = cloud_ice_water_mixing_ratio
  units = kg kg-1
  dimensions = (horizontal_loop_extent,vertical_layer_dimension)
  type = real | kind = kind_phys
  intent = in
[rain_cpl]
  standard_name = precipitation_rate
  units = kg m-2 s-1
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[cldf]
  standard_name = cloud_area_fraction
  units = fraction
  dimensions = (horizontal_loop_extent)
  type = real | kind = kind_phys
  intent = in
[dust_in]
  standard_name = dust_emission_flux
  units = kg m-2 s-1
  dimensions = (horizontal_loop_extent,12,5)
  type = real | kind = kind_phys
  intent = in
[constituent_props]
  standard_name = ccpp_constituent_properties
  units = none
  type = ccpp_constituent_prop_ptr_t
  dimensions = (number_of_ccpp_constituents)
  intent = in
[constituents]
  standard_name = ccpp_constituents
  units = none
  type = real | kind = kind_phys
  dimensions = (horizontal_loop_extent,vertical_layer_dimension,number_of_ccpp_constituents)
  intent = inout
[errmsg]
  standard_name = ccpp_error_message
  units = none
  type = character | kind = len=*
  dimensions = ()
  intent = out
[errflg]
  standard_name = ccpp_error_code
  units = 1
  type = integer
  dimensions = ()
  intent = out
```

- [ ] **Step 2: Commit metadata changes**

```bash
git add drivers/ccpp/ccpp_catchem_interface.meta
git commit -m "feat(ccpp): modernize metadata spec with dynamic constituent tables"
```

---

### Task 2: Implement CCPP Capgen Dynamic Driver Functions

We implement `ccpp_catchem_interface_register` for dynamic species metadata registration, `ccpp_catchem_interface_init` with pointer index resolution, and `ccpp_catchem_interface_run` with dynamic extraction.

**Files:**
- Modify: `drivers/ccpp/ccpp_catchem_interface.F90`

**Interfaces:**
- Consumes: `catchem_state_get_species_count`, `catchem_state_get_species_name_at`, `catchem_state_get_species_mw`, `catchem_state_get_species_is_advected` from `CATChem_API.F90`.
- Produces: Subroutines `ccpp_catchem_interface_register`, `ccpp_catchem_interface_init`, `ccpp_catchem_interface_run`, and `ccpp_catchem_interface_finalize`.

- [ ] **Step 1: Replace contents of `drivers/ccpp/ccpp_catchem_interface.F90`**

Write the complete implementation:

```fortran
!> \file ccpp_catchem_interface.F90
!! \brief CCPP interface module for CATChem integration with dynamic constituents
!!
module ccpp_catchem_interface

  use iso_c_binding,             only: c_loc, c_null_char, c_double, c_ptr, c_associated, c_int, c_char
  use machine,                   only: kind_phys
  use CATChem_API,               only: CATChem_Model
  use Error_Mod,                 only: CC_SUCCESS, CC_FAILURE
  use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t, ccpp_constituent_prop_ptr_t
  use precision_mod,             only: fp

  implicit none

  private

  public :: ccpp_catchem_interface_register, &
            ccpp_catchem_interface_init, &
            ccpp_catchem_interface_run, &
            ccpp_catchem_interface_finalize

  ! Module-scope wrapper mapping directly to the C++ core
  type(CATChem_Model), save :: cc_model

  ! Saved dynamic index mapper array to store host index for each species
  integer, allocatable, save :: catchem_indices_constituent_props(:)

  ! Linkable BIND(C) external declarations to avoid circular dependencies
  interface
     integer(c_int) function catchem_state_get_species_count(state_ptr) bind(C, name="catchem_state_get_species_count")
        import :: c_ptr, c_int
        type(c_ptr), value :: state_ptr
     end function

     subroutine catchem_state_get_species_name_at(state_ptr, index, name_out) bind(C, name="catchem_state_get_species_name_at")
        import :: c_ptr, c_int, c_char
        type(c_ptr), value :: state_ptr
        integer(c_int), value :: index
        character(kind=c_char), intent(out) :: name_out(*)
     end subroutine

     real(c_double) function catchem_state_get_species_mw(state_ptr, index) bind(C, name="catchem_state_get_species_mw")
        import :: c_ptr, c_int, c_double
        type(c_ptr), value :: state_ptr
        integer(c_int), value :: index
     end function

     integer(c_int) function catchem_state_get_species_is_advected(state_ptr, index) bind(C, name="catchem_state_get_species_is_advected")
        import :: c_ptr, c_int
        type(c_ptr), value :: state_ptr
        integer(c_int), value :: index
     end function
  end interface

contains

   ! Helper to convert null-terminated C buffers back to Fortran fixed strings
   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in) :: c_str(*)
      character(len=*), intent(out) :: f_str
      integer :: i
      f_str = ""
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

   !> \brief Register dynamic constituent properties with CCPP
   subroutine ccpp_catchem_interface_register(constituent_props, errmsg, errflg)
      implicit none

      type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
      character(len=*),                                 intent(out) :: errmsg
      integer,                                          intent(out) :: errflg

      ! Local variables
      character(len=512)     :: config_path
      character(kind=c_char) :: c_buf(128)
      character(len=64)      :: f_name
      real(kind_phys)        :: molar_mass_g_mol, molar_mass_kg_mol
      logical                :: is_advected
      integer                :: n_spec, i
      type(c_ptr)            :: state_mgr

      errmsg = ''
      errflg = 0

      ! 1. Locate YAML configuration path via environment variable
      call get_environment_variable("CATCHEM_CONFIG", config_path)
      if (trim(config_path) == "") then
         config_path = "./tests/CATChem_config.yml"
      end if

      ! 2. Lightweight core initialization to load species configuration
      call cc_model%initialize(config_path, 1, 1, 127, 3, 5, 20, errflg)
      if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Register: Failed to load species configuration.'
         return
      end if

      state_mgr = cc_model%get_state_manager()
      n_spec = int(catchem_state_get_species_count(state_mgr))

      ! 3. Allocate and populate constituent_props
      allocate(constituent_props(n_spec), stat=errflg)
      if (errflg /= 0) then
         errmsg = 'CATChem Register: Memory allocation failure.'
         return
      end if

      do i = 1, n_spec
         call catchem_state_get_species_name_at(state_mgr, int(i, c_int), c_buf)
         call c_to_f_string(c_buf, f_name)

         molar_mass_g_mol = real(catchem_state_get_species_mw(state_mgr, int(i, c_int)), kind_phys)
         molar_mass_kg_mol = molar_mass_g_mol * 1.0e-3_kind_phys
         is_advected = (catchem_state_get_species_is_advected(state_mgr, int(i, c_int)) == 1)

         call constituent_props(i)%instantiate( &
            std_name      = trim(f_name), &
            long_name     = trim(f_name), &
            diag_name     = trim(f_name), &
            units         = 'kg kg-1', &
            vertical_dim  = 'vertical_layer_dimension', &
            default_value = 0.0_kind_phys, &
            min_value     = 0.0_kind_phys, &
            molar_mass    = molar_mass_kg_mol, &
            advected      = is_advected, &
            errcode       = errflg, &
            errmsg        = errmsg )
         if (errflg /= 0) return
      end do

   end subroutine ccpp_catchem_interface_register

   !> \brief Initialize the CATChem CCPP interface
   subroutine ccpp_catchem_interface_init(im, do_catchem, catchem_configfile_in, &
                                          constituent_props_ptr, errmsg, errflg)
      use ccpp_const_utils, only: ccpp_const_get_idx
      implicit none

      integer,                           intent(in)  :: im
      logical,                           intent(in)  :: do_catchem
      character(len=*),                  intent(in)  :: catchem_configfile_in
      type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props_ptr(:)
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      character(kind=c_char) :: c_buf(128)
      character(len=64)      :: f_name
      integer                :: n_spec, i
      type(c_ptr)            :: state_mgr

      errmsg = ''
      errflg = 0

      if (.not. do_catchem) return

      ! 1. Fully initialize C++ Core manager
      call cc_model%initialize(catchem_configfile_in, im, 1, 127, 3, 5, 20, errflg)
      if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Init: Failed to initialize C++ Core via cc_model'
         return
      end if

      state_mgr = cc_model%get_state_manager()
      n_spec = int(catchem_state_get_species_count(state_mgr))

      ! 2. Resolve index mappings from CCPP global constituent properties pointer
      if (allocated(catchem_indices_constituent_props)) deallocate(catchem_indices_constituent_props)
      allocate(catchem_indices_constituent_props(n_spec))

      do i = 1, n_spec
         call catchem_state_get_species_name_at(state_mgr, int(i, c_int), c_buf)
         call c_to_f_string(c_buf, f_name)

         call ccpp_const_get_idx(constituent_props_ptr, trim(f_name), &
                                 catchem_indices_constituent_props(i), errmsg, errflg)
         if (errflg /= 0) then
            errmsg = "CATChem Init: Missing required tracer: " // trim(f_name)
            return
         end if
      end do

   end subroutine ccpp_catchem_interface_init

  !> \brief Finalize the CATChem CCPP interface
  subroutine ccpp_catchem_interface_finalize(do_catchem, errmsg, errflg)
     implicit none

     logical, intent(in) :: do_catchem

     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     errmsg = ''
     errflg = 0

     if (.not. do_catchem) return

     call cc_model%finalize(errflg)
     if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Finalize: Error finalising C++ Core via cc_model'
     end if

     if (allocated(catchem_indices_constituent_props)) deallocate(catchem_indices_constituent_props)

  end subroutine ccpp_catchem_interface_finalize

  !> \brief Execute CATChem chemistry calculations with dynamic tracer support
  subroutine ccpp_catchem_interface_run(im, kte, kme, garea, nsoil, nlndcat, nsoilcat, &
     lat, lon, &
     do_catchem, &
     dt, jdate, &
     xcosz, &
     lwi, frlanduse, gvf, seaicefrac, oceanfrac, lakefrac, landfrac, &
     stype, vtype, snowdepth, frsnow, lai, frsoil, pores, resid, &
     ustar, u10m, v10m, tskin, ts, hf2d, lf2d, znt, prsfc, pblh, &
     dswsfc, nirbmdi, nirdfdi, visbmdi, visdfdi, &
     sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, &
     soilmoist, pr3d, phl3d, prl3d, tk3d, q3d, us3d, vs3d, rh, &
     delp, airden, pfl_lsan, pfl_isan, &
     rain_cpl, cldf, &
     dust_in, &
     constituent_props, constituents, &
     errmsg, errflg)

     implicit none

     integer, intent(in) :: im
     integer, intent(in) :: kte
     integer, intent(in) :: kme
     integer, intent(in) :: nsoil
     integer, intent(in) :: nlndcat
     integer, intent(in) :: nsoilcat
     real(kind_phys), dimension(im), intent(in), target :: garea
     real(kind_phys), dimension(im), intent(in), target :: lat
     real(kind_phys), dimension(im), intent(in), target :: lon
     real(kind_phys), dimension(im), intent(in), target :: xcosz

     real(kind_phys), intent(in) :: dt
     integer, intent(in) :: jdate(8)

     logical, intent(in) :: do_catchem

     integer, dimension(im), intent(in), target                :: lwi
     integer, dimension(im), intent(in), target                :: stype
     integer, dimension(im), intent(in), target                :: vtype

     real(kind_phys), dimension(im, nlndcat), intent(in), target :: frlanduse
     real(kind_phys), dimension(im, nsoilcat), intent(in), target :: frsoil
     real(kind_phys), dimension(30), intent(in), target        :: pores
     real(kind_phys), dimension(30), intent(in), target        :: resid
     real(kind_phys), dimension(im), intent(in), target        :: seaicefrac
     real(kind_phys), dimension(im), intent(in), target        :: oceanfrac
     real(kind_phys), dimension(im), intent(in), target        :: frsnow
     real(kind_phys), dimension(im), intent(in), target        :: lakefrac
     real(kind_phys), dimension(im), intent(in), target        :: landfrac
     real(kind_phys), dimension(im), intent(in), target        :: gvf
     real(kind_phys), dimension(im), intent(in), target        :: lai

     real(kind_phys), dimension(im, nsoil), intent(in), target :: soilmoist
     real(kind_phys), dimension(im), intent(in), target        :: snowdepth
     real(kind_phys), dimension(im), intent(in), target        :: prsfc
     real(kind_phys), dimension(im), intent(in), target        :: pblh
     real(kind_phys), dimension(im), intent(in), target        :: u10m
     real(kind_phys), dimension(im), intent(in), target        :: v10m
     real(kind_phys), dimension(im), intent(in), target        :: ustar
     real(kind_phys), dimension(im), intent(in), target        :: tskin
     real(kind_phys), dimension(im), intent(in), target        :: ts
     real(kind_phys), dimension(im), intent(in), target        :: hf2d
     real(kind_phys), dimension(im), intent(in), target        :: lf2d
     real(kind_phys), dimension(im), intent(in), target        :: znt
     real(kind_phys), dimension(im), intent(in), target        :: dswsfc
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_nir_dir
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_nir_dif
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_uvvis_dir
     real(kind_phys), dimension(im), intent(in), target        :: sfc_alb_uvvis_dif
     real(kind_phys), dimension(im), intent(in), target        :: nirbmdi
     real(kind_phys), dimension(im), intent(in), target        :: nirdfdi
     real(kind_phys), dimension(im), intent(in), target        :: visbmdi
     real(kind_phys), dimension(im), intent(in), target        :: visdfdi

     real(kind_phys), dimension(im, kme), intent(in), target :: pr3d
     real(kind_phys), dimension(im, kte), intent(in), target :: prl3d
     real(kind_phys), dimension(im, kte), intent(in), target :: delp
     real(kind_phys), dimension(im, kte), intent(in), target :: phl3d
     real(kind_phys), dimension(im, kte), intent(in), target :: tk3d
     real(kind_phys), dimension(im, kte), intent(in), target :: us3d
     real(kind_phys), dimension(im, kte), intent(in), target :: vs3d
     real(kind_phys), dimension(im, kte), intent(in), target :: q3d
     real(kind_phys), dimension(im, kte), intent(in), target :: airden
     real(kind_phys), dimension(im, kte), intent(in), target :: rh
     real(kind_phys), dimension(im, kte), intent(in), target :: pfl_lsan
     real(kind_phys), dimension(im, kte), intent(in), target :: pfl_isan

     real(kind_phys), dimension(im), intent(in), target        :: rain_cpl
     real(kind_phys), dimension(im), intent(in), target        :: cldf
     real(kind_phys), dimension(im, 12, 5), intent(in), target :: dust_in

     type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
     real(kind_phys), target,           intent(inout) :: constituents(:,:,:)

     character(len=*), intent(out) :: errmsg
     integer,          intent(out) :: errflg

     ! Local variables
     real(kind_phys), target, allocatable :: subset_constituents(:,:,:)
     integer :: n_spec, i

     errmsg = ''
     errflg = 0

     if (.not. do_catchem) return

     ! 1. Extract non-contiguous indices into a contiguous local subset array
     n_spec = size(catchem_indices_constituent_props)
     allocate(subset_constituents(im, kte, n_spec))

     do i = 1, n_spec
        subset_constituents(:,:,i) = constituents(:,:,catchem_indices_constituent_props(i))
     end do

     ! 2. Meterological bindings (Direct unmanaged LayoutLeft C++ Views mappings)
     call cc_model%bind_met_3d("T"//c_null_char, c_loc(tk3d(1,1)))
     call cc_model%bind_met_3d("QV"//c_null_char, c_loc(q3d(1,1)))
     call cc_model%bind_met_3d("RH"//c_null_char, c_loc(rh(1,1)))
     call cc_model%bind_met_3d("PMID"//c_null_char, c_loc(prl3d(1,1)))
     call cc_model%bind_met_3d("PEDGE"//c_null_char, c_loc(pr3d(1,1)))
     call cc_model%bind_met_3d("DELP"//c_null_char, c_loc(delp(1,1)))
     call cc_model%bind_met_3d("AIRDEN"//c_null_char, c_loc(airden(1,1)))

     call cc_model%bind_met_2d("PS"//c_null_char, c_loc(prsfc(1)))
     call cc_model%bind_met_2d("TS"//c_null_char, c_loc(ts(1)))
     call cc_model%bind_met_2d("LAT"//c_null_char, c_loc(lat(1)))
     call cc_model%bind_met_2d("LON"//c_null_char, c_loc(lon(1)))
     call cc_model%bind_met_2d("PBLH"//c_null_char, c_loc(pblh(1)))
     call cc_model%bind_met_2d("USTAR"//c_null_char, c_loc(ustar(1)))
     call cc_model%bind_met_2d("HFLUX"//c_null_char, c_loc(hf2d(1)))
     call cc_model%bind_met_2d("AREA_M2"//c_null_char, c_loc(garea(1)))

     ! 3. Bind local unified chemistry concentrations subset
     call cc_model%bind_unified_chemistry(c_loc(subset_constituents(1,1,1)))

     ! 4. Central C++ Core timestep solver execution
     call cc_model%run_timestep(real(dt, fp), errflg)
     if (errflg /= CC_SUCCESS) then
         errmsg = 'CATChem Run Error: Scheduled process execution failed inside C++ Core.'
         deallocate(subset_constituents)
         return
     end if

     ! 5. Dynamic synchronisation back to the global host tracer array
     do i = 1, n_spec
        constituents(:,:,catchem_indices_constituent_props(i)) = subset_constituents(:,:,i)
     end do

     deallocate(subset_constituents)

   end subroutine ccpp_catchem_interface_run

end module ccpp_catchem_interface
```

- [ ] **Step 2: Commit driver implementation**

```bash
git add drivers/ccpp/ccpp_catchem_interface.F90
git commit -m "feat(ccpp): implement dynamic registration, mapping, and extraction logic"
```

---

### Task 3: Modernize Driver Interface Preservation Test

We modernize `tests/test_DriverInterfacePreservation.cmake` to verify `CATChem_API.F90` instead of `catchem.F90`, as the legacy modules and `catchem.F90` have been fully retired in favor of the C++ Core in the current stable branch.

**Files:**
- Modify: `tests/test_DriverInterfacePreservation.cmake`

**Interfaces:**
- Consumes: `drivers/ccpp/ccpp_catchem_interface.F90`, `src/api/CATChem_API.F90`
- Produces: Correctly executing driver preservation verification test script.

- [ ] **Step 1: Replace legacy checks in `tests/test_DriverInterfacePreservation.cmake`**

Modify `tests/test_DriverInterfacePreservation.cmake` to point to `CATChem_API.F90` and modern APIs.

```cmake
# test_DriverInterfacePreservation.cmake
# Verification that CCPP and NUOPC driver caps compile and conform
# to interface preservation requirements with the C++ Core API.

cmake_minimum_required(VERSION 3.10)

function(assert_file_contains filepath search_string description)
  if(NOT EXISTS "${filepath}")
    message(FATAL_ERROR "FAILED: ${description}\n  File not found: ${filepath}")
  endif()
  file(READ "${filepath}" file_content)
  string(FIND "${file_content}" "${search_string}" found_pos)
  if(found_pos EQUAL -1)
    message(
      FATAL_ERROR
      "FAILED: ${description}\n  File: ${filepath}\n  Missing: ${search_string}"
    )
  else()
    message(STATUS "PASSED: ${description}")
  endif()
endfunction()

set(SRC_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")

message(STATUS "=== Driver Interface Preservation Tests ===")

# Verify CATChem CCPP interface exposes the dynamic capgen routines
assert_file_contains(
  "${SRC_ROOT}/drivers/ccpp/ccpp_catchem_interface.F90"
  "subroutine ccpp_catchem_interface_register"
  "CCPP Driver exposes register subroutine"
)

assert_file_contains(
  "${SRC_ROOT}/drivers/ccpp/ccpp_catchem_interface.F90"
  "subroutine ccpp_catchem_interface_init"
  "CCPP Driver exposes init subroutine"
)

assert_file_contains(
  "${SRC_ROOT}/drivers/ccpp/ccpp_catchem_interface.F90"
  "subroutine ccpp_catchem_interface_run"
  "CCPP Driver exposes run subroutine"
)

# Verify Modern Public API
assert_file_contains(
  "${SRC_ROOT}/src/api/CATChem_API.F90"
  "type CATChem_Model"
  "Public API exports CATChem_Model"
)

assert_file_contains(
  "${SRC_ROOT}/src/api/CATChem_API.F90"
  "procedure :: initialize"
  "CATChem_Model has initialize binding"
)

assert_file_contains(
  "${SRC_ROOT}/src/api/CATChem_API.F90"
  "procedure :: run_timestep"
  "CATChem_Model has run_timestep binding"
)

message(STATUS "All Driver Interface Preservation checks PASSED!")
```

- [ ] **Step 2: Verify and run the modernized check**

Run: `cmake -P tests/test_DriverInterfacePreservation.cmake`
Expected output: All Driver Interface Preservation checks PASSED!

- [ ] **Step 3: Commit verification test changes**

```bash
git add tests/test_DriverInterfacePreservation.cmake
git commit -m "test(ccpp): modernize driver preservation tests for C++ Core API"
```
