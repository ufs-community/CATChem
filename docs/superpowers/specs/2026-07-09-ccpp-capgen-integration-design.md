# Specification: CCPP Capgen Dynamic Constituent Integration

* **Status:** Approved
* **Authors:** Gemini CLI Architect
* **Created:** July 9, 2026
* **Target Version:** 2.1.0
* **Pillars:** Decoupled Atmospheric Physics, Standard-Compliant CCPP Integration, Zero-Overhead Zero-Copy Bindings

## 1. Executive Summary & Architecture

To fully align the **CATChem** chemistry driver with modern, state-of-the-art CCPP (Common Community Physics Package) standards and support modern **`ccpp-capgen`** features, we are migrating the CCPP interface layer from static contiguous tracer assumptions to the dynamic **CCPP Constituent Registry** model.

Currently, host models must pre-allocate, define, and pass a contiguous chunk of 3D arrays to the chemistry scheme, creating a tight compile-time coupling between host tracers and chemistry solvers.

By incorporating modern `ccpp-capgen` capabilities:
1. **Tracer Decoupling:** CATChem declares its chemical species dynamically to the CCPP framework during a newly introduced **`_register`** phase.
2. **Dynamic Matching:** The host model maps and passes pointers to the unified CCPP tracer array `constituents` during the **`_init`** phase.
3. **Flexible Run-Time Layouts:** During the **`_run`** phase, CATChem dynamically maps, extracts, binds (via C++ Kokkos LayoutLeft bindings), and synchronizes the active constituents back in-place, eliminating the contiguity requirement and enabling complete platform compatibility.

---

## 2. Interface Subroutines Definition

### 2.1 Registration Phase (`ccpp_catchem_interface_register`)

This phase reads the species configuration, queries the species count from the C++ database, and instantiates the metadata using `ccpp_constituent_properties_t`.

```fortran
  subroutine ccpp_catchem_interface_register(constituent_props, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use machine,                   only: kind_phys
    use iso_c_binding,             only: c_null_char, c_ptr
    implicit none

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=*),                                 intent(out) :: errmsg
    integer,                                          intent(out) :: errflg

    ! Local variables
    character(len=512) :: config_path
    character(len=64)  :: c_name
    real(kind_phys)    :: molar_mass_g_mol, molar_mass_kg_mol
    logical            :: is_advected
    integer            :: n_spec, i
    type(c_ptr)        :: state_mgr

    errmsg = ''
    errflg = 0

    ! 1. Locate and load the YAML configuration path via environment variable
    call get_environment_variable("CATCHEM_CONFIG", config_path)
    if (trim(config_path) == "") then
       config_path = "./CATChem_config.yml"
    end if

    ! 2. Initialize a local state to parse species metadata
    call cc_model%initialize(config_path, 1, 1, 127, 3, 5, 20, errflg)
    if (errflg /= 0) then
       errmsg = 'CATChem Register Error: Failed to initialize model configuration.'
       return
    end if

    state_mgr = cc_model%get_state_manager()
    n_spec = int(catchem_state_get_species_count(state_mgr))

    ! 3. Allocate and instantiate constituent properties
    allocate(constituent_props(n_spec))

    do i = 1, n_spec
       call catchem_state_get_species_name_at(state_mgr, i, c_name)
       molar_mass_g_mol = real(catchem_state_get_species_mw(state_mgr, i), kind_phys)
       molar_mass_kg_mol = molar_mass_g_mol * 1.0e-3_kind_phys
       is_advected = (catchem_state_get_species_is_advected(state_mgr, i) == 1)

       call constituent_props(i)%instantiate( &
          std_name      = trim(c_name), &
          long_name     = trim(c_name), &
          diag_name     = trim(c_name), &
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
```

### 2.2 Initialization Phase (`ccpp_catchem_interface_init`)

The driver queries and stores the dynamic index mapping of each registered species from the CCPP environment using `ccpp_const_get_idx`.

```fortran
  subroutine ccpp_catchem_interface_init(im, do_catchem, catchem_configfile_in, &
                                         constituent_props_ptr, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx
    implicit none

    integer,                           intent(in)  :: im
    logical,                           intent(in)  :: do_catchem
    character(len=*),                  intent(in)  :: catchem_configfile_in
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props_ptr(:)
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    character(len=64) :: c_name
    integer           :: n_spec, i

    errmsg = ''
    errflg = 0

    if (.not. do_catchem) return

    ! 1. Initialize modern C++ Core orchestrator
    call cc_model%initialize(catchem_configfile_in, im, 1, 127, 3, 5, 20, errflg)
    if (errflg /= 0) then
       errmsg = 'CATChem Init Error: Failed to initialize C++ Core via cc_model'
       return
    end if

    ! 2. Parse and map each registered species to host tracer index
    n_spec = int(catchem_state_get_species_count(cc_model%get_state_manager()))

    if (allocated(catchem_indices_constituent_props)) deallocate(catchem_indices_constituent_props)
    allocate(catchem_indices_constituent_props(n_spec))

    do i = 1, n_spec
       call catchem_state_get_species_name_at(cc_model%get_state_manager(), i, c_name)

       call ccpp_const_get_idx(constituent_props_ptr, trim(c_name), &
                               catchem_indices_constituent_props(i), errmsg, errflg)
       if (errflg /= 0) then
          errmsg = "CATChem Init Error: Missing required tracer: " // trim(c_name)
          return
       end if
    end do

  end subroutine ccpp_catchem_interface_init
```

### 2.3 Execution Phase (`ccpp_catchem_interface_run`)

The execution routine extracts non-contiguous indices into a local contiguous buffer, binds it directly to Kokkos/C++ unmanaged views, runs the thread-safe calculations, and updates the global host array in-place.

```fortran
  subroutine ccpp_catchem_interface_run(im, kte, kme, garea, nsoil, nlndcat, nsoilcat, &
     lat, lon, do_catchem, dt, jdate, xcosz, &
     lwi, frlanduse, gvf, seaicefrac, oceanfrac, lakefrac, landfrac, &
     stype, vtype, snowdepth, frsnow, lai, frsoil, pores, resid, &
     ustar, u10m, v10m, tskin, ts, hf2d, lf2d, znt, prsfc, pblh, &
     dswsfc, nirbmdi, nirdfdi, visbmdi, visdfdi, &
     sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif, &
     soilmoist, pr3d, phl3d, prl3d, tk3d, q3d, us3d, vs3d, rh, &
     delp, airden, pfl_lsan, pfl_isan, &
     rain_cpl, cldf, dust_in, &
     constituent_props, constituents, &
     errmsg, errflg)

     ! All standard intent declarations as defined in the current API ...

     type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
     real(kind_phys), target,           intent(inout) :: constituents(:,:,:)

     real(kind_phys), target, allocatable :: subset_constituents(:,:,:)
     integer :: n_spec, i

     errmsg = ''
     errflg = 0

     if (.not. do_catchem) return

     ! 1. Direct dynamic subset extraction
     n_spec = size(catchem_indices_constituent_props)
     allocate(subset_constituents(im, kte, n_spec))

     do i = 1, n_spec
        subset_constituents(:,:,i) = constituents(:,:,catchem_indices_constituent_props(i))
     end do

     ! 2. Meterological Pointer Bindings (Exactly matches current zero-copy bindings)
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

     ! 3. Bind dynamic subset directly
     call cc_model%bind_unified_chemistry(c_loc(subset_constituents(1,1,1)))

     ! 4. Timestep Execution
     call cc_model%run_timestep(real(dt, fp), errflg)
     if (errflg /= 0) then
        errmsg = 'CATChem Run Error: Scheduled processes failed inside modern C++ core'
        deallocate(subset_constituents)
        return
     end if

     ! 5. Sync back changes
     do i = 1, n_spec
        constituents(:,:,catchem_indices_constituent_props(i)) = subset_constituents(:,:,i)
     end do

     deallocate(subset_constituents)

  end subroutine ccpp_catchem_interface_run
```

---

## 3. Metadata Configuration (`ccpp_catchem_interface.meta`)

The modernized `ccpp_catchem_interface.meta` declares the argument tables, replacing the static variables `ntrac`, `ntchs`, `ntchm`, `chemarr_phys`, and `chemarr` with dynamic constituents pointers:

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
