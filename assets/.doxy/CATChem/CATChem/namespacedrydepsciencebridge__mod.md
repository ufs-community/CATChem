

# Namespace drydepsciencebridge\_mod



[**Namespace List**](namespaces.md) **>** [**drydepsciencebridge\_mod**](namespacedrydepsciencebridge__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**run\_drydep\_science\_bridge**](#function-run_drydep_science_bridge) (integer(c\_int), value n\_cols, integer(c\_int), value n\_levels, integer(c\_int), value n\_species, real(c\_double), value dt, character(kind=c\_char), dimension(\*), intent(in) gas\_scheme, character(kind=c\_char), dimension(\*), intent(in) aero\_scheme, integer(c\_int), value diagnostics, real(c\_double), value wesely\_scale\_factor, integer(c\_int), value wesely\_co2\_effect, real(c\_double), value wesely\_co2\_level, real(c\_double), value wesely\_co2\_reference, real(c\_double), value gocart\_scale\_factor, integer(c\_int), value gocart\_resuspension, integer(c\_int), value gocart\_dust\_resusp\_only, real(c\_double), value zhang\_scale\_factor, integer(c\_int), value wesely\_skip\_so2, integer(c\_int), value gocart\_skip\_sulfate\_aero, type(c\_ptr), value c\_bxheight, type(c\_ptr), value c\_airden, type(c\_ptr), value c\_t\_air, type(c\_ptr), value c\_z\_edges, type(c\_ptr), value c\_rh, type(c\_ptr), value c\_cldfrc, type(c\_ptr), value c\_frlai, type(c\_ptr), value c\_frlanduse, type(c\_ptr), value c\_iland, type(c\_ptr), value c\_is\_ice, type(c\_ptr), value c\_is\_land, type(c\_ptr), value c\_is\_snow, type(c\_ptr), value c\_lat, type(c\_ptr), value c\_lon, type(c\_ptr), value c\_obk, type(c\_ptr), value c\_ps, type(c\_ptr), value c\_salinity, type(c\_ptr), value c\_suncosmid, type(c\_ptr), value c\_swgdn, type(c\_ptr), value c\_ts, type(c\_ptr), value c\_tskin, type(c\_ptr), value c\_ustar, type(c\_ptr), value c\_z0, type(c\_ptr), value c\_frlake, type(c\_ptr), value c\_gwettop, type(c\_ptr), value c\_hflux, type(c\_ptr), value c\_lwi, type(c\_ptr), value c\_pblh, type(c\_ptr), value c\_u10m, type(c\_ptr), value c\_v10m, type(c\_ptr), value c\_z0h, real(c\_double), dimension(n\_species), intent(in) species\_mw\_g, real(c\_double), dimension(n\_species), intent(in) species\_dd\_f0, real(c\_double), dimension(n\_species), intent(in) species\_dd\_hstar, real(c\_double), dimension(n\_species), intent(in) species\_dd\_dvzaersnow, real(c\_double), dimension(n\_species), intent(in) species\_dd\_dvzminval\_snow, real(c\_double), dimension(n\_species), intent(in) species\_dd\_dvzminval\_land, real(c\_double), dimension(n\_species), intent(in) species\_density, real(c\_double), dimension(n\_species), intent(in) species\_radius, logical(c\_bool), dimension(n\_species), intent(in) species\_is\_seasalt, logical(c\_bool), dimension(n\_species), intent(in) species\_is\_dust, real(c\_double), dimension(n\_species), intent(in) species\_lower\_radius, real(c\_double), dimension(n\_species), intent(in) species\_upper\_radius, logical(c\_bool), dimension(n\_species), intent(in) is\_gas\_arr, type(c\_ptr), value c\_conc, type(c\_ptr), value c\_tendency, character(kind=c\_char), dimension(32,n\_species), intent(in) species\_names, type(c\_ptr), value c\_diag\_con, type(c\_ptr), value c\_diag\_vel, integer(c\_int), dimension(n\_diag\_species), intent(in) diagnostic\_species\_id, integer(c\_int), value n\_diag\_species) <br> |




























## Public Functions Documentation




### function run\_drydep\_science\_bridge 

```Fortran
subroutine drydepsciencebridge_mod::run_drydep_science_bridge (
    integer(c_int), value n_cols,
    integer(c_int), value n_levels,
    integer(c_int), value n_species,
    real(c_double), value dt,
    character(kind=c_char), dimension(*), intent(in) gas_scheme,
    character(kind=c_char), dimension(*), intent(in) aero_scheme,
    integer(c_int), value diagnostics,
    real(c_double), value wesely_scale_factor,
    integer(c_int), value wesely_co2_effect,
    real(c_double), value wesely_co2_level,
    real(c_double), value wesely_co2_reference,
    real(c_double), value gocart_scale_factor,
    integer(c_int), value gocart_resuspension,
    integer(c_int), value gocart_dust_resusp_only,
    real(c_double), value zhang_scale_factor,
    integer(c_int), value wesely_skip_so2,
    integer(c_int), value gocart_skip_sulfate_aero,
    type(c_ptr), value c_bxheight,
    type(c_ptr), value c_airden,
    type(c_ptr), value c_t_air,
    type(c_ptr), value c_z_edges,
    type(c_ptr), value c_rh,
    type(c_ptr), value c_cldfrc,
    type(c_ptr), value c_frlai,
    type(c_ptr), value c_frlanduse,
    type(c_ptr), value c_iland,
    type(c_ptr), value c_is_ice,
    type(c_ptr), value c_is_land,
    type(c_ptr), value c_is_snow,
    type(c_ptr), value c_lat,
    type(c_ptr), value c_lon,
    type(c_ptr), value c_obk,
    type(c_ptr), value c_ps,
    type(c_ptr), value c_salinity,
    type(c_ptr), value c_suncosmid,
    type(c_ptr), value c_swgdn,
    type(c_ptr), value c_ts,
    type(c_ptr), value c_tskin,
    type(c_ptr), value c_ustar,
    type(c_ptr), value c_z0,
    type(c_ptr), value c_frlake,
    type(c_ptr), value c_gwettop,
    type(c_ptr), value c_hflux,
    type(c_ptr), value c_lwi,
    type(c_ptr), value c_pblh,
    type(c_ptr), value c_u10m,
    type(c_ptr), value c_v10m,
    type(c_ptr), value c_z0h,
    real(c_double), dimension(n_species), intent(in) species_mw_g,
    real(c_double), dimension(n_species), intent(in) species_dd_f0,
    real(c_double), dimension(n_species), intent(in) species_dd_hstar,
    real(c_double), dimension(n_species), intent(in) species_dd_dvzaersnow,
    real(c_double), dimension(n_species), intent(in) species_dd_dvzminval_snow,
    real(c_double), dimension(n_species), intent(in) species_dd_dvzminval_land,
    real(c_double), dimension(n_species), intent(in) species_density,
    real(c_double), dimension(n_species), intent(in) species_radius,
    logical(c_bool), dimension(n_species), intent(in) species_is_seasalt,
    logical(c_bool), dimension(n_species), intent(in) species_is_dust,
    real(c_double), dimension(n_species), intent(in) species_lower_radius,
    real(c_double), dimension(n_species), intent(in) species_upper_radius,
    logical(c_bool), dimension(n_species), intent(in) is_gas_arr,
    type(c_ptr), value c_conc,
    type(c_ptr), value c_tendency,
    character(kind=c_char), dimension(32,n_species), intent(in) species_names,
    type(c_ptr), value c_diag_con,
    type(c_ptr), value c_diag_vel,
    integer(c_int), dimension(n_diag_species), intent(in) diagnostic_species_id,
    integer(c_int), value n_diag_species
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/DryDepScienceBridge.F90`

