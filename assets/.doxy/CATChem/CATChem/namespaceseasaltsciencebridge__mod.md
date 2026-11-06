

# Namespace seasaltsciencebridge\_mod



[**Namespace List**](namespaces.md) **>** [**seasaltsciencebridge\_mod**](namespaceseasaltsciencebridge__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**run\_seasalt\_science\_bridge**](#function-run_seasalt_science_bridge) (integer(c\_int), value n\_cols, integer(c\_int), value n\_levels, integer(c\_int), value n\_species, integer(c\_int), value n\_total\_species, real(c\_double), value dt, character(kind=c\_char), dimension(\*), intent(in) active\_scheme, integer(c\_int), value diagnostics, real(c\_double), value gong97\_scale\_factor, integer(c\_int), value gong97\_weibull\_flag, real(c\_double), value gong03\_scale\_factor, integer(c\_int), value gong03\_weibull\_flag, real(c\_double), value geos12\_scale\_factor, integer(c\_int), value geos12\_weibull\_flag, type(c\_ptr), value c\_frocean, type(c\_ptr), value c\_frseaice, type(c\_ptr), value c\_lat, type(c\_ptr), value c\_lon, type(c\_ptr), value c\_sst, type(c\_ptr), value c\_u10m, type(c\_ptr), value c\_v10m, type(c\_ptr), value c\_ustar, type(c\_ptr), value c\_delp, real(c\_double), dimension(n\_species), intent(in) species\_density, real(c\_double), dimension(n\_species), intent(in) species\_radius, real(c\_double), dimension(n\_species), intent(in) species\_lower\_radius, real(c\_double), dimension(n\_species), intent(in) species\_upper\_radius, logical(c\_bool), dimension(n\_species), intent(in) is\_gas\_arr, real(c\_double), dimension(n\_species), intent(in) species\_mw\_g, character(kind=c\_char), dimension(32, n\_species), intent(in) bin\_species\_names, character(kind=c\_char), dimension(32, n\_total\_species), intent(in) species\_names, type(c\_ptr), value c\_conc, type(c\_ptr), value c\_tendency, type(c\_ptr), value c\_diag\_mass\_total, type(c\_ptr), value c\_diag\_num\_total, type(c\_ptr), value c\_diag\_mass\_bin, type(c\_ptr), value c\_diag\_num\_bin, integer(c\_int), dimension(max(n\_diag\_species,1)), intent(in) diagnostic\_species\_id, integer(c\_int), value n\_diag\_species) <br> |




























## Public Functions Documentation




### function run\_seasalt\_science\_bridge 

```Fortran
subroutine seasaltsciencebridge_mod::run_seasalt_science_bridge (
    integer(c_int), value n_cols,
    integer(c_int), value n_levels,
    integer(c_int), value n_species,
    integer(c_int), value n_total_species,
    real(c_double), value dt,
    character(kind=c_char), dimension(*), intent(in) active_scheme,
    integer(c_int), value diagnostics,
    real(c_double), value gong97_scale_factor,
    integer(c_int), value gong97_weibull_flag,
    real(c_double), value gong03_scale_factor,
    integer(c_int), value gong03_weibull_flag,
    real(c_double), value geos12_scale_factor,
    integer(c_int), value geos12_weibull_flag,
    type(c_ptr), value c_frocean,
    type(c_ptr), value c_frseaice,
    type(c_ptr), value c_lat,
    type(c_ptr), value c_lon,
    type(c_ptr), value c_sst,
    type(c_ptr), value c_u10m,
    type(c_ptr), value c_v10m,
    type(c_ptr), value c_ustar,
    type(c_ptr), value c_delp,
    real(c_double), dimension(n_species), intent(in) species_density,
    real(c_double), dimension(n_species), intent(in) species_radius,
    real(c_double), dimension(n_species), intent(in) species_lower_radius,
    real(c_double), dimension(n_species), intent(in) species_upper_radius,
    logical(c_bool), dimension(n_species), intent(in) is_gas_arr,
    real(c_double), dimension(n_species), intent(in) species_mw_g,
    character(kind=c_char), dimension(32, n_species), intent(in) bin_species_names,
    character(kind=c_char), dimension(32, n_total_species), intent(in) species_names,
    type(c_ptr), value c_conc,
    type(c_ptr), value c_tendency,
    type(c_ptr), value c_diag_mass_total,
    type(c_ptr), value c_diag_num_total,
    type(c_ptr), value c_diag_mass_bin,
    type(c_ptr), value c_diag_num_bin,
    integer(c_int), dimension(max(n_diag_species,1)), intent(in) diagnostic_species_id,
    integer(c_int), value n_diag_species
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/seasalt/SeaSaltScienceBridge.F90`

