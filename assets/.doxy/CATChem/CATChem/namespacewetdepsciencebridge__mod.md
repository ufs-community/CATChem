

# Namespace wetdepsciencebridge\_mod



[**Namespace List**](namespaces.md) **>** [**wetdepsciencebridge\_mod**](namespacewetdepsciencebridge__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**run\_wetdep\_gocart\_science\_bridge**](#function-run_wetdep_gocart_science_bridge) (integer(c\_int), value n\_cols, integer(c\_int), value n\_levels, integer(c\_int), value n\_species, real(c\_double), value dt, integer(c\_int), value diagnostics, real(c\_double), value gocart\_scale\_factor, real(c\_double), value gocart\_washout\_tuning, real(c\_double), value gocart\_radius\_threshold, type(c\_ptr), value c\_mairden, type(c\_ptr), value c\_pedge, type(c\_ptr), value c\_pfilsan, type(c\_ptr), value c\_pfllsan, type(c\_ptr), value c\_t\_air, type(c\_ptr), value c\_preccon, type(c\_ptr), value c\_preclsc, logical(c\_bool), dimension(n\_species), intent(in) species\_is\_aerosol, logical(c\_bool), dimension(n\_species), intent(in) species\_is\_wetdep, real(c\_double), dimension(n\_species, 3), intent(in) species\_wd\_rainouteff, real(c\_double), dimension(n\_species), intent(in) species\_radius, real(c\_double), dimension(n\_species), intent(in) species\_mw\_g, character(kind=c\_char), dimension(32, n\_species), intent(in) species\_names, type(c\_ptr), value c\_conc, type(c\_ptr), value c\_tendency, type(c\_ptr), value c\_diag\_mass, type(c\_ptr), value c\_diag\_flux, integer(c\_int), dimension(n\_diag\_species), intent(in) diagnostic\_species\_id, integer(c\_int), value n\_diag\_species) <br>_GOCART2G wet-removal entry point (scheme 'gocart')._  |
|  subroutine | [**run\_wetdep\_science\_bridge**](#function-run_wetdep_science_bridge) (integer(c\_int), value n\_cols, integer(c\_int), value n\_levels, integer(c\_int), value n\_species, real(c\_double), value dt, integer(c\_int), value diagnostics, real(c\_double), value jacob\_scale\_factor, real(c\_double), value jacob\_radius\_threshold, integer(c\_int), value jacob\_so4\_gocart\_resusp, real(c\_double), value jacob\_so4\_washout\_eff, type(c\_ptr), value c\_airden\_dry, type(c\_ptr), value c\_mairden, type(c\_ptr), value c\_pedge, type(c\_ptr), value c\_pfilsan, type(c\_ptr), value c\_pfllsan, type(c\_ptr), value c\_reevapls, type(c\_ptr), value c\_t\_air, logical(c\_bool), dimension(n\_species), intent(in) species\_is\_aerosol, logical(c\_bool), dimension(n\_species), intent(in) species\_is\_wetdep, real(c\_double), dimension(n\_species), intent(in) species\_henry\_cr, real(c\_double), dimension(n\_species), intent(in) species\_henry\_k0, real(c\_double), dimension(n\_species), intent(in) species\_henry\_pka, real(c\_double), dimension(n\_species), intent(in) species\_wd\_retfactor, logical(c\_bool), dimension(n\_species), intent(in) species\_wd\_liqandgas, real(c\_double), dimension(n\_species), intent(in) species\_wd\_convfaci2g, real(c\_double), dimension(n\_species, 3), intent(in) species\_wd\_rainouteff, real(c\_double), dimension(n\_species), intent(in) species\_wd\_reevap\_frac, real(c\_double), dimension(n\_species), intent(in) species\_radius, real(c\_double), dimension(n\_species), intent(in) species\_mw\_g, character(kind=c\_char), dimension(32, n\_species), intent(in) species\_names, type(c\_ptr), value c\_conc, type(c\_ptr), value c\_tendency, type(c\_ptr), value c\_diag\_mass, type(c\_ptr), value c\_diag\_flux, integer(c\_int), dimension(n\_diag\_species), intent(in) diagnostic\_species\_id, integer(c\_int), value n\_diag\_species) <br> |




























## Public Functions Documentation




### function run\_wetdep\_gocart\_science\_bridge 

_GOCART2G wet-removal entry point (scheme 'gocart')._ 
```Fortran
subroutine wetdepsciencebridge_mod::run_wetdep_gocart_science_bridge (
    integer(c_int), value n_cols,
    integer(c_int), value n_levels,
    integer(c_int), value n_species,
    real(c_double), value dt,
    integer(c_int), value diagnostics,
    real(c_double), value gocart_scale_factor,
    real(c_double), value gocart_washout_tuning,
    real(c_double), value gocart_radius_threshold,
    type(c_ptr), value c_mairden,
    type(c_ptr), value c_pedge,
    type(c_ptr), value c_pfilsan,
    type(c_ptr), value c_pfllsan,
    type(c_ptr), value c_t_air,
    type(c_ptr), value c_preccon,
    type(c_ptr), value c_preclsc,
    logical(c_bool), dimension(n_species), intent(in) species_is_aerosol,
    logical(c_bool), dimension(n_species), intent(in) species_is_wetdep,
    real(c_double), dimension(n_species, 3), intent(in) species_wd_rainouteff,
    real(c_double), dimension(n_species), intent(in) species_radius,
    real(c_double), dimension(n_species), intent(in) species_mw_g,
    character(kind=c_char), dimension(32, n_species), intent(in) species_names,
    type(c_ptr), value c_conc,
    type(c_ptr), value c_tendency,
    type(c_ptr), value c_diag_mass,
    type(c_ptr), value c_diag_flux,
    integer(c_int), dimension(n_diag_species), intent(in) diagnostic_species_id,
    integer(c_int), value n_diag_species
) 
```



Separate bind(C) entry from the Jacob bridge so the protected Jacob surface (run\_wetdep\_science\_bridge) stays byte-identical. Shares this module's species-array marshalling idiom; differs in the scheme configuration record (WetDepSchemeGOCARTConfig) and the two additional surface precipitation inputs PRECCON/PRECLSC required by the GOCART kernel (instantaneous rain amounts [m], divided by tstep inside the kernel to rates [m/s]). 


        

<hr>



### function run\_wetdep\_science\_bridge 

```Fortran
subroutine wetdepsciencebridge_mod::run_wetdep_science_bridge (
    integer(c_int), value n_cols,
    integer(c_int), value n_levels,
    integer(c_int), value n_species,
    real(c_double), value dt,
    integer(c_int), value diagnostics,
    real(c_double), value jacob_scale_factor,
    real(c_double), value jacob_radius_threshold,
    integer(c_int), value jacob_so4_gocart_resusp,
    real(c_double), value jacob_so4_washout_eff,
    type(c_ptr), value c_airden_dry,
    type(c_ptr), value c_mairden,
    type(c_ptr), value c_pedge,
    type(c_ptr), value c_pfilsan,
    type(c_ptr), value c_pfllsan,
    type(c_ptr), value c_reevapls,
    type(c_ptr), value c_t_air,
    logical(c_bool), dimension(n_species), intent(in) species_is_aerosol,
    logical(c_bool), dimension(n_species), intent(in) species_is_wetdep,
    real(c_double), dimension(n_species), intent(in) species_henry_cr,
    real(c_double), dimension(n_species), intent(in) species_henry_k0,
    real(c_double), dimension(n_species), intent(in) species_henry_pka,
    real(c_double), dimension(n_species), intent(in) species_wd_retfactor,
    logical(c_bool), dimension(n_species), intent(in) species_wd_liqandgas,
    real(c_double), dimension(n_species), intent(in) species_wd_convfaci2g,
    real(c_double), dimension(n_species, 3), intent(in) species_wd_rainouteff,
    real(c_double), dimension(n_species), intent(in) species_wd_reevap_frac,
    real(c_double), dimension(n_species), intent(in) species_radius,
    real(c_double), dimension(n_species), intent(in) species_mw_g,
    character(kind=c_char), dimension(32, n_species), intent(in) species_names,
    type(c_ptr), value c_conc,
    type(c_ptr), value c_tendency,
    type(c_ptr), value c_diag_mass,
    type(c_ptr), value c_diag_flux,
    integer(c_int), dimension(n_diag_species), intent(in) diagnostic_species_id,
    integer(c_int), value n_diag_species
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/wetdep/WetDepScienceBridge.F90`

