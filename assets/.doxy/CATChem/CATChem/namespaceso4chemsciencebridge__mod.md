

# Namespace so4chemsciencebridge\_mod



[**Namespace List**](namespaces.md) **>** [**so4chemsciencebridge\_mod**](namespaceso4chemsciencebridge__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**run\_so4chem\_science\_bridge**](#function-run_so4chem_science_bridge) (integer(c\_int), value n\_cols, integer(c\_int), value n\_levels, integer(c\_int), value n\_species, real(c\_double), value dt, integer(c\_int), value diagnostics, integer(c\_int), value gocart\_update\_so2, integer(c\_int), value gocart\_do\_drydep, integer(c\_int), value year, integer(c\_int), value month, integer(c\_int), value day, integer(c\_int), value hour, integer(c\_int), value minute, integer(c\_int), value second, type(c\_ptr), value c\_airden, type(c\_ptr), value c\_cldf, type(c\_ptr), value c\_delp, type(c\_ptr), value c\_pmid, type(c\_ptr), value c\_t\_air, type(c\_ptr), value c\_z\_edges, type(c\_ptr), value c\_hflux, type(c\_ptr), value c\_lat, type(c\_ptr), value c\_lon, type(c\_ptr), value c\_lwi, type(c\_ptr), value c\_pblh, type(c\_ptr), value c\_u10m, type(c\_ptr), value c\_ustar, type(c\_ptr), value c\_v10m, type(c\_ptr), value c\_z0h, real(c\_double), dimension(n\_species), intent(in) species\_mw\_g, character(kind=c\_char), dimension(32, n\_species), intent(in) species\_names, type(c\_ptr), value c\_conc, type(c\_ptr), value c\_tendency, type(c\_ptr), value c\_firsttime, type(c\_ptr), value c\_nymd\_last, type(c\_ptr), value c\_nhms\_last\_recycle, type(c\_ptr), value c\_xh2o2\_init, type(c\_ptr), value c\_pso4\_g\_so2, type(c\_ptr), value c\_pso4\_aq\_so2, type(c\_ptr), value c\_pso2\_dms, type(c\_ptr), value c\_dms\_flux, type(c\_ptr), value c\_diag\_prod\_rate, integer(c\_int), dimension(max(n\_diag\_species,1)), intent(in) diagnostic\_species\_id, integer(c\_int), value n\_diag\_species) <br> |




























## Public Functions Documentation




### function run\_so4chem\_science\_bridge 

```Fortran
subroutine so4chemsciencebridge_mod::run_so4chem_science_bridge (
    integer(c_int), value n_cols,
    integer(c_int), value n_levels,
    integer(c_int), value n_species,
    real(c_double), value dt,
    integer(c_int), value diagnostics,
    integer(c_int), value gocart_update_so2,
    integer(c_int), value gocart_do_drydep,
    integer(c_int), value year,
    integer(c_int), value month,
    integer(c_int), value day,
    integer(c_int), value hour,
    integer(c_int), value minute,
    integer(c_int), value second,
    type(c_ptr), value c_airden,
    type(c_ptr), value c_cldf,
    type(c_ptr), value c_delp,
    type(c_ptr), value c_pmid,
    type(c_ptr), value c_t_air,
    type(c_ptr), value c_z_edges,
    type(c_ptr), value c_hflux,
    type(c_ptr), value c_lat,
    type(c_ptr), value c_lon,
    type(c_ptr), value c_lwi,
    type(c_ptr), value c_pblh,
    type(c_ptr), value c_u10m,
    type(c_ptr), value c_ustar,
    type(c_ptr), value c_v10m,
    type(c_ptr), value c_z0h,
    real(c_double), dimension(n_species), intent(in) species_mw_g,
    character(kind=c_char), dimension(32, n_species), intent(in) species_names,
    type(c_ptr), value c_conc,
    type(c_ptr), value c_tendency,
    type(c_ptr), value c_firsttime,
    type(c_ptr), value c_nymd_last,
    type(c_ptr), value c_nhms_last_recycle,
    type(c_ptr), value c_xh2o2_init,
    type(c_ptr), value c_pso4_g_so2,
    type(c_ptr), value c_pso4_aq_so2,
    type(c_ptr), value c_pso2_dms,
    type(c_ptr), value c_dms_flux,
    type(c_ptr), value c_diag_prod_rate,
    integer(c_int), dimension(max(n_diag_species,1)), intent(in) diagnostic_species_id,
    integer(c_int), value n_diag_species
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/so4chem/SO4chemScienceBridge.F90`

