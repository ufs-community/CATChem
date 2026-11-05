

# Namespace so4chemscheme\_gocart\_mod



[**Namespace List**](namespaces.md) **>** [**so4chemscheme\_gocart\_mod**](namespaceso4chemscheme__gocart__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_gocart**](#function-compute_gocart) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**so4chemschemegocartconfig**](namespaceso4chemcommon__mod.md#none-so4chemschemegocartconfig)), intent(in) params, real(fp), intent(in) g0, real(fp), intent(in) cpd, real(fp), intent(in) avo, real(fp), intent(in) von\_karman, real(fp), intent(in) airmw, real(fp), intent(in) pi, integer, intent(in) year, integer, intent(in) month, integer, intent(in) day, integer, intent(in) hour, integer, intent(in) minute, integer, intent(in) second, real(fp), dimension(num\_layers), intent(in) airden, real(fp), dimension(num\_layers), intent(in) cldf, real(fp), dimension(num\_layers), intent(in) delp, real(fp), intent(in) hflux, real(fp), intent(in) lat, real(fp), intent(in) lon, integer, intent(in) lwi, real(fp), intent(in) pblh, real(fp), dimension(num\_layers), intent(in) pmid, real(fp), dimension(num\_layers), intent(in) t, real(fp), intent(in) tstep, real(fp), intent(in) u10m, real(fp), intent(in) ustar, real(fp), intent(in) v10m, real(fp), dimension(num\_layers+1), intent(in) z, real(fp), intent(in) z0h, real(fp), dimension(:), intent(in) species\_mw\_g, character(len=32), dimension(:), intent(in) species\_short\_name, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, logical, intent(inout) firsttime, integer, intent(inout) nymd\_last, integer, intent(inout) nhms\_last\_recycle, real(fp), dimension(:), intent(inout), allocatable xh2o2\_init, real(fp), dimension(:,:), intent(inout), optional production\_rate\_per\_species\_per\_level, real(fp), dimension(:), intent(inout), optional pso4\_from\_gaseous\_so2\_per\_level, real(fp), dimension(:), intent(inout), optional pso4\_from\_aqueous\_so2\_per\_level, real(fp), intent(inout), optional dms\_emission\_flux, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for gocart scheme._  |




























## Public Functions Documentation




### function compute\_gocart 

_Pure science computation for gocart scheme._ 
```Fortran
subroutine, public so4chemscheme_gocart_mod::compute_gocart (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( so4chemschemegocartconfig ), intent(in) params,
    real(fp), intent(in) g0,
    real(fp), intent(in) cpd,
    real(fp), intent(in) avo,
    real(fp), intent(in) von_karman,
    real(fp), intent(in) airmw,
    real(fp), intent(in) pi,
    integer, intent(in) year,
    integer, intent(in) month,
    integer, intent(in) day,
    integer, intent(in) hour,
    integer, intent(in) minute,
    integer, intent(in) second,
    real(fp), dimension(num_layers), intent(in) airden,
    real(fp), dimension(num_layers), intent(in) cldf,
    real(fp), dimension(num_layers), intent(in) delp,
    real(fp), intent(in) hflux,
    real(fp), intent(in) lat,
    real(fp), intent(in) lon,
    integer, intent(in) lwi,
    real(fp), intent(in) pblh,
    real(fp), dimension(num_layers), intent(in) pmid,
    real(fp), dimension(num_layers), intent(in) t,
    real(fp), intent(in) tstep,
    real(fp), intent(in) u10m,
    real(fp), intent(in) ustar,
    real(fp), intent(in) v10m,
    real(fp), dimension(num_layers+1), intent(in) z,
    real(fp), intent(in) z0h,
    real(fp), dimension(:), intent(in) species_mw_g,
    character(len=32), dimension(:), intent(in) species_short_name,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    logical, intent(inout) firsttime,
    integer, intent(inout) nymd_last,
    integer, intent(inout) nhms_last_recycle,
    real(fp), dimension(:), intent(inout), allocatable xh2o2_init,
    real(fp), dimension(:,:), intent(inout), optional production_rate_per_species_per_level,
    real(fp), dimension(:), intent(inout), optional pso4_from_gaseous_so2_per_level,
    real(fp), dimension(:), intent(inout), optional pso4_from_aqueous_so2_per_level,
    real(fp), intent(inout), optional dms_emission_flux,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing GOCART SO2 to SO4 production scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `g0` Required constant from Constants module 
* `Cpd` Required constant from Constants module 
* `AVO` Required constant from Constants module 
* `VON_KARMAN` Required constant from Constants module 
* `AIRMW` Required constant from Constants module 
* `PI` Required constant from Constants module 
* `year` Time parameter from TimeState (year) 
* `month` Time parameter from TimeState (month) 
* `day` Time parameter from TimeState (day) 
* `hour` Time parameter from TimeState (hour) 
* `minute` Time parameter from TimeState (minute) 
* `second` Time parameter from TimeState (second) 
* `airden` AIRDEN field [appropriate units] 
* `cldf` CLDF field [appropriate units] 
* `delp` DELP field [appropriate units] 
* `hflux` HFLUX field [appropriate units] 
* `lat` LAT field [appropriate units] 
* `lon` LON field [appropriate units] 
* `lwi` LWI field [appropriate units] 
* `pblh` PBLH field [appropriate units] 
* `pmid` PMID field [appropriate units] 
* `t` T field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `u10m` U10M field [appropriate units] 
* `ustar` USTAR field [appropriate units] 
* `v10m` V10M field [appropriate units] 
* `z` Z field [appropriate units] 
* `z0h` Z0H field [appropriate units] 
* `species_mw_g` Species mw\_g property 
* `species_short_name` Species short\_name property 
* `species_conc` Species concentrations [ppm or ug/kg] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) Persistent state variables (per-column): 
* `firsttime` flag for first time step 
* `nymd_last` last day of H2O2 update 
* `nhms_last_recycle` last time step of H2O2 recycle 
* `xh2o2_init` H2O2 column initialization 
* `PSO4_from_SO2_per_level` total sulfate production rate from SO2 per level [kg/kg/s] (num\_layers) 
* `PSO4_from_gaseous_SO2_per_level` sulfate production rate from gaseous SO2 per level [kg/kg/s] (num\_layers) 
* `PSO4_from_aqueous_SO2_per_level` sulfate production rate from aqueous SO2 per level [kg/kg/s] (num\_layers) 
* `PSO2_from_DMS_per_level` SO2 production rate from DMS per level [kg/kg/s] (num\_layers) 
* `DMS_emission_flux` DMS emission flux at the surface [kg/m2/s] 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/so4chem/schemes/SO4chemScheme_GOCART_Mod.F90`

