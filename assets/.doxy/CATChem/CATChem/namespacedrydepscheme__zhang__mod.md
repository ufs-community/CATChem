

# Namespace drydepscheme\_zhang\_mod



[**Namespace List**](namespaces.md) **>** [**drydepscheme\_zhang\_mod**](namespacedrydepscheme__zhang__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_zhang**](#function-compute_zhang) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**drydepschemezhangconfig**](namespacedrydepcommon__mod.md#none-drydepschemezhangconfig)), intent(in) params, real(fp), dimension(num\_layers), intent(in) bxheight, real(fp), dimension(:), intent(in) frlanduse, integer, dimension(:), intent(in) iland, logical, intent(in) isice, logical, intent(in) issnow, character(len=255), intent(in) lucname, real(fp), intent(in) obk, real(fp), intent(in) ps, real(fp), dimension(num\_layers), intent(in) rh, real(fp), intent(in) ts, real(fp), intent(in) tstep, real(fp), intent(in) u10m, real(fp), intent(in) ustar, real(fp), intent(in) v10m, real(fp), intent(in) z0, real(fp), dimension(num\_species), intent(in) species\_mw\_g, real(fp), dimension(num\_species), intent(in) species\_radius, real(fp), dimension(num\_species), intent(in) species\_density, character(len=32), dimension(num\_species), intent(in) species\_short\_name, real(fp), dimension(num\_species), intent(in) species\_dd\_hstar, real(fp), dimension(num\_species), intent(in) species\_dd\_dvzaersnow, real(fp), dimension(num\_species), intent(in) species\_dd\_dvzminval\_snow, real(fp), dimension(num\_species), intent(in) species\_dd\_dvzminval\_land, real(fp), dimension(num\_species), intent(in) species\_lower\_radius, real(fp), dimension(num\_species), intent(in) species\_upper\_radius, logical, dimension(num\_species), intent(in) species\_is\_dust, logical, dimension(num\_species), intent(in) species\_is\_seasalt, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, logical, dimension(num\_species), intent(in) is\_gas, real(fp), dimension(:), intent(inout), optional drydep\_con\_per\_species, real(fp), dimension(:), intent(inout), optional drydep\_velocity\_per\_species, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for zhang scheme._  |




























## Public Functions Documentation




### function compute\_zhang 

_Pure science computation for zhang scheme._ 
```Fortran
subroutine, public drydepscheme_zhang_mod::compute_zhang (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( drydepschemezhangconfig ), intent(in) params,
    real(fp), dimension(num_layers), intent(in) bxheight,
    real(fp), dimension(:), intent(in) frlanduse,
    integer, dimension(:), intent(in) iland,
    logical, intent(in) isice,
    logical, intent(in) issnow,
    character(len=255), intent(in) lucname,
    real(fp), intent(in) obk,
    real(fp), intent(in) ps,
    real(fp), dimension(num_layers), intent(in) rh,
    real(fp), intent(in) ts,
    real(fp), intent(in) tstep,
    real(fp), intent(in) u10m,
    real(fp), intent(in) ustar,
    real(fp), intent(in) v10m,
    real(fp), intent(in) z0,
    real(fp), dimension(num_species), intent(in) species_mw_g,
    real(fp), dimension(num_species), intent(in) species_radius,
    real(fp), dimension(num_species), intent(in) species_density,
    character(len=32), dimension(num_species), intent(in) species_short_name,
    real(fp), dimension(num_species), intent(in) species_dd_hstar,
    real(fp), dimension(num_species), intent(in) species_dd_dvzaersnow,
    real(fp), dimension(num_species), intent(in) species_dd_dvzminval_snow,
    real(fp), dimension(num_species), intent(in) species_dd_dvzminval_land,
    real(fp), dimension(num_species), intent(in) species_lower_radius,
    real(fp), dimension(num_species), intent(in) species_upper_radius,
    logical, dimension(num_species), intent(in) species_is_dust,
    logical, dimension(num_species), intent(in) species_is_seasalt,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    logical, dimension(num_species), intent(in) is_gas,
    real(fp), dimension(:), intent(inout), optional drydep_con_per_species,
    real(fp), dimension(:), intent(inout), optional drydep_velocity_per_species,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing Zhang et al. [2001] scheme with Emerson et al. [2020] updates. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `bxheight` BXHEIGHT field [appropriate units] 
* `frlanduse` FRLANDUSE field [appropriate units] 
* `iland` ILAND field [appropriate units] 
* `isice` IsIce field [appropriate units] 
* `issnow` IsSnow field [appropriate units] 
* `lucname` LUCNAME field [appropriate units] 
* `obk` OBK field [appropriate units] 
* `ps` PS field [appropriate units] 
* `rh` RH field [appropriate units] 
* `ts` TS field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `u10m` U10M field [appropriate units] 
* `ustar` USTAR field [appropriate units] 
* `v10m` V10M field [appropriate units] 
* `z0` Z0 field [appropriate units] 
* `species_mw_g` Species mw\_g property 
* `species_radius` Species radius property 
* `species_density` Species density property 
* `species_short_name` Species short\_name property 
* `species_dd_hstar` Species dd\_hstar property 
* `species_dd_DvzAerSnow` Species dd\_DvzAerSnow property 
* `species_dd_DvzMinVal_snow` Species dd\_DvzMinVal\_snow property 
* `species_dd_DvzMinVal_land` Species dd\_DvzMinVal\_land property 
* `species_lower_radius` Species lower\_radius property 
* `species_upper_radius` Species upper\_radius property 
* `species_is_dust` Species is\_dust property 
* `species_is_seasalt` Species is\_seasalt property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `drydep_con_per_species` Dry deposition concentration per species [ug/kg or ppm] (num\_species) 
* `drydep_velocity_per_species` Dry deposition velocity [m/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/schemes/DryDepScheme_ZHANG_Mod.F90`

