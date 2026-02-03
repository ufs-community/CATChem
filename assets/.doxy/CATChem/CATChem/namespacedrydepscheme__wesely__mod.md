

# Namespace drydepscheme\_wesely\_mod



[**Namespace List**](namespaces.md) **>** [**drydepscheme\_wesely\_mod**](namespacedrydepscheme__wesely__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_wesely**](#function-compute_wesely) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**drydepschemeweselyconfig**](namespacedrydepcommon__mod.md#none-drydepschemeweselyconfig)), intent(in) params, real(fp), dimension(num\_layers), intent(in) bxheight, real(fp), intent(in) cldfrc, real(fp), dimension(:), intent(in) frlai, real(fp), dimension(:), intent(in) frlanduse, integer, dimension(:), intent(in) iland, logical, intent(in) isice, logical, intent(in) island, logical, intent(in) issnow, real(fp), intent(in) lat, real(fp), intent(in) lon, character(len=255), intent(in) lucname, real(fp), intent(in) obk, real(fp), intent(in) ps, real(fp), intent(in) salinity, real(fp), intent(in) suncosmid, real(fp), intent(in) swgdn, real(fp), intent(in) ts, real(fp), intent(in) tskin, real(fp), intent(in) tstep, real(fp), intent(in) ustar, real(fp), intent(in) z0, real(fp), dimension(num\_species), intent(in) species\_mw\_g, real(fp), dimension(num\_species), intent(in) species\_dd\_f0, character(len=32), dimension(num\_species), intent(in) species\_short\_name, real(fp), dimension(num\_species), intent(in) species\_dd\_hstar, real(fp), dimension(num\_species), intent(in) species\_dd\_dvzaersnow, real(fp), dimension(num\_species), intent(in) species\_dd\_dvzminval\_snow, real(fp), dimension(num\_species), intent(in) species\_dd\_dvzminval\_land, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, logical, dimension(num\_species), intent(in) is\_gas, real(fp), dimension(:), intent(inout), optional drydep\_con\_per\_species, real(fp), dimension(:), intent(inout), optional drydep\_velocity\_per\_species, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for wesely scheme._  |




























## Public Functions Documentation




### function compute\_wesely 

_Pure science computation for wesely scheme._ 
```Fortran
subroutine, public drydepscheme_wesely_mod::compute_wesely (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( drydepschemeweselyconfig ), intent(in) params,
    real(fp), dimension(num_layers), intent(in) bxheight,
    real(fp), intent(in) cldfrc,
    real(fp), dimension(:), intent(in) frlai,
    real(fp), dimension(:), intent(in) frlanduse,
    integer, dimension(:), intent(in) iland,
    logical, intent(in) isice,
    logical, intent(in) island,
    logical, intent(in) issnow,
    real(fp), intent(in) lat,
    real(fp), intent(in) lon,
    character(len=255), intent(in) lucname,
    real(fp), intent(in) obk,
    real(fp), intent(in) ps,
    real(fp), intent(in) salinity,
    real(fp), intent(in) suncosmid,
    real(fp), intent(in) swgdn,
    real(fp), intent(in) ts,
    real(fp), intent(in) tskin,
    real(fp), intent(in) tstep,
    real(fp), intent(in) ustar,
    real(fp), intent(in) z0,
    real(fp), dimension(num_species), intent(in) species_mw_g,
    real(fp), dimension(num_species), intent(in) species_dd_f0,
    character(len=32), dimension(num_species), intent(in) species_short_name,
    real(fp), dimension(num_species), intent(in) species_dd_hstar,
    real(fp), dimension(num_species), intent(in) species_dd_dvzaersnow,
    real(fp), dimension(num_species), intent(in) species_dd_dvzminval_snow,
    real(fp), dimension(num_species), intent(in) species_dd_dvzminval_land,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    logical, dimension(num_species), intent(in) is_gas,
    real(fp), dimension(:), intent(inout), optional drydep_con_per_species,
    real(fp), dimension(:), intent(inout), optional drydep_velocity_per_species,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing Wesely 1989 gas dry deposition scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `bxheight` BXHEIGHT field [appropriate units] 
* `cldfrc` CLDFRC field [appropriate units] 
* `frlai` FRLAI field [appropriate units] 
* `frlanduse` FRLANDUSE field [appropriate units] 
* `iland` ILAND field [appropriate units] 
* `isice` IsIce field [appropriate units] 
* `island` IsLand field [appropriate units] 
* `issnow` IsSnow field [appropriate units] 
* `lat` LAT field [appropriate units] 
* `lon` LON field [appropriate units] 
* `lucname` LUCNAME field [appropriate units] 
* `obk` OBK field [appropriate units] 
* `ps` PS field [appropriate units] 
* `salinity` SALINITY field [appropriate units] 
* `suncosmid` SUNCOSmid field [appropriate units] 
* `swgdn` SWGDN field [appropriate units] 
* `ts` TS field [appropriate units] 
* `tskin` TSKIN field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `ustar` USTAR field [appropriate units] 
* `z0` Z0 field [appropriate units] 
* `species_mw_g` Species mw\_g property 
* `species_dd_f0` Species dd\_f0 property 
* `species_short_name` Species short\_name property 
* `species_dd_hstar` Species dd\_hstar property 
* `species_dd_DvzAerSnow` Species dd\_DvzAerSnow property 
* `species_dd_DvzMinVal_snow` Species dd\_DvzMinVal\_snow property 
* `species_dd_DvzMinVal_land` Species dd\_DvzMinVal\_land property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `drydep_con_per_species` Dry deposition concentration per species [ug/kg or ppm] (num\_species) 
* `drydep_velocity_per_species` Dry deposition velocity [m/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/schemes/DryDepScheme_WESELY_Mod.F90`

