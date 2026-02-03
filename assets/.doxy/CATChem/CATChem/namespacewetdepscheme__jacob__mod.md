

# Namespace wetdepscheme\_jacob\_mod



[**Namespace List**](namespaces.md) **>** [**wetdepscheme\_jacob\_mod**](namespacewetdepscheme__jacob__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_jacob**](#function-compute_jacob) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**wetdepschemejacobconfig**](namespacewetdepcommon__mod.md#none-wetdepschemejacobconfig)), intent(in) params, real(fp), dimension(num\_layers), intent(in) airden\_dry, real(fp), dimension(num\_layers), intent(in) mairden, real(fp), dimension(num\_layers+1), intent(in) pedge, real(fp), dimension(num\_layers+1), intent(in) pfilsan, real(fp), dimension(num\_layers+1), intent(in) pfllsan, real(fp), dimension(num\_layers), intent(in) reevapls, real(fp), dimension(num\_layers), intent(in) t, real(fp), intent(in) tstep, logical, dimension(:), intent(in) species\_is\_aerosol, character(len=32), dimension(:), intent(in) species\_short\_name, real(fp), dimension(:), intent(in) species\_henry\_cr, real(fp), dimension(:), intent(in) species\_henry\_k0, real(fp), dimension(:), intent(in) species\_henry\_pka, real(fp), dimension(:), intent(in) species\_wd\_retfactor, logical, dimension(:), intent(in) species\_wd\_liqandgas, real(fp), dimension(:), intent(in) species\_wd\_convfaci2g, real(fp), dimension(:,:), intent(in) species\_wd\_rainouteff, real(fp), dimension(:), intent(in) species\_radius, real(fp), dimension(:), intent(in) species\_mw\_g, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), dimension(:,:), intent(inout), optional wetdep\_mass\_per\_species\_per\_level, real(fp), dimension(:,:), intent(inout), optional wetdep\_flux\_per\_species\_per\_level, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for jacob scheme._  |




























## Public Functions Documentation




### function compute\_jacob 

_Pure science computation for jacob scheme._ 
```Fortran
subroutine, public wetdepscheme_jacob_mod::compute_jacob (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( wetdepschemejacobconfig ), intent(in) params,
    real(fp), dimension(num_layers), intent(in) airden_dry,
    real(fp), dimension(num_layers), intent(in) mairden,
    real(fp), dimension(num_layers+1), intent(in) pedge,
    real(fp), dimension(num_layers+1), intent(in) pfilsan,
    real(fp), dimension(num_layers+1), intent(in) pfllsan,
    real(fp), dimension(num_layers), intent(in) reevapls,
    real(fp), dimension(num_layers), intent(in) t,
    real(fp), intent(in) tstep,
    logical, dimension(:), intent(in) species_is_aerosol,
    character(len=32), dimension(:), intent(in) species_short_name,
    real(fp), dimension(:), intent(in) species_henry_cr,
    real(fp), dimension(:), intent(in) species_henry_k0,
    real(fp), dimension(:), intent(in) species_henry_pka,
    real(fp), dimension(:), intent(in) species_wd_retfactor,
    logical, dimension(:), intent(in) species_wd_liqandgas,
    real(fp), dimension(:), intent(in) species_wd_convfaci2g,
    real(fp), dimension(:,:), intent(in) species_wd_rainouteff,
    real(fp), dimension(:), intent(in) species_radius,
    real(fp), dimension(:), intent(in) species_mw_g,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), dimension(:,:), intent(inout), optional wetdep_mass_per_species_per_level,
    real(fp), dimension(:,:), intent(inout), optional wetdep_flux_per_species_per_level,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing Jacob et al. [2000] wet deposition scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `airden_dry` AIRDEN\_DRY field [appropriate units] 
* `mairden` MAIRDEN field [appropriate units] 
* `pedge` PEDGE field [appropriate units] 
* `pfilsan` PFILSAN field [appropriate units] 
* `pfllsan` PFLLSAN field [appropriate units] 
* `reevapls` REEVAPLS field [appropriate units] 
* `t` T field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `species_is_aerosol` Species is\_aerosol property 
* `species_short_name` Species short\_name property 
* `species_henry_cr` Species henry\_cr property 
* `species_henry_k0` Species henry\_k0 property 
* `species_henry_pKa` Species henry\_pKa property 
* `species_wd_retfactor` Species wd\_retfactor property 
* `species_wd_LiqAndGas` Species wd\_LiqAndGas property 
* `species_wd_convfacI2G` Species wd\_convfacI2G property 
* `species_wd_rainouteff` Species wd\_rainouteff property 
* `species_radius` Species radius property 
* `species_mw_g` Species mw\_g property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `wetdep_mass_per_species_per_level` Wet deposition mass loss per species per level [kg/m2] (num\_species) 
* `wetdep_flux_per_species_per_level` Wet deposition flux per species per level [kg/m2/s] (num\_species) 
* `wetdep_mass_per_species_per_level` Wet deposition mass loss per species per level [kg/m2] (num\_species) 
* `wetdep_flux_per_species_per_level` Wet deposition flux per species per level [kg/m2/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/wetdep/schemes/WetDepScheme_JACOB_Mod.F90`

