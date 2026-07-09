

# Namespace dustscheme\_ginoux\_mod



[**Namespace List**](namespaces.md) **>** [**dustscheme\_ginoux\_mod**](namespacedustscheme__ginoux__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  pure subroutine, public | [**compute\_ginoux**](#function-compute_ginoux) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**dustschemeginouxconfig**](namespacedustcommon__mod.md#none-dustschemeginouxconfig)), intent(in) params, real(fp), intent(in) g0, real(fp), dimension(num\_layers), intent(in) airden, real(fp), intent(in) frlake, real(fp), intent(in) frsno, real(fp), intent(in) gwettop, integer, intent(in) lwi, real(fp), intent(in) ssm, real(fp), intent(in) tskin, real(fp), intent(in) u10m, real(fp), intent(in) v10m, real(fp), dimension(:), intent(in) species\_density, real(fp), dimension(:), intent(in) species\_radius, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), intent(inout), optional dust\_emission\_total, real(fp), dimension(:), intent(inout), optional dust\_emission\_per\_bin, real(fp), dimension(:), intent(inout), optional utar\_threshold\_per\_bin, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for ginoux scheme._  |




























## Public Functions Documentation




### function compute\_ginoux 

_Pure science computation for ginoux scheme._ 
```Fortran
pure subroutine, public dustscheme_ginoux_mod::compute_ginoux (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( dustschemeginouxconfig ), intent(in) params,
    real(fp), intent(in) g0,
    real(fp), dimension(num_layers), intent(in) airden,
    real(fp), intent(in) frlake,
    real(fp), intent(in) frsno,
    real(fp), intent(in) gwettop,
    integer, intent(in) lwi,
    real(fp), intent(in) ssm,
    real(fp), intent(in) tskin,
    real(fp), intent(in) u10m,
    real(fp), intent(in) v10m,
    real(fp), dimension(:), intent(in) species_density,
    real(fp), dimension(:), intent(in) species_radius,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), intent(inout), optional dust_emission_total,
    real(fp), dimension(:), intent(inout), optional dust_emission_per_bin,
    real(fp), dimension(:), intent(inout), optional utar_threshold_per_bin,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing Ginoux dust emission scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `g0` Required constant from Constants module 
* `airden` AIRDEN field [appropriate units] 
* `frlake` FRLAKE field [appropriate units] 
* `frsno` FRSNO field [appropriate units] 
* `gwettop` GWETTOP field [appropriate units] 
* `lwi` LWI field [appropriate units] 
* `ssm` SSM field [appropriate units] 
* `tskin` TSKIN field [appropriate units] 
* `u10m` U10M field [appropriate units] 
* `v10m` V10M field [appropriate units] 
* `species_density` Species density property 
* `species_radius` Species radius property 
* `species_conc` Species concentrations [ppm or ug/kg] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `dust_emission_total` Total dust emissions for all bins [kg/m2/s] 
* `dust_emission_per_bin` Dust emission flux per bin [kg/m2/s] (num\_species) 
* `utar_threshold` Threshold friction velocity to have dust emission [m/s] 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/dust/schemes/DustScheme_GINOUX_Mod.F90`

