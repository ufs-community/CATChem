

# Namespace seasaltscheme\_geos12\_mod



[**Namespace List**](namespaces.md) **>** [**seasaltscheme\_geos12\_mod**](namespaceseasaltscheme__geos12__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  pure subroutine, public | [**compute\_geos12**](#function-compute_geos12) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**seasaltschemegeos12config**](namespaceseasaltcommon__mod.md#none-seasaltschemegeos12config)), intent(in) params, real(fp), intent(in) frocean, real(fp), intent(in) frseaice, real(fp), intent(in) sst, real(fp), intent(in) ustar, real(fp), dimension(num\_species), intent(in) species\_density, real(fp), dimension(num\_species), intent(in) species\_radius, real(fp), dimension(num\_species), intent(in) species\_lower\_radius, real(fp), dimension(num\_species), intent(in) species\_upper\_radius, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), intent(inout), optional seasalt\_mass\_emission\_total, real(fp), intent(inout), optional seasalt\_number\_emission\_total, real(fp), dimension(:), intent(inout), optional seasalt\_mass\_emission\_per\_bin, real(fp), dimension(:), intent(inout), optional seasalt\_number\_emission\_per\_bin, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for geos12 scheme._  |




























## Public Functions Documentation




### function compute\_geos12 

_Pure science computation for geos12 scheme._ 
```Fortran
pure subroutine, public seasaltscheme_geos12_mod::compute_geos12 (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( seasaltschemegeos12config ), intent(in) params,
    real(fp), intent(in) frocean,
    real(fp), intent(in) frseaice,
    real(fp), intent(in) sst,
    real(fp), intent(in) ustar,
    real(fp), dimension(num_species), intent(in) species_density,
    real(fp), dimension(num_species), intent(in) species_radius,
    real(fp), dimension(num_species), intent(in) species_lower_radius,
    real(fp), dimension(num_species), intent(in) species_upper_radius,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), intent(inout), optional seasalt_mass_emission_total,
    real(fp), intent(inout), optional seasalt_number_emission_total,
    real(fp), dimension(:), intent(inout), optional seasalt_mass_emission_per_bin,
    real(fp), dimension(:), intent(inout), optional seasalt_number_emission_per_bin,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing GEOS-Chem 2012 sea salt emission scheme with observational constraints. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `frocean` FROCEAN field [appropriate units] 
* `frseaice` FRSEAICE field [appropriate units] 
* `sst` SST field [appropriate units] 
* `ustar` USTAR field [appropriate units] 
* `species_density` Species density property 
* `species_radius` Species radius property 
* `species_lower_radius` Species lower\_radius property 
* `species_upper_radius` Species upper\_radius property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `seasalt_mass_emission_total` Total mass emission diagnostic [ug/m2/s] 
* `seasalt_number_emission_total` Total number emission diagnostic [#/m2/s] 
* `seasalt_mass_emission_per_bin` Mass emission per bin diagnostic [kg/m2/s] (num\_species) 
* `seasalt_number_emission_per_bin` Number emission per bin diagnostic [#/m2/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/seasalt/schemes/SeaSaltScheme_GEOS12_Mod.F90`

