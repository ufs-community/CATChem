

# Namespace dustscheme\_fengsha\_mod



[**Namespace List**](namespaces.md) **>** [**dustscheme\_fengsha\_mod**](namespacedustscheme__fengsha__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_fengsha**](#function-compute_fengsha) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**dustschemefengshaconfig**](namespacedustcommon__mod.md#none-dustschemefengshaconfig)), intent(in) params, real(fp), intent(in) g0, real(fp), dimension(num\_layers), intent(in) airden, real(fp), intent(in) clayfrac, real(fp), intent(in) frlake, real(fp), intent(in) frsno, real(fp), intent(in) gvf, real(fp), intent(in) lai, integer, intent(in) lwi, real(fp), intent(in) rdrag, real(fp), intent(in) sandfrac, real(fp), dimension(:), intent(in) soilm, real(fp), intent(in) ssm, real(fp), intent(in) tskin, real(fp), intent(in) ustar, real(fp), intent(in) ustar\_threshold, real(fp), intent(in) z0, real(fp), dimension(:), intent(in) species\_radius, real(fp), dimension(:), intent(in) species\_lower\_radius, real(fp), dimension(:), intent(in) species\_upper\_radius, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), intent(inout), optional dust\_emission\_total, real(fp), dimension(:), intent(inout), optional dust\_emission\_per\_bin, real(fp), intent(inout), optional dust\_horizontal\_flux, real(fp), intent(inout), optional dust\_moisture\_correction, real(fp), intent(inout), optional dust\_effective\_threshold, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for fengsha scheme._  |




























## Public Functions Documentation




### function compute\_fengsha 

_Pure science computation for fengsha scheme._ 
```Fortran
subroutine, public dustscheme_fengsha_mod::compute_fengsha (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( dustschemefengshaconfig ), intent(in) params,
    real(fp), intent(in) g0,
    real(fp), dimension(num_layers), intent(in) airden,
    real(fp), intent(in) clayfrac,
    real(fp), intent(in) frlake,
    real(fp), intent(in) frsno,
    real(fp), intent(in) gvf,
    real(fp), intent(in) lai,
    integer, intent(in) lwi,
    real(fp), intent(in) rdrag,
    real(fp), intent(in) sandfrac,
    real(fp), dimension(:), intent(in) soilm,
    real(fp), intent(in) ssm,
    real(fp), intent(in) tskin,
    real(fp), intent(in) ustar,
    real(fp), intent(in) ustar_threshold,
    real(fp), intent(in) z0,
    real(fp), dimension(:), intent(in) species_radius,
    real(fp), dimension(:), intent(in) species_lower_radius,
    real(fp), dimension(:), intent(in) species_upper_radius,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), intent(inout), optional dust_emission_total,
    real(fp), dimension(:), intent(inout), optional dust_emission_per_bin,
    real(fp), intent(inout), optional dust_horizontal_flux,
    real(fp), intent(inout), optional dust_moisture_correction,
    real(fp), intent(inout), optional dust_effective_threshold,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `g0` Required constant from Constants module 
* `airden` AIRDEN field [appropriate units] 
* `clayfrac` CLAYFRAC field [appropriate units] 
* `frlake` FRLAKE field [appropriate units] 
* `frsno` FRSNO field [appropriate units] 
* `gvf` GVF field [appropriate units] 
* `lai` LAI field [appropriate units] 
* `lwi` LWI field [appropriate units] 
* `rdrag` RDRAG field [appropriate units] 
* `sandfrac` SANDFRAC field [appropriate units] 
* `soilm` SOILM field [appropriate units] 
* `ssm` SSM field [appropriate units] 
* `tskin` TSKIN field [appropriate units] 
* `ustar` USTAR field [appropriate units] 
* `ustar_threshold` USTAR\_THRESHOLD field [appropriate units] 
* `z0` Z0 field [appropriate units] 
* `species_radius` Species radius property 
* `species_lower_radius` Species lower\_radius property 
* `species_upper_radius` Species upper\_radius property 
* `species_conc` Species concentrations [ppm or ug/kg] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `dust_emission_total` Total dust emissions for all bins [kg/m2/s] 
* `dust_emission_per_bin` Dust emission flux per bin [kg/m2/s] (num\_species) 
* `dust_horizontal_flux` Total horizontal flux - Q [kg/m2/s] 
* `dust_moisture_correction` Moisture Correction - H [1.0] 
* `dust_effective_threshold` Effective Dust threshold friction velocity: u\_thres \* H / R [m/s] 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/dust/schemes/DustScheme_FENGSHA_Mod.F90`

