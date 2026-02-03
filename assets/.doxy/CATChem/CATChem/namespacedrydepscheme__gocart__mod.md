

# Namespace drydepscheme\_gocart\_mod



[**Namespace List**](namespaces.md) **>** [**drydepscheme\_gocart\_mod**](namespacedrydepscheme__gocart__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_gocart**](#function-compute_gocart) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**drydepschemegocartconfig**](namespacedrydepcommon__mod.md#none-drydepschemegocartconfig)), intent(in) params, real(fp), dimension(num\_layers), intent(in) airden, real(fp), intent(in) frlake, real(fp), intent(in) gwettop, real(fp), intent(in) hflux, integer, intent(in) lwi, real(fp), intent(in) pblh, real(fp), dimension(num\_layers), intent(in) t, real(fp), intent(in) tstep, real(fp), intent(in) u10m, real(fp), intent(in) ustar, real(fp), intent(in) v10m, real(fp), dimension(num\_layers+1), intent(in) z, real(fp), intent(in) z0h, real(fp), dimension(num\_species), intent(in) species\_density, real(fp), dimension(num\_species), intent(in) species\_radius, logical, dimension(num\_species), intent(in) species\_is\_seasalt, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, logical, dimension(num\_species), intent(in) is\_gas, real(fp), dimension(:), intent(inout), optional drydep\_con\_per\_species, real(fp), dimension(:), intent(inout), optional drydep\_velocity\_per\_species, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for gocart scheme._  |




























## Public Functions Documentation




### function compute\_gocart 

_Pure science computation for gocart scheme._ 
```Fortran
subroutine, public drydepscheme_gocart_mod::compute_gocart (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( drydepschemegocartconfig ), intent(in) params,
    real(fp), dimension(num_layers), intent(in) airden,
    real(fp), intent(in) frlake,
    real(fp), intent(in) gwettop,
    real(fp), intent(in) hflux,
    integer, intent(in) lwi,
    real(fp), intent(in) pblh,
    real(fp), dimension(num_layers), intent(in) t,
    real(fp), intent(in) tstep,
    real(fp), intent(in) u10m,
    real(fp), intent(in) ustar,
    real(fp), intent(in) v10m,
    real(fp), dimension(num_layers+1), intent(in) z,
    real(fp), intent(in) z0h,
    real(fp), dimension(num_species), intent(in) species_density,
    real(fp), dimension(num_species), intent(in) species_radius,
    logical, dimension(num_species), intent(in) species_is_seasalt,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    logical, dimension(num_species), intent(in) is_gas,
    real(fp), dimension(:), intent(inout), optional drydep_con_per_species,
    real(fp), dimension(:), intent(inout), optional drydep_velocity_per_species,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing GOCART-2G aerosol dry deposition scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `airden` AIRDEN field [appropriate units] 
* `frlake` FRLAKE field [appropriate units] 
* `gwettop` GWETTOP field [appropriate units] 
* `hflux` HFLUX field [appropriate units] 
* `lwi` LWI field [appropriate units] 
* `pblh` PBLH field [appropriate units] 
* `t` T field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `u10m` U10M field [appropriate units] 
* `ustar` USTAR field [appropriate units] 
* `v10m` V10M field [appropriate units] 
* `z0h` Z0H field [appropriate units] 
* `z` Z field [appropriate units] 
* `species_density` Species density property 
* `species_radius` Species radius property 
* `species_is_seasalt` Species is\_seasalt property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `drydep_con_per_species` Dry deposition concentration per species [ug/kg or ppm] (num\_species) 
* `drydep_velocity_per_species` Dry deposition velocity [m/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/schemes/DryDepScheme_GOCART_Mod.F90`

