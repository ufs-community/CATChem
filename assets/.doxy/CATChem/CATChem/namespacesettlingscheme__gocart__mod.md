

# Namespace settlingscheme\_gocart\_mod



[**Namespace List**](namespaces.md) **>** [**settlingscheme\_gocart\_mod**](namespacesettlingscheme__gocart__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_gocart**](#function-compute_gocart) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**settlingschemegocartconfig**](namespacesettlingcommon__mod.md#none-settlingschemegocartconfig)), intent(in) params, real(fp), dimension(num\_layers), intent(in) airden, real(fp), dimension(num\_layers), intent(in) delp, real(fp), dimension(num\_layers), intent(in) pmid, real(fp), dimension(num\_layers), intent(in) rh, real(fp), dimension(num\_layers), intent(in) t, real(fp), intent(in) tstep, real(fp), dimension(num\_layers+1), intent(in) z, character(len=32), dimension(:), intent(in) species\_short\_name, type(gocart2g\_mie), dimension(:), intent(in) mie\_data, integer, dimension(num\_species), intent(in) species\_mie\_map, real(fp), dimension(:), intent(in) species\_radius, real(fp), dimension(:), intent(in) species\_density, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), dimension(:,:), intent(inout), optional settling\_velocity\_per\_species\_per\_level, real(fp), dimension(:), intent(inout), optional settling\_flux\_per\_species, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for gocart scheme._  |




























## Public Functions Documentation




### function compute\_gocart 

_Pure science computation for gocart scheme._ 
```Fortran
subroutine, public settlingscheme_gocart_mod::compute_gocart (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( settlingschemegocartconfig ), intent(in) params,
    real(fp), dimension(num_layers), intent(in) airden,
    real(fp), dimension(num_layers), intent(in) delp,
    real(fp), dimension(num_layers), intent(in) pmid,
    real(fp), dimension(num_layers), intent(in) rh,
    real(fp), dimension(num_layers), intent(in) t,
    real(fp), intent(in) tstep,
    real(fp), dimension(num_layers+1), intent(in) z,
    character(len=32), dimension(:), intent(in) species_short_name,
    type(gocart2g_mie), dimension(:), intent(in) mie_data,
    integer, dimension(num_species), intent(in) species_mie_map,
    real(fp), dimension(:), intent(in) species_radius,
    real(fp), dimension(:), intent(in) species_density,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), dimension(:,:), intent(inout), optional settling_velocity_per_species_per_level,
    real(fp), dimension(:), intent(inout), optional settling_flux_per_species,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing GOCART gravitational settling scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `airden` AIRDEN field [appropriate units] 
* `delp` DELP field [appropriate units] 
* `pmid` PMID field [appropriate units] 
* `rh` RH field [appropriate units] 
* `t` T field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `z` Z field [appropriate units] 
* `species_short_name` Species short\_name property 
* `mie_data` Complete Mie data array from ChemState 
* `species_mie_map` Mapping from process species to MieData indices 
* `species_radius` Species radius property 
* `species_density` Species density property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `settling_velocity_per_species_per_level` settling velocity per species per level [m/s] (num\_layers, num\_species) 
* `settling_flux_per_species` settling flux per species across column [kg/m2/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/settling/schemes/SettlingScheme_GOCART_Mod.F90`

