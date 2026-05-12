

# Namespace carbchemscheme\_gocart\_mod



[**Namespace List**](namespaces.md) **>** [**carbchemscheme\_gocart\_mod**](namespacecarbchemscheme__gocart__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_gocart**](#function-compute_gocart) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**carbchemschemegocartconfig**](namespacecarbchemcommon__mod.md#none-carbchemschemegocartconfig)), intent(in) params, real(fp), intent(in) g0, integer, intent(in) year, integer, intent(in) month, integer, intent(in) day, integer, intent(in) hour, integer, intent(in) minute, integer, intent(in) second, real(fp), dimension(num\_layers), intent(in) airden, real(fp), dimension(num\_layers), intent(in) delp, real(fp), dimension(num\_layers), intent(in) pmid, real(fp), intent(in) tstep, real(fp), dimension(:), intent(in) species\_t\_chem\_loss, character(len=32), dimension(:), intent(in) species\_short\_name, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), dimension(:,:), intent(inout), optional production\_mass\_per\_species\_per\_level, real(fp), dimension(:), intent(inout), optional loss\_flux\_per\_species, real(fp), dimension(:,:), intent(inout), optional phobictophilic\_mass\_per\_species\_per\_level, real(fp), dimension(:), intent(inout), optional phobictophilic\_flux\_per\_species, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Pure science computation for gocart scheme._  |




























## Public Functions Documentation




### function compute\_gocart 

_Pure science computation for gocart scheme._ 
```Fortran
subroutine, public carbchemscheme_gocart_mod::compute_gocart (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( carbchemschemegocartconfig ), intent(in) params,
    real(fp), intent(in) g0,
    integer, intent(in) year,
    integer, intent(in) month,
    integer, intent(in) day,
    integer, intent(in) hour,
    integer, intent(in) minute,
    integer, intent(in) second,
    real(fp), dimension(num_layers), intent(in) airden,
    real(fp), dimension(num_layers), intent(in) delp,
    real(fp), dimension(num_layers), intent(in) pmid,
    real(fp), intent(in) tstep,
    real(fp), dimension(:), intent(in) species_t_chem_loss,
    character(len=32), dimension(:), intent(in) species_short_name,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), dimension(:,:), intent(inout), optional production_mass_per_species_per_level,
    real(fp), dimension(:), intent(inout), optional loss_flux_per_species,
    real(fp), dimension(:,:), intent(inout), optional phobictophilic_mass_per_species_per_level,
    real(fp), dimension(:), intent(inout), optional phobictophilic_flux_per_species,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



This is a pure computational kernel implementing GOCART carbon species chemical production and loss scheme. NO error checking, validation, or infrastructure concerns. Host model must ensure all inputs are valid before calling.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `g0` Required constant from Constants module 
* `year` Time parameter from TimeState (year) 
* `month` Time parameter from TimeState (month) 
* `day` Time parameter from TimeState (day) 
* `hour` Time parameter from TimeState (hour) 
* `minute` Time parameter from TimeState (minute) 
* `second` Time parameter from TimeState (second) 
* `airden` AIRDEN field [appropriate units] 
* `delp` DELP field [appropriate units] 
* `pmid` PMID field [appropriate units] 
* `tstep` Time step [s] - retrieved from process interface 
* `species_t_chem_loss` Species t\_chem\_loss property 
* `species_short_name` Species short\_name property 
* `species_conc` Species concentrations [mol/mol] (num\_layers, num\_species) 
* `species_tendencies` Species tendency terms [mol/mol/s] (num\_layers, num\_species) 
* `Production_mass_per_species_per_level` Production mass (negative for loss) per species per level [kg/kg] (num\_layers, num\_species) 
* `loss_flux_per_species` chemical loss flux per species [kg/m2/s] (num\_species) 
* `PhobicToPhilic_flux_per_species` conversion flux from hydrophobic to hydrophilic per species [kg/m2/s] (num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/carbchem/schemes/CarbChemScheme_GOCART_Mod.F90`

