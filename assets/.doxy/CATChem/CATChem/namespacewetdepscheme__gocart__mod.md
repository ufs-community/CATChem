

# Namespace wetdepscheme\_gocart\_mod



[**Namespace List**](namespaces.md) **>** [**wetdepscheme\_gocart\_mod**](namespacewetdepscheme__gocart__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**compute\_gocart**](#function-compute_gocart) (integer, intent(in) num\_layers, integer, intent(in) num\_species, type([**wetdepschemegocartconfig**](namespacewetdepcommon__mod.md#none-wetdepschemegocartconfig)), intent(in) params, real(fp), dimension(num\_layers), intent(in) mairden, real(fp), dimension(num\_layers+1), intent(in) pedge, real(fp), dimension(num\_layers+1), intent(in) pfilsan, real(fp), dimension(num\_layers+1), intent(in) pfllsan, real(fp), intent(in) preccon, real(fp), intent(in) preclsc, real(fp), dimension(num\_layers), intent(in) t, real(fp), intent(in) tstep, logical, dimension(:), intent(in) species\_is\_aerosol, character(len=32), dimension(:), intent(in) species\_short\_name, real(fp), dimension(:), intent(in) species\_henry\_cr, real(fp), dimension(:), intent(in) species\_henry\_k0, real(fp), dimension(:), intent(in) species\_henry\_pka, real(fp), dimension(:), intent(in) species\_wd\_retfactor, logical, dimension(:), intent(in) species\_wd\_liqandgas, real(fp), dimension(:), intent(in) species\_wd\_convfaci2g, real(fp), dimension(:,:), intent(in) species\_wd\_rainouteff, real(fp), dimension(:), intent(in) species\_wd\_reevap\_frac, real(fp), dimension(:), intent(in) species\_radius, real(fp), dimension(:), intent(in) species\_mw\_g, real(fp), dimension(num\_layers, num\_species), intent(in) species\_conc, real(fp), dimension(num\_layers, num\_species), intent(inout) species\_tendencies, real(fp), dimension(:,:), intent(inout), optional wetdep\_mass\_per\_species\_per\_level, real(fp), dimension(:,:), intent(inout), optional wetdep\_flux\_per\_species\_per\_level, integer, dimension(:), intent(in), optional diagnostic\_species\_id) <br>_Science computation for gocart wet-removal scheme._  |




























## Public Functions Documentation




### function compute\_gocart 

_Science computation for gocart wet-removal scheme._ 
```Fortran
subroutine, public wetdepscheme_gocart_mod::compute_gocart (
    integer, intent(in) num_layers,
    integer, intent(in) num_species,
    type( wetdepschemegocartconfig ), intent(in) params,
    real(fp), dimension(num_layers), intent(in) mairden,
    real(fp), dimension(num_layers+1), intent(in) pedge,
    real(fp), dimension(num_layers+1), intent(in) pfilsan,
    real(fp), dimension(num_layers+1), intent(in) pfllsan,
    real(fp), intent(in) preccon,
    real(fp), intent(in) preclsc,
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
    real(fp), dimension(:), intent(in) species_wd_reevap_frac,
    real(fp), dimension(:), intent(in) species_radius,
    real(fp), dimension(:), intent(in) species_mw_g,
    real(fp), dimension(num_layers, num_species), intent(in) species_conc,
    real(fp), dimension(num_layers, num_species), intent(inout) species_tendencies,
    real(fp), dimension(:,:), intent(inout), optional wetdep_mass_per_species_per_level,
    real(fp), dimension(:,:), intent(inout), optional wetdep_flux_per_species_per_level,
    integer, dimension(:), intent(in), optional diagnostic_species_id
) 
```



Routes sulfate species (SO4/SO2/DMS/MSA) through GOCART's SU\_Wet\_Removal and all other aerosol species through GOCART's WetRemovalUFS. H2O2 is used only as an auxiliary input to SU\_Wet\_Removal (to cap soluble SO2) and is not itself wet-deposited.




**Parameters:**


* `num_layers` Number of vertical layers 
* `num_species` Number of chemical species 
* `params` Scheme parameters (pre-validated by host) 
* `mairden` MAIRDEN field (moist air density) [kg/m3] 
* `pedge` PEDGE field (level-edge air pressure) [Pa] (num\_layers+1) 
* `pfilsan` PFILSAN field (3D flux of ice nonconvective precip) [kg/m2/s] (num\_layers+1) 
* `pfllsan` PFLLSAN field (3D flux of liquid nonconvective precip) [kg/m2/s] (num\_layers+1) 
* `preccon` PRECCON field (convective precip at ground) [m] - divided by tstep -&gt; [m/s] 
* `preclsc` PRECLSC field (large-scale precip at ground) [m] - divided by tstep -&gt; [m/s] 
* `t` T field (temperature) [K] 
* `tstep` Time step [s] - retrieved from process interface 
* `species_is_aerosol` Species is\_aerosol property 
* `species_short_name` Species short\_name property 
* `species_henry_cr` Species henry\_cr property (unused; GOCART routines carry their own) 
* `species_henry_k0` Species henry\_k0 property (unused) 
* `species_henry_pKa` Species henry\_pKa property (unused) 
* `species_wd_retfactor` Species wd\_retfactor property (unused) 
* `species_wd_LiqAndGas` Species wd\_LiqAndGas property (unused) 
* `species_wd_convfacI2G` Species wd\_convfacI2G property (unused) 
* `species_wd_rainouteff` Species wd\_rainouteff property (3 temperature-dependent efficiencies) (num\_species,3) 
* `species_wd_reevap_frac` Species wd\_reevap\_frac property (unused) 
* `species_radius` Species radius property [um] 
* `species_mw_g` Species mw\_g property 
* `species_conc` Species concentrations [ppm for gases, ug/kg for aerosols] (num\_layers, num\_species) 
* `species_tendencies` Updated species concentrations (replacement mode) (num\_layers, num\_species) 
* `wetdep_mass_per_species_per_level` Wet deposition mass loss per species per level [kg/m2] (num\_layers, num\_species) 
* `wetdep_flux_per_species_per_level` Wet deposition flux per species per level [kg/m2/s] (num\_layers, num\_species) 
* `diagnostic_species_id` Indices mapping diagnostic species to species array (optional, for per-species diagnostics) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/wetdep/schemes/WetDepScheme_GOCART_Mod.F90`

