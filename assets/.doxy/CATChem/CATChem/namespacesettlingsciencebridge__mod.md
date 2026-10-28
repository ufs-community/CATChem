

# Namespace settlingsciencebridge\_mod



[**Namespace List**](namespaces.md) **>** [**settlingsciencebridge\_mod**](namespacesettlingsciencebridge__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  type(gocart2g\_mie), dimension(0), target, save | [**empty\_mie**](#variable-empty_mie)  <br> |
|  character(len=32), dimension(:), allocatable, target, save | [**mie\_names\_store**](#variable-mie_names_store)  <br> |
|  type(gocart2g\_mie), dimension(:), allocatable, target, save | [**mie\_store**](#variable-mie_store)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**run\_emission\_mie\_factor**](#function-run_emission_mie_factor) (integer(c\_int), value n\_columns, integer(c\_int), value n\_levels, character(kind=c\_char), dimension(64), intent(in) mie\_name\_c, real(c\_double), dimension(n\_columns,n\_levels), intent(in) emission\_flux, real(c\_double), dimension(n\_columns,n\_levels), intent(in) rh, real(c\_double), value scale\_factor, real(c\_double), value dt, real(c\_double), dimension(n\_columns), intent(out) factor, integer(c\_int), intent(out) limiter\_rc) <br> |
|  subroutine | [**run\_settling\_mie\_init**](#function-run_settling_mie_init) (integer(c\_int), value n\_files, character(kind=c\_char), dimension(32,n\_files), intent(in) type\_names, character(kind=c\_char), dimension(512,n\_files), intent(in) file\_paths, integer(c\_int), intent(out) init\_rc) <br> |
|  subroutine | [**run\_settling\_science\_bridge**](#function-run_settling_science_bridge) (integer(c\_int), value n\_columns, integer(c\_int), value n\_levels, integer(c\_int), value n\_aerosols, integer(c\_int), value n\_total\_species, real(c\_double), value dt, real(c\_double), value scale\_factor, integer(c\_int), value swelling\_method, integer(c\_int), value correction\_maring, integer(c\_int), value maring\_dust\_only, real(c\_double), dimension(n\_columns,n\_levels), intent(in) airden, real(c\_double), dimension(n\_columns,n\_levels), intent(in) delp, real(c\_double), dimension(n\_columns,n\_levels), intent(in) pmid, real(c\_double), dimension(n\_columns,n\_levels), intent(in) rh, real(c\_double), dimension(n\_columns,n\_levels), intent(in) temperature, real(c\_double), dimension(n\_columns,n\_levels+1), intent(in) z\_edge, character(kind=c\_char), dimension(32,n\_aerosols), intent(in) aerosol\_species\_names, character(kind=c\_char), dimension(32,n\_total\_species), intent(in) species\_names, integer(c\_int), dimension(n\_aerosols), intent(in) species\_is\_dust, real(c\_double), dimension(n\_aerosols), intent(in) radius, real(c\_double), dimension(n\_aerosols), intent(in) density, real(c\_double), dimension(n\_columns,n\_levels,n\_total\_species), intent(inout) concentration, integer(c\_int), value simple\_scheme, character(kind=c\_char), dimension(32,n\_aerosols), intent(in) aerosol\_mie\_names, type(c\_ptr), value diag\_velocity, type(c\_ptr), value diag\_flux, integer(c\_int), dimension(max(n\_diag\_species,1)), intent(in) diagnostic\_species\_id, integer(c\_int), value n\_diag\_species, integer(c\_int), intent(out) bridge\_rc) <br> |




























## Public Attributes Documentation




### variable empty\_mie 

```Fortran
type(gocart2g_mie), dimension(0), target, save settlingsciencebridge_mod::empty_mie;
```




<hr>



### variable mie\_names\_store 

```Fortran
character(len=32), dimension(:), allocatable, target, save settlingsciencebridge_mod::mie_names_store;
```




<hr>



### variable mie\_store 

```Fortran
type(gocart2g_mie), dimension(:), allocatable, target, save settlingsciencebridge_mod::mie_store;
```




<hr>
## Public Functions Documentation




### function run\_emission\_mie\_factor 

```Fortran
subroutine settlingsciencebridge_mod::run_emission_mie_factor (
    integer(c_int), value n_columns,
    integer(c_int), value n_levels,
    character(kind=c_char), dimension(64), intent(in) mie_name_c,
    real(c_double), dimension(n_columns,n_levels), intent(in) emission_flux,
    real(c_double), dimension(n_columns,n_levels), intent(in) rh,
    real(c_double), value scale_factor,
    real(c_double), value dt,
    real(c_double), dimension(n_columns), intent(out) factor,
    integer(c_int), intent(out) limiter_rc
) 
```




<hr>



### function run\_settling\_mie\_init 

```Fortran
subroutine settlingsciencebridge_mod::run_settling_mie_init (
    integer(c_int), value n_files,
    character(kind=c_char), dimension(32,n_files), intent(in) type_names,
    character(kind=c_char), dimension(512,n_files), intent(in) file_paths,
    integer(c_int), intent(out) init_rc
) 
```




<hr>



### function run\_settling\_science\_bridge 

```Fortran
subroutine settlingsciencebridge_mod::run_settling_science_bridge (
    integer(c_int), value n_columns,
    integer(c_int), value n_levels,
    integer(c_int), value n_aerosols,
    integer(c_int), value n_total_species,
    real(c_double), value dt,
    real(c_double), value scale_factor,
    integer(c_int), value swelling_method,
    integer(c_int), value correction_maring,
    integer(c_int), value maring_dust_only,
    real(c_double), dimension(n_columns,n_levels), intent(in) airden,
    real(c_double), dimension(n_columns,n_levels), intent(in) delp,
    real(c_double), dimension(n_columns,n_levels), intent(in) pmid,
    real(c_double), dimension(n_columns,n_levels), intent(in) rh,
    real(c_double), dimension(n_columns,n_levels), intent(in) temperature,
    real(c_double), dimension(n_columns,n_levels+1), intent(in) z_edge,
    character(kind=c_char), dimension(32,n_aerosols), intent(in) aerosol_species_names,
    character(kind=c_char), dimension(32,n_total_species), intent(in) species_names,
    integer(c_int), dimension(n_aerosols), intent(in) species_is_dust,
    real(c_double), dimension(n_aerosols), intent(in) radius,
    real(c_double), dimension(n_aerosols), intent(in) density,
    real(c_double), dimension(n_columns,n_levels,n_total_species), intent(inout) concentration,
    integer(c_int), value simple_scheme,
    character(kind=c_char), dimension(32,n_aerosols), intent(in) aerosol_mie_names,
    type(c_ptr), value diag_velocity,
    type(c_ptr), value diag_flux,
    integer(c_int), dimension(max(n_diag_species,1)), intent(in) diagnostic_species_id,
    integer(c_int), value n_diag_species,
    integer(c_int), intent(out) bridge_rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/settling/SettlingScienceBridge.F90`

