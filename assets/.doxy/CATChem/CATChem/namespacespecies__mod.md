

# Namespace species\_mod



[**Namespace List**](namespaces.md) **>** [**species\_mod**](namespacespecies__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(fp), parameter, public | [**missing\_vv**](#variable-missing_vv)   = `1.0e-20\_fp`<br>_Missing species concentration value._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**create\_species\_database**](#function-create_species_database) (type([**speciesmanagertype**](namespacespecies__mod.md#none-speciesmanagertype)), intent(out) species\_mgr, integer, intent(out) rc) <br>_Create basic species database._  |
|  integer function, public | [**find\_species\_by\_name**](#function-find_species_by_name) (type([**speciestype**](namespacespecies__mod.md#none-speciestype)), dimension(:), intent(in) species\_db, integer, intent(in) num\_species, character(len=\*), intent(in) species\_name, integer, intent(out) rc) <br>_Find species by name (standalone function)_  |
|  logical function, public | [**validate\_species**](#function-validate_species) (type([**speciestype**](namespacespecies__mod.md#none-speciestype)), intent(in) species, integer, intent(out) rc) <br>_Standalone species validation function._  |




























## Public Attributes Documentation




### variable missing\_vv 

_Missing species concentration value._ 
```Fortran
real(fp), parameter, public species_mod::missing_vv;
```




<hr>
## Public Functions Documentation




### function create\_species\_database 

_Create basic species database._ 
```Fortran
subroutine, public species_mod::create_species_database (
    type( speciesmanagertype ), intent(out) species_mgr,
    integer, intent(out) rc
) 
```



This subroutine creates a basic species database with common species.




**Parameters:**


* `species_mgr` Initialized species manager 
* `rc` Return code 




        

<hr>



### function find\_species\_by\_name 

_Find species by name (standalone function)_ 
```Fortran
integer function, public species_mod::find_species_by_name (
    type( speciestype ), dimension(:), intent(in) species_db,
    integer, intent(in) num_species,
    character(len=*), intent(in) species_name,
    integer, intent(out) rc
) 
```



This function provides species lookup functionality.




**Parameters:**


* `species_db` Array of species 
* `num_species` Number of species in database 
* `species_name` Name to search for 
* `rc` Return code 




        

<hr>



### function validate\_species 

_Standalone species validation function._ 
```Fortran
logical function, public species_mod::validate_species (
    type( speciestype ), intent(in) species,
    integer, intent(out) rc
) 
```



This function provides species validation outside of the class methods.




**Parameters:**


* `species` Species to validate 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/species_mod.F90`

