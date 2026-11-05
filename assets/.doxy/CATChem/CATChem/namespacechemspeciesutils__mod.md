

# Namespace chemspeciesutils\_mod



[**Namespace List**](namespaces.md) **>** [**chemspeciesutils\_mod**](namespacechemspeciesutils__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  logical function, public | [**check\_species\_exists**](#function-check_species_exists) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), intent(in) species\_name, integer, intent(out) rc) <br>_Check if a chemical species exists in the mechanism._  |
|  subroutine, public | [**create\_species\_mapping**](#function-create_species_mapping) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), dimension(:), intent(in) process\_species, integer, dimension(:), intent(out) species\_mapping, integer, intent(out) rc) <br>_Create mapping from process species to mechanism species._  |
|  subroutine, public | [**filter\_species\_by\_type**](#function-filter_species_by_type) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), intent(in) species\_type, character(len=32), dimension(:), intent(out), allocatable filtered\_species, integer, intent(out) rc) <br>_Filter species by type._  |
|  subroutine, public | [**get\_aerosol\_species\_list**](#function-get_aerosol_species_list) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=32), dimension(:), intent(out), allocatable aerosol\_species\_list, integer, intent(out) rc) <br>_Get list of aerosol species._  |
|  subroutine, public | [**get\_dust\_species\_list**](#function-get_dust_species_list) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=32), dimension(:), intent(out), allocatable dust\_species\_list, integer, intent(out) rc) <br>_Get list of dust species._  |
|  subroutine, public | [**get\_gas\_species\_list**](#function-get_gas_species_list) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=32), dimension(:), intent(out), allocatable gas\_species\_list, integer, intent(out) rc) <br>_Get list of gas species._  |
|  subroutine, public | [**get\_seasalt\_species\_list**](#function-get_seasalt_species_list) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=32), dimension(:), intent(out), allocatable seasalt\_species\_list, integer, intent(out) rc) <br>_Get list of sea salt species._  |
|  subroutine, public | [**get\_species\_concentration\_units**](#function-get_species_concentration_units) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), intent(in) species\_name, character(len=10), intent(out) units, integer, intent(out) rc) <br>_Get concentration units for a species._  |
|  integer function, public | [**get\_species\_index**](#function-get_species_index) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), intent(in) species\_name, integer, intent(out) rc) <br>_Get the index of a single chemical species._  |
|  subroutine, public | [**get\_species\_indices**](#function-get_species_indices) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), dimension(:), intent(in) species\_names, integer, dimension(:), intent(out) species\_indices, integer, intent(out) rc) <br>_Get indices for multiple chemical species._  |
|  subroutine, public | [**get\_species\_properties**](#function-get_species_properties) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=\*), intent(in) species\_name, type([**speciestype**](namespacespecies__mod.md#none-speciestype)), intent(out) species, integer, intent(out) rc) <br>_Get properties of a chemical species._  |
|  subroutine, public | [**get\_tracer\_species\_list**](#function-get_tracer_species_list) (type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, character(len=32), dimension(:), intent(out), allocatable tracer\_species\_list, integer, intent(out) rc) <br>_Get list of tracer species._  |
|  subroutine, public | [**parse\_species\_list**](#function-parse_species_list) (character(len=\*), intent(in) species\_string, character(len=32), dimension(:), intent(out) species\_array, integer, intent(out) n\_species) <br>_Parse a comma-separated species list._  |




























## Public Functions Documentation




### function check\_species\_exists 

_Check if a chemical species exists in the mechanism._ 
```Fortran
logical function, public chemspeciesutils_mod::check_species_exists (
    type( statemanagertype ), intent(inout) container,
    character(len=*), intent(in) species_name,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `species_name` Name of the species to check 
* `rc` Return code 



**Returns:**

True if species exists, false otherwise 





        

<hr>



### function create\_species\_mapping 

_Create mapping from process species to mechanism species._ 
```Fortran
subroutine, public chemspeciesutils_mod::create_species_mapping (
    type( statemanagertype ), intent(inout) container,
    character(len=*), dimension(:), intent(in) process_species,
    integer, dimension(:), intent(out) species_mapping,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `process_species` Array of process species names 
* `species_mapping` Array of mechanism species indices 
* `rc` Return code 




        

<hr>



### function filter\_species\_by\_type 

_Filter species by type._ 
```Fortran
subroutine, public chemspeciesutils_mod::filter_species_by_type (
    type( statemanagertype ), intent(inout) container,
    character(len=*), intent(in) species_type,
    character(len=32), dimension(:), intent(out), allocatable filtered_species,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `species_type` Type of species ('gas', 'aerosol', 'dust', 'seasalt', 'tracer') 
* `filtered_species` List of species names matching the type 
* `rc` Return code 




        

<hr>



### function get\_aerosol\_species\_list 

_Get list of aerosol species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_aerosol_species_list (
    type( statemanagertype ), intent(inout) container,
    character(len=32), dimension(:), intent(out), allocatable aerosol_species_list,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `aerosol_species_list` List of aerosol species names 
* `rc` Return code 




        

<hr>



### function get\_dust\_species\_list 

_Get list of dust species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_dust_species_list (
    type( statemanagertype ), intent(inout) container,
    character(len=32), dimension(:), intent(out), allocatable dust_species_list,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `dust_species_list` List of dust species names 
* `rc` Return code 




        

<hr>



### function get\_gas\_species\_list 

_Get list of gas species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_gas_species_list (
    type( statemanagertype ), intent(inout) container,
    character(len=32), dimension(:), intent(out), allocatable gas_species_list,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `gas_species_list` List of gas species names 
* `rc` Return code 




        

<hr>



### function get\_seasalt\_species\_list 

_Get list of sea salt species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_seasalt_species_list (
    type( statemanagertype ), intent(inout) container,
    character(len=32), dimension(:), intent(out), allocatable seasalt_species_list,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `seasalt_species_list` List of sea salt species names 
* `rc` Return code 




        

<hr>



### function get\_species\_concentration\_units 

_Get concentration units for a species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_species_concentration_units (
    type( statemanagertype ), intent(inout) container,
    character(len=*), intent(in) species_name,
    character(len=10), intent(out) units,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `species_name` Name of the species 
* `units` Units string ('v/v', 'kg/kg', etc.) 
* `rc` Return code 




        

<hr>



### function get\_species\_index 

_Get the index of a single chemical species._ 
```Fortran
integer function, public chemspeciesutils_mod::get_species_index (
    type( statemanagertype ), intent(inout) container,
    character(len=*), intent(in) species_name,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `species_name` Name of the species to find 
* `rc` Return code 



**Returns:**

Species index (&gt; 0 if found, 0 if not found) 





        

<hr>



### function get\_species\_indices 

_Get indices for multiple chemical species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_species_indices (
    type( statemanagertype ), intent(inout) container,
    character(len=*), dimension(:), intent(in) species_names,
    integer, dimension(:), intent(out) species_indices,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `species_names` Array of species names to find 
* `species_indices` Array of species indices (0 if not found) 
* `rc` Return code 




        

<hr>



### function get\_species\_properties 

_Get properties of a chemical species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_species_properties (
    type( statemanagertype ), intent(inout) container,
    character(len=*), intent(in) species_name,
    type( speciestype ), intent(out) species,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `species_name` Name of the species 
* `species` Properties of the species 
* `rc` Return code 




        

<hr>



### function get\_tracer\_species\_list 

_Get list of tracer species._ 
```Fortran
subroutine, public chemspeciesutils_mod::get_tracer_species_list (
    type( statemanagertype ), intent(inout) container,
    character(len=32), dimension(:), intent(out), allocatable tracer_species_list,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` StateManager containing chemical state 
* `tracer_species_list` List of tracer species names 
* `rc` Return code 




        

<hr>



### function parse\_species\_list 

_Parse a comma-separated species list._ 
```Fortran
subroutine, public chemspeciesutils_mod::parse_species_list (
    character(len=*), intent(in) species_string,
    character(len=32), dimension(:), intent(out) species_array,
    integer, intent(out) n_species
) 
```





**Parameters:**


* `species_string` Comma-separated species names 
* `species_array` Array of parsed species names 
* `n_species` Number of species parsed 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/ChemSpeciesUtils_Mod.F90`

