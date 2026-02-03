

# Namespace configmanager\_mod



[**Namespace List**](namespaces.md) **>** [**configmanager\_mod**](namespaceconfigmanager__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**config\_strategy\_fallback**](#variable-config_strategy_fallback)   = `3`<br>_Use defaults on errors._  |
|  integer, parameter, public | [**config\_strategy\_permissive**](#variable-config_strategy_permissive)   = `2`<br>_Warning on validation errors._  |
|  integer, parameter, public | [**config\_strategy\_strict**](#variable-config_strategy_strict)   = `1`<br>_Configuration loading strategies._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**discover\_nested\_yaml\_section\_items**](#function-discover_nested_yaml_section_items) (character(len=\*), intent(in) filename, character(len=\*), intent(in) section\_path, character(len=\*), dimension(:), intent(inout) item\_names, integer, intent(out) n\_items, integer, intent(out) rc, character(len=\*), intent(in), optional search\_mode) <br>_Discover items in a nested YAML section (supports arbitrary depth paths)_  |
|  subroutine, public | [**discover\_yaml\_section\_items**](#function-discover_yaml_section_items) (character(len=\*), intent(in) filename, character(len=\*), intent(in) section\_name, character(len=\*), intent(in) parse\_mode, character(len=64), dimension(:), intent(inout) item\_names, integer, intent(out) n\_items, integer, intent(out) rc) <br>_Generic YAML text parser for discovering items in a section This function reads the YAML file as plain text and extracts item names by parsing the structure line by line, looking for keys with ':' under the section._  |
|  subroutine | [**schema\_init**](#function-schema_init) (class([**configschematype**](namespaceconfigmanager__mod.md#none-configschematype)), intent(inout) this, character(len=\*), intent(in) name, character(len=\*), intent(in) description, logical, intent(in), optional strict) <br>_Initialize configuration schema._  |




























## Public Attributes Documentation




### variable config\_strategy\_fallback 

_Use defaults on errors._ 
```Fortran
integer, parameter, public configmanager_mod::config_strategy_fallback;
```




<hr>



### variable config\_strategy\_permissive 

_Warning on validation errors._ 
```Fortran
integer, parameter, public configmanager_mod::config_strategy_permissive;
```




<hr>



### variable config\_strategy\_strict 

_Configuration loading strategies._ 
```Fortran
integer, parameter, public configmanager_mod::config_strategy_strict;
```



Fail on any validation error 


        

<hr>
## Public Functions Documentation




### function discover\_nested\_yaml\_section\_items 

_Discover items in a nested YAML section (supports arbitrary depth paths)_ 
```Fortran
subroutine, public configmanager_mod::discover_nested_yaml_section_items (
    character(len=*), intent(in) filename,
    character(len=*), intent(in) section_path,
    character(len=*), dimension(:), intent(inout) item_names,
    integer, intent(out) n_items,
    integer, intent(out) rc,
    character(len=*), intent(in), optional search_mode
) 
```



This function discovers direct child items in any nested YAML section. It supports arbitrary nesting depth and flexible indentation.


Examples:
* "processes/extemis" -&gt; finds anthro, point, fire, fengsha
* "processes/extemis/anthro" -&gt; finds activate, scale\_factor, source\_file, etc.
* "mie/files" with mode 'key\_value\_pairs' -&gt; finds "SS: opticsBands\_SS.v3\_3.RRTMG.nc"






**Parameters:**


* `filename` YAML file to parse 
* `section_path` Nested path (e.g., "processes/extemis/anthro") 
* `item_names` Array to store discovered item names 
* `n_items` Number of items found 
* `rc` Return code 
* `search_mode` Optional: 'section\_headers' (default) or 'key\_value\_pairs' 




        

<hr>



### function discover\_yaml\_section\_items 

_Generic YAML text parser for discovering items in a section This function reads the YAML file as plain text and extracts item names by parsing the structure line by line, looking for keys with ':' under the section._ 
```Fortran
subroutine, public configmanager_mod::discover_yaml_section_items (
    character(len=*), intent(in) filename,
    character(len=*), intent(in) section_name,
    character(len=*), intent(in) parse_mode,
    character(len=64), dimension(:), intent(inout) item_names,
    integer, intent(out) n_items,
    integer, intent(out) rc
) 
```





**Parameters:**


* `filename` YAML configuration file path 
* `section_name` Section name to search for items 
* `parse_mode` Parsing mode: 'simple', 'emission\_fields' 
* `item_names` Array to store discovered item names 
* `n_items` Number of items discovered 
* `rc` Return code 




        

<hr>



### function schema\_init 

_Initialize configuration schema._ 
```Fortran
subroutine configmanager_mod::schema_init (
    class( configschematype ), intent(inout) this,
    character(len=*), intent(in) name,
    character(len=*), intent(in) description,
    logical, intent(in), optional strict
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/ConfigManager_Mod.F90`

