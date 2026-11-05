

# Namespace carbchemcommon\_mod



[**Namespace List**](namespaces.md) **>** [**carbchemcommon\_mod**](namespacecarbchemcommon__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  character(len=32) function, public | [**int\_to\_string**](#function-int_to_string) (integer, intent(in) int\_val) <br>_Convert integer to string (utility function)_  |
|  subroutine | [**validate\_carbchem\_config**](#function-validate_carbchem_config) (class([**carbchemconfig**](namespacecarbchemcommon__mod.md#none-carbchemconfig)), intent(inout) this, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout) error\_handler) <br>_Validate carbchem configuration._  |




























## Public Functions Documentation




### function int\_to\_string 

_Convert integer to string (utility function)_ 
```Fortran
character(len=32) function, public carbchemcommon_mod::int_to_string (
    integer, intent(in) int_val
) 
```




<hr>



### function validate\_carbchem\_config 

_Validate carbchem configuration._ 
```Fortran
subroutine carbchemcommon_mod::validate_carbchem_config (
    class( carbchemconfig ), intent(inout) this,
    type( errormanagertype ), intent(inout) error_handler
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/carbchem/CarbChemCommon_Mod.F90`

