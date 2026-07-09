

# Namespace dustcommon\_mod



[**Namespace List**](namespaces.md) **>** [**dustcommon\_mod**](namespacedustcommon__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  character(len=32) function, public | [**int\_to\_string**](#function-int_to_string) (integer, intent(in) int\_val) <br>_Convert integer to string (utility function)_  |
|  subroutine | [**validate\_dust\_config**](#function-validate_dust_config) (class([**dustconfig**](namespacedustcommon__mod.md#none-dustconfig)), intent(inout) this, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout) error\_handler) <br>_Validate dust configuration._  |




























## Public Functions Documentation




### function int\_to\_string 

_Convert integer to string (utility function)_ 
```Fortran
character(len=32) function, public dustcommon_mod::int_to_string (
    integer, intent(in) int_val
) 
```




<hr>



### function validate\_dust\_config 

_Validate dust configuration._ 
```Fortran
subroutine dustcommon_mod::validate_dust_config (
    class( dustconfig ), intent(inout) this,
    type( errormanagertype ), intent(inout) error_handler
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/dust/DustCommon_Mod.F90`

