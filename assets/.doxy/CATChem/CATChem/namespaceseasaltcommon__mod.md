

# Namespace seasaltcommon\_mod



[**Namespace List**](namespaces.md) **>** [**seasaltcommon\_mod**](namespaceseasaltcommon__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  character(len=32) function, public | [**int\_to\_string**](#function-int_to_string) (integer, intent(in) int\_val) <br>_Convert integer to string (utility function)_  |
|  subroutine | [**validate\_seasalt\_config**](#function-validate_seasalt_config) (class([**seasaltconfig**](namespaceseasaltcommon__mod.md#none-seasaltconfig)), intent(inout) this, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout) error\_handler) <br>_Validate seasalt configuration._  |




























## Public Functions Documentation




### function int\_to\_string 

_Convert integer to string (utility function)_ 
```Fortran
character(len=32) function, public seasaltcommon_mod::int_to_string (
    integer, intent(in) int_val
) 
```




<hr>



### function validate\_seasalt\_config 

_Validate seasalt configuration._ 
```Fortran
subroutine seasaltcommon_mod::validate_seasalt_config (
    class( seasaltconfig ), intent(inout) this,
    type( errormanagertype ), intent(inout) error_handler
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/seasalt/SeaSaltCommon_Mod.F90`

