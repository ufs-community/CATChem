

# Namespace settlingcommon\_mod



[**Namespace List**](namespaces.md) **>** [**settlingcommon\_mod**](namespacesettlingcommon__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  character(len=32) function, public | [**int\_to\_string**](#function-int_to_string) (integer, intent(in) int\_val) <br>_Convert integer to string (utility function)_  |
|  subroutine | [**validate\_settling\_config**](#function-validate_settling_config) (class([**settlingconfig**](namespacesettlingcommon__mod.md#none-settlingconfig)), intent(inout) this, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout) error\_handler) <br>_Validate settling configuration._  |




























## Public Functions Documentation




### function int\_to\_string 

_Convert integer to string (utility function)_ 
```Fortran
character(len=32) function, public settlingcommon_mod::int_to_string (
    integer, intent(in) int_val
) 
```




<hr>



### function validate\_settling\_config 

_Validate settling configuration._ 
```Fortran
subroutine settlingcommon_mod::validate_settling_config (
    class( settlingconfig ), intent(inout) this,
    type( errormanagertype ), intent(inout) error_handler
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/settling/SettlingCommon_Mod.F90`

