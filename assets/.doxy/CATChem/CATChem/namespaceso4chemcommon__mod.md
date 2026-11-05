

# Namespace so4chemcommon\_mod



[**Namespace List**](namespaces.md) **>** [**so4chemcommon\_mod**](namespaceso4chemcommon__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  character(len=32) function, public | [**int\_to\_string**](#function-int_to_string) (integer, intent(in) int\_val) <br>_Convert integer to string (utility function)_  |
|  subroutine | [**validate\_so4chem\_config**](#function-validate_so4chem_config) (class([**so4chemconfig**](namespaceso4chemcommon__mod.md#none-so4chemconfig)), intent(inout) this, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout) error\_handler) <br>_Validate so4chem configuration._  |




























## Public Functions Documentation




### function int\_to\_string 

_Convert integer to string (utility function)_ 
```Fortran
character(len=32) function, public so4chemcommon_mod::int_to_string (
    integer, intent(in) int_val
) 
```




<hr>



### function validate\_so4chem\_config 

_Validate so4chem configuration._ 
```Fortran
subroutine so4chemcommon_mod::validate_so4chem_config (
    class( so4chemconfig ), intent(inout) this,
    type( errormanagertype ), intent(inout) error_handler
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/so4chem/SO4chemCommon_Mod.F90`

