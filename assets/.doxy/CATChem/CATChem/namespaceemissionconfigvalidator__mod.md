

# Namespace emissionconfigvalidator\_mod



[**Namespace List**](namespaces.md) **>** [**emissionconfigvalidator\_mod**](namespaceemissionconfigvalidator__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**validation\_error**](#variable-validation_error)   = `2`<br> |
|  integer, parameter, public | [**validation\_success**](#variable-validation_success)   = `0`<br> |
|  integer, parameter, public | [**validation\_warning**](#variable-validation_warning)   = `1`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**validate\_emission\_config**](#function-validate_emission_config) (class([**emissionconfigvalidatortype**](namespaceemissionconfigvalidator__mod.md#none-emissionconfigvalidatortype)), intent(inout) this, character(len=\*), intent(in) config\_file, type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, type([**validationresulttype**](namespaceemissionconfigvalidator__mod.md#none-validationresulttype)), dimension(:), intent(out), allocatable results, integer, intent(out) rc) <br>_Validate complete emission configuration._  |




























## Public Attributes Documentation




### variable validation\_error 

```Fortran
integer, parameter, public emissionconfigvalidator_mod::validation_error;
```




<hr>



### variable validation\_success 

```Fortran
integer, parameter, public emissionconfigvalidator_mod::validation_success;
```




<hr>



### variable validation\_warning 

```Fortran
integer, parameter, public emissionconfigvalidator_mod::validation_warning;
```




<hr>
## Public Functions Documentation




### function validate\_emission\_config 

_Validate complete emission configuration._ 
```Fortran
subroutine emissionconfigvalidator_mod::validate_emission_config (
    class( emissionconfigvalidatortype ), intent(inout) this,
    character(len=*), intent(in) config_file,
    type( statemanagertype ), intent(inout) container,
    type( validationresulttype ), dimension(:), intent(out), allocatable results,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/EmissionConfigValidator_Mod.F90`

