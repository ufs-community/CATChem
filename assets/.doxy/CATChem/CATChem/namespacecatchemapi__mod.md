

# Namespace catchemapi\_mod



[**Namespace List**](namespaces.md) **>** [**catchemapi\_mod**](namespacecatchemapi__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**catchem\_failure**](#variable-catchem_failure)   = `-1`<br> |
|  integer, parameter, public | [**catchem\_success**](#variable-catchem_success)   = `0`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**instance\_init**](#function-instance_init) (class([**catcheminstancetype**](namespacecatchemapi__mod.md#none-catcheminstancetype)), intent(inout) this, type([**catchemconfigtype**](namespacecatchemapi__mod.md#none-catchemconfigtype)), intent(in) config, integer, intent(out) rc) <br>_Initialize CATChem with basic configuration._  |




























## Public Attributes Documentation




### variable catchem\_failure 

```Fortran
integer, parameter, public catchemapi_mod::catchem_failure;
```




<hr>



### variable catchem\_success 

```Fortran
integer, parameter, public catchemapi_mod::catchem_success;
```




<hr>
## Public Functions Documentation




### function instance\_init 

_Initialize CATChem with basic configuration._ 
```Fortran
subroutine catchemapi_mod::instance_init (
    class( catcheminstancetype ), intent(inout) this,
    type( catchemconfigtype ), intent(in) config,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/CATChemAPI_Mod.F90`

