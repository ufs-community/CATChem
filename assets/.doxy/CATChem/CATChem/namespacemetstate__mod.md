

# Namespace metstate\_mod



[**Namespace List**](namespaces.md) **>** [**metstate\_mod**](namespacemetstate__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**metstate\_init**](#function-metstate_init) (class([**metstatetype**](namespacemetstate__mod.md#none-metstatetype)), intent(inout) this, integer, intent(in) nx, integer, intent(in) ny, integer, intent(in) nlevs, integer, intent(in), optional nsoil, integer, intent(in), optional nsoiltype, integer, intent(in), optional nsurftype, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout), pointer error\_mgr, integer, intent(out) rc) <br>_Initialize a MetStateType object._  |




























## Public Functions Documentation




### function metstate\_init 

_Initialize a MetStateType object._ 
```Fortran
subroutine metstate_mod::metstate_init (
    class( metstatetype ), intent(inout) this,
    integer, intent(in) nx,
    integer, intent(in) ny,
    integer, intent(in) nlevs,
    integer, intent(in), optional nsoil,
    integer, intent(in), optional nsoiltype,
    integer, intent(in), optional nsurftype,
    type( errormanagertype ), intent(inout), pointer error_mgr,
    integer, intent(out) rc
) 
```



Initializes the meteorological state object, sets default values, and allocates required arrays.




**Parameters:**


* `this` MetStateType object to initialize 
* `nx` Number of grid points in x direction 
* `ny` Number of grid points in y direction 
* `nlevs` Number of vertical levels 
* `nsoil` Number of soil layers 
* `nsoiltype` Number of soil types 
* `nsurftype` Number of surface types 
* `error_mgr` Error manager for context and error reporting 
* `rc` Return code (CC\_SUCCESS or error code) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/metstate_mod.F90`

