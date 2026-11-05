

# Namespace catchemcore\_mod



[**Namespace List**](namespaces.md) **>** [**catchemcore\_mod**](namespacecatchemcore__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**core\_init**](#function-core_init) (class([**catchemcoretype**](namespacecatchemcore__mod.md#none-catchemcoretype)), intent(inout) this, character(len=\*), intent(in), optional name, integer, intent(out) rc) <br>_Initialize the CATChem core._  |




























## Public Functions Documentation




### function core\_init 

_Initialize the CATChem core._ 
```Fortran
subroutine catchemcore_mod::core_init (
    class( catchemcoretype ), intent(inout) this,
    character(len=*), intent(in), optional name,
    integer, intent(out) rc
) 
```



Sets up the core framework and initializes the error manager. Components are initialized later during configuration. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/CATChemCore_Mod.F90`

