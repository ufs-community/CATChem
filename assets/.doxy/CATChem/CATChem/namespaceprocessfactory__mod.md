

# Namespace processfactory\_mod



[**Namespace List**](namespaces.md) **>** [**processfactory\_mod**](namespaceprocessfactory__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  class([**processinterface**](namespaceprocessinterface__mod.md#none-processinterface)) function, allocatable, public | [**create\_process**](#function-create_process) (character(len=\*), intent(in) process\_name, type([**statemanagertype**](namespacestatemanager__mod.md#none-statemanagertype)), intent(inout) container, integer, intent(out) rc) <br>_Module-level convenience function._  |
|  subroutine | [**factory\_init**](#function-factory_init) (class([**processfactorytype**](namespaceprocessfactory__mod.md#none-processfactorytype)), intent(inout) this, integer, intent(out) rc) <br> |




























## Public Functions Documentation




### function create\_process 

_Module-level convenience function._ 
```Fortran
class( processinterface ) function, allocatable, public processfactory_mod::create_process (
    character(len=*), intent(in) process_name,
    type( statemanagertype ), intent(inout) container,
    integer, intent(out) rc
) 
```




<hr>



### function factory\_init 

```Fortran
subroutine processfactory_mod::factory_init (
    class( processfactorytype ), intent(inout) this,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/ProcessFactory_Mod.F90`

