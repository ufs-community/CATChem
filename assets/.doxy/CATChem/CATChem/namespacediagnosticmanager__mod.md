

# Namespace diagnosticmanager\_mod



[**Namespace List**](namespaces.md) **>** [**diagnosticmanager\_mod**](namespacediagnosticmanager__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**diagnostic\_manager\_init**](#function-diagnostic_manager_init) (class([**diagnosticmanagertype**](namespacediagnosticmanager__mod.md#none-diagnosticmanagertype)), intent(inout) this, type([**errormanagertype**](namespaceerror__mod.md#none-errormanagertype)), intent(inout), target error\_mgr, integer, intent(out) rc) <br>_Initialize diagnostic manager._  |




























## Public Functions Documentation




### function diagnostic\_manager\_init 

_Initialize diagnostic manager._ 
```Fortran
subroutine diagnosticmanager_mod::diagnostic_manager_init (
    class( diagnosticmanagertype ), intent(inout) this,
    type( errormanagertype ), intent(inout), target error_mgr,
    integer, intent(out) rc
) 
```





**Parameters:**


* `this` DiagnosticManagerType instance 
* `error_mgr` ErrorManager for error handling 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/DiagnosticManager_Mod.F90`

