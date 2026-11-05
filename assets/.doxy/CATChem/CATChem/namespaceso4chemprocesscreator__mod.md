

# Namespace so4chemprocesscreator\_mod



[**Namespace List**](namespaces.md) **>** [**so4chemprocesscreator\_mod**](namespaceso4chemprocesscreator__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**create\_so4chem\_process**](#function-create_so4chem_process) (class([**processinterface**](namespaceprocessinterface__mod.md#none-processinterface)), intent(out), allocatable process, integer, intent(out) rc) <br>_Create a new so4chem process instance._  |
|  subroutine, public | [**get\_so4chem\_default\_config**](#function-get_so4chem_default_config) (character(len=\*), intent(out) config\_data) <br>_Get default configuration for so4chem process._  |
|  subroutine, public | [**register\_so4chem\_process**](#function-register_so4chem_process) (type([**processmanagertype**](namespaceprocessmanager__mod.md#none-processmanagertype)), intent(inout) process\_mgr, integer, intent(out) rc) <br>_Register the so4chem process with a ProcessManager._  |




























## Public Functions Documentation




### function create\_so4chem\_process 

_Create a new so4chem process instance._ 
```Fortran
subroutine, public so4chemprocesscreator_mod::create_so4chem_process (
    class( processinterface ), intent(out), allocatable process,
    integer, intent(out) rc
) 
```



This factory function creates and returns a new instance of the so4chem process. The process is not initialized - the caller must call the init() method with appropriate configuration.




**Parameters:**


* `process` Allocated process instance 
* `rc` Return code 




        

<hr>



### function get\_so4chem\_default\_config 

_Get default configuration for so4chem process._ 
```Fortran
subroutine, public so4chemprocesscreator_mod::get_so4chem_default_config (
    character(len=*), intent(out) config_data
) 
```



This function returns a default configuration string that can be used to initialize the so4chem process with reasonable defaults.




**Parameters:**


* `config_data` Default configuration string 




        

<hr>



### function register\_so4chem\_process 

_Register the so4chem process with a ProcessManager._ 
```Fortran
subroutine, public so4chemprocesscreator_mod::register_so4chem_process (
    type( processmanagertype ), intent(inout) process_mgr,
    integer, intent(out) rc
) 
```



This subroutine registers the so4chem process with a ProcessManager's factory. This is the correct way to register processes for use in applications and integration tests.




**Parameters:**


* `process_mgr` The ProcessManager to register with 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/so4chem/SO4chemProcessCreator_Mod.F90`

