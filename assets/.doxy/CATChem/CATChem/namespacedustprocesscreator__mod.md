

# Namespace dustprocesscreator\_mod



[**Namespace List**](namespaces.md) **>** [**dustprocesscreator\_mod**](namespacedustprocesscreator__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**create\_dust\_process**](#function-create_dust_process) (class([**processinterface**](namespaceprocessinterface__mod.md#none-processinterface)), intent(out), allocatable process, integer, intent(out) rc) <br>_Create a new dust process instance._  |
|  subroutine, public | [**get\_dust\_default\_config**](#function-get_dust_default_config) (character(len=\*), intent(out) config\_data) <br>_Get default configuration for dust process._  |
|  subroutine, public | [**register\_dust\_process**](#function-register_dust_process) (type([**processmanagertype**](namespaceprocessmanager__mod.md#none-processmanagertype)), intent(inout) process\_mgr, integer, intent(out) rc) <br>_Register the dust process with a ProcessManager._  |




























## Public Functions Documentation




### function create\_dust\_process 

_Create a new dust process instance._ 
```Fortran
subroutine, public dustprocesscreator_mod::create_dust_process (
    class( processinterface ), intent(out), allocatable process,
    integer, intent(out) rc
) 
```



This factory function creates and returns a new instance of the dust process. The process is not initialized - the caller must call the init() method with appropriate configuration.




**Parameters:**


* `process` Allocated process instance 
* `rc` Return code 




        

<hr>



### function get\_dust\_default\_config 

_Get default configuration for dust process._ 
```Fortran
subroutine, public dustprocesscreator_mod::get_dust_default_config (
    character(len=*), intent(out) config_data
) 
```



This function returns a default configuration string that can be used to initialize the dust process with reasonable defaults.




**Parameters:**


* `config_data` Default configuration string 




        

<hr>



### function register\_dust\_process 

_Register the dust process with a ProcessManager._ 
```Fortran
subroutine, public dustprocesscreator_mod::register_dust_process (
    type( processmanagertype ), intent(inout) process_mgr,
    integer, intent(out) rc
) 
```



This subroutine registers the dust process with a ProcessManager's factory. This is the correct way to register processes for use in applications and integration tests.




**Parameters:**


* `process_mgr` The ProcessManager to register with 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/dust/DustProcessCreator_Mod.F90`

