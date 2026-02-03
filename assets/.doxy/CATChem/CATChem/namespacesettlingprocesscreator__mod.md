

# Namespace settlingprocesscreator\_mod



[**Namespace List**](namespaces.md) **>** [**settlingprocesscreator\_mod**](namespacesettlingprocesscreator__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**create\_settling\_process**](#function-create_settling_process) (class([**processinterface**](namespaceprocessinterface__mod.md#none-processinterface)), intent(out), allocatable process, integer, intent(out) rc) <br>_Create a new settling process instance._  |
|  subroutine, public | [**get\_settling\_default\_config**](#function-get_settling_default_config) (character(len=\*), intent(out) config\_data) <br>_Get default configuration for settling process._  |
|  subroutine, public | [**register\_settling\_process**](#function-register_settling_process) (type([**processmanagertype**](namespaceprocessmanager__mod.md#none-processmanagertype)), intent(inout) process\_mgr, integer, intent(out) rc) <br>_Register the settling process with a ProcessManager._  |




























## Public Functions Documentation




### function create\_settling\_process 

_Create a new settling process instance._ 
```Fortran
subroutine, public settlingprocesscreator_mod::create_settling_process (
    class( processinterface ), intent(out), allocatable process,
    integer, intent(out) rc
) 
```



This factory function creates and returns a new instance of the settling process. The process is not initialized - the caller must call the init() method with appropriate configuration.




**Parameters:**


* `process` Allocated process instance 
* `rc` Return code 




        

<hr>



### function get\_settling\_default\_config 

_Get default configuration for settling process._ 
```Fortran
subroutine, public settlingprocesscreator_mod::get_settling_default_config (
    character(len=*), intent(out) config_data
) 
```



This function returns a default configuration string that can be used to initialize the settling process with reasonable defaults.




**Parameters:**


* `config_data` Default configuration string 




        

<hr>



### function register\_settling\_process 

_Register the settling process with a ProcessManager._ 
```Fortran
subroutine, public settlingprocesscreator_mod::register_settling_process (
    type( processmanagertype ), intent(inout) process_mgr,
    integer, intent(out) rc
) 
```



This subroutine registers the settling process with a ProcessManager's factory. This is the correct way to register processes for use in applications and integration tests.




**Parameters:**


* `process_mgr` The ProcessManager to register with 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/settling/SettlingProcessCreator_Mod.F90`

