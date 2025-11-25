

# Namespace drydepprocesscreator\_mod



[**Namespace List**](namespaces.md) **>** [**drydepprocesscreator\_mod**](namespacedrydepprocesscreator__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**create\_drydep\_process**](#function-create_drydep_process) (class([**processinterface**](namespaceprocessinterface__mod.md#none-processinterface)), intent(out), allocatable process, integer, intent(out) rc) <br>_Create a new drydep process instance._  |
|  subroutine, public | [**get\_drydep\_default\_config**](#function-get_drydep_default_config) (character(len=\*), intent(out) config\_data) <br>_Get default configuration for drydep process._  |
|  subroutine, public | [**register\_drydep\_process**](#function-register_drydep_process) (type([**processmanagertype**](namespaceprocessmanager__mod.md#none-processmanagertype)), intent(inout) process\_mgr, integer, intent(out) rc) <br>_Register the drydep process with a ProcessManager._  |




























## Public Functions Documentation




### function create\_drydep\_process 

_Create a new drydep process instance._ 
```Fortran
subroutine, public drydepprocesscreator_mod::create_drydep_process (
    class( processinterface ), intent(out), allocatable process,
    integer, intent(out) rc
) 
```



This factory function creates and returns a new instance of the drydep process. The process is not initialized - the caller must call the init() method with appropriate configuration.




**Parameters:**


* `process` Allocated process instance 
* `rc` Return code 




        

<hr>



### function get\_drydep\_default\_config 

_Get default configuration for drydep process._ 
```Fortran
subroutine, public drydepprocesscreator_mod::get_drydep_default_config (
    character(len=*), intent(out) config_data
) 
```



This function returns a default configuration string that can be used to initialize the drydep process with reasonable defaults.




**Parameters:**


* `config_data` Default configuration string 




        

<hr>



### function register\_drydep\_process 

_Register the drydep process with a ProcessManager._ 
```Fortran
subroutine, public drydepprocesscreator_mod::register_drydep_process (
    type( processmanagertype ), intent(inout) process_mgr,
    integer, intent(out) rc
) 
```



This subroutine registers the drydep process with a ProcessManager's factory. This is the correct way to register processes for use in applications and integration tests.




**Parameters:**


* `process_mgr` The ProcessManager to register with 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/drydep/DryDepProcessCreator_Mod.F90`

