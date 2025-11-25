

# Namespace init\_mod



[**Namespace List**](namespaces.md) **>** [**init\_mod**](namespaceinit__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**finalize\_catchem**](#function-finalize_catchem) (type(statecontainertype), intent(inout) container, integer, intent(out) rc) <br>_Finalize CATChem and cleanup resources._  |
|  subroutine, public | [**initialize\_catchem**](#function-initialize_catchem) (type(statecontainertype), intent(out) container, character(len=\*), intent(in) config\_file, integer, intent(out) rc, character(len=\*), intent(in), optional container\_name) <br>_Main CATChem initialization routine._  |
|  subroutine, public | [**initialize\_core\_states**](#function-initialize_core_states) (type(statecontainertype), intent(inout) container, integer, intent(out) rc) <br>_Initialize core state objects._  |
|  subroutine, public | [**initialize\_processes**](#function-initialize_processes) (type(statecontainertype), intent(inout) container, integer, intent(out) rc) <br>_Initialize atmospheric chemistry processes._  |
|  subroutine, public | [**validate\_initialization**](#function-validate_initialization) (type(statecontainertype), intent(inout) container, integer, intent(out) rc) <br>_Validate complete initialization._  |




























## Public Functions Documentation




### function finalize\_catchem 

_Finalize CATChem and cleanup resources._ 
```Fortran
subroutine, public init_mod::finalize_catchem (
    type(statecontainertype), intent(inout) container,
    integer, intent(out) rc
) 
```



This routine performs cleanup of all CATChem components and resources.




**Parameters:**


* `container` StateContainer to finalize 
* `rc` Return code 




        

<hr>



### function initialize\_catchem 

_Main CATChem initialization routine._ 
```Fortran
subroutine, public init_mod::initialize_catchem (
    type(statecontainertype), intent(out) container,
    character(len=*), intent(in) config_file,
    integer, intent(out) rc,
    character(len=*), intent(in), optional container_name
) 
```



This is the primary initialization routine that sets up a complete CATChem instance using the StateContainer pattern with comprehensive error handling.




**Parameters:**


* `container` Initialized StateContainer with all components 
* `config_file` Path to configuration file 
* `container_name` Optional name for the container 
* `rc` Return code 




        

<hr>



### function initialize\_core\_states 

_Initialize core state objects._ 
```Fortran
subroutine, public init_mod::initialize_core_states (
    type(statecontainertype), intent(inout) container,
    integer, intent(out) rc
) 
```



This routine initializes all core state objects (Met, Chem, Emis, Diag) within the StateContainer.




**Parameters:**


* `container` StateContainer to initialize 
* `rc` Return code 




        

<hr>



### function initialize\_processes 

_Initialize atmospheric chemistry processes._ 
```Fortran
subroutine, public init_mod::initialize_processes (
    type(statecontainertype), intent(inout) container,
    integer, intent(out) rc
) 
```



This routine initializes all atmospheric chemistry processes (dust, seasalt, dry deposition, etc.) based on the loaded configuration.




**Parameters:**


* `container` StateContainer with processes to initialize 
* `rc` Return code 




        

<hr>



### function validate\_initialization 

_Validate complete initialization._ 
```Fortran
subroutine, public init_mod::validate_initialization (
    type(statecontainertype), intent(inout) container,
    integer, intent(out) rc
) 
```



This routine performs comprehensive validation of the initialized StateContainer to ensure all components are properly set up and consistent.




**Parameters:**


* `container` StateContainer to validate 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/init_mod.F90`

