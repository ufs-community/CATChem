

# Namespace run\_mod



[**Namespace List**](namespaces.md) **>** [**run\_mod**](namespacerun__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**finalize\_process**](#function-finalize_process) (type(duststatetype), intent(inout) duststate, type(seasaltstatetype), intent(inout) seasaltstate, type(drydepstatetype), intent(inout) drydepstate, integer, intent(out) rc) <br>_Finalize all CATChem atmospheric chemistry processes._  |
|  subroutine, public | [**init\_process**](#function-init_process) (type(configtype), intent(in) config, type([**chemstatetype**](namespacechemstate__mod.md#none-chemstatetype)), intent(inout) chemstate, type(emisstatetype), intent(inout) emisstate, type(duststatetype), intent(inout) duststate, type(seasaltstatetype), intent(inout) seasaltstate, type(drydepstatetype), intent(inout) drydepstate, integer, intent(inout) rc) <br>_Initialize all CATChem atmospheric chemistry processes._  |
|  subroutine, public | [**run\_process**](#function-run_process) (type([**metstatetype**](namespacemetstate__mod.md#none-metstatetype)), intent(in) metstate, type(diagstatetype), intent(inout) diagstate, type([**chemstatetype**](namespacechemstate__mod.md#none-chemstatetype)), intent(inout) chemstate, type(emisstatetype), intent(inout) emisstate, type(duststatetype), intent(inout) duststate, type(seasaltstatetype), intent(inout) seasaltstate, type(drydepstatetype), intent(inout) drydepstate, integer, intent(out) rc) <br>_Run all CATChem atmospheric chemistry processes._  |




























## Public Functions Documentation




### function finalize\_process 

_Finalize all CATChem atmospheric chemistry processes._ 
```Fortran
subroutine, public run_mod::finalize_process (
    type(duststatetype), intent(inout) duststate,
    type(seasaltstatetype), intent(inout) seasaltstate,
    type(drydepstatetype), intent(inout) drydepstate,
    integer, intent(out) rc
) 
```



This subroutine performs cleanup and finalization for all atmospheric chemistry processes, deallocating memory and closing any open resources.




**Parameters:**


* `DustState` The dust process state to finalize 
* `SeaSaltState` The sea salt process state to finalize 
* `DryDepState` The dry deposition process state to finalize 
* `RC` The return code indicating success (CC\_SUCCESS) or failure

The finalization includes:
* Cleanup of dust emission process state
* Cleanup of sea salt emission process state
* Cleanup of dry deposition process state
* Cleanup of wet deposition process state (planned)






**Note:**

This routine should be called at the end of simulation 




**Warning:**

RC should be checked after calling this routine 





        

<hr>



### function init\_process 

_Initialize all CATChem atmospheric chemistry processes._ 
```Fortran
subroutine, public run_mod::init_process (
    type(configtype), intent(in) config,
    type( chemstatetype ), intent(inout) chemstate,
    type(emisstatetype), intent(inout) emisstate,
    type(duststatetype), intent(inout) duststate,
    type(seasaltstatetype), intent(inout) seasaltstate,
    type(drydepstatetype), intent(inout) drydepstate,
    integer, intent(inout) rc
) 
```



This subroutine performs initialization of all atmospheric chemistry processes including dust emission, sea salt emission, and dry deposition. It must be called before any process execution routines.




**Parameters:**


* `config` The CATChem configuration object containing setup parameters 
* `ChemState` The chemical state to be initialized with species data 
* `EmisState` The emission state to be initialized for emission processes 
* `DustState` The dust process state to be initialized 
* `SeaSaltState` The sea salt process state to be initialized 
* `DryDepState` The dry deposition process state to be initialized 
* `RC` The return code indicating success (CC\_SUCCESS) or failure

The initialization includes:
* Setting up dust emission process parameters and lookup tables
* Configuring sea salt emission parameterizations
* Preparing dry deposition velocity calculations
* Validating process configuration consistency






**Note:**

This routine must be called after state allocation but before process execution 




**Warning:**

RC should be checked after calling this routine 





        

<hr>



### function run\_process 

_Run all CATChem atmospheric chemistry processes._ 
```Fortran
subroutine, public run_mod::run_process (
    type( metstatetype ), intent(in) metstate,
    type(diagstatetype), intent(inout) diagstate,
    type( chemstatetype ), intent(inout) chemstate,
    type(emisstatetype), intent(inout) emisstate,
    type(duststatetype), intent(inout) duststate,
    type(seasaltstatetype), intent(inout) seasaltstate,
    type(drydepstatetype), intent(inout) drydepstate,
    integer, intent(out) rc
) 
```



This subroutine executes all enabled atmospheric chemistry processes including emissions, dry deposition, and other chemical processes for a single time step.




**Parameters:**


* `MetState` The meteorological state containing atmospheric conditions 
* `DiagState` The diagnostic state for storing process outputs 
* `ChemState` The chemical state containing species concentrations 
* `EmisState` The emission state containing emission data 
* `DustState` The dust process state 
* `SeaSaltState` The sea salt process state 
* `DryDepState` The dry deposition process state 
* `RC` The return code indicating success (CC\_SUCCESS) or failure

The processes are executed in the following order:
* Emission processes (dust, sea salt, etc.)
* Dry deposition calculations
* Wet deposition calculations (planned)






**Note:**

All processes update the ChemState and DiagState objects 




**Warning:**

RC should be checked after calling this routine 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/run_mod.F90`

