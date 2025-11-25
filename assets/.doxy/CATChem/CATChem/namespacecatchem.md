

# Namespace catchem



[**Namespace List**](namespaces.md) **>** [**catchem**](namespacecatchem.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**cc\_create\_virtual\_column**](#function-cc_create_virtual_column) (type(statecontainertype), intent(inout) container, integer, intent(in) col\_idx, type(virtualcolumntype), intent(out) virtual\_col, integer, intent(out) rc) <br>_Create a virtual column for processing._  |
|  subroutine, public | [**cc\_init\_grid\_manager**](#function-cc_init_grid_manager) (type(statecontainertype), intent(inout) container, integer, intent(in) nx, integer, intent(in) ny, integer, intent(in) nz, integer, intent(out) rc) <br>_Initialize grid manager with column virtualization._  |
|  subroutine, public | [**cc\_process\_all\_columns**](#function-cc_process_all_columns) (type([**processmanagertype**](namespaceprocessmanager__mod.md#none-processmanagertype)), intent(inout) proc\_mgr, character(len=\*), intent(in) process\_name, type(statecontainertype), intent(inout) container, integer, intent(out) rc) <br>_Process all columns with a specific process._  |
|  subroutine, public | [**cc\_run\_column\_processes**](#function-cc_run_column_processes) (type([**processmanagertype**](namespaceprocessmanager__mod.md#none-processmanagertype)), intent(inout) proc\_mgr, type(statecontainertype), intent(inout) container, integer, intent(out) rc) <br>_Run all column-based processes using process manager._  |




























## Public Functions Documentation




### function cc\_create\_virtual\_column 

_Create a virtual column for processing._ 
```Fortran
subroutine, public catchem::cc_create_virtual_column (
    type(statecontainertype), intent(inout) container,
    integer, intent(in) col_idx,
    type(virtualcolumntype), intent(out) virtual_col,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` State container 
* `col_idx` Column index 
* `virtual_col` Virtual column object 
* `rc` Return code 




        

<hr>



### function cc\_init\_grid\_manager 

_Initialize grid manager with column virtualization._ 
```Fortran
subroutine, public catchem::cc_init_grid_manager (
    type(statecontainertype), intent(inout) container,
    integer, intent(in) nx,
    integer, intent(in) ny,
    integer, intent(in) nz,
    integer, intent(out) rc
) 
```





**Parameters:**


* `container` State container 
* `nx` Number of grid points in x direction 
* `ny` Number of grid points in y direction 
* `nz` Number of grid points in z direction 
* `rc` Return code 




        

<hr>



### function cc\_process\_all\_columns 

_Process all columns with a specific process._ 
```Fortran
subroutine, public catchem::cc_process_all_columns (
    type( processmanagertype ), intent(inout) proc_mgr,
    character(len=*), intent(in) process_name,
    type(statecontainertype), intent(inout) container,
    integer, intent(out) rc
) 
```





**Parameters:**


* `proc_mgr` Process manager 
* `process_name` Name of process to run 
* `container` State container 
* `rc` Return code 




        

<hr>



### function cc\_run\_column\_processes 

_Run all column-based processes using process manager._ 
```Fortran
subroutine, public catchem::cc_run_column_processes (
    type( processmanagertype ), intent(inout) proc_mgr,
    type(statecontainertype), intent(inout) container,
    integer, intent(out) rc
) 
```





**Parameters:**


* `proc_mgr` Process manager 
* `container` State container 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/catchem.F90`

