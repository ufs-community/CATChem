

# Namespace catchem\_bridge\_error



[**Namespace List**](namespaces.md) **>** [**catchem\_bridge\_error**](namespacecatchem__bridge__error.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**category\_computation**](#variable-category_computation)   = `2`<br> |
|  integer, parameter, public | [**category\_general**](#variable-category_general)   = `0`<br> |
|  integer, parameter, public | [**category\_input**](#variable-category_input)   = `1`<br> |
|  integer, parameter, public | [**category\_io**](#variable-category_io)   = `4`<br> |
|  integer, parameter, public | [**category\_memory**](#variable-category_memory)   = `3`<br> |
|  integer, parameter, public | [**category\_mpi**](#variable-category_mpi)   = `5`<br> |
|  integer, parameter, public | [**category\_process**](#variable-category_process)   = `6`<br> |
|  integer, parameter, public | [**cc\_failure**](#variable-cc_failure)   = `-1`<br> |
|  integer, parameter, public | [**cc\_success**](#variable-cc_success)   = `0`<br> |
|  integer, parameter, public | [**error\_bounds\_check**](#variable-error_bounds_check)   = `1010`<br> |
|  integer, parameter, public | [**error\_convergence**](#variable-error_convergence)   = `1011`<br> |
|  integer, parameter, public | [**error\_dimension\_mismatch**](#variable-error_dimension_mismatch)   = `1009`<br> |
|  integer, parameter, public | [**error\_duplicate\_entry**](#variable-error_duplicate_entry)   = `1017`<br> |
|  integer, parameter, public | [**error\_file\_not\_found**](#variable-error_file_not_found)   = `1004`<br> |
|  integer, parameter, public | [**error\_file\_read**](#variable-error_file_read)   = `1005`<br> |
|  integer, parameter, public | [**error\_file\_write**](#variable-error_file_write)   = `1006`<br> |
|  integer, parameter, public | [**error\_invalid\_config**](#variable-error_invalid_config)   = `1002`<br> |
|  integer, parameter, public | [**error\_invalid\_input**](#variable-error_invalid_input)   = `1001`<br> |
|  integer, parameter, public | [**error\_invalid\_state**](#variable-error_invalid_state)   = `1003`<br> |
|  integer, parameter, public | [**error\_memory\_allocation**](#variable-error_memory_allocation)   = `1007`<br> |
|  integer, parameter, public | [**error\_memory\_deallocation**](#variable-error_memory_deallocation)   = `1008`<br> |
|  integer, parameter, public | [**error\_mpi\_communication**](#variable-error_mpi_communication)   = `1013`<br> |
|  integer, parameter, public | [**error\_none**](#variable-error_none)   = `0`<br> |
|  integer, parameter, public | [**error\_not\_found**](#variable-error_not_found)   = `1018`<br> |
|  integer, parameter, public | [**error\_numerical\_instability**](#variable-error_numerical_instability)   = `1012`<br> |
|  integer, parameter, public | [**error\_process\_initialization**](#variable-error_process_initialization)   = `1014`<br> |
|  integer, parameter, public | [**error\_state\_inconsistency**](#variable-error_state_inconsistency)   = `1015`<br> |
|  integer, parameter, public | [**error\_unsupported\_operation**](#variable-error_unsupported_operation)   = `1016`<br> |
|  integer, parameter, public | [**severity\_critical**](#variable-severity_critical)   = `3`<br> |
|  integer, parameter, public | [**severity\_error**](#variable-severity_error)   = `2`<br> |
|  integer, parameter, public | [**severity\_fatal**](#variable-severity_fatal)   = `4`<br> |
|  integer, parameter, public | [**severity\_info**](#variable-severity_info)   = `0`<br> |
|  integer, parameter, public | [**severity\_warning**](#variable-severity_warning)   = `1`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**cc\_checkvar**](#function-cc_checkvar) (character(len=\*), intent(in) variable, integer, intent(in) operation, integer, intent(inout) rc) <br>_CC\_CheckVar._  |
|  subroutine, public | [**cc\_error**](#function-cc_error) (character(len=\*), intent(in) errmsg, integer, intent(inout) rc, character(len=\*), intent(in), optional thisloc, character(len=\*), intent(in), optional instr) <br> |
|  subroutine, public | [**cc\_fatal**](#function-cc_fatal) (character(len=\*), intent(in) errmsg, character(len=\*), intent(in) thisloc) <br> |
|  subroutine, public | [**cc\_warning**](#function-cc_warning) (character(len=\*), intent(in) warnmsg, integer, intent(inout) rc, character(len=\*), intent(in), optional thisloc, character(len=\*), intent(in), optional instr) <br>_CC\_Warning._  |




























## Public Attributes Documentation




### variable category\_computation 

```Fortran
integer, parameter, public catchem_bridge_error::category_computation;
```




<hr>



### variable category\_general 

```Fortran
integer, parameter, public catchem_bridge_error::category_general;
```




<hr>



### variable category\_input 

```Fortran
integer, parameter, public catchem_bridge_error::category_input;
```




<hr>



### variable category\_io 

```Fortran
integer, parameter, public catchem_bridge_error::category_io;
```




<hr>



### variable category\_memory 

```Fortran
integer, parameter, public catchem_bridge_error::category_memory;
```




<hr>



### variable category\_mpi 

```Fortran
integer, parameter, public catchem_bridge_error::category_mpi;
```




<hr>



### variable category\_process 

```Fortran
integer, parameter, public catchem_bridge_error::category_process;
```




<hr>



### variable cc\_failure 

```Fortran
integer, parameter, public catchem_bridge_error::cc_failure;
```




<hr>



### variable cc\_success 

```Fortran
integer, parameter, public catchem_bridge_error::cc_success;
```




<hr>



### variable error\_bounds\_check 

```Fortran
integer, parameter, public catchem_bridge_error::error_bounds_check;
```




<hr>



### variable error\_convergence 

```Fortran
integer, parameter, public catchem_bridge_error::error_convergence;
```




<hr>



### variable error\_dimension\_mismatch 

```Fortran
integer, parameter, public catchem_bridge_error::error_dimension_mismatch;
```




<hr>



### variable error\_duplicate\_entry 

```Fortran
integer, parameter, public catchem_bridge_error::error_duplicate_entry;
```




<hr>



### variable error\_file\_not\_found 

```Fortran
integer, parameter, public catchem_bridge_error::error_file_not_found;
```




<hr>



### variable error\_file\_read 

```Fortran
integer, parameter, public catchem_bridge_error::error_file_read;
```




<hr>



### variable error\_file\_write 

```Fortran
integer, parameter, public catchem_bridge_error::error_file_write;
```




<hr>



### variable error\_invalid\_config 

```Fortran
integer, parameter, public catchem_bridge_error::error_invalid_config;
```




<hr>



### variable error\_invalid\_input 

```Fortran
integer, parameter, public catchem_bridge_error::error_invalid_input;
```




<hr>



### variable error\_invalid\_state 

```Fortran
integer, parameter, public catchem_bridge_error::error_invalid_state;
```




<hr>



### variable error\_memory\_allocation 

```Fortran
integer, parameter, public catchem_bridge_error::error_memory_allocation;
```




<hr>



### variable error\_memory\_deallocation 

```Fortran
integer, parameter, public catchem_bridge_error::error_memory_deallocation;
```




<hr>



### variable error\_mpi\_communication 

```Fortran
integer, parameter, public catchem_bridge_error::error_mpi_communication;
```




<hr>



### variable error\_none 

```Fortran
integer, parameter, public catchem_bridge_error::error_none;
```




<hr>



### variable error\_not\_found 

```Fortran
integer, parameter, public catchem_bridge_error::error_not_found;
```




<hr>



### variable error\_numerical\_instability 

```Fortran
integer, parameter, public catchem_bridge_error::error_numerical_instability;
```




<hr>



### variable error\_process\_initialization 

```Fortran
integer, parameter, public catchem_bridge_error::error_process_initialization;
```




<hr>



### variable error\_state\_inconsistency 

```Fortran
integer, parameter, public catchem_bridge_error::error_state_inconsistency;
```




<hr>



### variable error\_unsupported\_operation 

```Fortran
integer, parameter, public catchem_bridge_error::error_unsupported_operation;
```




<hr>



### variable severity\_critical 

```Fortran
integer, parameter, public catchem_bridge_error::severity_critical;
```




<hr>



### variable severity\_error 

```Fortran
integer, parameter, public catchem_bridge_error::severity_error;
```




<hr>



### variable severity\_fatal 

```Fortran
integer, parameter, public catchem_bridge_error::severity_fatal;
```




<hr>



### variable severity\_info 

```Fortran
integer, parameter, public catchem_bridge_error::severity_info;
```




<hr>



### variable severity\_warning 

```Fortran
integer, parameter, public catchem_bridge_error::severity_warning;
```




<hr>
## Public Functions Documentation




### function cc\_checkvar 

_CC\_CheckVar._ 
```Fortran
subroutine, public catchem_bridge_error::cc_checkvar (
    character(len=*), intent(in) variable,
    integer, intent(in) operation,
    integer, intent(inout) rc
) 
```



This subroutine checks if a variable is allocated.




**Parameters:**


* `Variable` The variable to check 
* `Operation` 0=Allocate 1=Register 2=Deallocate 
* `RC` The return code


>  



        

<hr>



### function cc\_error 

```Fortran
subroutine, public catchem_bridge_error::cc_error (
    character(len=*), intent(in) errmsg,
    integer, intent(inout) rc,
    character(len=*), intent(in), optional thisloc,
    character(len=*), intent(in), optional instr
) 
```




<hr>



### function cc\_fatal 

```Fortran
subroutine, public catchem_bridge_error::cc_fatal (
    character(len=*), intent(in) errmsg,
    character(len=*), intent(in) thisloc
) 
```




<hr>



### function cc\_warning 

_CC\_Warning._ 
```Fortran
subroutine, public catchem_bridge_error::cc_warning (
    character(len=*), intent(in) warnmsg,
    integer, intent(inout) rc,
    character(len=*), intent(in), optional thisloc,
    character(len=*), intent(in), optional instr
) 
```



This subroutine prints a warning message and sets RC to CC\_SUCCESS.




**Parameters:**


* `WarnMsg` The warning message 
* `RC` The return code 
* `ThisLoc` The location of the warning 
* `Instr` Other instructions 
>  





        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_bridge_utilities.F90`

