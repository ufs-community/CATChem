

# Namespace virtualcolumn\_mod



[**Namespace List**](namespaces.md) **>** [**virtualcolumn\_mod**](namespacevirtualcolumn__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**virtual\_met\_cleanup**](#function-virtual_met_cleanup) (class(virtualmettype), intent(inout) this) <br>_Clean up virtual met pointers._  |




























## Public Functions Documentation




### function virtual\_met\_cleanup 

_Clean up virtual met pointers._ 
```Fortran
subroutine virtualcolumn_mod::virtual_met_cleanup (
    class(virtualmettype), intent(inout) this
) 
```



Nullifies all pointers - does not deallocate since pointers point to MetState data which is managed elsewhere 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/VirtualColumn_Mod.F90`

