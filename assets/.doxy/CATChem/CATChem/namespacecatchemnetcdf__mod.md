

# Namespace catchemnetcdf\_mod



[**Namespace List**](namespaces.md) **>** [**catchemnetcdf\_mod**](namespacecatchemnetcdf__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**netcdf\_append**](#variable-netcdf_append)   = `3`<br> |
|  integer, parameter, public | [**netcdf\_read**](#variable-netcdf_read)   = `1`<br> |
|  integer, parameter, public | [**netcdf\_write**](#variable-netcdf_write)   = `2`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  logical function, public | [**mpi\_netcdf\_available**](#function-mpi_netcdf_available) () <br>_Check if parallel NetCDF with MPI support is available._  |
|  logical function, public | [**netcdf\_available**](#function-netcdf_available) () <br>_Check if NetCDF support is available._  |




























## Public Attributes Documentation




### variable netcdf\_append 

```Fortran
integer, parameter, public catchemnetcdf_mod::netcdf_append;
```




<hr>



### variable netcdf\_read 

```Fortran
integer, parameter, public catchemnetcdf_mod::netcdf_read;
```




<hr>



### variable netcdf\_write 

```Fortran
integer, parameter, public catchemnetcdf_mod::netcdf_write;
```




<hr>
## Public Functions Documentation




### function mpi\_netcdf\_available 

_Check if parallel NetCDF with MPI support is available._ 
```Fortran
logical function, public catchemnetcdf_mod::mpi_netcdf_available () 
```




<hr>



### function netcdf\_available 

_Check if NetCDF support is available._ 
```Fortran
logical function, public catchemnetcdf_mod::netcdf_available () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/api/CATChemNetCDF_Mod.F90`

