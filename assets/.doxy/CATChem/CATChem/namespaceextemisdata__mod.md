

# Namespace extemisdata\_mod



[**Namespace List**](namespaces.md) **>** [**extemisdata\_mod**](namespaceextemisdata__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**extemifield\_init**](#function-extemifield_init) (class([**extemisfieldtype**](namespaceextemisdata__mod.md#none-extemisfieldtype)), intent(inout) this, character(len=\*), intent(in) field\_name, integer, intent(in) nx, integer, intent(in) ny, integer, intent(in), optional nz, integer, intent(in), optional n\_times, character(len=\*), intent(in), optional units, integer, intent(out) rc) <br>_Initialize an emission field with metadata._  |




























## Public Functions Documentation




### function extemifield\_init 

_Initialize an emission field with metadata._ 
```Fortran
subroutine extemisdata_mod::extemifield_init (
    class( extemisfieldtype ), intent(inout) this,
    character(len=*), intent(in) field_name,
    integer, intent(in) nx,
    integer, intent(in) ny,
    integer, intent(in), optional nz,
    integer, intent(in), optional n_times,
    character(len=*), intent(in), optional units,
    integer, intent(out) rc
) 
```



Sets up field metadata including dimensions and coordinates




**Parameters:**


* `this` The ExtEmisFieldType object 
* `field_name` Field name 
* `nx` Number of longitude points 
* `ny` Number of latitude points 
* `nz` Number of vertical levels (optional, default=1) 
* `n_times` Number of time steps (optional, default=1) 
* `units` Units string (optional, default='kg/m2/s') 
* `rc` Return code 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/ExtEmisData_Mod.F90`

