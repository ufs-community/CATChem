

# Namespace gridmanager\_mod



[**Namespace List**](namespaces.md) **>** [**gridmanager\_mod**](namespacegridmanager__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  integer, parameter, public | [**coord\_cartesian**](#variable-coord_cartesian)   = `1`<br>_Cartesian coordinates._  |
|  integer, parameter, public | [**coord\_lonlat**](#variable-coord_lonlat)   = `2`<br>_Lon/lat coordinates._  |
|  integer, parameter, public | [**coord\_projected**](#variable-coord_projected)   = `3`<br>_Projected coordinates._  |
|  integer, parameter, public | [**grid\_type\_2d**](#variable-grid_type_2d)   = `2`<br>_2D (x-z) model_  |
|  integer, parameter, public | [**grid\_type\_3d**](#variable-grid_type_3d)   = `3`<br>_3D (x-y-z) model_  |
|  integer, parameter, public | [**grid\_type\_column**](#variable-grid_type_column)   = `1`<br>_Pure column model._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**geometry\_init**](#function-geometry_init) (class([**gridgeometrytype**](namespacegridmanager__mod.md#none-gridgeometrytype)), intent(inout) this, integer, intent(in) nx, integer, intent(in) ny, integer, intent(in) nz, integer, intent(in), optional grid\_type, integer, intent(in), optional coord\_system, integer, intent(out) rc) <br>_Initialize grid geometry._  |




























## Public Attributes Documentation




### variable coord\_cartesian 

_Cartesian coordinates._ 
```Fortran
integer, parameter, public gridmanager_mod::coord_cartesian;
```




<hr>



### variable coord\_lonlat 

_Lon/lat coordinates._ 
```Fortran
integer, parameter, public gridmanager_mod::coord_lonlat;
```




<hr>



### variable coord\_projected 

_Projected coordinates._ 
```Fortran
integer, parameter, public gridmanager_mod::coord_projected;
```




<hr>



### variable grid\_type\_2d 

_2D (x-z) model_ 
```Fortran
integer, parameter, public gridmanager_mod::grid_type_2d;
```




<hr>



### variable grid\_type\_3d 

_3D (x-y-z) model_ 
```Fortran
integer, parameter, public gridmanager_mod::grid_type_3d;
```




<hr>



### variable grid\_type\_column 

_Pure column model._ 
```Fortran
integer, parameter, public gridmanager_mod::grid_type_column;
```




<hr>
## Public Functions Documentation




### function geometry\_init 

_Initialize grid geometry._ 
```Fortran
subroutine gridmanager_mod::geometry_init (
    class( gridgeometrytype ), intent(inout) this,
    integer, intent(in) nx,
    integer, intent(in) ny,
    integer, intent(in) nz,
    integer, intent(in), optional grid_type,
    integer, intent(in), optional coord_system,
    integer, intent(out) rc
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/GridManager_Mod.F90`

