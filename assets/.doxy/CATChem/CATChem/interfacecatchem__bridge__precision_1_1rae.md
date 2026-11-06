

# Interface catchem\_bridge\_precision::rae



[**ClassList**](annotated.md) **>** [**catchem\_bridge\_precision**](namespacecatchem__bridge__precision.md) **>** [**rae**](interfacecatchem__bridge__precision_1_1rae.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  logical function | [**rae\_f4**](#function-rae_f4) (real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) a, real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) b) <br>_Real approximately equal:_ `abs(a - b) < tiny(a)` __ |
|  logical function | [**rae\_f8**](#function-rae_f8) (real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) a, real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) b) <br>_Real approximately equal:_ `abs(a - b) < tiny(a)` __ |




























## Public Functions Documentation




### function rae\_f4 

_Real approximately equal:_ `abs(a - b) < tiny(a)` __
```Fortran
logical function catchem_bridge_precision::rae::rae_f4 (
    real( f4 ), intent(in) a,
    real( f4 ), intent(in) b
) 
```




<hr>



### function rae\_f8 

_Real approximately equal:_ `abs(a - b) < tiny(a)` __
```Fortran
logical function catchem_bridge_precision::rae::rae_f8 (
    real( f8 ), intent(in) a,
    real( f8 ), intent(in) b
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_bridge_utilities.F90`

