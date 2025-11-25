

# Interface precision\_mod::rae



[**ClassList**](annotated.md) **>** [**precision\_mod**](namespaceprecision__mod.md) **>** [**rae**](interfaceprecision__mod_1_1rae.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  logical function | [**rae\_f4**](#function-rae_f4) (real([**f4**](namespaceprecision__mod.md#variable-f4)), intent(in) a, real([**f4**](namespaceprecision__mod.md#variable-f4)), intent(in) b) <br>_Real approximately equal:_ `abs(a - b) < tiny(a)` __ |
|  logical function | [**rae\_f8**](#function-rae_f8) (real([**f8**](namespaceprecision__mod.md#variable-f8)), intent(in) a, real([**f8**](namespaceprecision__mod.md#variable-f8)), intent(in) b) <br>_Real approximately equal:_ `abs(a - b) < tiny(a)` __ |




























## Public Functions Documentation




### function rae\_f4 

_Real approximately equal:_ `abs(a - b) < tiny(a)` __
```Fortran
logical function precision_mod::rae::rae_f4 (
    real( f4 ), intent(in) a,
    real( f4 ), intent(in) b
) 
```




<hr>



### function rae\_f8 

_Real approximately equal:_ `abs(a - b) < tiny(a)` __
```Fortran
logical function precision_mod::rae::rae_f8 (
    real( f8 ), intent(in) a,
    real( f8 ), intent(in) b
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/core/Precision_Mod.F90`

