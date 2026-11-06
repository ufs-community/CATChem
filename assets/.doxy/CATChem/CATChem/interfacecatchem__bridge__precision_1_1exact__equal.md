

# Interface catchem\_bridge\_precision::exact\_equal



[**ClassList**](annotated.md) **>** [**catchem\_bridge\_precision**](namespacecatchem__bridge__precision.md) **>** [**exact\_equal**](interfacecatchem__bridge__precision_1_1exact__equal.md)



_Exact IEEE equality, spelled without == so gfortran's -Wcompare-reals stays quiet at deliberate exact comparisons (sentinels, literal zeros)._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  elemental logical function | [**exact\_equal\_f4**](#function-exact_equal_f4) (real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) a, real([**f4**](namespacecatchem__bridge__precision.md#variable-f4)), intent(in) b) <br>_Exact equality (see the_ `exact_equal` _interface)_ |
|  elemental logical function | [**exact\_equal\_f8**](#function-exact_equal_f8) (real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) a, real([**f8**](namespacecatchem__bridge__precision.md#variable-f8)), intent(in) b) <br>_Exact equality (see the_ `exact_equal` _interface)_ |




























## Detailed Description


(a &lt;= b .and. a &gt;= b) is precisely a == b: false for NaN, true for equal infinities and for +0/-0. Not a tolerance test; see `rae` for that. 


    
## Public Functions Documentation




### function exact\_equal\_f4 

_Exact equality (see the_ `exact_equal` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::exact_equal::exact_equal_f4 (
    real( f4 ), intent(in) a,
    real( f4 ), intent(in) b
) 
```




<hr>



### function exact\_equal\_f8 

_Exact equality (see the_ `exact_equal` _interface)_
```Fortran
elemental logical function catchem_bridge_precision::exact_equal::exact_equal_f8 (
    real( f8 ), intent(in) a,
    real( f8 ), intent(in) b
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/process/catchem_bridge_utilities.F90`

