

# Group core\_modules



[**Modules**](modules.md) **>** [**core\_modules**](group__core__modules.md)



_Core modules and data types for CATChem._ [More...](#detailed-description)








## Files

| Type | Name |
| ---: | :--- |
| file | [**CATChemCore\_Mod.F90**](_c_a_t_chem_core___mod_8_f90.md) <br>_Central CATChem core framework with unified component management._  |
| file | [**ChemSpeciesUtils\_Mod.F90**](_chem_species_utils___mod_8_f90.md) <br>_Utility functions for chemical species access and manipulation._  |
| file | [**ConfigManager\_Mod.F90**](_config_manager___mod_8_f90.md) <br>_Enhanced configuration management for CATChem._  |
| file | [**DiagnosticInterface\_Mod.F90**](_diagnostic_interface___mod_8_f90.md) <br>_Dynamic diagnostic system interfaces and types._  |
| file | [**DiagnosticManager\_Mod.F90**](_diagnostic_manager___mod_8_f90.md) <br>_Central diagnostic manager integrating with CATChem framework._  |
| file | [**ExtEmisData\_Mod.F90**](_ext_emis_data___mod_8_f90.md) <br>_Module for external emission data storage._  |
| file | [**GridManager\_Mod.F90**](_grid_manager___mod_8_f90.md) <br>_Advanced grid management with column virtualization support._  |
| file | [**Precision\_Mod.F90**](_precision___mod_8_f90.md) <br>_Module PRECISION\_MOD is used to change the precision of many variables throughout catchem at compile-time._  |
| file | [**StateManager\_Mod.F90**](_state_manager___mod_8_f90.md) <br>_Unified state management module for CATChem._  |
| file | [**UnitConversion\_Mod.F90**](_unit_conversion___mod_8_f90.md) <br>_Comprehensive unit conversion utilities for atmospheric chemistry._  |
| file | [**VirtualColumn\_Mod.F90**](_virtual_column___mod_8_f90.md) <br>_Virtual column data container for CATChem processes with macro-generated meteorological fields._  |
| file | [**chemstate\_mod.F90**](chemstate__mod_8_f90.md) <br>_Contains the_ `ChemStateType` _data type and related subroutines and functions._ |
| file | [**constants.F90**](constants_8_f90.md) <br>_Physical and mathematical constants for CATChem._  |
| file | [**met\_utilities\_mod.F90**](met__utilities__mod_8_f90.md) <br>_Meteorological utility functions for CATChem._  |
| file | [**run\_mod.F90**](run__mod_8_f90.md) <br>_Run module for CATChem atmospheric chemistry processes._  |
| file | [**species\_mod.F90**](species__mod_8_f90.md) <br>_Modern species definition and management for CATChem._  |
| file | [**state\_interface\_mod.F90**](state__interface__mod_8_f90.md) <br>_State interface module for CATChem._  |
| file | [**utilities\_mod.F90**](utilities__mod_8_f90.md) <br>_General utility functions for CATChem._  |






























## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**cc\_checkvar**](#function-cc_checkvar) (character(len=\*), intent(in) variable, integer, intent(in) operation, integer, intent(inout) rc) <br>_CC\_CheckVar._  |
|  subroutine, public | [**cc\_warning**](#function-cc_warning) (character(len=\*), intent(in) warnmsg, integer, intent(inout) rc, character(len=\*), intent(in), optional thisloc, character(len=\*), intent(in), optional instr) <br>_CC\_Warning._  |
|  subroutine, public | [**find\_number\_of\_species**](#function-find_number_of_species) (type(chemstatetype), intent(inout) chemstate, integer, intent(out) rc) <br>_Find the number of species._  |




























## Detailed Description


This group contains all core modules, state definitions, and fundamental data types used throughout the CATChem system. 


    
## Public Functions Documentation




### function cc\_checkvar 

_CC\_CheckVar._ 
```
subroutine, public cc_checkvar (
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



### function cc\_warning 

_CC\_Warning._ 
```
subroutine, public cc_warning (
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



### function find\_number\_of\_species 

_Find the number of species._ 
```
subroutine, public find_number_of_species (
    type(chemstatetype), intent(inout) chemstate,
    integer, intent(out) rc
) 
```



This subroutine finds the number of species




**Parameters:**


* `ChemState` The ChemState object 
* `RC` The return code


>  



        

<hr>

------------------------------


