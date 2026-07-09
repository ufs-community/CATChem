
# File List

Here is a list of all files with brief descriptions:


* **dir** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md)     
    * **dir** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md)     
        * **file** [**CATChemAPI\_Mod.F90**](_c_a_t_chem_a_p_i___mod_8_f90.md) _High-level CATChem API for easy integration._     
        * **file** [**CATChemNetCDF\_Mod.F90**](_c_a_t_chem_net_c_d_f___mod_8_f90.md) _High-level NetCDF I/O interface for CATChem with MPI support._     
        * **file** [**CATChem\_API.F90**](_c_a_t_chem___a_p_i_8_f90.md) _Streamlined CATChem API for host model integration._     
        * **file** [**CATChem\_HighLevel\_API.F90**](_c_a_t_chem___high_level___a_p_i_8_f90.md) _High-Level CATChem API for easy integration into modeling systems._     
        * **file** [**FieldMapping\_Mod.F90**](_field_mapping___mod_8_f90.md) _Field mapping system for CATChem high-level API._     
        * **file** [**catchem.F90**](catchem_8_f90.md) _CATChem core data types and routines._     
        * **file** [**run\_mod.F90**](run__mod_8_f90.md) _Run module for CATChem atmospheric chemistry processes._     
    * **dir** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md)     
        * **file** [**CATChemCore\_Mod.F90**](_c_a_t_chem_core___mod_8_f90.md) _Central CATChem core framework with unified component management._     
        * **file** [**ChemSpeciesUtils\_Mod.F90**](_chem_species_utils___mod_8_f90.md) _Utility functions for chemical species access and manipulation._     
        * **file** [**ColumnInterface\_Mod.F90**](_column_interface___mod_8_f90.md)     
        * **file** [**ConfigManager\_Mod.F90**](_config_manager___mod_8_f90.md) _Enhanced configuration management for CATChem._     
        * **file** [**DiagnosticInterface\_Mod.F90**](_diagnostic_interface___mod_8_f90.md) _Dynamic diagnostic system interfaces and types._     
        * **file** [**DiagnosticManager\_Mod.F90**](_diagnostic_manager___mod_8_f90.md) _Central diagnostic manager integrating with CATChem framework._     
        * **file** [**EmissionConfigValidator\_Mod.F90**](_emission_config_validator___mod_8_f90.md) _Configuration validation for emission species mapping._     
        * **file** [**ExtEmisData\_Mod.F90**](_ext_emis_data___mod_8_f90.md) _Module for external emission data storage._     
        * **file** [**GridGeometry\_Mod.F90**](_grid_geometry___mod_8_f90.md)     
        * **file** [**GridManager\_Mod.F90**](_grid_manager___mod_8_f90.md) _Advanced grid management with column virtualization support._     
        * **file** [**Precision\_Mod.F90**](_precision___mod_8_f90.md) _Module PRECISION\_MOD is used to change the precision of many variables throughout catchem at compile-time._     
        * **file** [**ProcessFactory\_Mod.F90**](_process_factory___mod_8_f90.md) _Process factory for dynamic process creation following architecture guide._     
        * **file** [**ProcessInterface\_Mod.F90**](_process_interface___mod_8_f90.md) _Abstract base class interface for all atmospheric processes._     
        * **file** [**ProcessManager\_Mod.F90**](_process_manager___mod_8_f90.md) _High-level process management following the architecture guide._     
        * **file** [**ProcessRegistry\_Mod.F90**](_process_registry___mod_8_f90.md) _Process registration system for dynamic process discovery._     
        * **file** [**StateManager\_Mod.F90**](_state_manager___mod_8_f90.md) _Unified state management module for CATChem._     
        * **file** [**TimeState\_Mod.F90**](_time_state___mod_8_f90.md) _Time state and common time/solar functions for atmospheric chemistry._     
        * **file** [**UnitConversion\_Mod.F90**](_unit_conversion___mod_8_f90.md) _Comprehensive unit conversion utilities for atmospheric chemistry._     
        * **file** [**VirtualColumn\_Mod.F90**](_virtual_column___mod_8_f90.md) _Virtual column data container for CATChem processes with macro-generated meteorological fields._     
        * **file** [**chemstate\_mod.F90**](chemstate__mod_8_f90.md) _Contains the_ `ChemStateType` _data type and related subroutines and functions._    
        * **file** [**constants.F90**](constants_8_f90.md) _Physical and mathematical constants for CATChem._     
        * **file** [**error\_mod.F90**](error__mod_8_f90.md)     
        * **file** [**init\_mod.F90**](init__mod_8_f90.md)     
        * **file** [**met\_utilities\_mod.F90**](met__utilities__mod_8_f90.md) _Meteorological utility functions for CATChem._     
        * **file** [**metstate\_mod.F90**](metstate__mod_8_f90.md)     
        * **file** [**species\_mod.F90**](species__mod_8_f90.md) _Modern species definition and management for CATChem._     
        * **file** [**state\_interface\_mod.F90**](state__interface__mod_8_f90.md) _State interface module for CATChem._     
        * **file** [**utilities\_mod.F90**](utilities__mod_8_f90.md) _General utility functions for CATChem._     
    * **dir** [**external**](dir_805a0af995e93a362739e98abd740eb2.md)     
        * **dir** [**yaml\_interface**](dir_d0b1a67acd809cff502adc02c61e9ebd.md)     
            * **file** [**yaml\_interface\_mod.F90**](yaml__interface__mod_8_f90.md) _High-level Fortran interface for yaml-cpp._     
    * **dir** [**process**](dir_c0cd66d8ddae4fc5bc5dc2f24e29763b.md)     
        * **dir** [**carbchem**](dir_5dbdd03f815becc4c35f94e0692b4e09.md)     
            * **file** [**CarbChemCommon\_Mod.F90**](_carb_chem_common___mod_8_f90.md) _Common types and utilities for carbchem process._     
            * **file** [**CarbChemProcessCreator\_Mod.F90**](_carb_chem_process_creator___mod_8_f90.md) _Factory for creating carbchem process instances._     
            * **file** [**ProcessCarbChemInterface\_Mod.F90**](_process_carb_chem_interface___mod_8_f90.md)     
            * **dir** [**schemes**](dir_2aa0506a6ee0350ff25dc7952f19a30f.md)     
                * **file** [**CarbChemScheme\_GOCART\_Mod.F90**](_carb_chem_scheme___g_o_c_a_r_t___mod_8_f90.md) _GOCART carbon species chemical production and loss scheme._     
        * **dir** [**drydep**](dir_57fb5aa14ddb2cd518a6d90b65ffd000.md)     
            * **file** [**DryDepCommon\_Mod.F90**](_dry_dep_common___mod_8_f90.md) _Common types and utilities for drydep process._     
            * **file** [**DryDepProcessCreator\_Mod.F90**](_dry_dep_process_creator___mod_8_f90.md) _Factory for creating drydep process instances._     
            * **file** [**ProcessDryDepInterface\_Mod.F90**](_process_dry_dep_interface___mod_8_f90.md)     
            * **dir** [**schemes**](dir_5a3c86e36f17958630366ebc2b7ca21b.md)     
                * **file** [**DryDepScheme\_GOCART\_Mod.F90**](_dry_dep_scheme___g_o_c_a_r_t___mod_8_f90.md) _GOCART-2G aerosol dry deposition scheme._     
                * **file** [**DryDepScheme\_WESELY\_Mod.F90**](_dry_dep_scheme___w_e_s_e_l_y___mod_8_f90.md) _Wesely 1989 gas dry deposition scheme._     
                * **file** [**DryDepScheme\_ZHANG\_Mod.F90**](_dry_dep_scheme___z_h_a_n_g___mod_8_f90.md) _Zhang et al._     
        * **dir** [**dust**](dir_1c14dfbaca1e3f4c2e26e74290119ebd.md)     
            * **file** [**DustCommon\_Mod.F90**](_dust_common___mod_8_f90.md) _Common types and utilities for dust process._     
            * **file** [**DustProcessCreator\_Mod.F90**](_dust_process_creator___mod_8_f90.md) _Factory for creating dust process instances._     
            * **file** [**ProcessDustInterface\_Mod.F90**](_process_dust_interface___mod_8_f90.md)     
            * **dir** [**schemes**](dir_11b2254edcf6ee5df673de29b129f986.md)     
                * **file** [**DustScheme\_FENGSHA\_Mod.F90**](_dust_scheme___f_e_n_g_s_h_a___mod_8_f90.md) _Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS._     
                * **file** [**DustScheme\_GINOUX\_Mod.F90**](_dust_scheme___g_i_n_o_u_x___mod_8_f90.md) _Ginoux dust emission scheme._     
        * **dir** [**seasalt**](dir_3e6fd2e121e43ca7d4114b6c0b4e05b6.md)     
            * **file** [**ProcessSeaSaltInterface\_Mod.F90**](_process_sea_salt_interface___mod_8_f90.md)     
            * **file** [**SeaSaltCommon\_Mod.F90**](_sea_salt_common___mod_8_f90.md) _Common types and utilities for seasalt process._     
            * **file** [**SeaSaltProcessCreator\_Mod.F90**](_sea_salt_process_creator___mod_8_f90.md) _Factory for creating seasalt process instances._     
            * **dir** [**schemes**](dir_ec083b49fedbd640552af85049fd7226.md)     
                * **file** [**SeaSaltScheme\_GEOS12\_Mod.F90**](_sea_salt_scheme___g_e_o_s12___mod_8_f90.md) _GEOS-Chem 2012 sea salt emission scheme with observational constraints._     
                * **file** [**SeaSaltScheme\_GONG03\_Mod.F90**](_sea_salt_scheme___g_o_n_g03___mod_8_f90.md) _Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment._     
                * **file** [**SeaSaltScheme\_GONG97\_Mod.F90**](_sea_salt_scheme___g_o_n_g97___mod_8_f90.md) _Gong 1997 sea salt emission scheme._     
        * **dir** [**settling**](dir_1a0bba2ffdf6e6637fcb76856471cb75.md)     
            * **file** [**ProcessSettlingInterface\_Mod.F90**](_process_settling_interface___mod_8_f90.md)     
            * **file** [**SettlingCommon\_Mod.F90**](_settling_common___mod_8_f90.md) _Common types and utilities for settling process._     
            * **file** [**SettlingProcessCreator\_Mod.F90**](_settling_process_creator___mod_8_f90.md) _Factory for creating settling process instances._     
            * **dir** [**schemes**](dir_34df91cc26d24067840a7381fe21b817.md)     
                * **file** [**SettlingScheme\_GOCART\_Mod.F90**](_settling_scheme___g_o_c_a_r_t___mod_8_f90.md) _GOCART gravitational settling scheme._     
        * **dir** [**so4chem**](dir_fb8fc0df5ebe1b02f5e46b98d91cbc63.md)     
            * **file** [**ProcessSO4chemInterface\_Mod.F90**](_process_s_o4chem_interface___mod_8_f90.md)     
            * **file** [**SO4chemCommon\_Mod.F90**](_s_o4chem_common___mod_8_f90.md) _Common types and utilities for so4chem process._     
            * **file** [**SO4chemProcessCreator\_Mod.F90**](_s_o4chem_process_creator___mod_8_f90.md) _Factory for creating so4chem process instances._     
            * **dir** [**schemes**](dir_429bccfa51a729cf5e11bef5cf73cbc7.md)     
                * **file** [**SO4chemScheme\_GOCART\_Mod.F90**](_s_o4chem_scheme___g_o_c_a_r_t___mod_8_f90.md) _GOCART SO2 to SO4 production scheme._     
        * **dir** [**wetdep**](dir_8b9a0ce556ea4a65f6920dfb49dcd69d.md)     
            * **file** [**ProcessWetDepInterface\_Mod.F90**](_process_wet_dep_interface___mod_8_f90.md)     
            * **file** [**WetDepCommon\_Mod.F90**](_wet_dep_common___mod_8_f90.md) _Common types and utilities for wetdep process._     
            * **file** [**WetDepProcessCreator\_Mod.F90**](_wet_dep_process_creator___mod_8_f90.md) _Factory for creating wetdep process instances._     
            * **dir** [**schemes**](dir_8ca87c5e2f5cf830ab1a41055168a46b.md)     
                * **file** [**WetDepScheme\_JACOB\_Mod.F90**](_wet_dep_scheme___j_a_c_o_b___mod_8_f90.md) _Jacob et al._     

