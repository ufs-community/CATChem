!> \file
!! \brief CATChem core data types and routines.
!!
!! \defgroup catchem_api
!! CATChem core data types and routines.
!!
!!!>
module CATChem

   !---------------
   ! CATChem States
   !---------------
   use ChemState_Mod,  only: ChemStateType    !< Chemical State
   use GridState_Mod,  only: GridStateType    !< Grid State
   use DiagState_Mod,  only: DiagStateType    !< Diagnostic State
   use EmisState_Mod,  only: EmisStateType    !< Emission State
   use MetState_Mod,   only: MetStateType     !< Meteorology State
   use species_mod,    only: SpeciesType      !< Species State
   use Config_Opt_Mod, only: ConfigType

   !----------------
   ! Core routines
   !----------------
   ! chemstate
   use ChemState_Mod, only: cc_find_species_by_name => FindSpecByName
   use ChemState_Mod, only: cc_get_species_conc => GetSpecConc
   use ChemState_Mod, only: cc_get_species_conc_by_name => GetSpecConcByName
   use ChemState_Mod, only: cc_get_species_conc_by_index => GetSpecConcByIndex
   use ChemState_Mod, only: cc_allocate_chemstate => Chem_Allocate
   ! metstate
   use MetState_Mod, only: cc_allocate_metstate => Met_Allocate
   ! diagstate
   use DiagState_Mod, only: cc_allocate_diagstate => Diag_Allocate
   ! emisstate
   use EmisState_Mod, only: cc_allocate_emisstate => Emis_Allocate
   use EmisState_Mod, only: cc_deallocate_emisstate => EmisState_CleanUp
   use EmisState_Mod, only: cc_emis_to_chem_map => Emis_Find_Chem_Map_Index
   use EmisState_Mod, only: cc_apply_emis_to_chem => Apply_Emis_to_Chem

   !-------------------
   ! Configuration Read
   !-------------------
   use Config_Mod, only: cc_read_config => Read_Input_File  !< Method for reading the configuration file

   !---------------
   ! Error Handling
   !---------------
   use Error_Mod, only: cc_check_var => CC_CheckVar     !< Method for checking variables
   use Error_Mod, only: cc_emit_error => CC_Error       !< Method for emitting errors
   use Error_Mod, only: CC_FAILURE                      !< CATCHem Failure return code
   use Error_Mod, only: CC_SUCCESS                      !< CATChem Successful return code
   use Error_Mod, only: cc_emit_warning => CC_Warning   !< Method for emitting warnings
   !
   use init_mod, only: cc_init_diag => Init_Diag        !< Method for initializing the diag state
   use init_mod, only: cc_init_met => Init_Met          !< Method for initializing the met state

   !------------------
   ! CATChem Precision
   !------------------
   use precision_mod, only: cc_rk => fp                 !< Real Precision

   !------------------
   ! CATChem Processes
   !------------------
   ! Dust
   use CCPr_Dust_Common_Mod, only: DustStateType                                !< Dust State
   use CCPr_Dust_mod, only: cc_dust_init => CCPr_Dust_Init                      !< Dust Process Initialization Routine
   use CCPr_Dust_mod, only: cc_dust_run => CCPr_Dust_Run                        !< Dust Process Run Routine
   use CCPr_Dust_mod, only: cc_dust_finalize => CCPr_Dust_Finalize              !< Dust Process Finalization Routine
   ! Seasalt
   use CCPr_SeaSalt_Common_Mod, only: SeaSaltStateType                          !< SeaSalt State
   use CCPr_SeaSalt_mod, only: cc_seasalt_init => CCPr_SeaSalt_Init             !< SeaSalt Process Initialization Routine
   use CCPr_SeaSalt_mod, only: cc_seasalt_run => CCPr_SeaSalt_Run               !< SeaSalt Process Run Routine
   use CCPr_SeaSalt_mod, only: cc_seasalt_finalize => CCPr_SeaSalt_Finalize     !< SeaSalt Process Finalization Routine
   ! Plumerise
   use CCPr_Plumerise_mod, only: PlumeRiseStateType                               !< Plumerise State
   use CCPr_Plumerise_mod, only: cc_plumerise_init => CCPr_Plumerise_Init         !< Plumerise Process Initialization Routine
   use CCPr_Plumerise_mod, only: cc_plumerise_run => CCPr_Plumerise_Run           !< Plumerise Process Run Routine
   use CCPr_Plumerise_mod, only: cc_plumerise_finalize => CCPr_Plumerise_Finalize !< Plumerise Process Finalization Routine
   ! Dry Dep
   use CCPr_DryDep_mod, only: DryDepStateType                            !< DryDep State
   use CCPr_DryDep_mod, only: cc_drydep_init => CCPr_DryDep_Init         !< DryDep Process Initialization Routine
   use CCPr_DryDep_mod, only: cc_drydep_run => CCPr_DryDep_Run           !< DryDep Process Run Routine
   use CCPr_DryDep_mod, only: cc_drydep_finalize => CCPr_DryDep_Finalize !< DryDep Process Finalization Routine
   ! SUVolcanicEmissions
   use CCPr_SUVolcanicEmissions_mod, only: SUVolcanicStateType                                         !< SUVolcanic State
   use CCPr_SUVolcanicEmissions_mod, only: cc_suvolcanic_init => CCPr_SUVolcanicEmissions_Init         !< SUVolcanicEmissions Process Initialization Routine
   use CCPr_SUVolcanicEmissions_mod, only: cc_suvolcanic_run => CCPr_SUVolcanicEmissions_Run           !< SUVolcanicEmissions Process Run Routine
   use CCPr_SUVolcanicEmissions_mod, only: cc_suvolcanic_finalize => CCPr_SUVolcanicEmissions_Finalize !< SUVolcanicEmissions Process Finalization Routine
   ! DMS
   use CCPr_DMS_mod, only: DMSStateType                         !< DMS State
   use CCPr_DMS_mod, only: cc_dms_init => CCPr_DMS_Init         !< DMS Process Initialization Routine
   use CCPr_DMS_mod, only: cc_dms_run => CCPr_DMS_Run           !< DMS Process Run Routine
   use CCPr_DMS_mod, only: cc_dms_finalize => CCPr_DMS_Finalize !< DMS Process Finalization Routine

   implicit none

   public

end module CATChem
