

# File catchem.F90

[**File List**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**catchem.F90**](catchem_8_f90.md)

[Go to the documentation of this file](catchem_8_f90.md)


```Fortran

module catchem

   !---------------
   ! CATChem States
   !---------------
   use chemstate_mod,  only: chemstatetype    ! Chemical State
   ! use GridState_Mod,  only: GridStateType    ! Grid State (DEPRECATED - removed)
   use diagstate_mod,  only: diagstatetype    ! Diagnostic State
   use emisstate_mod,  only: emisstatetype    ! Emission State
   use metstate_mod,   only: metstatetype     ! Meteorology State
   use species_mod,    only: speciestype      ! Species State
   use config_opt_mod, only: configtype
   use state_mod,      only: statecontainertype, statebuildertype  ! Modern state management

   !------------------
   ! Grid Management
   !------------------
   use gridmanager_mod, only: gridmanagertype, gridgeometrytype, columniteratortype
   use columninterface_mod, only: virtualcolumntype, columnprocessortype

   !--------------------
   ! Process Management
   !--------------------
   use processinterface_mod, only: processinterface, columnprocessinterface
   use processmanager_mod, only: processmanagertype
   use processfactory_mod, only: processfactorytype

   !----------------
   ! Core routines
   !----------------
   ! chemstate
   use chemstate_mod, only: cc_find_species_by_name => findspecbyname
   use chemstate_mod, only: cc_get_species_conc => getspecconc
   use chemstate_mod, only: cc_get_species_conc_by_name => getspecconcbyname
   use chemstate_mod, only: cc_get_species_conc_by_index => getspecconcbyindex
   ! use ChemState_Mod, only: cc_allocate_chemstate => Chem_Allocate  ! DEPRECATED
   ! metstate
   ! use MetState_Mod, only: cc_allocate_metstate => Met_Allocate  ! DEPRECATED
   ! diagstate
   ! diagstate - Already modernized, use type-bound procedures:
   ! - DiagState%init() for initialization
   ! - DiagState%cleanup() for cleanup
   ! - DiagState%validate() for validation
   ! emisstate - All legacy functions removed. Use modern type-bound procedures:
   ! - EmisState%init() for initialization
   ! - EmisState%cleanup() for cleanup
   ! - EmisState%find_species() for species mapping
   ! - EmisState%apply_emissions() for applying emissions

   !-------------------
   ! Configuration Read
   !-------------------
   ! Legacy config_mod disabled - use modern ConfigManager instead:
   ! - ConfigManager%load_config() for configuration loading
   ! - ConfigManager%get_value() for accessing values

   !---------------
   ! Error Handling
   !---------------
   use error_mod, only: cc_check_var => cc_checkvar     ! Method for checking variables
   use error_mod, only: cc_emit_error => cc_error       ! Method for emitting errors
   use error_mod, only: cc_failure                      ! CATCHem Failure return code
   use error_mod, only: cc_success                      ! CATChem Successful return code
   use error_mod, only: cc_emit_warning => cc_warning   ! Method for emitting warnings
   !
   use init_mod, only: cc_init_diag => init_diag        ! Method for initializing the diag state
   use init_mod, only: cc_init_met => init_met          ! Method for initializing the met state
   use init_mod, only: cc_init_emis => init_emis        ! Method for initializing the emission state
   use run_mod, only: cc_init_process => init_process  ! Method for initializing the process state
   use run_mod,  only: cc_run_process => run_process    ! Method for running the processes
   use run_mod,  only: cc_finalize_process => finalize_process    ! Method for finalizing the processes

   !-----------------------------
   ! Column Virtualization API
   !-----------------------------
   public :: cc_init_grid_manager       ! Initialize grid manager with column virtualization
   public :: cc_create_virtual_column   ! Create a virtual column for processing
   public :: cc_run_column_processes    ! Run all column-based processes
   public :: cc_process_all_columns     ! Process all columns with specified process

   !------------------
   ! CATChem Precision
   !------------------
   use precision_mod, only: cc_rk => fp                 ! Real Precision

   !------------------
   ! CATChem Processes
   !------------------
   ! Dust
   use ccpr_dust_common_mod, only: duststatetype                                ! Dust State
   use ccpr_dust_mod, only: cc_dust_init => ccpr_dust_init                      ! Dust Process Initialization Routine
   use ccpr_dust_mod, only: cc_dust_run => ccpr_dust_run                        ! Dust Process Run Routine
   use ccpr_dust_mod, only: cc_dust_finalize => ccpr_dust_finalize              ! Dust Process Finalization Routine
   ! Seasalt
   use ccpr_seasalt_common_mod, only: seasaltstatetype                          ! SeaSalt State
   use ccpr_seasalt_mod, only: cc_seasalt_init => ccpr_seasalt_init             ! SeaSalt Process Initialization Routine
   use ccpr_seasalt_mod, only: cc_seasalt_run => ccpr_seasalt_run               ! SeaSalt Process Run Routine
   use ccpr_seasalt_mod, only: cc_seasalt_finalize => ccpr_seasalt_finalize     ! SeaSalt Process Finalization Routine
   ! Plumerise
   use ccpr_plumerise_mod, only: plumerisestatetype                               ! Plumerise State
   use ccpr_plumerise_mod, only: cc_plumerise_init => ccpr_plumerise_init         ! Plumerise Process Initialization Routine
   use ccpr_plumerise_mod, only: cc_plumerise_run => ccpr_plumerise_run           ! Plumerise Process Run Routine
   use ccpr_plumerise_mod, only: cc_plumerise_finalize => ccpr_plumerise_finalize ! Plumerise Process Finalization Routine
   ! Dry Dep
   use ccpr_drydep_mod, only: drydepstatetype                            ! DryDep State
   use ccpr_drydep_mod, only: cc_drydep_init => ccpr_drydep_init         ! DryDep Process Initialization Routine
   use ccpr_drydep_mod, only: cc_drydep_run => ccpr_drydep_run           ! DryDep Process Run Routine
   use ccpr_drydep_mod, only: cc_drydep_finalize => ccpr_drydep_finalize ! DryDep Process Finalization Routine
   ! Chemical mechanism solver (TODO: turn off now)
   !use CCPr_Chem_mod, only: cc_get_micm_version => get_micm_version

   implicit none

   public

contains

   subroutine cc_init_grid_manager(container, nx, ny, nz, rc)
      type(StateContainerType), intent(inout) :: container
      integer, intent(in) :: nx, ny, nz
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr
      type(GridGeometryType) :: geometry

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = cc_failure
         return
      endif

      ! Initialize geometry
      call geometry%init(nx, ny, nz, rc)
      if (rc /= cc_success) return

      ! Initialize grid manager
      call grid_mgr%init(geometry, rc)
      if (rc /= cc_success) return

      rc = cc_success
   end subroutine cc_init_grid_manager

   subroutine cc_create_virtual_column(container, col_idx, virtual_col, rc)
      type(StateContainerType), intent(inout) :: container
      integer, intent(in) :: col_idx
      type(VirtualColumnType), intent(out) :: virtual_col
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = cc_failure
         return
      endif

      ! Create virtual column
      call grid_mgr%create_virtual_column(col_idx, virtual_col, rc)
      if (rc /= cc_success) return

      ! Extract data from container
      call virtual_col%extract_from_container(container, rc)
   end subroutine cc_create_virtual_column

   subroutine cc_run_column_processes(proc_mgr, container, rc)
      type(ProcessManagerType), intent(inout) :: proc_mgr
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      call proc_mgr%run_column_processes(container, rc)
   end subroutine cc_run_column_processes

   subroutine cc_process_all_columns(proc_mgr, process_name, container, rc)
      type(ProcessManagerType), intent(inout) :: proc_mgr
      character(len=*), intent(in) :: process_name
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      call proc_mgr%run_process(process_name, container, rc)
   end subroutine cc_process_all_columns

end module catchem
```


