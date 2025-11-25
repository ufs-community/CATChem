

# File run\_mod.F90

[**File List**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**run\_mod.F90**](run__mod_8_f90.md)

[Go to the documentation of this file](run__mod_8_f90.md)


```Fortran

module run_mod
   USE metstate_mod, ONLY : metstatetype
   USE chemstate_mod, ONLY : chemstatetype
   USE diagstate_mod, ONLY : diagstatetype
   USE emisstate_mod, ONLY : emisstatetype
   use config_opt_mod, only: configtype
   USE ccpr_dust_common_mod, ONLY : duststatetype
   USE ccpr_seasalt_common_mod, ONLY : seasaltstatetype
   USE ccpr_drydep_mod, ONLY : drydepstatetype
   USE error_mod

   implicit none

   PRIVATE

   PUBLIC :: run_process
   PUBLIC :: finalize_process
   PUBLIC :: init_process

contains


   subroutine run_process(MetState, DiagState, ChemState, EmisState, DustState, SeaSaltState, DryDepState, RC)
      !
      use ccpr_drydep_mod, Only : ccpr_drydep_run

      implicit none

      ! Arguments
      type(MetStateType), INTENT(IN) :: MetState
      type(DiagStateType), INTENT(INOUT) :: DiagState
      type(ChemStateType), INTENT(INOUT) :: ChemState
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      type(DustStateType), INTENT(INOUT) :: DustState
      type(SeaSaltStateType), INTENT(INOUT) :: SeaSaltState
      type(DryDepStateType), INTENT(INOUT) :: DryDepState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      rc = 0
      errmsg = ''
      thisloc = ' -> at Run_Process (in core/run_mod.F90)'

      !run all the emission processes
      call run_emis(metstate, diagstate, chemstate, emisstate, duststate, seasaltstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error running emission processes.'
         call cc_error(errmsg, rc , thisloc)
      endif

      !run dry deposition
      call ccpr_drydep_run( metstate, diagstate, drydepstate, chemstate, rc )
      if (rc /= cc_success) then
         errmsg = 'Error running dry deposition.'
         call cc_error(errmsg, rc , thisloc)
      endif

      !run wet deposition

   end subroutine run_process


   subroutine run_emis(MetState, DiagState, ChemState, EmisState, DustState, SeaSaltState, RC)
      use ccpr_dust_mod, ONLY : ccpr_dust_run
      use ccpr_seasalt_mod, ONLY : ccpr_seasalt_run
      use emisstate_mod, only: apply_emis_to_chem
      implicit none
      ! Arguments
      type(ChemStateType), intent(inout) :: ChemState
      type(DiagStateType), intent(inout) :: DiagState
      type(MetStateType), intent(in)  :: MetState
      type(EmisStateType), intent(inout) :: EmisState
      type(DustStateType), intent(inout) :: DustState
      type(SeaSaltStateType), intent(inout) :: SeaSaltState
      integer, intent(inout) :: RC

      !local variables
      !integer :: i, k

      ! Error handling
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc
      thisloc = ' -> at Run_Emis (in core/run_mo.F90)'


      !TODO: make sure each emission process has filled the EmisState flux since cc_apply_emis_to_chem is based on the EmisState
      ! Call dust emission
      !TODO: dust per bin emission only works on the default 5 biins of the Ginux scheme despite what defined in the emission namelist
      call ccpr_dust_run(metstate, diagstate, duststate, emisstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error running dust emissions'
         call cc_error(errmsg, rc , thisloc)
         return
      end if

      ! Call sea salt emission
      call ccpr_seasalt_run(metstate, seasaltstate, emisstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error running seal salt emissions'
         call cc_error(errmsg, rc , thisloc)
         return
      end if

      !apply emissions to chemistry state
      call apply_emis_to_chem(emisstate, metstate, chemstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error applying emissions to chemistry state'
         call cc_error(errmsg, rc , thisloc)
         return
      end if


   end subroutine run_emis


   subroutine finalize_process(DustState, SeaSaltState, DryDepState, RC)
      !
      use ccpr_dust_mod, Only : ccpr_dust_finalize
      use ccpr_seasalt_mod, Only : ccpr_seasalt_finalize
      use ccpr_drydep_mod, Only : ccpr_drydep_finalize

      implicit none

      ! Arguments
      type(DustStateType), INTENT(INOUT) :: DustState
      type(SeaSaltStateType), INTENT(INOUT) :: SeaSaltState
      type(DryDepStateType), INTENT(INOUT) :: DryDepState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      rc = 0
      errmsg = ''
      thisloc = ' -> at Finalize_Process (in core/run_mod.F90)'

      !finalize all the emission processes
      call ccpr_dust_finalize(duststate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error finalizing dust emissions.'
         call cc_error(errmsg, rc , thisloc)
      endif
      call ccpr_seasalt_finalize(seasaltstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error finalizing sea salt emissions.'
         call cc_error(errmsg, rc , thisloc)
      endif

      !finalize dry deposition
      call ccpr_drydep_finalize(drydepstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error finalizing dry deposition.'
         call cc_error(errmsg, rc , thisloc)
      endif

      !finalize wet deposition

   end subroutine finalize_process

   subroutine init_process (config, ChemState, EmisState, DustState, SeaSaltState, DryDepState, RC)
      use ccpr_dust_mod, only: ccpr_dust_init
      use ccpr_seasalt_mod, only: ccpr_seasalt_init
      use ccpr_drydep_mod, only: ccpr_drydep_init
      implicit none

      type(ConfigType), intent(in) :: config
      type(ChemStateType), intent(inout) :: ChemState
      type(EmisStateType), intent(inout) :: EmisState
      type(DustStateType), intent(inout) :: DustState
      type(SeaSaltStateType), intent(inout) :: SeaSaltState
      type(DryDepStateType), intent(inout) :: DryDepState
      integer, intent(inout) :: RC

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc
      thisloc = ' -> at init_process (in core/init_mod.F90)'

      ! Initialize DustState
      call ccpr_dust_init(config, duststate, chemstate, emisstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error in CCPR_Dust_Init'
         call cc_error(errmsg, rc, thisloc )
         return
      end if

      ! Initialize SealSaltState
      call ccpr_seasalt_init(config, seasaltstate, chemstate, emisstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error in CCPr_SeaSalt_Init'
         call cc_error(errmsg, rc, thisloc )
         return
      end if

      ! Initialize DryDepState
      call ccpr_drydep_init(config, drydepstate, chemstate, rc)
      if (rc /= cc_success) then
         errmsg = 'Error in CCPr_DryDep_Init'
         call cc_error(errmsg, rc, thisloc )
         return
      end if

   end subroutine init_process

end module run_mod
```


