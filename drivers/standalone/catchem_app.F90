!> \file catchem_app.F90
!! \brief Main application for standalone CATChem execution
!!
!! \details
!! This is the top-level program for running CATChem as a standalone ESMF/NUOPC
!! application. It mirrors the structure of the UFS main program (UFS.F90) and
!! the CECE standalone main application (mainApp.F90):
!!
!!   1. Initialize the ESMF framework.
!!   2. Determine the driver configure file (command-line argument or default).
!!   3. Create the top-level driver gridded component.
!!   4. Register the driver's Initialize/Run/Finalize via SetServices.
!!   5. Execute Initialize -> Run -> Finalize.
!!   6. Finalize the ESMF framework.
!!
!! The driver itself (catchem_driver) owns the clock and adds the CATChem
!! cap (and, in the future, any additional components) as children. This keeps
!! the main program minimal and reusable.
!!
!! \par Usage
!! \code
!!   mpirun -np 1 catchem_app [driver.configure]
!! \endcode
!! If no configure file is given, "catchem_standalone.configure" is used.
!!
!! \author CATChem standalone driver
!! \ingroup catchem_nuopc_group
program catchem_app

   use ESMF
   use NUOPC
   use catchem_driver, only: driver_SS => SetServices, set_driver_config_file

   use Constants, only: MAX_LEN_PATH

   implicit none

   character(len=*), parameter :: u_FILE_u = __FILE__

   integer             :: rc, urc
   type(ESMF_VM)       :: vm
   type(ESMF_GridComp) :: driver
   integer             :: localPet
   character(len=MAX_LEN_PATH)  :: driver_cfg_file

   ! ------------------------------------------------------------------
   ! Initialize ESMF
   ! ------------------------------------------------------------------
   call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
      logKindFlag=ESMF_LOGKIND_MULTI, vm=vm, rc=rc)
   if (rc /= ESMF_SUCCESS) then
      write(*,'(A,I0)') "ERROR: ESMF_Initialize failed, rc=", rc
      error stop 1
   end if

   call ESMF_LogWrite("CATChem standalone application starting", &
      ESMF_LOGMSG_INFO, rc=rc)

   call ESMF_VMGet(vm, localPet=localPet, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ! ------------------------------------------------------------------
   ! Resolve the driver configure file
   ! ------------------------------------------------------------------
   call get_command_argument(1, driver_cfg_file)
   if (len_trim(driver_cfg_file) == 0) then
      driver_cfg_file = "catchem_standalone.configure"
   end if
   call set_driver_config_file(trim(driver_cfg_file))

   if (localPet == 0) then
      call ESMF_LogWrite("[catchem_app] Driver configure file: "// &
         trim(driver_cfg_file), ESMF_LOGMSG_INFO, rc=rc)
   end if

   ! ------------------------------------------------------------------
   ! Create and register the driver component
   ! ------------------------------------------------------------------
   driver = ESMF_GridCompCreate(name="catchem_driver", rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetServices(driver, driver_SS, userRc=urc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ! ------------------------------------------------------------------
   ! Initialize
   ! ------------------------------------------------------------------
   call ESMF_GridCompInitialize(driver, userRc=urc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ! ------------------------------------------------------------------
   ! Run (NUOPC_Driver advances its internal clock to the stop time)
   ! ------------------------------------------------------------------
   call ESMF_GridCompRun(driver, userRc=urc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ! ------------------------------------------------------------------
   ! Finalize
   ! ------------------------------------------------------------------
   call ESMF_GridCompFinalize(driver, userRc=urc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridCompDestroy(driver, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=u_FILE_u)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (localPet == 0) then
      write(*,'(A)') "INFO: [catchem_app] CATChem standalone run completed"
   end if

   call ESMF_LogWrite("CATChem standalone application finished", &
      ESMF_LOGMSG_INFO, rc=rc)

   ! ------------------------------------------------------------------
   ! Finalize ESMF
   ! ------------------------------------------------------------------
   call ESMF_Finalize()

end program catchem_app
