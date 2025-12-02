!> \file catchem_nuopc_driver.F90
!! \brief Simple NUOPC driver for testing CATChem cap in standalone mode
!!
!! \details
!! This is a standalone NUOPC driver application that can be used to test
!! the CATChem NUOPC cap in isolation from larger Earth system models.
!! It sets up a minimal NUOPC application with only the CATChem component,
!! enabling testing and validation of the cap functionality.
!!
!! The driver follows standard NUOPC application patterns:
!! - Initializes ESMF framework
!! - Creates and configures the driver component
!! - Adds the CATChem component to the driver
!! - Executes the standard initialize-run-finalize sequence
!! - Provides logging and error handling throughout
!!
!! This standalone capability is essential for:
!! - Component testing and validation
!! - Performance benchmarking
!! - Debugging and development
!! - Integration testing before coupling with larger systems
!!
!! \note This driver requires ESMF/NUOPC libraries and appropriate
!!       configuration files for CATChem
!!
!! \author Barry Baker, NOAA/OAR/ARL
!! \date November 2024
!! \ingroup catchem_nuopc_group

program catchem_nuopc_driver

   use ESMF
   use NUOPC
   use NUOPC_Driver, driverSS => SetServices

   implicit none

   integer :: rc, userRc
   type(ESMF_GridComp) :: driver

   ! Initialize ESMF
   call ESMF_Initialize(rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   call ESMF_LogWrite("CATChem NUOPC Driver: Starting", ESMF_LOGMSG_INFO, rc=rc)

   ! Create the driver component
   driver = ESMF_GridCompCreate(name="CATChem_Driver", rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   ! Set driver services
   call ESMF_GridCompSetServices(driver, userRoutine=DriverSetServices, &
      userRc=userRc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop
   if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   ! Initialize the driver
   call ESMF_GridCompInitialize(driver, userRc=userRc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop
   if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   ! Run the driver
   call ESMF_GridCompRun(driver, userRc=userRc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop
   if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   ! Finalize the driver
   call ESMF_GridCompFinalize(driver, userRc=userRc, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop
   if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   ! Destroy the driver
   call ESMF_GridCompDestroy(driver, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) stop

   call ESMF_LogWrite("CATChem NUOPC Driver: Completed successfully", &
      ESMF_LOGMSG_INFO, rc=rc)

   ! Finalize ESMF
   call ESMF_Finalize(rc=rc)

contains

   ! Set services for the driver
   subroutine DriverSetServices(driver, rc)
      type(ESMF_GridComp) :: driver
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS

      ! Use NUOPC generic driver
      call NUOPC_CompDerive(driver, driverSS, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Set the driver config
      call NUOPC_CompAttributeSet(driver, name="Verbosity", value="high", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

      ! Add the CATChem component
      call NUOPC_DriverAddComp(driver, compLabel="CATChem", &
         compSetServicesRoutine=catchem_comp_set_services, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return

   end subroutine DriverSetServices

   ! Set services for CATChem component
   subroutine catchem_comp_set_services(comp, rc)
      use catchem_nuopc_cap, only: SetServices

      type(ESMF_GridComp) :: comp
      integer, intent(out) :: rc

      call SetServices(comp, rc)

   end subroutine catchem_comp_set_services

end program catchem_nuopc_driver
