program test_suvolcanic
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   type(ConfigType) :: Config
   type(ChemStateType) :: ChemState
   type(MetStateType) :: MetState
   type(DiagStateType) :: DiagState
   type(SUVolcanicEmissionsStateType) :: SUVolcanicEmissionsState
   type(GridStateType) :: GridState
   type(EmisStateType) :: EmisState

   ! Integers
   INTEGER:: rc          ! Success or failure

   character(len=:), allocatable :: title

   ! Error handling
   CHARACTER(LEN=512) :: errMsg
   CHARACTER(LEN=255) :: thisLoc
   CHARACTER(LEN=255), PARAMETER :: configFile ='Configs/Default/CATChem_config.yml'


   thisLoc = 'test_suvolcanic -> at read CATChem_Config.yml'
   errMsg = ''
   rc = CC_SUCCESS

   write(*,*) '   CCCCC      A     TTTTTTT   CCCCC  H'
   write(*,*) '  C          A A       T     C       H        EEEE   M       M'
   write(*,*) '  C         AAAAA      T     C       HHHHH   E    E  M M   M M'
   write(*,*) '  C        A     A     T     C       H   H   E EE    M   M   M'
   write(*,*) '   CCCCC  A       A    T      CCCCC  H   H    EEEEE  M       M'
   write(*,*) ''
   write(*,*) ''

   !----------------------------
   ! Test 1
   !----------------------------

   ! Read input file and initialize grid
   call cc_read_config(Config, GridState, EmisState, ChemState, rc, configFile)
   if (rc /= CC_success) then
      errMsg = 'Error reading configuration file: ' // TRIM( configFile )
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   endif


   title = 'Volcanic Test 1 | Read Config'
   SUVolcanicEmissionsState%Activate = .false.
   call print_info(Config, SUVolcanicEmissionsState, MetState, title)
   write (*,*) '-- '
   write (*,*) 'Completed ', title
   write (*,*) '--'

   !----------------------------
   ! Test 2
   !----------------------------
   ! Set number of Volcanic species

   ChemState%nSpeciesSUVolcanic = 2
   SUVolcanicEmissionsState%Activate = .true.

   ! Meteorological State
   MetState%NLEVS = 2
   allocate(MetState%BXHEIGHT(MetState%NLEVS))
   allocate(MetState%DELP(MetState%NLEVS))

   MetState%DELP(1:MetState%NLEVS)= 10000      ! Need to change to something more reasonable and check units.
   MetState%BXHEIGHT(1:MetState%NLEVS) = 100  ! temporary, change to something more reasonable and check units

   SUVolcanicEmissionsState%SchemeOpt = 1

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "SUVolcanic Test 2 | Test GOCART SUVolcanic defaults"

   call cc_suvolcanic_init(Config, SUVolcanicEmissionsState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_suvolcanic_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_suvolcanic_run(MetState, DiagState, &
      SUVolcanicEmissionsState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _suvolcanicemissions_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call print_info(Config, SUVolcanicEmissionsState, MetState, title)
   call cc_suvolcanic_finalize( SUVolcanicEmissionsState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _suvolcanic_finalize'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

contains

   subroutine print_info(Config_, SUVolcanicEmissionsState_, MetState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(SUVolcanicEmissionsStateType), intent(in) :: SUVolcanicEmissionsState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%suvolcanicemissions_activate = ', Config_%suvolcanicemissions_activate
      write(*,*) 'Config%suvolcanicemissions_scheme = ', Config_%suvolcanicemissions_scheme


      if (SUVolcanicEmissionsState_%Activate) then

         write(*,*) 'SUVolcanicEmissionsState%Activate = ', SUVolcanicEmissionsState_%Activate
         write(*,*) 'SUVolcanicEmissionsState%SchemeOpt = ', SUVolcanicEmissionsState_%SchemeOpt
         write(*,*) 'MetState%DELP =', MetState_%DELP
         write(*,*) 'MetState%BXHEIGHT = ', MetState_%BXHEIGHT

      end if

   end subroutine print_info


end program test_suvolcanic
