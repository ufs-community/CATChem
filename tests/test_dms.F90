program test_DMS
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   !type(ConfigType) :: Config
   !type(ChemStateType) :: ChemState
   !type(DMSStateType) :: DMSState
   TYPE(ConfigType), POINTER       :: Config    ! Module options
   type(ChemStateType), POINTER :: ChemState
   type(DMSStateType), POINTER :: DMSState
   type(MetStateType) :: MetState
   type(DiagStateType) :: DiagState
   type(GridStateType) :: GridState
   type(EmisStateType) :: EmisState

   ! Integers
   INTEGER:: rc          ! Success or failure
   INTEGER:: i           ! indexer

   character(len=:), allocatable :: title

   ! Error handling
   CHARACTER(LEN=512) :: errMsg
   CHARACTER(LEN=255) :: thisLoc
   CHARACTER(LEN=255), PARAMETER :: configFile ='Configs/Default/CATChem_config.yml'


   thisLoc = 'test_DMS -> at read CATChem_Config.yml'
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


   title = 'DMS Test 1 | Read Config'
   DMSState%Activate = .false.
   call print_info(Config, DMSState, MetState, title)
   write (*,*) '-- '
   write (*,*) 'Completed ', title
   write (*,*) '--'

   !----------------------------
   ! Test 2
   !----------------------------

   DMSState%Activate = .true.

   ! Meteorological State
   allocate(MetState%T(MetState%NLEVS))
   allocate(MetState%DELP(MetState%NLEVS))
   MetState%NLEVS = 2
   MetState%DELP(1:MetState%NLEVS)= 10000      ! Need to change to something more reasonable and check units.
   MetState%T(1:MetState%NLEVS) = 100  ! temporary, change to something more reasonable and check units
   MetState%U10M = 1.0_fp
   MetState%V10M = 1.0_fp
   MetState%LWI = 0


   do i = 1, MetState%NLEVS
      MetState%T(i)=273.15 + I   ! K -
      MetState%DELP(i) = 5000  ! check units (Pa), this is shallow for near sfc!
   end do

   DMSState%SchemeOpt = 1

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "DMS Test 2 | Test GOCART DMS defaults"

   call cc_dms_init(Config, DMSState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_dms_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_dms_run(MetState, DiagState, &
      DMSState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _dms_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call print_info(Config, DMSState, MetState, title)
   call cc_dms_finalize( DMSState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _dms_finalize'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

contains

   subroutine print_info(Config_, DMSState_, MetState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(DMSStateType), intent(in) :: DMSState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%dms_activate = ', Config_%dms_activate
      write(*,*) 'Config%dms_scheme = ', Config_%dms_scheme


      if (DMSState_%Activate) then

         write(*,*) 'DMSState%Activate = ', DMSState_%Activate
         write(*,*) 'DMSState%SchemeOpt = ', DMSState_%SchemeOpt
         write(*,*) 'MetState%DELP =', MetState_%DELP
         write(*,*) 'MetState%T = ', MetState_%T
         write(*,*) 'MetState%U10M = ', MetState_%U10M
         write(*,*) 'MetState%V10M = ', MetState_%V10M
         write(*,*) 'MetState%LWI = ', MetState_%LWI

      end if

   end subroutine print_info


end program test_dms
