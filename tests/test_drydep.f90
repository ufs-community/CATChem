program test_drydep
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   type(ConfigType) :: Config
   type(ChemStateType) :: ChemState
   type(MetStateType) :: MetState
   type(DiagStateType) :: DiagState
   type(DryDepStateType) :: DryDepState
   type(GridStateType) :: GridState
   type(EmisStateType) :: EmisState

   ! Integers
   INTEGER:: rc          ! Success or failure

   character(len=:), allocatable :: title
   integer :: i ! loop counter

   ! Error handling
   CHARACTER(LEN=512) :: errMsg
   CHARACTER(LEN=255) :: thisLoc
   CHARACTER(LEN=255), PARAMETER :: configFile ='Configs/Default/CATChem_config.yml'


   thisLoc = 'test_drydep -> at read CATChem_Config.yml'
   errMsg = ''
   rc = CC_SUCCESS

   write(*,*) '   CCCCC      A     TTTTTTT   CCCCC  H'
   write(*,*) '  C          A A       T     C       H       CCCC   EEEE   M       M'
   write(*,*) '  C         AAAAA      T     C       HHHHH  C      E    E  M M   M M'
   write(*,*) '  C        A     A     T     C       H   H  C      E EE    M   M   M'
   write(*,*) '   CCCCC  A       A    T      CCCCC  H   H   CCCC   EEEEE  M       M'
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


   title = 'drydep Test 1 | Read Config'
   !DryDepState%SchemeOpt = 1
   DryDepState%Activate = .false.
   call print_info(Config, DryDepState, MetState, title)
   write (*,*) '-- '
   write (*,*) 'Completed ', title
   write (*,*) '--'

   !----------------------------
   ! Test 2
   !----------------------------
   ! Set number of drydep species

   ChemState%nSpeciesAerodrydep = 2
   DryDepState%Activate = .true.

   ! Meteorological State
   MetState%LWI = 1.0_fp
   MetState%USTAR = 0.1_fp
   MetState%PBLH = 1000.0_fp
   MetState%HFLUX = 0.5_fp
   MetState%Z0H = 0.1_fp
   Metstate%NLEVS = 5
   Metstate%TSTEP = 60
   Metstate%U10M = 3.0
   Metstate%V10M = 3.0
   Metstate%FRLAKE = 0.0
   Metstate%GWETTOP = 0.00001
   allocate(MetState%AIRDEN(MetState%NLEVS))
   allocate(MetState%T(MetState%NLEVS))
   allocate(MetState%ZMID(MetState%NLEVS))
   gridstate%number_of_levels = MetState%NLEVS

   do i = 1, MetState%NLEVS
      MetState%T(i)=273.15 + I        ! K, roughly adiabatic
      MetState%AIRDEN(i) = 1.2   ! kg/m3
      MetState%ZMID(i) = (MetState%NLEVS*100 - I*100)   ! m
   end do

   DryDepState%SchemeOpt = 1
   ! Turn off resuspension
   DryDepState%Resuspension = .FALSE.

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "DryDep Test 2 | Test GOCART DryDep defaults"

   call cc_drydep_init(Config, DryDepState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_drydep_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   ! commenting out for now
   call cc_drydep_run(MetState, DiagState, DryDepState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_drydep_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call print_info(Config, DryDepState, MetState, title)
   call assert(DiagState%drydep_frequency(1) > 0.0_fp, "Test GOCART DryDep Scheme (no resuspension)")


   !----------------------------
   ! Test 3
   !----------------------------
   title = "drydep Test 3 | resuspension is .TRUE. "
   ChemState%nSpeciesAerodrydep = 1
   ! Turn on resuspension
   DryDepState%Resuspension = .TRUE.
   DryDepState%particleradius = 0.000001   ! [m]
   DryDepState%particledensity = 2500.   !  [kg/m3]

   call cc_drydep_run(MetState, DiagState, DryDepState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_drydep_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   ! Please revisit statements below - confirm only 1 valid value is being returned
   call print_info(Config, DryDepState, MetState, title)
   call assert(DiagState%drydep_frequency(1) > 0.0_fp, "Test 2 GOCART drydep Scheme (resuspension activated)")


contains

   subroutine print_info(Config_, DryDepState_, MetState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(DryDepStateType), intent(in) :: DryDepState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%drydep_activate = ', Config_%drydep_activate
      write(*,*) 'Config%drydep_scheme = ', Config_%drydep_scheme
      write(*,*) 'Config%drydep_resuspension = ', Config_%drydep_resuspension

      if (DryDepState_%Activate) then

         write(*,*) 'DryDepState%Activate = ', DryDepState_%Activate
         write(*,*) 'DryDepState%SchemeOpt = ', DryDepState_%SchemeOpt
         write(*,*) 'DryDepState%Resuspension = ', DryDepState_%Resuspension

         if (DryDepState_%Resuspension) then
            write(*,*) 'MetState%GWETTOP =', MetState_%GWETTOP
            write(*,*) 'MetState%USTAR =', MetState_%USTAR
         end if

         write(*,*) 'MetState%AIRDEN =', MetState_%AIRDEN
         write(*,*) 'DryDepState%drydepf = ', DryDepState_%drydep_frequency

      end if

   end subroutine print_info


end program test_drydep
