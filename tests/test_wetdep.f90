program test_wetdep
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   type(ConfigType) :: Config
   type(ChemStateType) :: ChemState
   type(MetStateType) :: MetState
   !type(DiagStateType) :: DiagState
   type(WetDepStateType) :: WetDepState
   type(GridStateType) :: GridState
   type(EmisStateType) :: EmisState

   ! Integers
   INTEGER:: rc          ! Success or failure

   character(len=:), allocatable :: title
   integer :: n ! loop counter
   !integer :: H2O2_id=0, SO4_id=0

   ! Error handling
   CHARACTER(LEN=512) :: errMsg
   CHARACTER(LEN=255) :: thisLoc
   CHARACTER(LEN=255), PARAMETER :: configFile ='Configs/Default/CATChem_config.yml'


   thisLoc = 'test_wetdep -> at read CATChem_Config.yml'
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
   Metstate%NLEVS = 8
   GridState%number_of_levels = MetState%NLEVS

   call cc_read_config(Config, GridState, EmisState, ChemState, rc, configFile)
   if (rc /= CC_success) then
      errMsg = 'Error reading configuration file: ' // TRIM( configFile )
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   endif

   title = 'wetdep Test 1 | Read Config'
   !WetDepState%AeroSchemeOpt = 1
   WetDepState%Activate = .false.
   call print_info(Config, WetDepState, MetState, ChemState, title)
   write (*,*) '-- '
   write (*,*) 'Completed ', title
   write (*,*) '--'


   !dummy MET variables used for Jacob scheme
   MetState%TSTEP = 300
   allocate(MetState%PEDGE_DRY(MetState%NLEVS), MetState%T(MetState%NLEVS), &
      MetState%MAIRDEN(MetState%NLEVS), MetState%PFLLSAN(MetState%NLEVS), &
      MetState%PFILSAN(MetState%NLEVS), MetState%REEVAPLS(MetState%NLEVS), &
      MetState%AIRDEN(MetState%NLEVS))

   !TODO: is this the right variable for ple?
   MetState%PEDGE_DRY = (/101325, 80000, 70000, 50000, 30000, 10000, 1000, 500/)  !Pa
   MetState%T = (/300, 285, 280, 275, 270, 265, 260, 255/) !K
   MetState%MAIRDEN = (/1.16, 1.06, 0.96, 0.76, 0.56, 0.36, 0.16, 0.06/) !kg/m3
   MetState%PFLLSAN = (/0.0, 0.0, 0.001, 0.001, 0.0008, 0.0006, 0.0004, 0.0/)  !kg/m2/s
   MetState%PFILSAN = (/0.0016, 0.0018, 0.0020, 0.0022, 0.0018, 0.0006, 0.0004, 0.0/) !kg/m2/s
   !TODO: why this is not used in GOCART??
   MetState%REEVAPLS = (/0.0_fp, 2.0e-6_fp, 5.0e-6_fp, 7.0e-6_fp, 6.0e-6_fp, 4.0e-6_fp, 0.0_fp, 0.0_fp/) !kg/kg/s
   !TODO: this is related to REEVAPLS; not sure if MAIRDEN can be used here to replace AIRDEN
   MetState%AIRDEN =(/1.20, 1.10, 1.00, 0.80, 0.60, 0.40, 0.20, 0.10/)

   ! do n = 1, ChemState%nSpecies
   !    if ('H2O2' == TRIM(ChemState%SpeciesNames(n))) then
   !       H2O2_id = n
   !    else if ('SO4' == TRIM(ChemState%SpeciesNames(n))) then
   !       SO4_id = n
   !    endif
   ! enddo

   !give all species concentration the same value, including H2O2 and SO4
   do n = 1, ChemState%nSpecies
      ChemState%ChemSpecies(n)%conc(:) = 5.0e-7_fp ! kg/kg
   enddo
   !ChemState%ChemSpecies(H2O2_id)%conc=  ! kg/kg
   !ChemState%ChemSpecies(SO4_id)%conc    !kg/kg

   !----------------------------
   ! Test 2
   !----------------------------
   title = "WetDep Test 2 | Test Jacob WetDep scheme"

   Config%wetdep_scheme = 1

   ! Allocate DiagState
   !call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   !if (rc /= CC_SUCCESS) then
   !   errMsg = 'Error in cc_allocate_diagstate'
   !   stop 1
   !endif

   call cc_wetdep_init(Config, WetDepState, MetState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_wetdep_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   ! commenting out for now
   call cc_wetdep_run(MetState, WetDepState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_wetdep_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call print_info(Config, WetDepState, MetState, ChemState, title)
   call assert(SUM(WetDepState%wetdep_flux(:,:)) > 0.0_fp, "Test Jacob WetDep Scheme")

   !clean up the test above for a different scheme test if any
   call cc_wetdep_finalize( WetDepState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_wetdep_finalize'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if


contains

   subroutine print_info(Config_, WetDepState_, MetState_, ChemState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(WetDepStateType), intent(in) :: WetDepState_
      type(ChemStateType), intent(in) :: ChemState_
      character(len=*), intent(in) :: title_
      integer :: i ! loop counter

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%wetdep_activate = ', Config_%wetdep_activate
      write(*,*) 'Config%wetdep_scheme = ', Config_%wetdep_scheme


      if (WetDepState_%Activate) then

         write(*,*) 'WetDepState%Activate = ', WetDepState_%Activate
         write(*,*) 'WetDepState%SchemeOpt = ', WetDepState_%SchemeOpt
         write(*,*) 'MetState%AIRDEN =', MetState_%AIRDEN
         write(*,*) 'MetState%MAIRDEN =', MetState_%MAIRDEN
         write(*,*) 'ChemState%nSpeciesWetdep = ', ChemState_%nSpeciesWetdep
         do i = 1, ChemState%nSpeciesWetDep
            write(*,*) 'ChemState%chemSpecies%name =', ChemState_%chemSpecies(ChemState%WetDepIndex(i))%short_name
            write(*,*) 'WetDepState_%wetdep_flux =', WetDepState_%wetdep_flux(i,:)
         enddo

      end if

   end subroutine print_info


end program test_wetdep
