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
   call print_info(Config, DryDepState, MetState, ChemState, title)
   write (*,*) '-- '
   write (*,*) 'Completed ', title
   write (*,*) '--'

   !----------------------------
   ! Test 2
   !----------------------------
   ! Set number of drydep species

   !ChemState%nSpeciesAerodrydep = 2
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

   call print_info(Config, DryDepState, MetState, ChemState, title)
   call assert(DiagState%drydep_frequency(1) > 0.0_fp, "Test GOCART DryDep Scheme (no resuspension)")


   !----------------------------
   ! Test 3
   !----------------------------
   title = "drydep Test 3 | resuspension is .TRUE. "
   !ChemState%nSpeciesAerodrydep = 1
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
   call print_info(Config, DryDepState, MetState, ChemState, title)
   call assert(DiagState%drydep_frequency(1) > 0.0_fp, "Test 2 GOCART drydep Scheme (resuspension activated)")


   !----------------------------
   ! Test 4
   !----------------------------
   title = "drydep Test 4 | scheme_opt=2 "

   !clean up the test above for a different scheme test (Test 4)
   call cc_drydep_finalize( DryDepState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_drydep_finalize'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   !deallocate the diagstate arrays
   if (allocated(DiagState%drydep_frequency)) deallocate(DiagState%drydep_frequency)
   if (allocated(DiagState%drydep_vel)) deallocate(DiagState%drydep_vel)

   !assign MetState values
   MetState%SWGDN = 500.0_fp
   MetState%TS = 301.0_fp
   MetState%SUNCOSmid=  0.97_fp
   MetState%USTAR = 0.1_fp
   MetState%OBK = -100
   MetState%CLDFRC = 0.1
   MetState%PBLH = 1000.0_fp
   allocate(MetState%BXHEIGHT(MetState%NLEVS))
   MetState%BXHEIGHT = 40.0_fp
   MetState%Z0 = 1000_fp
   allocate(MetState%RH(MetState%NLEVS))
   MetState%RH = 0.4661  !unitless (low values[<=0.466 in this test] will lead to the error of DEN for seasalt species)
   MetState%PS = 1000.0_fp ! hPa
   Metstate%U10M = 3.0
   Metstate%V10M = 3.0
   MetState%FRLAI = (/  3.0, 3.0, 1.0, 0.0/)      !TODO: whether LAI is separated to each land type?
   MetState%ILAND = (/   5,   6,  18,   1 /)
   MetState%FRLANDUSE = (/ 0.4, 0.4, 0.1, 0.1 /)
   MetState%SALINITY=10     ! greater than 20 (in ppt; part per thousand) is considered as ocean
   MetState%TSKIN = 305.0_fp
   MetState%IODIDE = 100_fp !in [nM; nanoMolar]
   MetState%LON = -92.0_fp
   MetState%LAT = 38.0_fp
   MetState%LNLPBL = .true.
   MetState%IsSnow = .false.
   MetState%IsIce = .false.
   MetState%IsLand = .true.

   Config%drydep_scheme = 2
   !ChemState%nSpeciesDrydep = 34

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   call cc_drydep_init(Config, DryDepState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_drydep_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_drydep_run(MetState, DiagState, DryDepState, ChemState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_drydep_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   ! Please revisit statements below - confirm only 1 valid value is being returned
   call print_info(Config, DryDepState, MetState, ChemState, title)
   call assert(DiagState%drydep_frequency(1) > 0.0_fp, "Test 4 Wesely drydep Scheme")


contains

   subroutine print_info(Config_, DryDepState_, MetState_, ChemState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(DryDepStateType), intent(in) :: DryDepState_
      type(ChemStateType), intent(in) :: ChemState_
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

         if (DryDepState_%SchemeOpt == 1) then
            write(*,*) 'ChemState%nSpeciesAerodrydep = ', ChemState_%nSpeciesAerodrydep
            write(*,*) 'ChemState%chemSpecies%name =', ChemState_%chemSpecies(ChemState%AeroDryDepIndex(:))%short_name
            write(*,*) 'DryDepState_%drydep_vel =', DryDepState_%drydep_vel
            write(*,*) 'DryDepState%drydepf = ', DryDepState_%drydep_frequency
         else if (DryDepState_%SchemeOpt == 2) then
            write(*,*) 'ChemState%nSpeciesDrydep = ', ChemState_%nSpeciesDrydep
            write(*,*) 'ChemState%chemSpecies%name =', ChemState_%chemSpecies(ChemState%DryDepIndex(:))%short_name
            write(*,*) 'DryDepState_%drydep_vel =', DryDepState_%drydep_vel
            write(*,*) 'DryDepState%drydepf = ', DryDepState_%drydep_frequency
         end if

      end if

   end subroutine print_info


end program test_drydep
