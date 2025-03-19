program test_suvolcanic
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   type(ConfigType) :: Config
   type(ChemStateType) :: ChemState
   type(MetStateType) :: MetState
   type(DiagStateType) :: DiagState
   type(SUVolcanicStateType) :: SUVolcanicState
   type(GridStateType) :: GridState
   type(EmisStateType) :: EmisState

   ! Integers
   INTEGER:: rc          ! Success or failure

   character(len=:), allocatable :: title
   integer :: c ,s  ! Loop counter for emission state

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
   SUVolcanicState%Activate = .false.
   call print_info(Config, SUVolcanicState, MetState, title)
   write (*,*) '-- '
   write (*,*) 'Completed ', title
   write (*,*) '--'

   !allocate emission state
   if (EmisState%nCats > 0) then
      do c = 1, EmisState%nCats
         do s = 1, EmisState%Cats(c)%nSpecies
            ALLOCATE(EmisState%Cats(c)%Species(s)%Flux(GridState%number_of_levels), STAT=RC)
            if (RC /= CC_SUCCESS) then
               ErrMsg = 'Error allocating "EmisState%Cats%Species%Flux"!'
               call cc_emit_error(ErrMsg, RC, ThisLoc)
               stop 1  !!Note here is not 'return'
            endif
         end do
      end do
   end if

   !----------------------------
   ! Test 2
   !----------------------------
   ! Set number of Volcanic species

   !ChemState%nSpeciesSUVolcanic = 2
   !SUVolcanicState%Activate = .true.

   ! Meteorological State
   MetState%YMD = 20220101
   MetState%HMS = 120000
   MetState%TSTEP = 300
   MetState%AREA_M2 = 10000.0_fp
   MetState%NLEVS = 8
   allocate(MetState%BXHEIGHT(MetState%NLEVS))

   !TODO: is the layer index reversed in the GOCART?
   allocate(MetState%DELP(MetState%NLEVS))
   MetState%DELP(1:MetState%NLEVS)= (/25000, 20000, 15000, 110000, 10000, 9000,5000,4000/)! Need to change to something more reasonable and check units.
   MetState%BXHEIGHT(1:MetState%NLEVS) =(/10000, 8000, 7000, 5000, 3000, 1000, 100, 50/)  ! temporary, change to something more reasonable and check units

   !SUVolcanicState%SchemeOpt = 1

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "SUVolcanic Test 2 | Test GOCART SUVolcanic defaults"
   Config%suvolcanic_activate = .TRUE.
   Config%suvolcanic_scheme = 1

   call cc_suvolcanic_init(Config, SUVolcanicState, EmisState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_suvolcanic_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_suvolcanic_run(MetState, SUVolcanicState, EmisState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _suvolcanicemissions_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call assert(sum(SUVolcanicState%TotalEmission) > 0.0_fp, "Test Sulfur Volcanic Emissions")
   call print_info(Config, SUVolcanicState, MetState, title)

   call cc_suvolcanic_finalize( SUVolcanicState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _suvolcanic_finalize'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

contains

   subroutine print_info(Config_, SUVolcanicState_, MetState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(SUVolcanicStateType), intent(in) :: SUVolcanicState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%suvolcanic_activate = ', Config_%suvolcanic_activate
      write(*,*) 'Config%suvolcanic_scheme = ', Config_%suvolcanic_scheme


      if (SUVolcanicState_%Activate) then

         write(*,*) 'SUVolcanicState%Activate = ', SUVolcanicState_%Activate
         write(*,*) 'SUVolcanicState%SchemeOpt = ', SUVolcanicState_%SchemeOpt
         write(*,*) 'MetState%DELP =', MetState_%DELP
         write(*,*) 'MetState%BXHEIGHT = ', MetState_%BXHEIGHT
         write(*,*) 'SUVolcanicState%CatIndex = ', SUVolcanicState_%CatIndex
         write(*,*) 'SUVolcanicState%nSUVolcanicSpecies = ', SUVolcanicState_%nSUVolcanicSpecies
         write(*,*) 'SUVolcanicState%SUVolcanicSpeciesName = ', SUVolcanicState_%SUVolcanicSpeciesName
         write(*,*) 'SUVolcanicState%EmissionPerSpecies = ', SUVolcanicState_%EmissionPerSpecies
         write(*,*) 'SUVolcanicState%TotalEmission = ', SUVolcanicState_%TotalEmission

      end if

   end subroutine print_info


end program test_suvolcanic
