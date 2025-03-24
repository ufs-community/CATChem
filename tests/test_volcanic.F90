program test_volcanic
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   type(ConfigType) :: Config
   type(ChemStateType) :: ChemState
   type(MetStateType) :: MetState
   type(DiagStateType) :: DiagState
   type(VolcanicStateType) :: VolcanicState
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


   thisLoc = 'test_volcanic -> at read CATChem_Config.yml'
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
   VolcanicState%Activate = .false.
   call print_info(Config, VolcanicState, MetState, title)
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

   !ChemState%nSpeciesVolcanic = 2
   !VolcanicState%Activate = .true.

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

   !VolcanicState%SchemeOpt = 1

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "Volcanic Test 2 | Test GOCART Volcanic defaults"
   Config%volcanic_activate = .TRUE.
   Config%volcanic_scheme = 1

   call cc_volcanic_init(Config, VolcanicState, EmisState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_volcanic_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_volcanic_run(MetState, VolcanicState, EmisState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _volcanicemissions_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   !TODO: This is specific to the inputs in this test only. Change it when inputs are changed.
   call assert( sum(VolcanicState%TotalEmission) > 1.0e-2_fp, "Test non-zero Sulfur Volcanic Emissions")
   call assert( rae(sum(VolcanicState%TotalEmission(1:3)), 0.0_fp) , "Test1 zero Sulfur Volcanic Emissions")
   call assert( rae(sum(VolcanicState%TotalEmission(8:28)),  0.0_fp),  "Test2 zero Sulfur Volcanic Emissions")
   call print_info(Config, VolcanicState, MetState, title)

   call cc_volcanic_finalize( VolcanicState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _volcanic_finalize'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

contains

   subroutine print_info(Config_, VolcanicState_, MetState_, title_)

      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(VolcanicStateType), intent(in) :: VolcanicState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%volcanic_activate = ', Config_%volcanic_activate
      write(*,*) 'Config%volcanic_scheme = ', Config_%volcanic_scheme


      if (VolcanicState_%Activate) then

         write(*,*) 'VolcanicState%Activate = ', VolcanicState_%Activate
         write(*,*) 'VolcanicState%SchemeOpt = ', VolcanicState_%SchemeOpt
         write(*,*) 'MetState%DELP =', MetState_%DELP
         write(*,*) 'MetState%BXHEIGHT = ', MetState_%BXHEIGHT
         write(*,*) 'VolcanicState%CatIndex = ', VolcanicState_%CatIndex
         write(*,*) 'VolcanicState%nVolcanicSpecies = ', VolcanicState_%nVolcanicSpecies
         write(*,*) 'VolcanicState%VolcanicSpeciesName = ', VolcanicState_%VolcanicSpeciesName
         write(*,*) 'VolcanicState%EmissionPerSpecies = ', VolcanicState_%EmissionPerSpecies
         write(*,*) 'VolcanicState%TotalEmission = ', VolcanicState_%TotalEmission

      end if

   end subroutine print_info


end program test_volcanic
