program test_DMS
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use precision_mod, only: rae
   implicit none

   type(ConfigType) :: Config
   type(ChemStateType) :: ChemState
   type(DMSStateType) :: DMSState
   type(MetStateType) :: MetState
   type(DiagStateType) :: DiagState
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

   ! Meteorological State
   MetState%TSTEP = 300
   MetState%NLEVS = 1
   allocate(MetState%T(MetState%NLEVS))
   allocate(MetState%DELP(MetState%NLEVS))
   MetState%DELP(1:MetState%NLEVS)= 5000   ! Need to change to something more reasonable and check units.
   MetState%T(1:MetState%NLEVS) = 300      ! temporary, change to something more reasonable and check units
   MetState%U10M = 5.0_fp
   MetState%V10M = 5.0_fp
   MetState%LWI = 0   !gocart OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   MetState%DMSO_CONC = 3.25_fp  !DMS ocean concentration [nmol/L];TODO: may read from ChemState in the future

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "DMS Test 2 | Test GOCART DMS defaults"
   !---------------------------------------------
   DMSState%Activate = .true.
   DMSState%SchemeOpt = 1

   call cc_dms_init(Config, DMSState, EmisState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_dms_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call cc_dms_run(MetState, DMSState, EmisState, rc)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in _dms_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call assert( DMSState%TotalEmission > 0.0_fp, "Test DMS Emissions")
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
         write(*,*) 'DMSState%CatIndex = ', DMSState_%CatIndex
         write(*,*) 'DMSState%nDMSSpecies = ', DMSState_%nDMSSpecies
         write(*,*) 'DMSState%DMSSpeciesName = ', DMSState_%DMSSpeciesName
         write(*,*) 'DMSState%EmissionPerSpecies = ', DMSState_%EmissionPerSpecies
         write(*,*) 'DMSState%TotalEmission = ', DMSState_%TotalEmission

      end if

   end subroutine print_info


end program test_dms
