program test_bvoc
   use CATChem, fp => cc_rk
   use testing_mod, only: assert
   use state_mod
   use ccpr_bvoc_common_mod, only : BvocStateType
   !use EmisState_Mod, only: Emis_Allocate

   implicit none


   type(BvocStateType) :: BvocState

   ! Integers
   INTEGER:: rc          ! Success or failure

   character(len=:), allocatable :: title

   integer :: c ! Loop counter for emission Cats
   integer :: s ! Loop counter for emitted species

   ! Error handling
   CHARACTER(LEN=512) :: errMsg
   CHARACTER(LEN=255) :: thisLoc
   CHARACTER(LEN=255), PARAMETER :: configFile ='Configs/Default/CATChem_config.yml'

   thisLoc = 'test_bvoc -> at read CATChem_Config.yml'
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
   call cc_read_config(Config, GridState, EmisState, ChemState, rc,configFile)
   if (rc /= CC_success) then
      errMsg = 'Error reading configuration file: ' // TRIM( configFile )
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   endif
   title = 'BVOC Test 1 | Read Config'
   write(*,*) 'title = ', title
   write(*,*) 'Config%bvoc_activate = ', Config%bvoc_activate
   write(*,*) 'Config%bvoc_scheme = ', Config%bvoc_scheme
   write(*,*) 'Config%megan_co2_inhib = ', Config%megan_co2_inhib
   write(*,*) 'Config%megan_co2_conc_ppm = ', Config%megan_co2_conc_ppm
   write(*,*) 'EmisState%nCats = ', EmisState%nCats
   !write(*,*) 'EmisState%Cats = ', EmisState%Cats  !cannot write allocatable variables here; need to put in the subroutin at the bottom

   !call Emis_Allocate(GridState, EmisState, RC) !Not sure why cannot call it even if I made it public and used EmisStateMod in this module
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

   ! Meteorological State (get from 20190620 18:00:00 HEMCO log out)
   MetState%TSTEP= 3600
   MetState%LAI= 3.9094312191009521
   allocate(MetState%PFT_16(16))
   MetState%PFT_16= (/0.00, 0.11120668053627014, 0.00, 0.00, 0.00, 0.00, 0.00, 0.35108909010887146, &
      0.00, 0.00, 0.00, 0.00, 0.00,0.18369837105274200,6.9862455129623413E-002,0.26875913143157959/)
   MetState%PMISOLAI= 3.90559316  ! needs to multiply PFTSUM
   MetState%Q_DIR_2= 368.02691650390625
   MetState%Q_DIFF_2= 55.688854217529297
   MetState%PARDR_LASTXDAYS= 69.8014755
   MetState%PARDF_LASTXDAYS= 37.4157791
   MetState%TS= 300.42892456054688
   MetState%T_LASTXDAYS= 294.497833
   MetState%T_LAST24H= 294.919861
   MetState%GWETROOT= 0.76981580257415771
   MetState%SUNCOS= 0.96700188446067026
   MetState%LAT=  38.00
   MetState%DOY= 171
   MetState%LocalHour= 12.00
   MetState%D_BTW_M=  1.00
   !These are read from file. Keep them here for now
   !MetState%AEF_ISOP= 1.8055753025901623E-009
   !MetState%AEF_MBOX= 4.9856540616165277E-013
   !MetState%AEF_BPIN= 1.3127712426530056E-011
   !MetState%AEF_CARE= 3.6563258621377127E-012
   !MetState%AEF_LIMO= 6.3985420662872528E-012
   !MetState%AEF_OCIM= 3.0705341449874024E-011
   !MetState%AEF_SABI= 1.0054971341792413E-011

   ! Allocate DiagState
   call cc_allocate_diagstate(Config, DiagState, ChemState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_allocate_diagstate'
      stop 1
   endif

   title = "BVOC Test 2 | Test each species"
   Config%bvoc_activate = .TRUE.

   !call cc_bvoc_init(Config, ChemState, EmisState, BvocState, RC)
   call cc_bvoc_init(Config, EmisState, BvocState, RC)
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_bvoc_init'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   !call cc_bvoc_run(MetState, EmisState, DiagState, BvocState, ChemState, RC )
   call cc_bvoc_run(MetState, EmisState,  DiagState, BvocState, RC )
   if (rc /= CC_SUCCESS) then
      errMsg = 'Error in cc_bvoc_run'
      call cc_emit_error(errMsg, rc, thisLoc)
      stop 1
   end if

   call print_info(Config, BvocState, MetState, DiagState, title)
   call assert(BvocState%TotalEmission > 1.0e-9_fp, "Test BVOC species")
   !add check to each species
   do s = 1, BvocState%nBvocSpecies
      call assert(BvocState%EmissionPerSpecies(s) > 1.0e-12_fp, &
         "Test lower emission limit for species " // TRIM(BvocState%BvocSpeciesName(s)))
      call assert(BvocState%EmissionPerSpecies(s) < 1.0e-9_fp, &
         "Test upper emission limit for species " // TRIM(BvocState%BvocSpeciesName(s)))
   end do


contains

   subroutine print_info(Config_, BvocState_, MetState_, DiagState_, title_)
      type(ConfigType), intent(in) :: Config_
      type(MetStateType), intent(in) :: MetState_
      type(BvocStateType), intent(in) :: BvocState_
      type(DiagStatetype), INTENT(in) :: DiagState_
      character(len=*), intent(in) :: title_

      write(*,*) '======================================='
      write(*,*) title_
      write(*,*) '======================================='
      write(*,*) '*************'
      write(*,*) 'Configuration '
      write(*,*) '*************'
      write(*,*) 'Config%bvoc_activate = ', Config_%bvoc_activate
      write(*,*) 'BvocState%activate = ', BvocState_%activate
      write(*,*) 'BvocState%CatIndex = ', BvocState_%CatIndex
      write(*,*) 'BvocState%CO2Inhib = ', BvocState_%CO2Inhib
      write(*,*) 'BvocState%CO2conc  = ', BvocState_%CO2conc
      write(*,*) 'MetState%LAI =', MetState_%LAI
      write(*,*) 'MetState%DOY =', MetState_%DOY
      !write(*,*) 'MetState%AEF_ISOP =', MetState_%AEF_ISOP
      write(*,*) 'MetState%PFT_16 =', MetState_%PFT_16
      write(*,*) 'BvocState%BvocSpeciesName=', BvocState_%BvocSpeciesName
      write(*,*) 'BvocState%EmissionPerSpecies=', BvocState_%EmissionPerSpecies
      write(*,*) 'BvocState%TotalEmission = ', BvocState_%TotalEmission
      write(*,*) 'DiagState%PARDR_LASTXDAYS = ', DiagState_%PARDR_LASTXDAYS

   end subroutine print_info

end program test_bvoc
