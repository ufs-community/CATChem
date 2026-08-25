!> \file test_UnitConversion.f90
!! \brief Property test for unit conversion relocation equivalence
!!
!! Tests that standalone convert_process_concentration_units and
!! convert_process_flux_units in UnitConversion_Mod produce correct
!! results for various unit pairs, molecular weights, temperatures,
!! and pressures.
!!
!! **Validates: Requirements 3.5**
!!
program test_UnitConversion
   use precision_mod, only: fp
   use testing_mod, only: assert, assert_close
   use UnitConversion_Mod, only: convert_process_concentration_units, &
      convert_process_flux_units
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   use Constants, only: RSTARG, AVO

   implicit none

   ! Relative tolerance appropriate for fp precision
   real(fp), parameter :: REL_TOL = 1.0e-5_fp

   write(*,*) 'Testing UnitConversion_Mod relocated procedures...'
   write(*,*) ''

   call test_concentration_identity_conversion()
   call test_flux_identity_conversion()
   call test_ppbv_to_ppmv()
   call test_ppbv_to_molec_cm3()
   call test_ppbv_to_ug_m3()
   call test_molec_cm3_to_ppbv()
   call test_molec_cm3_to_ug_m3()
   call test_flux_kg_to_molec()
   call test_flux_molec_to_kg()
   call test_various_molecular_weights()
   call test_various_temperatures_pressures()
   call test_unsupported_units_return_failure()

   write(*,*) 'All UnitConversion tests passed!'

contains

   ! ================================================================
   ! Test: Identity conversion (same units -> no change)
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_concentration_identity_conversion()
      real(fp) :: vals(3)
      integer :: rc

      write(*,*) 'Test: Concentration identity conversion (same units)'

      vals = (/ 1.0_fp, 50.0_fp, 1000.0_fp /)
      call convert_process_concentration_units(vals, 'ppbv', 'ppbv', &
         molecular_weight=48.0_fp, temperature=298.0_fp, pressure=101325.0_fp, rc=rc)
      call assert(rc == CC_SUCCESS, "Identity ppbv->ppbv should succeed")
      call assert_close(vals(1), 1.0_fp, 1.0e-12_fp, "ppbv identity val 1")
      call assert_close(vals(2), 50.0_fp, 1.0e-10_fp, "ppbv identity val 2")
      call assert_close(vals(3), 1000.0_fp, 1.0e-9_fp, "ppbv identity val 3")

      vals = (/ 1.0e10_fp, 2.5e12_fp, 7.0e8_fp /)
      call convert_process_concentration_units(vals, 'molec/cm3', 'molec/cm3', &
         molecular_weight=48.0_fp, temperature=298.0_fp, pressure=101325.0_fp, rc=rc)
      call assert(rc == CC_SUCCESS, "Identity molec/cm3->molec/cm3 should succeed")
      call assert_close(vals(1), 1.0e10_fp, 1.0_fp, "molec/cm3 identity val 1")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_concentration_identity_conversion

   ! ================================================================
   ! Test: Flux identity conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_flux_identity_conversion()
      real(fp) :: vals(2)
      integer :: rc

      write(*,*) 'Test: Flux identity conversion (same units)'

      vals = (/ 1.0e-8_fp, 5.0e-6_fp /)
      call convert_process_flux_units(vals, 'kg/m2/s', 'kg/m2/s', &
         molecular_weight=48.0_fp, rc=rc)
      call assert(rc == CC_SUCCESS, "Identity kg/m2/s->kg/m2/s should succeed")
      call assert_close(vals(1), 1.0e-8_fp, 1.0e-20_fp, "flux identity val 1")
      call assert_close(vals(2), 5.0e-6_fp, 1.0e-18_fp, "flux identity val 2")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_flux_identity_conversion

   ! ================================================================
   ! Test: ppbv -> ppmv conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_ppbv_to_ppmv()
      real(fp) :: vals(3)
      integer :: rc

      write(*,*) 'Test: ppbv -> ppmv conversion'

      vals = (/ 1000.0_fp, 500.0_fp, 1.0_fp /)
      call convert_process_concentration_units(vals, 'ppbv', 'ppmv', &
         molecular_weight=48.0_fp, temperature=298.0_fp, pressure=101325.0_fp, rc=rc)
      call assert(rc == CC_SUCCESS, "ppbv->ppmv should succeed")
      ! 1000 ppbv = 1 ppmv
      call assert_close(vals(1), 1.0_fp, REL_TOL, "1000 ppbv = 1 ppmv")
      call assert_close(vals(2), 0.5_fp, REL_TOL, "500 ppbv = 0.5 ppmv")
      call assert_close(vals(3), 0.001_fp, REL_TOL, "1 ppbv = 0.001 ppmv")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_ppbv_to_ppmv

   ! ================================================================
   ! Test: ppbv -> molec/cm3 conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_ppbv_to_molec_cm3()
      real(fp) :: vals(1), expected
      real(fp) :: temp, pres
      integer :: rc

      write(*,*) 'Test: ppbv -> molec/cm3 conversion'

      temp = 298.0_fp
      pres = 101325.0_fp
      vals = (/ 1.0_fp /)

      call convert_process_concentration_units(vals, 'ppbv', 'molec/cm3', &
         molecular_weight=48.0_fp, temperature=temp, pressure=pres, rc=rc)
      call assert(rc == CC_SUCCESS, "ppbv->molec/cm3 should succeed")

      ! Expected: 1 ppbv * (P/RT) * 1e-9 * NA * 1e-6
      expected = (pres / (RSTARG * temp)) * 1.0e-9_fp * AVO * 1.0e-6_fp
      call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
         "ppbv->molec/cm3 value check")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_ppbv_to_molec_cm3

   ! ================================================================
   ! Test: ppbv -> ug/m3 conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_ppbv_to_ug_m3()
      real(fp) :: vals(1), expected
      real(fp) :: temp, pres, mw
      integer :: rc

      write(*,*) 'Test: ppbv -> ug/m3 conversion'

      temp = 298.0_fp
      pres = 101325.0_fp
      mw = 48.0_fp  ! O3
      vals = (/ 100.0_fp /)

      call convert_process_concentration_units(vals, 'ppbv', 'ug/m3', &
         molecular_weight=mw, temperature=temp, pressure=pres, rc=rc)
      call assert(rc == CC_SUCCESS, "ppbv->ug/m3 should succeed")

      ! Expected: 100 * (P/RT) * MW * 1e-3
      expected = 100.0_fp * (pres / (RSTARG * temp)) * mw * 1.0e-3_fp
      call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
         "ppbv->ug/m3 value check")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_ppbv_to_ug_m3

   ! ================================================================
   ! Test: molec/cm3 -> ppbv conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_molec_cm3_to_ppbv()
      real(fp) :: vals(1), expected
      real(fp) :: temp, pres
      integer :: rc

      write(*,*) 'Test: molec/cm3 -> ppbv conversion'

      temp = 298.0_fp
      pres = 101325.0_fp
      vals = (/ 2.46e10_fp /)

      call convert_process_concentration_units(vals, 'molec/cm3', 'ppbv', &
         molecular_weight=48.0_fp, temperature=temp, pressure=pres, rc=rc)
      call assert(rc == CC_SUCCESS, "molec/cm3->ppbv should succeed")

      ! Expected: val * (RT/P) * 1e9 / NA * 1e6
      expected = 2.46e10_fp * (RSTARG * temp / pres) * 1.0e9_fp / AVO * 1.0e6_fp
      call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
         "molec/cm3->ppbv value check")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_molec_cm3_to_ppbv

   ! ================================================================
   ! Test: molec/cm3 -> ug/m3 conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_molec_cm3_to_ug_m3()
      real(fp) :: vals(1), expected
      real(fp) :: mw
      integer :: rc

      write(*,*) 'Test: molec/cm3 -> ug/m3 conversion'

      mw = 48.0_fp
      vals = (/ 2.46e10_fp /)

      call convert_process_concentration_units(vals, 'molec/cm3', 'ug/m3', &
         molecular_weight=mw, temperature=298.0_fp, pressure=101325.0_fp, rc=rc)
      call assert(rc == CC_SUCCESS, "molec/cm3->ug/m3 should succeed")

      ! Expected: val * MW / NA * 1e12
      expected = 2.46e10_fp * mw / AVO * 1.0e12_fp
      call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
         "molec/cm3->ug/m3 value check")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_molec_cm3_to_ug_m3

   ! ================================================================
   ! Test: kg/m2/s -> molec/cm2/s flux conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_flux_kg_to_molec()
      real(fp) :: vals(1), expected
      real(fp) :: mw
      integer :: rc

      write(*,*) 'Test: kg/m2/s -> molec/cm2/s flux conversion'

      mw = 48.0_fp
      vals = (/ 1.0e-8_fp /)

      call convert_process_flux_units(vals, 'kg/m2/s', 'molec/cm2/s', &
         molecular_weight=mw, rc=rc)
      call assert(rc == CC_SUCCESS, "kg/m2/s->molec/cm2/s should succeed")

      ! Expected: val * 1000 * (1/MW) * NA * 1e-4
      expected = 1.0e-8_fp * 1000.0_fp * (1.0_fp / mw) * AVO * 1.0e-4_fp
      call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
         "kg/m2/s->molec/cm2/s value check")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_flux_kg_to_molec

   ! ================================================================
   ! Test: molec/cm2/s -> kg/m2/s flux conversion
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_flux_molec_to_kg()
      real(fp) :: vals(1), expected
      real(fp) :: mw
      integer :: rc

      write(*,*) 'Test: molec/cm2/s -> kg/m2/s flux conversion'

      mw = 48.0_fp
      vals = (/ 1.0e12_fp /)

      call convert_process_flux_units(vals, 'molec/cm2/s', 'kg/m2/s', &
         molecular_weight=mw, rc=rc)
      call assert(rc == CC_SUCCESS, "molec/cm2/s->kg/m2/s should succeed")

      ! Expected: val * (1/NA) * MW * 1e-3 * 1e4
      expected = 1.0e12_fp * (1.0_fp / AVO) * mw * 1.0e-3_fp * 1.0e4_fp
      call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
         "molec/cm2/s->kg/m2/s value check")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_flux_molec_to_kg

   ! ================================================================
   ! Test: Various molecular weights
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_various_molecular_weights()
      real(fp) :: vals(1), expected
      real(fp) :: mws(4)
      real(fp) :: temp, pres
      integer :: rc, i

      write(*,*) 'Test: Various molecular weights for ppbv->ug/m3'

      temp = 298.0_fp
      pres = 101325.0_fp
      ! O3=48, NO2=46, SO2=64.1, NH3=17
      mws = (/ 48.0_fp, 46.0_fp, 64.1_fp, 17.0_fp /)

      do i = 1, 4
         vals = (/ 10.0_fp /)
         call convert_process_concentration_units(vals, 'ppbv', 'ug/m3', &
            molecular_weight=mws(i), temperature=temp, pressure=pres, rc=rc)
         call assert(rc == CC_SUCCESS, "ppbv->ug/m3 should succeed for various MW")

         expected = 10.0_fp * (pres / (RSTARG * temp)) * mws(i) * 1.0e-3_fp
         call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
            "ppbv->ug/m3 with varying MW")
      end do

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_various_molecular_weights

   ! ================================================================
   ! Test: Various temperatures and pressures
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_various_temperatures_pressures()
      real(fp) :: vals(1), expected
      real(fp) :: temps(3), press(3)
      integer :: rc, i

      write(*,*) 'Test: Various temperatures and pressures for ppbv->molec/cm3'

      temps = (/ 220.0_fp, 298.0_fp, 350.0_fp /)
      press = (/ 50000.0_fp, 101325.0_fp, 110000.0_fp /)

      do i = 1, 3
         vals = (/ 5.0_fp /)
         call convert_process_concentration_units(vals, 'ppbv', 'molec/cm3', &
            molecular_weight=48.0_fp, temperature=temps(i), pressure=press(i), rc=rc)
         call assert(rc == CC_SUCCESS, "ppbv->molec/cm3 should succeed for various T,P")

         expected = 5.0_fp * (press(i) / (RSTARG * temps(i))) * &
            1.0e-9_fp * AVO * 1.0e-6_fp
         call assert_close(vals(1), expected, abs(expected) * REL_TOL, &
            "ppbv->molec/cm3 with varying T,P")
      end do

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_various_temperatures_pressures

   ! ================================================================
   ! Test: Unsupported unit pairs return CC_FAILURE
   ! **Validates: Requirements 3.5**
   ! ================================================================
   subroutine test_unsupported_units_return_failure()
      real(fp) :: vals(1)
      integer :: rc

      write(*,*) 'Test: Unsupported unit pairs return CC_FAILURE'

      vals = (/ 1.0_fp /)
      call convert_process_concentration_units(vals, 'ppbv', 'kg/m3', &
         molecular_weight=48.0_fp, temperature=298.0_fp, pressure=101325.0_fp, rc=rc)
      call assert(rc == CC_FAILURE, "ppbv->kg/m3 should return CC_FAILURE")

      vals = (/ 1.0_fp /)
      call convert_process_flux_units(vals, 'kg/m2/s', 'g/m2/s', &
         molecular_weight=48.0_fp, rc=rc)
      call assert(rc == CC_FAILURE, "kg/m2/s->g/m2/s should return CC_FAILURE")

      write(*,*) 'Test passed!'
      write(*,*) ''
   end subroutine test_unsupported_units_return_failure

end program test_UnitConversion
