!> \file test_Precision.f90
!! \brief Test program for Precision module
!!
!!!>
program test_Precision
   use testing_mod, only: assert, assert_close
   use Precision_Mod, only: fp, f4, f8, MISSING, ZERO, ONE, rae

   implicit none

   real(fp) :: test_val
   real(f4) :: test_f4
   real(f8) :: test_f8
   logical :: result

   write(*,*) 'Testing Precision module...'
   write(*,*) ''

   ! Test 1: Precision definitions
   write(*,*) 'Test 1: Precision definitions'
   call assert(fp == f4 .or. fp == f8, "fp should be either f4 or f8")
   call assert(f4 == 4, "f4 should be 4 (kind parameter)")
   call assert(f8 == 8, "f8 should be 8 (kind parameter)")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Special values
   write(*,*) 'Test 2: Special values'
   call assert(ZERO == 0.0_fp, "ZERO should be 0.0")
   call assert(ONE == 1.0_fp, "ONE should be 1.0")
   call assert(MISSING == -999.0_fp, "MISSING should be -999.0")

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Real approximately equal function
   write(*,*) 'Test 3: Real approximately equal function'

   ! Test with f4 values
   test_f4 = 1.0_f4
   result = rae(test_f4, 1.0_f4)
   call assert(result, "RAE should return true for equal f4 values")

   ! Test with f8 values
   test_f8 = 1.0_f8
   result = rae(test_f8, 1.0_f8)
   call assert(result, "RAE should return true for equal f8 values")

   ! Test with slightly different values (should still be true due to tiny tolerance)
   test_f4 = 1.0_f4
   result = rae(test_f4, 1.0_f4 + 1.0e-20_f4)
   call assert(result, "RAE should return true for very close f4 values")

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   write(*,*) 'All Precision module tests passed!'

end program test_Precision
