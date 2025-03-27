program test_chem
   use CATChem, only: cc_get_micm_version
   use testing_mod, only: assert
   implicit none

   character(len=:), allocatable :: micm_version
   character(len=*), parameter :: expected_micm_version = "3.7.0"

   micm_version = trim(adjustl(cc_get_micm_version()))
   print "('MICM version:', 1x, '''', a, '''')", micm_version
   call assert(micm_version == expected_micm_version, &
      "MICM version should be "//expected_micm_version)

end program test_chem
