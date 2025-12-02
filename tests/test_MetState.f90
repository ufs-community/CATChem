program test_MetState
   !> Test program for MetState dimensional field assignment functionality
   !!
   !! This test verifies that the MetState field assignment interface correctly
   !! handles scalar, 2D, and 3D field assignments for REAL, INTEGER, and LOGICAL types.
   !!
   !! Tests:
   !! - Scalar field assignment (broadcasting to arrays)
   !! - 2D field assignment with proper dimensions
   !! - 3D field assignment with proper dimensions
   !! - Error handling for invalid field names
   !! - Type safety for dimensional interfaces

   use testing_mod, only: assert, assert_close
   use MetState_Mod, only: MetStateType
   use Error_Mod, only: ErrorManagerType, CC_SUCCESS
   use Precision_Mod, only: fp
   implicit none

   type(MetStateType) :: metstate
   type(ErrorManagerType), target :: error_mgr_target
   type(ErrorManagerType), pointer :: error_manager
   integer :: rc

   ! Scalar test variables
   real(fp) :: test_scalar_real
   logical :: test_scalar_logical
   integer :: test_scalar_integer

   ! 2D test arrays
   real(fp), allocatable :: test_2d_real(:,:)
   logical, allocatable :: test_2d_logical(:,:)
   integer, allocatable :: test_2d_integer(:,:)

   ! 3D test arrays
   real(fp), allocatable :: test_3d_real(:,:,:)
   integer, allocatable :: test_3d_integer(:,:,:)
   logical, allocatable :: test_3d_logical(:,:,:)

   ! Verification pointers
   real(fp), pointer :: area_ptr(:,:)
   logical, pointer :: land_ptr(:,:)
   integer, pointer :: lwi_ptr(:,:)
   real(fp), pointer :: temp_ptr(:)  ! For column pointer operations

   ! Test variables for multiple fields interface
   character(len=10), allocatable :: empty_field_names(:)

   integer :: nx, ny, nz, i, j, k

   print *, "Starting MetState dimensional field assignment tests..."

   ! Set up error manager
   error_manager => error_mgr_target

   ! Initialize metstate with basic dimensions
   nx = 4; ny = 6; nz = 12
   call metstate%init(nx, ny, nz, 3, 5, 10, error_manager, rc)
   call assert(rc == CC_SUCCESS, "MetState initialization failed")

   print *, "=== Testing Scalar Field Assignment (Broadcasting) ==="

   ! Test setting a REAL scalar field value (broadcasts to all array elements)
   test_scalar_real = 15.5_fp
   call metstate%set_field("AREA_M2", test_scalar_real, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set AREA_M2 scalar")

   ! Verify the scalar was broadcast to all elements
   test_scalar_real = metstate%get_2Dto0D_value("AREA_M2", 1, 1)
   call assert_close(test_scalar_real, 15.5_fp, 1e-6_fp, "AREA_M2[1,1] scalar broadcast")
   test_scalar_real = metstate%get_2Dto0D_value("AREA_M2", nx, ny)
   call assert_close(test_scalar_real, 15.5_fp, 1e-6_fp, "AREA_M2[nx,ny] scalar broadcast")
   print *, "SUCCESS: Scalar REAL field assignment verified"

   ! Test setting a LOGICAL scalar field value
   test_scalar_logical = .true.
   call metstate%set_field("IsLand", test_scalar_logical, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set IsLand scalar")

   test_scalar_logical = metstate%get_2Dto0D_value_logical("IsLand", 1, 1)
   call assert(test_scalar_logical .eqv. .true., "IsLand[1,1] scalar broadcast")
   test_scalar_logical = metstate%get_2Dto0D_value_logical("IsLand", nx, ny)
   call assert(test_scalar_logical .eqv. .true., "IsLand[nx,ny] scalar broadcast")
   print *, "SUCCESS: Scalar LOGICAL field assignment verified"

   ! Test setting an INTEGER scalar field value
   test_scalar_integer = 2
   call metstate%set_field("nSOIL", test_scalar_integer, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set nSOIL scalar")

   test_scalar_integer = metstate%get_scalar_value_int("nSOIL")
   call assert(test_scalar_integer == 2, "nSOIL scalar field retrieval")
   print *, "SUCCESS: Scalar INTEGER field assignment verified"

   print *, "=== Testing 2D Field Assignment ==="

   ! Allocate and populate 2D test arrays
   allocate(test_2d_real(nx, ny))
   allocate(test_2d_logical(nx, ny))
   allocate(test_2d_integer(nx, ny))

   ! Fill 2D arrays with test patterns
   do j = 1, ny
      do i = 1, nx
         test_2d_real(i, j) = real(i*10 + j, fp)
         test_2d_logical(i, j) = mod(i + j, 2) == 0
         test_2d_integer(i, j) = i * 100 + j
      end do
   end do

   ! Test setting 2D REAL field
   call metstate%set_field("FRVEG", test_2d_real, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set FRVEG 2D array")

   ! Test setting 2D LOGICAL field
   call metstate%set_field("IsWater", test_2d_logical, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set IsWater 2D array")

   ! Test setting 2D INTEGER field (overwrites previous scalar)
   call metstate%set_field("LWI", test_2d_integer, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set LWI 2D array")

   ! Verify 2D assignments
   test_scalar_integer = metstate%get_2Dto0D_value_int("LWI", 2, 3)
   call assert(test_scalar_integer == 203, "LWI[2,3] 2D assignment verification")
   test_scalar_integer = metstate%get_2Dto0D_value_int("LWI", nx, ny)
   call assert(test_scalar_integer == nx*100 + ny, "LWI[nx,ny] 2D assignment verification")

   print *, "SUCCESS: 2D field assignments verified"

   print *, "=== Testing 3D Field Assignment ==="

   ! Allocate and populate 3D test arrays
   allocate(test_3d_real(nx, ny, nz))
   allocate(test_3d_integer(nx, ny, nz))
   allocate(test_3d_logical(nx, ny, nz))

   ! Fill 3D arrays with test patterns
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            test_3d_real(i, j, k) = real(i + j*0.1_fp + k*0.01_fp, fp)
            test_3d_integer(i, j, k) = i * 10000 + j * 100 + k
            test_3d_logical(i, j, k) = mod(i + j + k, 3) == 0
         end do
      end do
   end do

   ! Test setting 3D REAL field
   call metstate%set_field("T", test_3d_real, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set T 3D array")

   ! Verify 3D assignment
   temp_ptr => metstate%get_column_ptr_func("T", 1, 1)
   call assert(associated(temp_ptr), "Failed to get T column pointer after 3D assignment")
   call assert_close(temp_ptr(1), 1.11_fp, 1e-5_fp, "T[1,1,1] 3D assignment verification")
   temp_ptr => metstate%get_column_ptr_func("T", 2, 3)
   call assert_close(temp_ptr(4), 2.34_fp, 1e-5_fp, "T[2,3,4] 3D assignment verification")

   print *, "SUCCESS: 3D field assignments verified"

   print *, "=== Testing Error Handling ==="

   ! Test setting an invalid field
   call metstate%set_field("nonexistent_field", test_scalar_real, error_manager, rc)
   call assert(rc /= CC_SUCCESS, "Should have rejected nonexistent_field")

   ! Test wrong dimension assignment (try to assign 3D array to 2D field)
   call metstate%set_field("AREA_M2", test_3d_real, error_manager, rc)
   call assert(rc /= CC_SUCCESS, "Should have rejected 3D array for 2D field")

   print *, "SUCCESS: Error handling verified"

   print *, "=== Testing Type Safety ==="

   ! Test type mismatch (try to assign REAL to INTEGER field)
   call metstate%set_field("LWI", test_scalar_real, error_manager, rc)
   call assert(rc /= CC_SUCCESS, "Should have rejected REAL value for INTEGER field")

   ! Test type mismatch (try to assign INTEGER to LOGICAL field)
   call metstate%set_field("IsLand", test_scalar_integer, error_manager, rc)
   call assert(rc /= CC_SUCCESS, "Should have rejected INTEGER value for LOGICAL field")

   print *, "SUCCESS: Type safety verified"

   print *, "=== Testing Field Access and Retrieval ==="

   ! Test accessing 3D fields through different methods
   temp_ptr => metstate%get_column_ptr_func("T", 3, 4)
   call assert(associated(temp_ptr), "Failed to get T column pointer at [3,4]")
   call assert_close(temp_ptr(5), 3.45_fp, 1e-5_fp, "T[3,4,5] column access verification")

   ! Test 2D field access
   test_scalar_real = metstate%get_2Dto0D_value("AREA_M2", 2, 3)
   call assert_close(test_scalar_real, 15.5_fp, 1e-6_fp, "AREA_M2[2,3] 2D field access")

   test_scalar_integer = metstate%get_2Dto0D_value_int("LWI", 3, 4)
   call assert(test_scalar_integer == 304, "LWI[3,4] 2D integer field access")

   ! Test scalar field access
   test_scalar_integer = metstate%get_scalar_value_int("nSOIL")
   call assert(test_scalar_integer == 2, "nSOIL scalar field access")

   print *, "SUCCESS: Field access and retrieval verified"

   print *, "=== Testing Field Allocation and Memory Management ==="

   ! Test individual field allocation
   call metstate%allocate_field("PMID", rc)
   call assert(rc == CC_SUCCESS, "Failed to allocate PMID field")

   ! Test memory usage calculation
   if (metstate%get_memory_usage() > 0) then
      print *, "SUCCESS: Memory usage calculation working"
   else
      call assert(.false., "Memory usage should be positive")
   endif

   ! Test field deallocation
   call metstate%deallocate_field("PMID", rc)
   call assert(rc == CC_SUCCESS, "Failed to deallocate PMID field")

   print *, "SUCCESS: Field allocation and memory management verified"

   print *, "=== Testing VirtualMet Integration ==="

   ! Test that VirtualMet-relevant fields are accessible
   ! Allocate and test a key atmospheric field
   call metstate%allocate_field("U", rc)
   call assert(rc == CC_SUCCESS, "Failed to allocate U field")

   temp_ptr => metstate%get_column_ptr_func("U", 1, 1)
   call assert(associated(temp_ptr), "U field should be accessible for VirtualMet")

   ! Test surface field access
   test_scalar_real = metstate%get_2Dto0D_value("AREA_M2", 1, 1)
   call assert_close(test_scalar_real, 15.5_fp, 1e-6_fp, "Surface field accessible for VirtualMet")

   print *, "SUCCESS: VirtualMet integration verified"

   print *, "=== Testing Field Validation and Edge Cases ==="

   ! Test accessing unallocated field (deallocate first to ensure it's unallocated)
   call metstate%deallocate_field("OMEGA", rc)
   call assert(rc == CC_SUCCESS, "Failed to deallocate OMEGA field")
   temp_ptr => metstate%get_column_ptr_func("OMEGA", 1, 1)
   call assert(.not. associated(temp_ptr), "Unallocated field should return null pointer")

   ! Test boundary conditions
   temp_ptr => metstate%get_column_ptr_func("T", nx, ny)
   call assert(associated(temp_ptr), "Boundary access [nx,ny] should work")

   ! Test out-of-bounds access (should be clamped)
   temp_ptr => metstate%get_column_ptr_func("T", nx+10, ny+10)
   call assert(associated(temp_ptr), "Out-of-bounds access should be clamped")

   print *, "SUCCESS: Field validation and edge cases verified"

   print *, "=== Testing Multi-Type Field Access ==="

   ! Test INTEGER 2D field access (TropLev is a 2D field)
   if (.not. allocated(test_2d_integer)) then
      allocate(test_2d_integer(nx, ny))
   endif
   do j = 1, ny
      do i = 1, nx
         test_2d_integer(i, j) = i*100 + j
      end do
   end do

   call metstate%set_field("TropLev", test_2d_integer, error_manager, rc)
   call assert(rc == CC_SUCCESS, "Failed to set TropLev 2D integer field")

   ! Test LOGICAL field access
   test_scalar_logical = metstate%get_2Dto0D_value_logical("IsLand", 2, 2)
   call assert(test_scalar_logical .eqv. .true., "IsLand logical field access")

   print *, "SUCCESS: Multi-type field access verified"

   print *, "=== Testing Multiple Fields Interface ==="

   ! Test setting multiple 2D REAL fields at once
   call metstate%set_multiple_fields( &
      field_names=['AREA_M2   ', 'FRVEG     ', 'LAT       '], &
      error_mgr=error_manager, &
      rc=rc, &
      AREA_M2_data=test_2d_real, &
      FRVEG_data=test_2d_real, &
      LAT_data=test_2d_real)
   call assert(rc == CC_SUCCESS, "Failed to set multiple 2D REAL fields")

   ! Verify the multiple field assignments
   test_scalar_real = metstate%get_2Dto0D_value("AREA_M2", 2, 3)
   call assert_close(test_scalar_real, 23.0_fp, 1e-5_fp, "AREA_M2[2,3] from multiple fields")
   test_scalar_real = metstate%get_2Dto0D_value("FRVEG", 3, 4)
   call assert_close(test_scalar_real, 34.0_fp, 1e-5_fp, "FRVEG[3,4] from multiple fields")
   test_scalar_real = metstate%get_2Dto0D_value("LAT", 4, 5)
   call assert_close(test_scalar_real, 45.0_fp, 1e-5_fp, "LAT[4,5] from multiple fields")
   print *, "SUCCESS: Multiple 2D REAL fields set correctly"

   ! Test setting multiple 3D REAL fields at once
   call metstate%set_multiple_fields( &
      field_names=['T         ', 'U         ', 'V         '], &
      error_mgr=error_manager, &
      rc=rc, &
      T_data=test_3d_real, &
      U_data=test_3d_real, &
      V_data=test_3d_real)
   call assert(rc == CC_SUCCESS, "Failed to set multiple 3D REAL fields")

   ! Verify 3D field assignments
   temp_ptr => metstate%get_column_ptr_func("T", 2, 3)
   call assert_close(temp_ptr(4), 2.34_fp, 1e-5_fp, "T[2,3,4] from multiple fields")
   temp_ptr => metstate%get_column_ptr_func("U", 1, 2)
   call assert_close(temp_ptr(3), 1.23_fp, 1e-5_fp, "U[1,2,3] from multiple fields")
   temp_ptr => metstate%get_column_ptr_func("V", 3, 4)
   call assert_close(temp_ptr(5), 3.45_fp, 1e-5_fp, "V[3,4,5] from multiple fields")
   print *, "SUCCESS: Multiple 3D REAL fields set correctly"

   ! Test setting multiple LOGICAL fields at once
   call metstate%set_multiple_fields( &
      field_names=['IsLand    ', 'IsWater   '], &
      error_mgr=error_manager, &
      rc=rc, &
      IsLand_data=test_2d_logical, &
      IsWater_data=test_2d_logical)
   call assert(rc == CC_SUCCESS, "Failed to set multiple LOGICAL fields")

   ! Verify LOGICAL field assignments
   test_scalar_logical = metstate%get_2Dto0D_value_logical("IsLand", 2, 2)
   call assert(test_scalar_logical .eqv. .true., "IsLand[2,2] from multiple fields")
   test_scalar_logical = metstate%get_2Dto0D_value_logical("IsWater", 3, 3)
   call assert(test_scalar_logical .eqv. .true., "IsWater[3,3] from multiple fields")
   print *, "SUCCESS: Multiple LOGICAL fields set correctly"

   ! Test setting multiple INTEGER fields at once (scalar and 2D)
   call metstate%set_multiple_fields( &
      field_names=['NLEVS     ', 'LWI       '], &
      error_mgr=error_manager, &
      rc=rc, &
      NLEVS_data=25, &
      LWI_data=test_2d_integer)
   call assert(rc == CC_SUCCESS, "Failed to set multiple INTEGER fields")

   ! Verify INTEGER field assignments
   test_scalar_integer = metstate%get_scalar_value_int("NLEVS")
   call assert(test_scalar_integer == 25, "NLEVS scalar from multiple fields")
   test_scalar_integer = metstate%get_2Dto0D_value_int("LWI", 2, 3)
   call assert(test_scalar_integer == 203, "LWI[2,3] from multiple fields")
   print *, "SUCCESS: Multiple INTEGER fields set correctly"

   ! Test mixed field types in single call
   call metstate%set_multiple_fields( &
      field_names=['AREA_M2   ', 'IsLand    ', 'NLEVS     ', 'T         '], &
      error_mgr=error_manager, &
      rc=rc, &
      AREA_M2_data=test_2d_real, &
      IsLand_data=test_2d_logical, &
      NLEVS_data=30, &
      T_data=test_3d_real)
   call assert(rc == CC_SUCCESS, "Failed to set mixed field types")

   ! Verify mixed field assignments
   test_scalar_real = metstate%get_2Dto0D_value("AREA_M2", 1, 1)
   call assert_close(test_scalar_real, 11.0_fp, 1e-5_fp, "AREA_M2 from mixed fields")
   test_scalar_logical = metstate%get_2Dto0D_value_logical("IsLand", 2, 4)
   call assert(test_scalar_logical .eqv. .true., "IsLand from mixed fields")
   test_scalar_integer = metstate%get_scalar_value_int("NLEVS")
   call assert(test_scalar_integer == 30, "NLEVS from mixed fields")
   temp_ptr => metstate%get_column_ptr_func("T", 1, 1)
   call assert_close(temp_ptr(1), 1.11_fp, 1e-5_fp, "T from mixed fields")
   print *, "SUCCESS: Mixed field types set correctly"

   ! Test error handling for multiple fields interface
   call metstate%set_multiple_fields( &
      field_names=['AREA_M2   ', 'BadField  '], &
      error_mgr=error_manager, &
      rc=rc, &
      AREA_M2_data=test_2d_real)
   call assert(rc /= CC_SUCCESS, "Should have failed with invalid field name")
   print *, "SUCCESS: Multiple fields error handling verified"

   ! Test partial field setting (only some fields provided)
   call metstate%set_multiple_fields( &
      field_names=['FRVEG     ', 'FRLAND    '], &
      error_mgr=error_manager, &
      rc=rc, &
      FRVEG_data=test_2d_real)
   call assert(rc /= CC_SUCCESS, "Should have failed with missing FRLAND data")
   print *, "SUCCESS: Partial field setting error handling verified"

   ! Test successful partial setting when only requesting fields we provide
   call metstate%set_multiple_fields( &
      field_names=['FRVEG     '], &
      error_mgr=error_manager, &
      rc=rc, &
      FRVEG_data=test_2d_real)
   call assert(rc == CC_SUCCESS, "Should succeed with single field")
   print *, "SUCCESS: Single field via multiple interface verified"

   ! Test bulk transfer scenario (like host model integration)
   call metstate%set_multiple_fields( &
      field_names=['T         ', 'U         ', 'V         ', 'QV        ', 'RH        ', &
      'LAT       ', 'LON       ', 'IsLand    ', 'IsWater   ', 'NLEVS     '], &
      error_mgr=error_manager, &
      rc=rc, &
      T_data=test_3d_real, &
      U_data=test_3d_real, &
      V_data=test_3d_real, &
      QV_data=test_3d_real, &
      RH_data=test_3d_real, &
      LAT_data=test_2d_real, &
      LON_data=test_2d_real, &
      IsLand_data=test_2d_logical, &
      IsWater_data=test_2d_logical, &
      NLEVS_data=25)
   call assert(rc == CC_SUCCESS, "Failed bulk transfer of 10 fields")

   ! Verify a few key fields from bulk transfer
   temp_ptr => metstate%get_column_ptr_func("QV", 2, 2)
   call assert_close(temp_ptr(2), 2.22_fp, 1e-5_fp, "QV from bulk transfer")
   test_scalar_real = metstate%get_2Dto0D_value("LON", 3, 2)
   call assert_close(test_scalar_real, 32.0_fp, 1e-5_fp, "LON from bulk transfer")
   test_scalar_integer = metstate%get_scalar_value_int("NLEVS")
   call assert(test_scalar_integer == 25, "NLEVS from bulk transfer")
   print *, "SUCCESS: Bulk transfer of 10 fields verified"

   ! Test 3D LOGICAL field assignment
   call metstate%set_multiple_fields( &
      field_names=['InPbl     '], &
      error_mgr=error_manager, &
      rc=rc, &
      InPbl_data=test_3d_logical)
   call assert(rc == CC_SUCCESS, "Failed to set 3D LOGICAL field")
   print *, "SUCCESS: 3D LOGICAL field assignment verified"

   ! Test scalar REAL and LOGICAL fields
   call metstate%set_multiple_fields( &
      field_names=['nSOIL     ', 'nSOILTYPE '], &
      error_mgr=error_manager, &
      rc=rc, &
      nSOIL_data=5, &
      nSOILTYPE_data=15)
   call assert(rc == CC_SUCCESS, "Failed to set scalar INTEGER fields")

   test_scalar_integer = metstate%get_scalar_value_int("nSOIL")
   call assert(test_scalar_integer == 5, "nSOIL scalar from multiple interface")
   test_scalar_integer = metstate%get_scalar_value_int("nSOILTYPE")
   call assert(test_scalar_integer == 15, "nSOILTYPE scalar from multiple interface")
   print *, "SUCCESS: Scalar field assignments verified"

   ! Test field name case insensitivity
   call metstate%set_multiple_fields( &
      field_names=['area_m2   ', 'frveg     '], &
      error_mgr=error_manager, &
      rc=rc, &
      AREA_M2_data=test_2d_real, &
      FRVEG_data=test_2d_real)
   call assert(rc == CC_SUCCESS, "Should handle lowercase field names")
   print *, "SUCCESS: Case insensitive field names verified"

   ! Test empty field list (edge case)
   allocate(character(len=10) :: empty_field_names(0))
   call metstate%set_multiple_fields( &
      field_names=empty_field_names, &
      error_mgr=error_manager, &
      rc=rc)
   call assert(rc == CC_SUCCESS, "Should handle empty field list")
   deallocate(empty_field_names)
   print *, "SUCCESS: Empty field list handling verified"

   print *, "=== Multiple Fields Interface Tests Complete ==="

   print *, "=== All MetState Field Assignment and Access Tests Passed! ==="

   ! Cleanup
   if (allocated(test_2d_real)) deallocate(test_2d_real)
   if (allocated(test_2d_logical)) deallocate(test_2d_logical)
   if (allocated(test_2d_integer)) deallocate(test_2d_integer)
   if (allocated(test_3d_real)) deallocate(test_3d_real)
   if (allocated(test_3d_integer)) deallocate(test_3d_integer)
   if (allocated(test_3d_logical)) deallocate(test_3d_logical)

   call metstate%cleanup("ALL", rc)

end program test_MetState
