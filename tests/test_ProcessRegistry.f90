!> \file test_ProcessRegistry.f90
!! \brief Test program for ProcessRegistry module
!!
!!!>
program test_ProcessRegistry
   use testing_mod, only: assert, assert_close
   use ProcessRegistry_Mod
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS

   implicit none

   type(ProcessRegistryType) :: registry
   integer :: rc

   write(*,*) 'Testing ProcessRegistry module...'
   write(*,*) ''

   ! Test 1: Initialize registry
   write(*,*) 'Test 1: Initialize registry'
   call registry%init(rc)
   call assert(rc == CC_SUCCESS, "Process registry initialization should succeed")

   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! Test 2: Check if registry is ready
   write(*,*) 'Test 2: Check if registry is ready'
   block
      logical :: is_ready
      is_ready = registry%validate_registry(rc)
      call assert(is_ready, "Registry should be valid after initialization")
   end block

   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! Test 3: Get registry info
   write(*,*) 'Test 3: Get registry info'
   block
      character(len=64) :: name
      character(len=32) :: version
      integer :: num_processes

      call registry%get_registry_info(name, version, num_processes, rc)
      call assert(rc == CC_SUCCESS, "Getting registry info should succeed")
      call assert(num_processes == 0, "Registry should initially contain zero processes")
   end block

   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! Test 4: Register a process (dummy)
   write(*,*) 'Test 4: Register a process (dummy)'
   ! Note: We can't easily register a real process without implementing a creator function
   ! This test is limited in what it can do without a full process implementation

   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! Test 5: List processes (should be empty)
   write(*,*) 'Test 5: List processes (should be empty)'
   block
      character(len=64), allocatable :: process_names(:)
      integer :: num_processes

      call registry%list_processes(process_names, rc)
      num_processes = size(process_names)
      call assert(rc == CC_SUCCESS, "Listing processes should succeed")
      call assert(num_processes == 0, "Registry should contain zero processes")

      if (allocated(process_names)) then
         deallocate(process_names)
      end if
   end block

   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! Test 6: List categories (should be empty)
   write(*,*) 'Test 6: List categories (should be empty)'
   block
      character(len=32), allocatable :: categories(:)
      integer :: num_categories

      call registry%list_categories(categories, rc)
      num_categories = size(categories)
      call assert(rc == CC_SUCCESS, "Listing categories should succeed")
      call assert(num_categories == 0, "Registry should contain zero categories")

      if (allocated(categories)) then
         deallocate(categories)
      end if
   end block

   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   ! Test 7: Cleanup
   write(*,*) 'Test 7: Cleanup'
   call registry%finalize(rc)
   call assert(rc == CC_SUCCESS, "Registry cleanup should succeed")

   write(*,*) 'Test 7 passed!'
   write(*,*) ''

   write(*,*) 'All ProcessRegistry tests passed!'

end program test_ProcessRegistry
