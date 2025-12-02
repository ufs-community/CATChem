!> \file ProcessRegistry_Mod.F90
!! \brief Process registration system for dynamic process discovery
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides a registry system for atmospheric processes, enabling
!! dynamic discovery and instantiation of processes at runtime. It supports
!! the plugin architecture outlined in the refactoring plan.
!!
!! \details
!! **Key Features:**
!! - Dynamic process registration and discovery
!! - Category-based process organization
!! - Runtime metadata querying
!! - Process availability checking
!! - Plugin-style architecture support
!!
!! **Usage Example:**
!! \code{.f90}
!! use ProcessRegistry_Mod
!! type(ProcessRegistryType) :: registry
!! class(ProcessInterface), allocatable :: process
!!
!! call registry%init()
!! call registry%register_process('dust_fengsha', create_dust_fengsha)
!! call registry%create_process('dust_fengsha', process, rc)
!! \endcode
!!
module ProcessRegistry_Mod
   use precision_mod
   use error_mod, only : ErrorManagerType, CC_SUCCESS, CC_FAILURE
   use ProcessInterface_Mod, only : ProcessInterface

   implicit none
   private

   public :: ProcessRegistryType
   public :: ProcessCreatorInterface
   public :: get_global_registry

   !> \brief Function pointer interface for process creators
   abstract interface
      subroutine ProcessCreatorInterface(process, rc)
         import :: ProcessInterface
         class(ProcessInterface), allocatable, intent(out) :: process
         integer, intent(out) :: rc
      end subroutine ProcessCreatorInterface
   end interface

   !> \brief Registry entry for a single process
   type :: ProcessRegistryEntry
      character(len=64) :: name = ''                          !< Process name
      character(len=32) :: category = ''                      !< Process category
      character(len=256) :: description = ''                  !< Process description
      procedure(ProcessCreatorInterface), pointer, nopass :: creator => null() !< Creator function
      logical :: is_available = .false.                       !< Availability status
   end type ProcessRegistryEntry

   !> \brief Main process registry for managing available processes
   !!
   !! This type provides centralized management of all available atmospheric
   !! processes in the system. It supports dynamic registration, discovery,
   !! and instantiation following the factory pattern.
   !!
   type :: ProcessRegistryType
      private

      ! Registry storage
      type(ProcessRegistryEntry), allocatable :: entries(:)
      integer :: num_entries = 0
      integer :: max_entries = 100                            !< Maximum registry entries

      ! Registry metadata
      character(len=64) :: registry_name = 'CATChem Process Registry'
      character(len=32) :: version = '2.0'
      logical :: is_initialized = .false.

   contains
      ! Core registry operations
      procedure :: init => registry_init
      procedure :: finalize => registry_finalize
      procedure :: register_process => registry_register_process
      procedure :: unregister_process => registry_unregister_process
      procedure :: create_process => registry_create_process

      ! Query operations
      procedure :: is_process_available => registry_is_process_available
      procedure :: get_process_description => registry_get_process_description
      procedure :: list_processes => registry_list_processes
      procedure :: list_categories => registry_list_categories
      procedure :: get_processes_by_category => registry_get_processes_by_category

      ! Utility operations
      procedure :: get_registry_info => registry_get_registry_info
      procedure :: print_summary => registry_print_summary
      procedure :: validate_registry => registry_validate_registry

      ! Private utilities
      procedure, private :: find_process_index => registry_find_process_index
      procedure, private :: ensure_capacity => registry_ensure_capacity
   end type ProcessRegistryType

   ! Global registry instance
   type(ProcessRegistryType), save, target :: global_registry
   logical, save :: global_registry_initialized = .false.

contains

   !> \brief Get the global process registry instance
   function get_global_registry() result(registry)
      type(ProcessRegistryType), pointer :: registry

      if (.not. global_registry_initialized) then
         call global_registry%init()
         global_registry_initialized = .true.
      endif

      registry => global_registry
   end function get_global_registry

   !> \brief Initialize the process registry
   subroutine registry_init(this, rc)
      implicit none
      class(ProcessRegistryType), intent(inout) :: this
      integer, optional, intent(out) :: rc

      integer :: local_rc

      local_rc = CC_SUCCESS

      if (this%is_initialized) then
         if (present(rc)) rc = CC_SUCCESS
         return
      endif

      ! Allocate storage for registry entries
      if (.not. allocated(this%entries)) then
         allocate(this%entries(this%max_entries), stat=local_rc)
         if (local_rc /= 0) then
            local_rc = CC_FAILURE
         endif
      endif

      this%num_entries = 0
      this%is_initialized = .true.

      if (present(rc)) rc = local_rc
   end subroutine registry_init

   !> \brief Finalize and cleanup the registry
   subroutine registry_finalize(this, rc)
      implicit none
      class(ProcessRegistryType), intent(inout) :: this
      integer, optional, intent(out) :: rc

      if (allocated(this%entries)) then
         deallocate(this%entries)
      endif

      this%num_entries = 0
      this%is_initialized = .false.

      if (present(rc)) rc = CC_SUCCESS
   end subroutine registry_finalize

   !> \brief Register a new process in the registry
   subroutine registry_register_process(this, name, category, description, creator, rc)
      implicit none
      class(ProcessRegistryType), intent(inout) :: this
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: category
      character(len=*), intent(in) :: description
      procedure(ProcessCreatorInterface) :: creator
      integer, optional, intent(out) :: rc

      integer :: local_rc, idx

      local_rc = CC_SUCCESS

      if (.not. this%is_initialized) then
         call this%init(local_rc)
         if (local_rc /= CC_SUCCESS) then
            if (present(rc)) rc = local_rc
            return
         endif
      endif

      ! Check if process already exists
      call this%find_process_index(name, idx)
      if (idx > 0) then
         ! Update existing entry
         this%entries(idx)%category = trim(category)
         this%entries(idx)%description = trim(description)
         this%entries(idx)%creator => creator
         this%entries(idx)%is_available = .true.
      else
         ! Add new entry
         call this%ensure_capacity(local_rc)
         if (local_rc /= CC_SUCCESS) then
            if (present(rc)) rc = local_rc
            return
         endif

         this%num_entries = this%num_entries + 1
         idx = this%num_entries

         this%entries(idx)%name = trim(name)
         this%entries(idx)%category = trim(category)
         this%entries(idx)%description = trim(description)
         this%entries(idx)%creator => creator
         this%entries(idx)%is_available = .true.
      endif

      if (present(rc)) rc = local_rc
   end subroutine registry_register_process

   !> \brief Unregister a process from the registry
   subroutine registry_unregister_process(this, name, rc)
      implicit none
      class(ProcessRegistryType), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, optional, intent(out) :: rc

      integer :: idx, i

      call this%find_process_index(name, idx)
      if (idx > 0) then
         ! Remove entry by shifting remaining entries
         do i = idx, this%num_entries - 1
            this%entries(i) = this%entries(i + 1)
         enddo
         this%num_entries = this%num_entries - 1
         if (present(rc)) rc = CC_SUCCESS
      else
         if (present(rc)) rc = CC_FAILURE
      endif
   end subroutine registry_unregister_process

   !> \brief Create a process instance by name
   subroutine registry_create_process(this, name, process, rc)
      implicit none
      class(ProcessRegistryType), intent(inout) :: this
      character(len=*), intent(in) :: name
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      integer :: idx

      call this%find_process_index(name, idx)
      if (idx > 0) then
         if (this%entries(idx)%is_available) then
            call this%entries(idx)%creator(process, rc)
         else
            rc = CC_FAILURE
         endif
      else
         rc = CC_FAILURE
      endif
   end subroutine registry_create_process

   !> \brief Check if a process is available
   function registry_is_process_available(this, name) result(available)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=*), intent(in) :: name
      logical :: available

      integer :: idx

      call this%find_process_index(name, idx)
      available = (idx > 0 .and. this%entries(idx)%is_available)
   end function registry_is_process_available

   !> \brief Get description for a registered process
   function registry_get_process_description(this, name, rc) result(description)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=*), intent(in) :: name
      integer, optional, intent(out) :: rc
      character(len=256) :: description

      integer :: idx

      call this%find_process_index(name, idx)
      if (idx > 0) then
         description = this%entries(idx)%description
         if (present(rc)) rc = CC_SUCCESS
      else
         description = ''
         if (present(rc)) rc = CC_FAILURE
      endif
   end function registry_get_process_description

   !> \brief List all registered process names
   subroutine registry_list_processes(this, process_names, rc)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: process_names(:)
      integer, optional, intent(out) :: rc

      integer :: i, local_rc

      local_rc = CC_SUCCESS

      if (this%num_entries > 0) then
         allocate(process_names(this%num_entries), stat=local_rc)
         if (local_rc == 0) then
            do i = 1, this%num_entries
               process_names(i) = this%entries(i)%name
            enddo
         else
            local_rc = CC_FAILURE
         endif
      else
         allocate(process_names(0), stat=local_rc)
      endif

      if (present(rc)) rc = local_rc
   end subroutine registry_list_processes

   !> \brief List all process categories
   subroutine registry_list_categories(this, categories, rc)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=32), allocatable, intent(out) :: categories(:)
      integer, optional, intent(out) :: rc

      character(len=32) :: unique_categories(this%num_entries)
      integer :: i, j, num_unique, local_rc
      logical :: found

      local_rc = CC_SUCCESS
      num_unique = 0

      ! Find unique categories
      do i = 1, this%num_entries
         found = .false.
         do j = 1, num_unique
            if (trim(this%entries(i)%category) == trim(unique_categories(j))) then
               found = .true.
               exit
            endif
         enddo
         if (.not. found) then
            num_unique = num_unique + 1
            unique_categories(num_unique) = this%entries(i)%category
         endif
      enddo

      ! Allocate and fill result
      allocate(categories(num_unique), stat=local_rc)
      if (local_rc == 0) then
         categories(1:num_unique) = unique_categories(1:num_unique)
      else
         local_rc = CC_FAILURE
      endif

      if (present(rc)) rc = local_rc
   end subroutine registry_list_categories

   !> \brief Get processes in a specific category
   subroutine registry_get_processes_by_category(this, category, process_names, rc)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=*), intent(in) :: category
      character(len=64), allocatable, intent(out) :: process_names(:)
      integer, optional, intent(out) :: rc

      character(len=64) :: temp_names(this%num_entries)
      integer :: i, count, local_rc

      local_rc = CC_SUCCESS
      count = 0

      ! Find processes in category
      do i = 1, this%num_entries
         if (trim(this%entries(i)%category) == trim(category)) then
            count = count + 1
            temp_names(count) = this%entries(i)%name
         endif
      enddo

      ! Allocate and fill result
      allocate(process_names(count), stat=local_rc)
      if (local_rc == 0 .and. count > 0) then
         process_names(1:count) = temp_names(1:count)
      endif

      if (present(rc)) rc = local_rc
   end subroutine registry_get_processes_by_category

   !> \brief Get registry information
   subroutine registry_get_registry_info(this, name, version, num_processes, rc)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=*), optional, intent(out) :: name
      character(len=*), optional, intent(out) :: version
      integer, optional, intent(out) :: num_processes
      integer, optional, intent(out) :: rc

      if (present(name)) name = this%registry_name
      if (present(version)) version = this%version
      if (present(num_processes)) num_processes = this%num_entries
      if (present(rc)) rc = CC_SUCCESS
   end subroutine registry_get_registry_info

   !> \brief Print registry summary
   subroutine registry_print_summary(this)
      implicit none
      class(ProcessRegistryType), intent(in) :: this

      integer :: i

      write(*,'(A)') '=== Process Registry Summary ==='
      write(*,'(A,A)') 'Registry: ', trim(this%registry_name)
      write(*,'(A,A)') 'Version: ', trim(this%version)
      write(*,'(A,I0)') 'Registered Processes: ', this%num_entries

      if (this%num_entries > 0) then
         write(*,'(A)') ''
         write(*,'(A)') 'Available Processes:'
         do i = 1, this%num_entries
            write(*,'(A,A,A,A,A)') '  ', trim(this%entries(i)%name), &
               ' [', trim(this%entries(i)%category), ']'
         enddo
      endif
      write(*,'(A)') '================================='
   end subroutine registry_print_summary

   !> \brief Validate registry consistency
   function registry_validate_registry(this, rc) result(is_valid)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      integer, optional, intent(out) :: rc
      logical :: is_valid

      integer :: i, j, local_rc

      local_rc = CC_SUCCESS
      is_valid = .true.

      ! Check for duplicate names
      do i = 1, this%num_entries - 1
         do j = i + 1, this%num_entries
            if (trim(this%entries(i)%name) == trim(this%entries(j)%name)) then
               is_valid = .false.
               local_rc = CC_FAILURE
               exit
            endif
         enddo
         if (.not. is_valid) exit
      enddo

      if (present(rc)) rc = local_rc
   end function registry_validate_registry

   !> \brief Find index of process by name
   subroutine registry_find_process_index(this, name, index)
      implicit none
      class(ProcessRegistryType), intent(in) :: this
      character(len=*), intent(in) :: name
      integer, intent(out) :: index

      integer :: i

      index = 0
      do i = 1, this%num_entries
         if (trim(this%entries(i)%name) == trim(name)) then
            index = i
            exit
         endif
      enddo
   end subroutine registry_find_process_index

   !> \brief Ensure registry has capacity for new entries
   subroutine registry_ensure_capacity(this, rc)
      implicit none
      class(ProcessRegistryType), intent(inout) :: this
      integer, intent(out) :: rc

      type(ProcessRegistryEntry), allocatable :: temp_entries(:)
      integer :: new_size

      rc = CC_SUCCESS

      if (this%num_entries >= this%max_entries) then
         ! Grow the registry
         new_size = this%max_entries * 2
         allocate(temp_entries(new_size), stat=rc)
         if (rc == 0) then
            temp_entries(1:this%num_entries) = this%entries(1:this%num_entries)
            call move_alloc(temp_entries, this%entries)
            this%max_entries = new_size
         else
            rc = CC_FAILURE
         endif
      endif
   end subroutine registry_ensure_capacity

end module ProcessRegistry_Mod
