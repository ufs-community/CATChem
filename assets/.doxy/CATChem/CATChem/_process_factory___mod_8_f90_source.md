

# File ProcessFactory\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ProcessFactory\_Mod.F90**](_process_factory___mod_8_f90.md)

[Go to the documentation of this file](_process_factory___mod_8_f90.md)


```Fortran

module processfactory_mod
   use precision_mod
   use statemanager_mod, only : statemanagertype
   use metstate_mod, only : metstatetype
   use error_mod
   use processinterface_mod
   use processregistry_mod

   implicit none
   private

   public :: processfactorytype
   public :: create_process

   type :: processfactorytype
      private
      type(ProcessRegistryType) :: registry
   contains
      procedure :: init => factory_init
      procedure :: create_process => factory_create_process
      procedure :: list_available => factory_list_available
      procedure :: register_process => factory_register_process
      procedure, private :: register_builtin_processes
   end type processfactorytype

contains

   subroutine factory_init(this, rc)
      class(ProcessFactoryType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%registry%init(rc)
      if (rc /= cc_success) return

      ! Register built-in processes
      call this%register_builtin_processes(rc)
   end subroutine factory_init

   subroutine factory_create_process(this, process_name, container, process, rc)
      class(ProcessFactoryType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      type(StateManagerType), intent(inout) :: container
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=32), allocatable :: met_fields(:)
      integer :: i, alloc_rc
      type(MetStateType), pointer :: met_state

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('factory_create_process', &
         'creating process: ' // trim(process_name))

      ! Create process from registry (scheme is read from process configuration)
      call this%registry%create_process(process_name, process, rc)
      if (rc /= cc_success) then
         call error_mgr%report_error(error_invalid_config, &
            'Unknown process: ' // trim(process_name), rc, &
            'factory_create_process', &
            'Check available processes with list_available()')
         call error_mgr%pop_context()
         return
      endif

      ! Allocate only required met fields for this process
      met_fields = process%get_required_met_fields()
      met_state => container%get_met_state_ptr()
      if (associated(met_state) .and. allocated(met_fields)) then
         do i = 1, size(met_fields)
            call met_state%allocate_field(met_fields(i), alloc_rc)
            ! Optionally handle alloc_rc errors here
         end do
      endif

      call error_mgr%pop_context()
   end subroutine factory_create_process

   subroutine factory_list_available(this, process_names, rc)
      class(ProcessFactoryType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: process_names(:)
      integer, intent(out) :: rc

      call this%registry%list_processes(process_names, rc)
   end subroutine factory_list_available

   subroutine factory_register_process(this, name, category, description, creator, rc)
      class(ProcessFactoryType), intent(inout) :: this
      character(len=*), intent(in) :: name, category, description
      procedure(ProcessCreatorInterface) :: creator
      integer, intent(out) :: rc

      call this%registry%register_process(name, category, description, creator, rc)
   end subroutine factory_register_process

   subroutine register_builtin_processes(this, rc)
      class(ProcessFactoryType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      ! Core processes would be registered here, but in the current
      ! architecture, processes register themselves during module
      ! initialization. This method is kept for future core processes
      ! that might be built directly into the factory.

      ! Note: Process modules (dust, seasalt, etc.) should register
      ! themselves using the public register_process method during
      ! their initialization phase.

   end subroutine register_builtin_processes

   function create_process(process_name, container, rc) result(process)
      character(len=*), intent(in) :: process_name
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc
      class(ProcessInterface), allocatable :: process

      type(ProcessFactoryType) :: factory

      call factory%init(rc)
      if (rc /= cc_success) return

      call factory%create_process(process_name, container, process, rc)
   end function create_process

end module processfactory_mod
```


