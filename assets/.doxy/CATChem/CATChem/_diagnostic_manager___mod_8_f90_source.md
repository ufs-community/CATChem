

# File DiagnosticManager\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**DiagnosticManager\_Mod.F90**](_diagnostic_manager___mod_8_f90.md)

[Go to the documentation of this file](_diagnostic_manager___mod_8_f90.md)


```Fortran

module diagnosticmanager_mod
   use precision_mod, only: fp
   use error_mod, only: errormanagertype, cc_success, cc_failure, &
      error_invalid_input, error_not_found, error_memory_allocation
   use diagnosticinterface_mod, only: diagnosticregistrytype, diagnosticfieldtype, &
      diagnosticdatatype, diag_real_scalar, diag_real_1d, &
      diag_real_2d, diag_real_3d
   ! Removed StateManager_Mod import to break circular dependency

   implicit none
   private

   public :: diagnosticmanagertype

   type :: diagnosticmanagertype
      private

      ! Error management
      type(ErrorManagerType), pointer :: error_mgr => null()

      ! Process registries management
      type(DiagnosticRegistryType), allocatable :: process_registries(:)
      character(len=64), allocatable :: process_names(:)
      integer :: num_processes = 0
      integer :: max_processes = 50

      ! Configuration and output management
      logical :: output_enabled = .true.
      character(len=256) :: output_prefix = 'catchem_diag'
      character(len=256) :: output_directory = './'
      integer :: output_frequency = 1  ! timesteps

      ! Collection state
      logical :: collection_enabled = .true.
      integer :: current_timestep = 0

   contains
      ! Initialization and management
      procedure :: init => diagnostic_manager_init
      procedure :: finalize => diagnostic_manager_finalize
      procedure :: reset => diagnostic_manager_reset

      ! Process registry management
      procedure :: register_process => diagnostic_manager_register_process
      procedure :: get_process_registry => diagnostic_manager_get_process_registry
      procedure :: remove_process => diagnostic_manager_remove_process
      procedure :: list_processes => diagnostic_manager_list_processes

      ! Configuration management
      procedure :: configure_output => diagnostic_manager_configure_output
      procedure :: set_output_frequency => diagnostic_manager_set_output_frequency
      procedure :: enable_collection => diagnostic_manager_enable_collection
      procedure :: disable_collection => diagnostic_manager_disable_collection

      ! Diagnostic collection and output
      procedure :: collect_all_diagnostics => diagnostic_manager_collect_all
      procedure :: collect_process_diagnostics => diagnostic_manager_collect_process
      procedure :: get_field_value => diagnostic_manager_get_field_value
      procedure :: write_output => diagnostic_manager_write_output
      procedure :: advance_timestep => diagnostic_manager_advance_timestep

      ! Utility methods
      procedure :: get_total_diagnostics => diagnostic_manager_get_total_diagnostics
      procedure :: print_summary => diagnostic_manager_print_summary
      procedure :: validate_state => diagnostic_manager_validate_state

   end type diagnosticmanagertype

contains

   subroutine diagnostic_manager_init(this, error_mgr, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      type(ErrorManagerType), intent(inout), target :: error_mgr
      integer, intent(out) :: rc

      rc = cc_success

      ! Store ErrorManager reference for internal use
      this%error_mgr => error_mgr

      ! Initialize diagnostic manager
      call this%error_mgr%push_context('diagnostic_manager_init', 'Initializing diagnostic manager')

      ! Allocate process arrays
      if (.not. allocated(this%process_registries)) then
         allocate(this%process_registries(this%max_processes), stat=rc)
         if (rc /= 0) then
            call this%error_mgr%report_error(error_memory_allocation, &
               'Failed to allocate process registries', rc)
            call this%error_mgr%pop_context()
            return
         endif
      endif

      if (.not. allocated(this%process_names)) then
         allocate(this%process_names(this%max_processes), stat=rc)
         if (rc /= 0) then
            call this%error_mgr%report_error(error_memory_allocation, &
               'Failed to allocate process names', rc)
            call this%error_mgr%pop_context()
            return
         endif
      endif

      ! Initialize counters
      this%num_processes = 0
      this%current_timestep = 0

      ! TODO: Read configuration from container's ConfigManager
      ! For now, use defaults

      call this%error_mgr%pop_context()

   end subroutine diagnostic_manager_init

   subroutine diagnostic_manager_finalize(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = cc_success

      ! Finalize all process registries
      do i = 1, this%num_processes
         call this%process_registries(i)%finalize(local_rc)
         if (local_rc /= cc_success) then
            rc = local_rc
         endif
      enddo

      ! Deallocate arrays
      if (allocated(this%process_registries)) then
         deallocate(this%process_registries)
      endif

      if (allocated(this%process_names)) then
         deallocate(this%process_names)
      endif

      this%num_processes = 0

   end subroutine diagnostic_manager_finalize

   subroutine diagnostic_manager_reset(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = cc_success

      ! Reset all process registries
      do i = 1, this%num_processes
         call this%process_registries(i)%reset(local_rc)
         if (local_rc /= cc_success) then
            rc = local_rc
         endif
      enddo

      this%current_timestep = 0

   end subroutine diagnostic_manager_reset

   subroutine diagnostic_manager_register_process(this, process_name, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      integer, intent(out) :: rc

      integer :: i

      rc = cc_success

      ! Check if process already registered
      do i = 1, this%num_processes
         if (trim(this%process_names(i)) == trim(process_name)) then
            call this%error_mgr%report_error(error_invalid_input, &
               'Process already registered: ' // trim(process_name), rc)
            return
         endif
      enddo

      ! Check capacity
      if (this%num_processes >= this%max_processes) then
         call this%error_mgr%report_error(error_memory_allocation, &
            'Maximum number of processes reached', rc)
         return
      endif

      ! Add new process
      this%num_processes = this%num_processes + 1
      this%process_names(this%num_processes) = trim(process_name)

      ! Initialize process registry
      call this%process_registries(this%num_processes)%init(process_name, this%error_mgr, rc)

   end subroutine diagnostic_manager_register_process

   subroutine diagnostic_manager_get_process_registry(this, process_name, registry, rc)
      class(DiagnosticManagerType), intent(inout), target :: this
      character(len=*), intent(in) :: process_name
      type(DiagnosticRegistryType), pointer, intent(out) :: registry
      integer, intent(out) :: rc

      integer :: i

      rc = cc_failure
      nullify(registry)

      ! Find process
      do i = 1, this%num_processes
         if (trim(this%process_names(i)) == trim(process_name)) then
            registry => this%process_registries(i)
            rc = cc_success
            return
         endif
      enddo

   end subroutine diagnostic_manager_get_process_registry

   subroutine diagnostic_manager_remove_process(this, process_name, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      integer, intent(out) :: rc

      integer :: i, j, local_rc

      rc = cc_failure

      ! Find and remove process
      do i = 1, this%num_processes
         if (trim(this%process_names(i)) == trim(process_name)) then
            ! Finalize registry
            call this%process_registries(i)%finalize(local_rc)

            ! Shift remaining processes
            do j = i, this%num_processes - 1
               this%process_names(j) = this%process_names(j + 1)
               this%process_registries(j) = this%process_registries(j + 1)
            enddo

            this%num_processes = this%num_processes - 1
            rc = cc_success
            return
         endif
      enddo

   end subroutine diagnostic_manager_remove_process

   subroutine diagnostic_manager_list_processes(this, process_list, num_processes, rc)
      class(DiagnosticManagerType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: process_list(:)
      integer, intent(out) :: num_processes
      integer, intent(out) :: rc

      rc = cc_success
      num_processes = this%num_processes

      if (num_processes > 0) then
         allocate(process_list(num_processes), stat=rc)
         if (rc /= 0) then
            rc = error_memory_allocation
            return
         endif

         process_list(1:num_processes) = this%process_names(1:num_processes)
      endif

   end subroutine diagnostic_manager_list_processes

   subroutine diagnostic_manager_configure_output(this, rc, output_prefix, &
      output_directory, output_frequency)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc
      character(len=*), optional, intent(in) :: output_prefix
      character(len=*), optional, intent(in) :: output_directory
      integer, optional, intent(in) :: output_frequency

      rc = cc_success

      if (present(output_prefix)) then
         this%output_prefix = trim(output_prefix)
      endif

      if (present(output_directory)) then
         this%output_directory = trim(output_directory)
      endif

      if (present(output_frequency)) then
         this%output_frequency = output_frequency
      endif

   end subroutine diagnostic_manager_configure_output

   subroutine diagnostic_manager_set_output_frequency(this, frequency, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(in) :: frequency
      integer, intent(out) :: rc

      rc = cc_success

      if (frequency > 0) then
         this%output_frequency = frequency
      else
         rc = error_invalid_input
      endif

   end subroutine diagnostic_manager_set_output_frequency

   subroutine diagnostic_manager_enable_collection(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success
      this%collection_enabled = .true.

   end subroutine diagnostic_manager_enable_collection

   subroutine diagnostic_manager_disable_collection(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success
      this%collection_enabled = .false.

   end subroutine diagnostic_manager_disable_collection

   subroutine diagnostic_manager_collect_all(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc, total_fields
      character(len=64) :: current_process

      rc = cc_success

      if (.not. this%collection_enabled) return

      call this%error_mgr%push_context('diagnostic_manager_collect_all', &
         'Collecting diagnostics from all processes')

      total_fields = 0

      ! Collect from each registered process
      do i = 1, this%num_processes
         current_process = trim(this%process_names(i))

         call this%collect_process_diagnostics(current_process, local_rc)
         if (local_rc /= cc_success) then
            call this%error_mgr%report_error(local_rc, &
               'Failed to collect diagnostics from: ' // &
               trim(current_process), local_rc)
            ! Continue with other processes - don't let one failure stop collection
         else
            ! Count fields collected from this process
            total_fields = total_fields + this%process_registries(i)%get_field_count()
         endif
      enddo

      ! Log collection summary
      if (this%num_processes > 0) then
         write(*,'(A,I0,A,I0,A)') 'DiagnosticManager: Collected diagnostics from ', &
            this%num_processes, ' processes (', total_fields, ' total fields)'
      else
         write(*,'(A)') 'DiagnosticManager: No processes registered for diagnostic collection'
      end if

      call this%error_mgr%pop_context()

   end subroutine diagnostic_manager_collect_all

   subroutine diagnostic_manager_collect_process(this, process_name, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      integer, intent(out) :: rc

      type(DiagnosticRegistryType), pointer :: registry
      character(len=64), allocatable :: field_names(:)
      integer :: num_fields, i, local_rc, data_type
      real(fp) :: scalar_value
      real(fp), pointer :: array_1d_ptr(:) => null()
      real(fp), pointer :: array_2d_ptr(:,:) => null()
      real(fp), pointer :: array_3d_ptr(:,:,:) => null()
      character(len=64) :: field_name

      rc = cc_success

      call this%error_mgr%push_context('diagnostic_manager_collect_process', &
         'Collecting diagnostics from process: ' // trim(process_name))

      ! Get process registry
      call this%get_process_registry(process_name, registry, rc)
      if (rc /= cc_success) then
         call this%error_mgr%pop_context()
         return
      end if

      ! Get the number of registered diagnostic fields
      num_fields = registry%get_field_count()
      if (num_fields == 0) then
         call this%error_mgr%pop_context()
         return  ! No fields to collect
      end if

      ! Allocate array for field names
      allocate(field_names(num_fields), stat=local_rc)
      if (local_rc /= 0) then
         call this%error_mgr%report_error(error_memory_allocation, &
            'Failed to allocate field names array', rc)
         call this%error_mgr%pop_context()
         return
      end if

      ! Get all field names
      call registry%list_fields(field_names, num_fields)

      ! Iterate through all diagnostic fields and collect their current values
      do i = 1, num_fields
         field_name = trim(field_names(i))

         ! Get field value using existing get_field_value method
         call this%get_field_value(process_name, field_name, &
            scalar_value, array_1d_ptr, array_2d_ptr, array_3d_ptr, &
            data_type, rc=local_rc)

         if (local_rc /= cc_success) then
            ! Log warning but continue with other fields
            call this%error_mgr%report_error(local_rc, &
               'Failed to collect field: ' // trim(field_name), local_rc)
            ! Don't fail the entire collection for one field
         else
            ! TODO: Write diagnostic fields to output file
            ! This is where we would write each field's values to a diagnostic output file
            ! The file format could be NetCDF, CSV, or custom binary format
            ! Example: call write_field_to_file(process_name, field_name, data_type, values)
            ! For now, we just verify the field is accessible and report successful collection
         end if

         ! Clean up pointers
         if (associated(array_1d_ptr)) nullify(array_1d_ptr)
         if (associated(array_2d_ptr)) nullify(array_2d_ptr)
         if (associated(array_3d_ptr)) nullify(array_3d_ptr)
      end do

      ! Clean up
      deallocate(field_names)

      call this%error_mgr%pop_context()

   end subroutine diagnostic_manager_collect_process

   subroutine diagnostic_manager_get_field_value(this, process_name, field_name, &
      scalar_value, array_1d_ptr, array_2d_ptr, array_3d_ptr, &
      data_type, description, units, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: field_name
      real(fp), intent(out), optional :: scalar_value
      real(fp), pointer, intent(out), optional :: array_1d_ptr(:)
      real(fp), pointer, intent(out), optional :: array_2d_ptr(:,:)
      real(fp), pointer, intent(out), optional :: array_3d_ptr(:,:,:)
      integer, intent(out), optional :: data_type
      character(len=*), intent(out), optional :: description
      character(len=*), intent(out), optional :: units
      integer, intent(out) :: rc

      type(DiagnosticRegistryType), pointer :: registry
      type(DiagnosticFieldType), pointer :: field_ptr
      type(DiagnosticDataType), pointer :: data_ptr
      integer :: local_data_type

      rc = cc_success

      ! Initialize outputs
      if (present(scalar_value)) scalar_value = 0.0_fp
      if (present(array_1d_ptr)) nullify(array_1d_ptr)
      if (present(array_2d_ptr)) nullify(array_2d_ptr)
      if (present(array_3d_ptr)) nullify(array_3d_ptr)
      if (present(data_type)) data_type = 0
      if (present(description)) description = ''
      if (present(units)) units = ''

      ! Get process registry
      call this%get_process_registry(process_name, registry, rc)
      if (rc /= cc_success) return

      ! Get diagnostic field
      field_ptr => registry%get_field_ptr(field_name)
      if (.not. associated(field_ptr)) then
         rc = error_not_found
         call this%error_mgr%report_error(error_not_found, &
            'Diagnostic field not found: ' // trim(field_name), rc)
         return
      end if

      ! Check if field is ready and enabled
      if (.not. field_ptr%is_ready() .or. .not. field_ptr%get_is_enabled()) then
         rc = error_invalid_input
         call this%error_mgr%report_error(error_invalid_input, &
            'Diagnostic field not ready or disabled: ' // trim(field_name), rc)
         return
      end if

      ! Get data pointer
      data_ptr => field_ptr%get_data_ptr()
      if (.not. associated(data_ptr)) then
         rc = error_invalid_input
         return
      end if

      ! Get data type and extract values
      local_data_type = data_ptr%get_data_type()
      if (present(data_type)) data_type = local_data_type

      ! Get field metadata
      if (present(description)) then
         description = field_ptr%get_description()
      end if

      if (present(units)) then
         units = field_ptr%get_units()
      end if

      select case (local_data_type)
       case (diag_real_scalar)
         if (present(scalar_value)) then
            scalar_value = data_ptr%get_real_scalar()
         end if

       case (diag_real_1d)
         if (present(array_1d_ptr)) then
            array_1d_ptr => data_ptr%get_real_1d_ptr()
         end if

       case (diag_real_2d)
         if (present(array_2d_ptr)) then
            array_2d_ptr => data_ptr%get_real_2d_ptr()
         end if

       case (diag_real_3d)
         if (present(array_3d_ptr)) then
            array_3d_ptr => data_ptr%get_real_3d_ptr()
         end if

       case default
         rc = error_invalid_input
         call this%error_mgr%report_error(error_invalid_input, &
            'Unsupported diagnostic data type for field: ' // trim(field_name), rc)
      end select

   end subroutine diagnostic_manager_get_field_value

   subroutine diagnostic_manager_write_output(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success

      if (.not. this%output_enabled) return

      ! Check if it's time to output
      if (mod(this%current_timestep, this%output_frequency) /= 0) return

      call this%error_mgr%push_context('diagnostic_manager_write_output', &
         'Writing diagnostic output')

      ! TODO: Implement actual file output
      ! This would write diagnostic data to NetCDF files or other formats
      ! For now, placeholder

      call this%error_mgr%pop_context()

   end subroutine diagnostic_manager_write_output

   subroutine diagnostic_manager_advance_timestep(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = cc_success
      this%current_timestep = this%current_timestep + 1

   end subroutine diagnostic_manager_advance_timestep

   function diagnostic_manager_get_total_diagnostics(this) result(total_diagnostics)
      class(DiagnosticManagerType), intent(in) :: this
      integer :: total_diagnostics

      integer :: i

      total_diagnostics = 0

      do i = 1, this%num_processes
         total_diagnostics = total_diagnostics + this%process_registries(i)%get_num_diagnostics()
      enddo

   end function diagnostic_manager_get_total_diagnostics

   subroutine diagnostic_manager_print_summary(this)
      class(DiagnosticManagerType), intent(in) :: this

      integer :: i

      write(*,'(A)') '=== DiagnosticManager Summary ==='
      write(*,'(A,I0)') 'Number of processes: ', this%num_processes
      write(*,'(A,I0)') 'Total diagnostics: ', this%get_total_diagnostics()
      write(*,'(A,L1)') 'Collection enabled: ', this%collection_enabled
      write(*,'(A,L1)') 'Output enabled: ', this%output_enabled
      write(*,'(A,I0)') 'Output frequency: ', this%output_frequency
      write(*,'(A,I0)') 'Current timestep: ', this%current_timestep

      write(*,'(A)') 'Registered processes:'
      do i = 1, this%num_processes
         write(*,'(A,I0,A,A,A,I0,A)') '  ', i, ': ', trim(this%process_names(i)), &
            ' (', this%process_registries(i)%get_num_diagnostics(), ' diagnostics)'
      enddo

      write(*,'(A)') '================================='

   end subroutine diagnostic_manager_print_summary

   subroutine diagnostic_manager_validate_state(this, rc)
      class(DiagnosticManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = cc_success

      ! Validate each process registry
      do i = 1, this%num_processes
         call this%process_registries(i)%validate(this%error_mgr, local_rc)
         if (local_rc /= cc_success) then
            call this%error_mgr%report_error(local_rc, &
               'Validation failed for process: ' // &
               trim(this%process_names(i)), rc)
            return
         endif
      enddo

   end subroutine diagnostic_manager_validate_state

end module diagnosticmanager_mod
```


