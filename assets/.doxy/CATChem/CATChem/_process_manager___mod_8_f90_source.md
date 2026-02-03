

# File ProcessManager\_Mod.F90

[**File List**](files.md) **>** [**core**](dir_aebb8dcc11953d78e620bbef0b9e2183.md) **>** [**ProcessManager\_Mod.F90**](_process_manager___mod_8_f90.md)

[Go to the documentation of this file](_process_manager___mod_8_f90.md)


```Fortran

module processmanager_mod
   use precision_mod
   use statemanager_mod, only : statemanagertype
   use error_mod, only : cc_success, cc_failure
   use processinterface_mod, only : processinterface, columnprocessinterface
   use processfactory_mod, only : processfactorytype
   use gridmanager_mod, only : gridmanagertype, columniteratortype
   use columninterface_mod, only : columnprocessortype, columnviewtype
   use virtualcolumn_mod, only : virtualcolumntype
   use configmanager_mod, only : configdatatype, runphasetype, processconfigtype

   implicit none
   private

   public :: processmanagertype

   ! 1. Define a wrapper otherwise the polymorphic array allocation fails
   type :: processcontainertype
      ! The 'allocatable' keyword here is the magic sauce
      class(ProcessInterface), allocatable :: item
   end type processcontainertype

   type :: processmanagertype
      private
      class(ProcessContainerType), allocatable, public  :: processes(:)
      integer :: num_processes = 0
      integer :: max_processes = 50
      type(ProcessFactoryType) :: factory
      type(ColumnProcessorType) :: column_processor
      character(len=64), public, allocatable :: required_met_fields(:)
   contains
      procedure :: init => manager_init
      procedure :: add_process => manager_add_process
      procedure :: run_all => manager_run_all
      procedure :: run_process => manager_run_process
      procedure :: run_column_processes => manager_run_column_processes
      procedure :: run_process_on_columns => manager_run_process_on_columns
      procedure :: finalize => manager_finalize
      procedure :: list_processes => manager_list_processes
      procedure :: get_column_processes => manager_get_column_processes
      procedure :: run_phase => manager_run_phase
      procedure :: run_all_phases => manager_run_all_phases
      procedure :: run_all_processes => manager_run_all_processes
      procedure :: set_max_processes => manager_set_max_processes
      procedure :: enable_column_batching => manager_enable_column_batching
      procedure :: register_process => manager_register_process
      procedure, private :: add_met_fields_from_process => manager_add_met_fields_from_process
   end type processmanagertype

contains

   subroutine manager_init(this, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%factory%init(rc)
      if (rc /= cc_success) return

      ! Initialize counters - allocatable array will be allocated on first assignment
      this%num_processes = 0

      ! Initialize empty required met fields array
      if (allocated(this%required_met_fields)) deallocate(this%required_met_fields)
      allocate(this%required_met_fields(0))

      ! Initialize column processor with default max columns
      call this%column_processor%init(100, rc)  ! Default max columns
      if (rc /= cc_success) return

      rc = cc_success
   end subroutine manager_init

   subroutine manager_add_process(this, process_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i
      class(ProcessInterface), allocatable :: new_process
      class(ProcessContainerType), allocatable :: tmp(:)

      if (this%num_processes >= this%max_processes) then
         rc = cc_failure
         return
      endif

      ! Create the process (scheme is read from configuration)
      call this%factory%create_process(process_name, container, new_process, rc)
      if (rc /= cc_success) return

      ! Initialize the process
      call new_process%init(container, rc)
      if (rc /= cc_success) return

      ! Set timestep for each process to the same value.
      call new_process%set_timestep(container%tstep)

      ! Collect required met fields from this process
      call this%add_met_fields_from_process(new_process, rc)
      if (rc /= cc_success) return

      ! Add to manager
      this%num_processes = this%num_processes + 1

      ! Check bounds first
      if (this%num_processes > this%max_processes) then
         rc = cc_failure
         return
      endif

      ! ! Handle polymorphic array allocation on first use
      ! if (.not. allocated(this%processes)) then
      !    ! For polymorphic arrays, we need to allocate with proper bounds
      !    ! We'll allocate the whole array when adding the first process
      !    block
      !       class(ProcessInterface), allocatable :: temp_array(:)
      !       allocate(temp_array(this%max_processes), source=new_process)
      !       call move_alloc(temp_array, this%processes)
      !    end block
      ! else
      !    ! Subsequent assignments - just copy into the allocated slot
      !    allocate(this%processes(this%num_processes), source=new_process)
      ! endif

      ! Allocate a temporary array of WRAPPERS with the new size
      allocate(tmp(this%num_processes))

      ! Move existing items into the new array
      ! We check if there was previous data to move
      if (allocated(this%processes)) then
         do i = 1, this%num_processes - 1
            ! We can move the internal allocatable item safely!
            call move_alloc(from=this%processes(i)%item, to=tmp(i)%item)
         end do
         ! Optional: Deallocate the old empty shell (Fortran does this auto, but safe to be explicit)
         deallocate(this%processes)
      endif

      ! Move the new process into the last slot
      call move_alloc(from=new_process, to=tmp(this%num_processes)%item)

      ! Move the wrapper array back to the manager
      call move_alloc(from=tmp, to=this%processes)

      rc = cc_success
   end subroutine manager_add_process

   subroutine manager_run_all(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, local_rc
      type(GridManagerType), pointer :: grid_mgr

      ! If there are no processes, succeed immediately
      if (this%num_processes == 0) then
         rc = cc_success
         return
      endif

      rc = cc_success

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = cc_failure
         return
      endif

      do i = 1, this%num_processes
         if (this%processes(i)%item%is_ready()) then
            ! Check if this is a column process
            select type(proc => this%processes(i)%item)
             class is (columnprocessinterface)
               ! Run column-based process
               call this%run_process_on_columns(i, container, local_rc)
             class default
               ! Run traditional 3D process
               call this%processes(i)%item%run(container, local_rc)
            end select

            if (local_rc /= cc_success) then
               rc = local_rc
               return
            endif
         endif
      enddo
   end subroutine manager_run_all

   subroutine manager_run_column_processes(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, local_rc, col_i, col_j
      type(GridManagerType), pointer :: grid_mgr
      type(ColumnIteratorType) :: col_iter
      type(VirtualColumnType) :: virtual_col

      rc = cc_success

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = cc_failure
         return
      endif

      ! Initialize column iterator using create_column_iterator
      col_iter = grid_mgr%create_column_iterator()

      ! Process each column
      do while (col_iter%has_next())
         call col_iter%next(rc)
         if (rc /= cc_success) return

         ! Get current column indices (i, j)
         call col_iter%get_current_indices(col_i, col_j)

         ! Create and populate virtual column for this (i, j)
         call container%create_virtual_column(col_i, col_j, virtual_col, rc)
         if (rc /= cc_success) return

         ! Run all column processes on this column
         do i = 1, this%num_processes
            select type(proc => this%processes(i)%item)
             class is (columnprocessinterface)
               if (proc%is_ready()) then
                  call proc%run_column(virtual_col, container, local_rc)
                  if (local_rc /= cc_success) then
                     rc = local_rc
                     return
                  endif
               endif
            end select
         enddo

         ! Apply virtual column changes back to container
         call container%apply_virtual_column(virtual_col, rc)
         if (rc /= cc_success) return
      enddo
      ! Clean up virtual column
      if (virtual_col%is_valid)  call virtual_col%cleanup()

   end subroutine manager_run_column_processes

   subroutine manager_run_process_on_columns(this, process_index, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(in) :: process_index
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr
      type(ColumnIteratorType) :: col_iter
      type(VirtualColumnType) :: virtual_col
      integer :: col_i, col_j

      rc = cc_success

      if (process_index < 1 .or. process_index > this%num_processes) then
         rc = cc_failure
         return
      endif

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = cc_failure
         return
      endif

      select type(proc => this%processes(process_index)%item)
       class is (columnprocessinterface)
         ! Initialize column iterator using create_column_iterator
         col_iter = grid_mgr%create_column_iterator()

         ! Process each column
         do while (col_iter%has_next())
            call col_iter%next(rc)
            if (rc /= cc_success) return

            ! Get current column indices
            call col_iter%get_current_indices(col_i, col_j)

            call container%create_virtual_column(col_i, col_j, virtual_col, rc)
            if (rc /= cc_success) return

            call proc%run_column(virtual_col, container, rc)
            if (rc /= cc_success) return

            ! Apply virtual column changes back to container for each column
            call container%apply_virtual_column(virtual_col, rc)
            if (rc /= cc_success) return
         enddo
         ! Clean up virtual column
         if (virtual_col%is_valid)  call virtual_col%cleanup()

      end select
   end subroutine manager_run_process_on_columns

   subroutine manager_run_process(this, process_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i

      rc = cc_failure
      do i = 1, this%num_processes
         if (trim(this%processes(i)%item%get_name()) == trim(process_name)) then
            ! Check if this is a column process and run appropriately
            select type(proc => this%processes(i)%item)
             class is (columnprocessinterface)
               call this%run_process_on_columns(i, container, rc)
             class default
               call this%processes(i)%item%run(container, rc)
            end select
            return
         endif
      enddo
   end subroutine manager_run_process

   subroutine manager_get_column_processes(this, column_indices, count)
      class(ProcessManagerType), intent(in) :: this
      integer, intent(out) :: column_indices(:)
      integer, intent(out) :: count

      integer :: i, max_count

      count = 0
      max_count = min(this%num_processes, size(column_indices))

      do i = 1, this%num_processes
         select type(proc => this%processes(i)%item)
          class is (columnprocessinterface)
            count = count + 1
            if (count <= max_count) then
               column_indices(count) = i
            endif
         end select
      enddo
   end subroutine manager_get_column_processes

   subroutine manager_run_phase(this, phase_name, config_data, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: phase_name
      type(ConfigDataType), intent(in) :: config_data
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, j, local_rc, phase_idx, process_idx
      type(RunPhaseType) :: current_phase
      type(ProcessConfigType) :: process_config
      logical :: phase_found

      rc = cc_success
      phase_found = .false.

      ! Check if run phases are configured
      if (.not. allocated(config_data%run_phases)) then
         write(*,*) 'WARNING: No run phases configured, skipping phase execution'
         return
      endif

      ! Find the requested phase
      do phase_idx = 1, size(config_data%run_phases)
         if (trim(config_data%run_phases(phase_idx)%name) == trim(phase_name)) then
            current_phase = config_data%run_phases(phase_idx)
            phase_found = .true.
            exit
         endif
      enddo

      if (.not. phase_found) then
         write(*,*) 'WARNING: Phase "', trim(phase_name), '" not found in configuration'
         rc = cc_failure
         return
      endif

      write(*,*) 'INFO: Running phase "', trim(phase_name), '" with ', current_phase%num_processes, ' processes'

      ! Loop through each process in this phase
      do j = 1, current_phase%num_processes
         process_config = current_phase%processes(j)

         ! Check if process is enabled
         if (.not. process_config%enabled) then
            write(*,*) 'INFO: Skipping disabled process: ', trim(process_config%name)
            cycle
         endif

         ! Map process_index to actual process in manager's array
         process_idx = process_config%process_index

         ! Validate process index bounds
         if (process_idx < 1 .or. process_idx > this%num_processes) then
            write(*,*) 'WARNING: Process index ', process_idx, ' out of bounds for process: ', &
               trim(process_config%name)
            cycle
         endif

         write(*,*) 'INFO: Running process: ', trim(process_config%name), ' (index=', process_idx, ')'

         ! Run the process based on its type
         if (this%processes(process_idx)%item%is_ready()) then
            select type(proc => this%processes(process_idx)%item)
             class is (columnprocessinterface)
               !!write(*,*) 'Test phase process', process_idx, trim(proc%name) !debug only
               call this%run_process_on_columns(process_idx, container, local_rc)
             class default
               call this%processes(process_idx)%item%run(container, local_rc)
            end select

            if (local_rc /= cc_success) then
               write(*,*) 'ERROR: Process ', trim(process_config%name), ' failed with code: ', local_rc
               rc = local_rc
               return
            endif
         else
            write(*,*) 'WARNING: Process ', trim(process_config%name), ' is not ready'
         endif
      enddo

      write(*,*) 'INFO: Phase "', trim(phase_name), '" completed successfully'

   end subroutine manager_run_phase

   subroutine manager_run_all_phases(this, config_data, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(ConfigDataType), intent(in) :: config_data
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: phase_idx, local_rc
      character(len=64) :: phase_name

      rc = cc_success

      ! Check if run phases are configured
      if (.not. allocated(config_data%run_phases)) then
         write(*,*) 'WARNING: No run phases configured, skipping all phases execution'
         return
      endif

      if (size(config_data%run_phases) == 0) then
         write(*,*) 'WARNING: Empty run phases array, skipping all phases execution'
         return
      endif

      write(*,*) 'INFO: Running all phases - total phases: ', size(config_data%run_phases)

      ! Loop through all phases and run each one
      do phase_idx = 1, size(config_data%run_phases)
         phase_name = config_data%run_phases(phase_idx)%name

         write(*,*) 'INFO: Starting phase ', phase_idx, ' of ', size(config_data%run_phases), &
            ': "', trim(phase_name), '"'

         ! Run this phase
         call this%run_phase(phase_name, config_data, container, local_rc)
         if (local_rc /= cc_success) then
            write(*,*) 'ERROR: Phase "', trim(phase_name), '" failed with code: ', local_rc
            rc = local_rc
            return
         endif

         write(*,*) 'INFO: Phase "', trim(phase_name), '" completed successfully'
      enddo

      write(*,*) 'INFO: All ', size(config_data%run_phases), ' phases completed successfully'

   end subroutine manager_run_all_phases

   subroutine manager_run_all_processes(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      call this%run_all(container, rc)

   end subroutine manager_run_all_processes

   subroutine manager_set_max_processes(this, max_processes, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(in) :: max_processes
      integer, intent(out) :: rc

      ! TODO: Reallocate process array if needed
      this%max_processes = max_processes
      rc = cc_success

   end subroutine manager_set_max_processes

   subroutine manager_enable_column_batching(this, enable, rc)
      class(ProcessManagerType), intent(inout) :: this
      logical, intent(in) :: enable
      integer, intent(out) :: rc

      ! TODO: Configure column processor batching mode
      rc = cc_success

   end subroutine manager_enable_column_batching

   subroutine manager_finalize(this, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = cc_success
      do i = 1, this%num_processes
         call this%processes(i)%item%finalize(local_rc)
         if (local_rc /= cc_success) then
            rc = local_rc
         endif
      enddo

      ! Clean up column processor
      call this%column_processor%cleanup()

      if (allocated(this%processes)) deallocate(this%processes)
      if (allocated(this%required_met_fields)) deallocate(this%required_met_fields)
      this%num_processes = 0
   end subroutine manager_finalize

   subroutine manager_list_processes(this, process_names, count)
      class(ProcessManagerType), intent(in) :: this
      character(len=64), intent(out) :: process_names(:)
      integer, intent(out) :: count

      integer :: i, max_count

      max_count = min(this%num_processes, size(process_names))
      do i = 1, max_count
         process_names(i) = this%processes(i)%item%get_name()
      enddo
      count = max_count
   end subroutine manager_list_processes

   subroutine manager_register_process(this, name, category, description, creator, rc)
      use processregistry_mod, only: processcreatorinterface

      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: name, category, description
      procedure(ProcessCreatorInterface) :: creator
      integer, intent(out) :: rc

      call this%factory%register_process(name, category, description, creator, rc)
   end subroutine manager_register_process

   subroutine manager_add_met_fields_from_process(this, process, rc)
      class(ProcessManagerType), intent(inout) :: this
      class(ProcessInterface), intent(in) :: process
      integer, intent(out) :: rc

      character(len=64), allocatable :: new_fields(:), merged_fields(:)
      character(len=64), allocatable :: current_fields(:)
      integer :: i, j, current_size, new_size, merged_size
      logical :: field_exists

      rc = cc_success

      ! Get required fields from the new process
      new_fields = process%get_required_met_fields()
      new_size = size(new_fields)

      ! If no new fields, nothing to do
      if (new_size == 0) return

      ! Get current fields (make a copy for merging)
      if (allocated(this%required_met_fields)) then
         current_size = size(this%required_met_fields)
         allocate(current_fields(current_size))
         current_fields = this%required_met_fields
      else
         current_size = 0
         allocate(current_fields(0))
      endif

      ! Merge fields, avoiding duplicates and filtering out TSTEP
      ! Worst case: all new fields are unique, so allocate maximum possible size
      allocate(merged_fields(current_size + new_size))
      merged_size = 0

      ! Start with current fields, but filter out TSTEP
      do i = 1, current_size
         if (trim(adjustl(current_fields(i))) /= 'TSTEP' .and. &
            trim(adjustl(current_fields(i))) /= 'LON'   .and. &
            trim(adjustl(current_fields(i))) /= 'LAT') then
            merged_size = merged_size + 1
            merged_fields(merged_size) = current_fields(i)
         endif
      end do

      ! Add new fields if they're not already present and not TSTEP
      do i = 1, new_size
         field_exists = .false.

         ! Skip TSTEP field (case insensitive)
         if (trim(adjustl(new_fields(i))) == 'TSTEP' .or. &
            trim(adjustl(new_fields(i))) == 'LON'   .or. &
            trim(adjustl(new_fields(i))) == 'LAT') then
            cycle
         endif

         ! Check if this field already exists (case insensitive)
         do j = 1, merged_size
            if (trim(adjustl(new_fields(i))) == trim(adjustl(merged_fields(j)))) then
               field_exists = .true.
               exit
            endif
         end do

         ! Add if it's a new field
         if (.not. field_exists) then
            merged_size = merged_size + 1
            merged_fields(merged_size) = new_fields(i)
         endif
      end do

      ! Update the manager's required fields list with proper size
      if (allocated(this%required_met_fields)) deallocate(this%required_met_fields)
      allocate(this%required_met_fields(merged_size))
      this%required_met_fields(1:merged_size) = merged_fields(1:merged_size)

      ! Cleanup
      if (allocated(current_fields)) deallocate(current_fields)
      if (allocated(new_fields)) deallocate(new_fields)
      if (allocated(merged_fields)) deallocate(merged_fields)

   end subroutine manager_add_met_fields_from_process

end module processmanager_mod
```


