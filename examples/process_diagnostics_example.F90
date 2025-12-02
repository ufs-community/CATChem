!> \file process_diagnostics_example.F90
!! \brief Example of accessing process-level diagnostics with CATChem high-level API
!! \author CATChem Development Team
!! \date 2025

program process_diagnostics_example
   use CATChemAPI_Mod
   use precision_mod
   implicit none

   type(CATChemInstanceType) :: catchem
   type(CATChemConfigType) :: config
   type(CATChemDataType) :: input_data, output_data
   type(CATChemDiagnosticType) :: diagnostic
   type(CATChemDiagnosticType), allocatable :: all_diagnostics(:)

   character(len=64), allocatable :: diagnostic_names(:)
   character(len=32), allocatable :: process_names(:)
   character(len=256), allocatable :: diagnostic_info(:)

   integer :: rc, i, j
   character(len=256) :: error_msg

   ! ============================================================================
   ! Example 1: Discover available process diagnostics
   ! ============================================================================

   write(*,'(A)') '=== CATChem Process Diagnostics Example ==='
   write(*,*)

   ! Setup basic configuration
   config%nx = 10
   config%ny = 10
   config%nz = 20
   config%nspecies = 50
   config%enable_dust = .true.
   config%enable_seasalt = .true.
   config%enable_drydep = .true.
   config%enable_external_emis = .true.
   config%enable_diagnostics = .true.
   config%use_column_processing = .true.

   ! Initialize CATChem
   call catchem%init(config, rc)
   if (rc /= CATCHEM_SUCCESS) then
      call catchem%get_error_message(error_msg)
      write(*,'(A,A)') 'ERROR: Failed to initialize CATChem: ', trim(error_msg)
      stop 1
   endif

   ! Add processes
   call catchem%add_process('Dust', rc=rc)
   call catchem%add_process('SeaSalt', rc=rc)
   call catchem%add_process('DryDeposition', rc=rc)
   call catchem%add_process('ExternalEmission', rc=rc)

   write(*,'(A)') '1. Discovering available process diagnostics...'

   ! List all process diagnostics
   call catchem%list_process_diagnostics(process_names, diagnostic_info, rc)
   if (rc == CATCHEM_SUCCESS) then
      write(*,'(A,I0,A)') '   Found ', size(process_names), ' processes with diagnostics:'
      do i = 1, size(diagnostic_info)
         write(*,'(A,A)') '   - ', trim(diagnostic_info(i))
      end do
   else
      write(*,'(A)') '   No process diagnostics available'
   endif
   write(*,*)

   ! ============================================================================
   ! Example 2: Access specific process diagnostics
   ! ============================================================================

   write(*,'(A)') '2. Accessing specific process diagnostics...'

   ! Get dust process diagnostics
   if (any(process_names == 'Dust')) then
      write(*,'(A)') '   Dust Process Diagnostics:'

      call catchem%get_available_diagnostics('Dust', diagnostic_names, rc)
      if (rc == CATCHEM_SUCCESS) then
         do i = 1, size(diagnostic_names)
            call catchem%get_process_diagnostic('Dust', diagnostic_names(i), diagnostic, rc)
            if (rc == CATCHEM_SUCCESS .and. diagnostic%is_available) then
               write(*,'(A,A,A,A,A,A,A,I0,A)') &
                  '     - ', trim(diagnostic%field_name), ' (', trim(diagnostic%units), &
                  '): ', trim(diagnostic%description), ' [', diagnostic%dimensions, 'D]'

               ! Show some sample data
               if (diagnostic%dimensions == 2 .and. allocated(diagnostic%data_2d)) then
                  write(*,'(A,E12.4)') '       Sample value: ', diagnostic%data_2d(1,1)
               elseif (diagnostic%dimensions == 3 .and. allocated(diagnostic%data_3d)) then
                  write(*,'(A,E12.4)') '       Sample value: ', diagnostic%data_3d(1,1,1)
               elseif (diagnostic%dimensions == 0) then
                  write(*,'(A,E12.4)') '       Scalar value: ', diagnostic%scalar_data
               endif
            endif
         end do
      endif
   endif
   write(*,*)

   ! ============================================================================
   ! Example 3: Get all diagnostics for a process at once
   ! ============================================================================

   write(*,'(A)') '3. Getting all diagnostics for External Emission process...'

   call catchem%get_all_process_diagnostics('ExternalEmission', all_diagnostics, rc)
   if (rc == CATCHEM_SUCCESS) then
      write(*,'(A,I0,A)') '   Found ', size(all_diagnostics), ' diagnostic fields:'

      do i = 1, size(all_diagnostics)
         if (all_diagnostics(i)%is_available) then
            write(*,'(A,A,A,A,A)') &
               '     - ', trim(all_diagnostics(i)%field_name), &
               ' (', trim(all_diagnostics(i)%units), ')'
            write(*,'(A,A)') '       Description: ', trim(all_diagnostics(i)%description)
         else
            write(*,'(A,A,A)') '     - ', trim(all_diagnostics(i)%field_name), ' (not available)'
         endif
      end do
   else
      write(*,'(A)') '   No diagnostics available for External Emission process'
   endif
   write(*,*)

   ! ============================================================================
   ! Example 4: Access diagnostics during simulation
   ! ============================================================================

   write(*,'(A)') '4. Accessing diagnostics during simulation...'

   ! Setup input data (simplified for example)
   allocate(input_data%temperature(config%nx, config%ny, config%nz))
   allocate(input_data%pressure(config%nx, config%ny, config%nz))
   allocate(input_data%concentrations(config%nx, config%ny, config%nz, config%nspecies))

   ! Fill with example data
   input_data%temperature = 288.0_fp
   input_data%pressure = 101325.0_fp
   input_data%concentrations = 1.0e-9_fp

   ! Run a time step
   call catchem%run_timestep(input_data, output_data, rc)
   if (rc == CATCHEM_SUCCESS) then
      write(*,'(A)') '   Simulation step completed successfully'

      ! Now access updated diagnostics
      call catchem%get_process_diagnostic('Dust', 'surface_emission_flux', diagnostic, rc)
      if (rc == CATCHEM_SUCCESS .and. diagnostic%is_available) then
         if (allocated(diagnostic%data_2d)) then
            write(*,'(A,E12.4,A,A)') '   Dust surface emission flux: ', &
               maxval(diagnostic%data_2d), ' ', trim(diagnostic%units)
         endif
      endif

      call catchem%get_process_diagnostic('SeaSalt', 'whitecap_fraction', diagnostic, rc)
      if (rc == CATCHEM_SUCCESS .and. diagnostic%is_available) then
         if (allocated(diagnostic%data_2d)) then
            write(*,'(A,F8.4)') '   Sea salt whitecap fraction: ', &
               maxval(diagnostic%data_2d)
         endif
      endif

   else
      call catchem%get_error_message(error_msg)
      write(*,'(A,A)') '   Simulation failed: ', trim(error_msg)
   endif
   write(*,*)

   ! ============================================================================
   ! Example 5: Monitoring specific diagnostic over time
   ! ============================================================================

   write(*,'(A)') '5. Monitoring diagnostics over time...'

   do i = 1, 5
      ! Run time step
      call catchem%run_timestep(input_data, output_data, rc)
      if (rc /= CATCHEM_SUCCESS) exit

      ! Get dust emission diagnostic
      call catchem%get_process_diagnostic('Dust', 'total_column_emission', diagnostic, rc)
      if (rc == CATCHEM_SUCCESS .and. diagnostic%is_available) then
         if (allocated(diagnostic%data_2d)) then
            write(*,'(A,I0,A,E12.4,A,A)') '   Step ', i, ': Total dust emission = ', &
               sum(diagnostic%data_2d), ' ', trim(diagnostic%units)
         endif
      endif
   end do
   write(*,*)

   ! ============================================================================
   ! Example 6: Process-specific diagnostic analysis
   ! ============================================================================

   write(*,'(A)') '6. Process-specific diagnostic analysis...'

   ! Analyze external emission process diagnostics
   call catchem%get_available_diagnostics('ExternalEmission', diagnostic_names, rc)
   if (rc == CATCHEM_SUCCESS) then
      write(*,'(A)') '   External Emission Process Analysis:'

      do i = 1, size(diagnostic_names)
         call catchem%get_process_diagnostic('ExternalEmission', diagnostic_names(i), diagnostic, rc)
         if (rc == CATCHEM_SUCCESS .and. diagnostic%is_available) then

            select case (trim(diagnostic%field_name))
             case ('species_tendency')
               if (allocated(diagnostic%data_3d)) then
                  write(*,'(A,E12.4,A)') '     Species tendency range: ', &
                     maxval(diagnostic%data_3d) - minval(diagnostic%data_3d), ' mol/mol/s'
               endif

             case ('emission_rate_applied')
               if (allocated(diagnostic%data_2d)) then
                  write(*,'(A,E12.4,A)') '     Applied emission rate: ', &
                     sum(diagnostic%data_2d), ' mol/m²/s'
               endif

             case ('species_validation_status')
               write(*,'(A,F4.1)') '     Species validation status: ', diagnostic%scalar_data

            end select
         endif
      end do
   endif
   write(*,*)

   ! Clean up
   call catchem%finalize(rc)

   write(*,'(A)') '=== Process Diagnostics Example Complete ==='

end program process_diagnostics_example
