!> \file test_yaml_parsing.F90
!! \brief Test program for YAML field configuration parsing
!!
!! This program tests the new parse_field_section functionality
!! to verify that YAML field arrays are read correctly.
!!
!! \author Assistant
!! \date October 2025

module field_mapping_module
   implicit none

   ! Field mapping type definition
   type :: field_mapping_type
      character(len=128) :: standard_name !< NUOPC/CF standard field name
      character(len=128) :: catchem_var   !< Corresponding CATChem variable path
      integer :: dimensions               !< Number of spatial dimensions (2D/3D)
      character(len=64) :: units          !< Physical units for conversion
      logical :: optional = .false.       !< Whether field is required or optional
   end type field_mapping_type

   ! Field configuration type definition
   type :: field_config_type
      integer :: n_import_fields = 0                         !< Number of import fields
      integer :: n_export_fields = 0                         !< Number of export fields
      type(field_mapping_type), allocatable :: import_fields(:) !< Import field mapping array
      type(field_mapping_type), allocatable :: export_fields(:) !< Export field mapping array
   end type field_config_type

   ! Test constants
   integer, parameter :: CC_SUCCESS = 0
   integer, parameter :: CC_FAILURE = 1

contains

   !> Parse a field section (import_fields or export_fields) from YAML file
   !!
   !! \param filename YAML configuration file
   !! \param section_name Section name ('import_fields' or 'export_fields')
   !! \param fields Array to store parsed fields
   !! \param n_fields Number of fields found
   !! \param errflg Error flag
   !! \param errmsg Error message
   !!
   subroutine parse_field_section(filename, section_name, fields, n_fields, errflg, errmsg)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: section_name
      type(field_mapping_type), allocatable, intent(out) :: fields(:)
      integer, intent(out) :: n_fields
      integer, intent(out) :: errflg
      character(len=*), intent(out) :: errmsg

      integer :: unit_num, io_stat, indent_level, section_indent
      character(len=256) :: line, trimmed_line, field_name, field_value
      logical :: in_section, found_section
      integer :: line_number, field_idx, colon_pos
      type(field_mapping_type), allocatable :: temp_fields(:)
      type(field_mapping_type) :: current_field
      logical :: in_field_item

      errflg = CC_SUCCESS
      errmsg = ''
      n_fields = 0
      in_section = .false.
      found_section = .false.
      section_indent = -1
      line_number = 0
      field_idx = 0
      in_field_item = .false.

      ! Initialize current field
      current_field%standard_name = ''
      current_field%catchem_var = ''
      current_field%dimensions = 0
      current_field%units = ''
      current_field%optional = .false.

      ! Open file for reading
      open(newunit=unit_num, file=trim(filename), status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         write(errmsg, '(A,A)') 'Cannot open configuration file: ', trim(filename)
         errflg = CC_FAILURE
         return
      endif

      ! Allocate temporary storage for up to 50 fields
      allocate(temp_fields(50))

      ! Read file line by line
      do
         read(unit_num, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit  ! End of file or error

         line_number = line_number + 1
         trimmed_line = trim(adjustl(line))

         ! Skip empty lines and comments
         if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '#') cycle

         ! Calculate indentation level
         do indent_level = 1, len_trim(line)
            if (line(indent_level:indent_level) /= ' ') exit
         end do
         indent_level = indent_level - 1

         ! Look for section header
         if (index(trimmed_line, ':') > 0) then
            colon_pos = index(trimmed_line, ':')
            field_name = trim(adjustl(trimmed_line(1:colon_pos-1)))
            field_value = trim(adjustl(trimmed_line(colon_pos+1:)))

            ! Check if this is the section we're looking for
            if (trim(field_name) == trim(section_name)) then
               found_section = .true.
               in_section = .true.
               section_indent = indent_level
               cycle
            endif

            ! If we're in the section and find a same-level or higher section, we're done
            if (in_section .and. indent_level <= section_indent) then
               exit
            endif

            ! Process field items within the section
            if (in_section) then
               if (trimmed_line(1:2) == '- ') then
                  ! Save previous field if we were processing one
                  if (in_field_item .and. current_field%standard_name /= '') then
                     n_fields = n_fields + 1
                     if (n_fields <= size(temp_fields)) then
                        temp_fields(n_fields) = current_field
                     endif
                  endif

                  ! Start new field item
                  in_field_item = .true.
                  current_field%standard_name = ''
                  current_field%catchem_var = ''
                  current_field%dimensions = 0
                  current_field%units = ''
                  current_field%optional = .false.

                  ! Check if there's a property on the same line as the dash
                  if (colon_pos > 2) then
                     field_name = trim(adjustl(trimmed_line(3:colon_pos-1)))
                     field_value = trim(adjustl(trimmed_line(colon_pos+1:)))
                     call parse_field_property(field_name, field_value, current_field)
                  endif

               else if (in_field_item .and. indent_level > section_indent) then
                  ! This is a property of the current field item
                  call parse_field_property(field_name, field_value, current_field)
               endif
            endif
         endif
      end do

      ! Save the last field if we're still processing one
      if (in_field_item .and. current_field%standard_name /= '') then
         n_fields = n_fields + 1
         if (n_fields <= size(temp_fields)) then
            temp_fields(n_fields) = current_field
         endif
      endif

      close(unit_num)

      ! Check if we found the section
      if (.not. found_section) then
         write(errmsg, '(A,A,A)') 'Section "', trim(section_name), '" not found in configuration file'
         errflg = CC_FAILURE
         deallocate(temp_fields)
         return
      endif

      ! Allocate final array and copy data
      if (n_fields > 0) then
         allocate(fields(n_fields))
         fields(1:n_fields) = temp_fields(1:n_fields)
      endif

      deallocate(temp_fields)

   end subroutine parse_field_section

   !> Parse a field property and set it in the field structure
   !!
   !! \param property_name Name of the property
   !! \param property_value Value of the property
   !! \param field Field structure to update
   !!
   subroutine parse_field_property(property_name, property_value, field)
      character(len=*), intent(in) :: property_name
      character(len=*), intent(in) :: property_value
      type(field_mapping_type), intent(inout) :: field

      character(len=256) :: clean_value
      integer :: read_stat

      ! Remove quotes from string values
      clean_value = property_value
      if (len_trim(clean_value) >= 2) then
         if ((clean_value(1:1) == '"' .and. clean_value(len_trim(clean_value):len_trim(clean_value)) == '"') .or. &
            (clean_value(1:1) == "'" .and. clean_value(len_trim(clean_value):len_trim(clean_value)) == "'")) then
            clean_value = clean_value(2:len_trim(clean_value)-1)
         endif
      endif

      select case (trim(property_name))
       case ('standard_name')
         field%standard_name = trim(clean_value)
       case ('catchem_var')
         field%catchem_var = trim(clean_value)
       case ('dimensions')
         read(clean_value, *, iostat=read_stat) field%dimensions
         if (read_stat /= 0) field%dimensions = 0
       case ('units')
         field%units = trim(clean_value)
       case ('optional')
         select case (trim(clean_value))
          case ('true', 'True', 'TRUE', '.true.')
            field%optional = .true.
          case ('false', 'False', 'FALSE', '.false.')
            field%optional = .false.
          case default
            field%optional = .false.
         end select
      end select

   end subroutine parse_field_property

end module field_mapping_module

program test_yaml_parsing
   use field_mapping_module
   implicit none

   ! Test variables
   type(field_config_type) :: field_config
   integer :: errflg
   character(len=512) :: errmsg
   integer :: i

   write(*,*) '=================================================='
   write(*,*) '=== YAML FIELD PARSING TEST ===================='
   write(*,*) '=================================================='
   write(*,*) ''

   write(*,*) 'Testing YAML field parsing from file...'
   write(*,*) 'File: CATChem_field_mapping.yml'
   write(*,*) ''

   ! Test import fields parsing
   write(*,*) 'Parsing import_fields section...'
   call parse_field_section('CATChem_field_mapping.yml', 'import_fields', &
      field_config%import_fields, field_config%n_import_fields, errflg, errmsg)

   if (errflg == CC_SUCCESS) then
      write(*,'(A,I0,A)') '✓ Found ', field_config%n_import_fields, ' import fields:'
      do i = 1, field_config%n_import_fields
         write(*,'(A,I0,A)') '  Field ', i, ':'
         write(*,'(A,A)') '    Standard Name: ', trim(field_config%import_fields(i)%standard_name)
         write(*,'(A,A)') '    CATChem Var:   ', trim(field_config%import_fields(i)%catchem_var)
         write(*,'(A,I0)') '    Dimensions:    ', field_config%import_fields(i)%dimensions
         write(*,'(A,A)') '    Units:         ', trim(field_config%import_fields(i)%units)
         write(*,'(A,L1)') '    Optional:      ', field_config%import_fields(i)%optional
         write(*,*) ''
      end do
   else
      write(*,*) '✗ Import fields parsing failed!'
      write(*,*) 'Error: ', trim(errmsg)
      stop 1
   endif

   ! Test export fields parsing
   write(*,*) 'Parsing export_fields section...'
   call parse_field_section('CATChem_field_mapping.yml', 'export_fields', &
      field_config%export_fields, field_config%n_export_fields, errflg, errmsg)

   if (errflg == CC_SUCCESS) then
      write(*,'(A,I0,A)') '✓ Found ', field_config%n_export_fields, ' export fields:'
      do i = 1, field_config%n_export_fields
         write(*,'(A,I0,A)') '  Field ', i, ':'
         write(*,'(A,A)') '    Standard Name: ', trim(field_config%export_fields(i)%standard_name)
         write(*,'(A,A)') '    CATChem Var:   ', trim(field_config%export_fields(i)%catchem_var)
         write(*,'(A,I0)') '    Dimensions:    ', field_config%export_fields(i)%dimensions
         write(*,'(A,A)') '    Units:         ', trim(field_config%export_fields(i)%units)
         write(*,'(A,L1)') '    Optional:      ', field_config%export_fields(i)%optional
         write(*,*) ''
      end do
   else
      write(*,*) '✗ Export fields parsing failed!'
      write(*,*) 'Error: ', trim(errmsg)
      stop 1
   endif

   write(*,*) '=================================================='
   write(*,*) '=== YAML PARSING TEST PASSED ==================='
   write(*,*) 'Successfully read field configuration from YAML!'
   write(*,*) '=================================================='

   ! Cleanup
   if (allocated(field_config%import_fields)) deallocate(field_config%import_fields)
   if (allocated(field_config%export_fields)) deallocate(field_config%export_fields)

end program test_yaml_parsing
