

# File yaml\_interface\_mod.F90

[**File List**](files.md) **>** [**external**](dir_805a0af995e93a362739e98abd740eb2.md) **>** [**yaml\_interface**](dir_d0b1a67acd809cff502adc02c61e9ebd.md) **>** [**yaml\_interface\_mod.F90**](yaml__interface__mod_8_f90.md)

[Go to the documentation of this file](yaml__interface__mod_8_f90.md)


```Fortran

module yaml_interface_mod
   use iso_c_binding
   use iso_fortran_env, only: real32, real64
   use precision_mod, only: fp  ! Use project-wide precision definition
   implicit none
   private

   ! No longer define fp here - use the one from Precision_Mod

   ! Public types
   public :: yaml_node_t

   ! Public procedures - Low level
   public :: yaml_load_file, yaml_load_string, yaml_destroy_node
   public :: yaml_get_string, yaml_get_integer, yaml_get_real, yaml_get_logical
   public :: yaml_get_real_array, yaml_get_integer_array, yaml_get_string_array
   public :: yaml_has_key, yaml_get_size, yaml_is_sequence, yaml_is_map
   public :: yaml_set_string, yaml_set_integer, yaml_set_real, yaml_set_logical
   public :: yaml_save_file, yaml_get_all_keys

   ! Safe conversion functions (no yaml-cpp error messages)
   public :: safe_yaml_get_real, safe_yaml_get_logical, safe_yaml_get_integer

   ! High-level generic interfaces
   public :: yaml_get, yaml_set, yaml_get_array

   type :: yaml_node_t
      type(c_ptr) :: ptr = c_null_ptr
   end type yaml_node_t

   interface yaml_get
      module procedure yaml_get_string_generic
      module procedure yaml_get_integer_generic
      module procedure yaml_get_logical_generic
      module procedure yaml_get_real_sp_generic
      module procedure yaml_get_real_dp_generic
   end interface yaml_get

   interface yaml_set
      module procedure yaml_set_string_generic
      module procedure yaml_set_integer_generic
      module procedure yaml_set_logical_generic
      module procedure yaml_set_real_dp_generic
   end interface yaml_set

   interface yaml_get_array
      module procedure yaml_get_string_array_generic
      module procedure yaml_get_integer_array_generic
      module procedure yaml_get_real_dp_array_generic
   end interface yaml_get_array

   ! C interface declarations
   interface
      ! Node management
      function c_yaml_load_file(filename) bind(C, name="yaml_load_file")
         import :: c_ptr, c_char
         character(kind=c_char), intent(in) :: filename(*)
         type(c_ptr) :: c_yaml_load_file
      end function c_yaml_load_file

      function c_yaml_load_string(yaml_string) bind(C, name="yaml_load_string")
         import :: c_ptr, c_char
         character(kind=c_char), intent(in) :: yaml_string(*)
         type(c_ptr) :: c_yaml_load_string
      end function c_yaml_load_string

      subroutine c_yaml_destroy_node(node) bind(C, name="yaml_destroy_node")
         import :: c_ptr
         type(c_ptr), value :: node
      end subroutine c_yaml_destroy_node

      ! Getters
      function c_yaml_get_string(node, key, value, max_len) bind(C, name="yaml_get_string")
         import :: c_ptr, c_char, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         character(kind=c_char), intent(out) :: value(*)
         integer(c_int), value :: max_len
         logical(c_bool) :: c_yaml_get_string
      end function c_yaml_get_string

      function c_yaml_get_integer(node, key, value) bind(C, name="yaml_get_integer")
         import :: c_ptr, c_char, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         integer(c_int), intent(out) :: value
         logical(c_bool) :: c_yaml_get_integer
      end function c_yaml_get_integer

      function c_yaml_get_real(node, key, value) bind(C, name="yaml_get_real")
         import :: c_ptr, c_char, c_double, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         real(c_double), intent(out) :: value
         logical(c_bool) :: c_yaml_get_real
      end function c_yaml_get_real

      function c_yaml_get_logical(node, key, value) bind(C, name="yaml_get_logical")
         import :: c_ptr, c_char, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         logical(c_bool), intent(out) :: value
         logical(c_bool) :: c_yaml_get_logical
      end function c_yaml_get_logical

      ! Array getters
      function c_yaml_get_real_array(node, key, values, max_size, actual_size) bind(C, name="yaml_get_real_array")
         import :: c_ptr, c_char, c_double, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         real(c_double), intent(out) :: values(*)
         integer(c_int), value :: max_size
         integer(c_int), intent(out) :: actual_size
         logical(c_bool) :: c_yaml_get_real_array
      end function c_yaml_get_real_array

      function c_yaml_get_integer_array(node, key, values, max_size, actual_size) bind(C, name="yaml_get_integer_array")
         import :: c_ptr, c_char, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         integer(c_int), intent(out) :: values(*)
         integer(c_int), value :: max_size
         integer(c_int), intent(out) :: actual_size
         logical(c_bool) :: c_yaml_get_integer_array
      end function c_yaml_get_integer_array

      function c_yaml_get_string_array(node, key, values, max_strings, max_len, actual_size) bind(C, name="yaml_get_string_array")
         import :: c_ptr, c_char, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         character(kind=c_char), intent(out) :: values(*)
         integer(c_int), value :: max_strings
         integer(c_int), value :: max_len
         integer(c_int), intent(out) :: actual_size
         logical(c_bool) :: c_yaml_get_string_array
      end function c_yaml_get_string_array

      ! Utility functions
      function c_yaml_has_key(node, key) bind(C, name="yaml_has_key")
         import :: c_ptr, c_char, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         logical(c_bool) :: c_yaml_has_key
      end function c_yaml_has_key

      function c_yaml_get_size(node) bind(C, name="yaml_get_size")
         import :: c_ptr, c_int
         type(c_ptr), value :: node
         integer(c_int) :: c_yaml_get_size
      end function c_yaml_get_size

      function c_yaml_is_sequence(node) bind(C, name="yaml_is_sequence")
         import :: c_ptr, c_bool
         type(c_ptr), value :: node
         logical(c_bool) :: c_yaml_is_sequence
      end function c_yaml_is_sequence

      function c_yaml_is_map(node) bind(C, name="yaml_is_map")
         import :: c_ptr, c_bool
         type(c_ptr), value :: node
         logical(c_bool) :: c_yaml_is_map
      end function c_yaml_is_map

      ! Setters
      function c_yaml_set_string(node, key, value) bind(C, name="yaml_set_string")
         import :: c_ptr, c_char, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         character(kind=c_char), intent(in) :: value(*)
         logical(c_bool) :: c_yaml_set_string
      end function c_yaml_set_string

      function c_yaml_set_integer(node, key, value) bind(C, name="yaml_set_integer")
         import :: c_ptr, c_char, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         integer(c_int), value :: value
         logical(c_bool) :: c_yaml_set_integer
      end function c_yaml_set_integer

      function c_yaml_set_real(node, key, value) bind(C, name="yaml_set_real")
         import :: c_ptr, c_char, c_double, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         real(c_double), value :: value
         logical(c_bool) :: c_yaml_set_real
      end function c_yaml_set_real

      function c_yaml_set_logical(node, key, value) bind(C, name="yaml_set_logical")
         import :: c_ptr, c_char, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: key(*)
         logical(c_bool), value :: value
         logical(c_bool) :: c_yaml_set_logical
      end function c_yaml_set_logical

      ! File operations
      function c_yaml_save_file(node, filename) bind(C, name="yaml_save_file")
         import :: c_ptr, c_char, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(in) :: filename(*)
         logical(c_bool) :: c_yaml_save_file
      end function c_yaml_save_file

      function c_yaml_get_all_keys(node, keys, max_keys, max_key_len, actual_count) bind(C, name="yaml_get_all_keys")
         import :: c_ptr, c_char, c_int, c_bool
         type(c_ptr), value :: node
         character(kind=c_char), intent(out) :: keys(*)
         integer(c_int), value :: max_keys
         integer(c_int), value :: max_key_len
         integer(c_int), intent(out) :: actual_count
         logical(c_bool) :: c_yaml_get_all_keys
      end function c_yaml_get_all_keys

   end interface

contains

   !========================================================================
   ! Low-level wrapper functions
   !========================================================================

   function yaml_load_file(filename) result(node)
      character(len=*), intent(in) :: filename
      type(yaml_node_t) :: node

      node%ptr = c_yaml_load_file(trim(filename)//c_null_char)
   end function yaml_load_file

   function yaml_load_string(yaml_string) result(node)
      character(len=*), intent(in) :: yaml_string
      type(yaml_node_t) :: node

      node%ptr = c_yaml_load_string(trim(yaml_string)//c_null_char)
   end function yaml_load_string

   subroutine yaml_destroy_node(node)
      type(yaml_node_t), intent(inout) :: node

      if (c_associated(node%ptr)) then
         call c_yaml_destroy_node(node%ptr)
         node%ptr = c_null_ptr
      endif
   end subroutine yaml_destroy_node

   function yaml_get_string(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: value
      logical :: success

      character(len=len(value)) :: c_value
      integer :: i

      success = c_yaml_get_string(node%ptr, trim(key)//c_null_char, c_value, len(value))
      if (success) then
         ! Clean null characters from C string before returning to Fortran
         do i = 1, len(c_value)
            if (ichar(c_value(i:i)) == 0) then
               c_value(i:i) = ' '
            endif
         end do
         value = c_value
      endif
   end function yaml_get_string

   function yaml_get_integer(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      integer, intent(out) :: value
      logical :: success

      integer(c_int) :: c_value
      success = c_yaml_get_integer(node%ptr, trim(key)//c_null_char, c_value)
      if (success) then
         value = c_value
      endif
   end function yaml_get_integer

   function yaml_get_real(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(fp), intent(out) :: value
      logical :: success

      real(c_double) :: c_value
      success = c_yaml_get_real(node%ptr, trim(key)//c_null_char, c_value)
      if (success) then
         value = real(c_value, fp)
      endif
   end function yaml_get_real

   function yaml_get_logical(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      logical, intent(out) :: value
      logical :: success

      logical(c_bool) :: c_value
      success = c_yaml_get_logical(node%ptr, trim(key)//c_null_char, c_value)
      if (success) then
         value = c_value
      endif
   end function yaml_get_logical

   function yaml_get_real_array(node, key, values, actual_size) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(fp), intent(out) :: values(:)
      integer, intent(out) :: actual_size
      logical :: success

      real(c_double) :: c_values(size(values))
      integer(c_int) :: c_actual_size

      success = c_yaml_get_real_array(node%ptr, trim(key)//c_null_char, &
         c_values, size(values), c_actual_size)
      if (success) then
         actual_size = c_actual_size
         values(1:actual_size) = real(c_values(1:actual_size), fp)
      endif
   end function yaml_get_real_array

   function yaml_get_integer_array(node, key, values, actual_size) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      integer, intent(out) :: values(:)
      integer, intent(out) :: actual_size
      logical :: success

      integer(c_int) :: c_values(size(values))
      integer(c_int) :: c_actual_size

      success = c_yaml_get_integer_array(node%ptr, trim(key)//c_null_char, &
         c_values, size(values), c_actual_size)
      if (success) then
         actual_size = c_actual_size
         values(1:actual_size) = c_values(1:actual_size)
      endif
   end function yaml_get_integer_array

   function yaml_get_string_array(node, key, values, actual_size) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: values(:)
      integer, intent(out) :: actual_size
      logical :: success

      character(len=len(values)) :: c_values(size(values))
      integer(c_int) :: c_actual_size
      integer :: i, j

      success = c_yaml_get_string_array(node%ptr, trim(key)//c_null_char, &
         c_values, size(values), len(values), c_actual_size)
      if (success) then
         actual_size = c_actual_size
         do i = 1, actual_size
            ! Clean null characters from each C string before returning to Fortran
            do j = 1, len(c_values(i))
               if (ichar(c_values(i)(j:j)) == 0) then
                  c_values(i)(j:j) = ' '
               endif
            end do
            values(i) = c_values(i)
         end do
      endif
   end function yaml_get_string_array

   function yaml_has_key(node, key) result(exists)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      logical :: exists

      exists = c_yaml_has_key(node%ptr, trim(key)//c_null_char)
   end function yaml_has_key

   function yaml_get_size(node) result(size)
      type(yaml_node_t), intent(in) :: node
      integer :: size

      size = c_yaml_get_size(node%ptr)
   end function yaml_get_size

   function yaml_is_sequence(node) result(is_seq)
      type(yaml_node_t), intent(in) :: node
      logical :: is_seq

      is_seq = c_yaml_is_sequence(node%ptr)
   end function yaml_is_sequence

   function yaml_is_map(node) result(is_map)
      type(yaml_node_t), intent(in) :: node
      logical :: is_map

      is_map = c_yaml_is_map(node%ptr)
   end function yaml_is_map

   function yaml_set_string(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      logical :: success

      success = c_yaml_set_string(node%ptr, trim(key)//c_null_char, trim(value)//c_null_char)
   end function yaml_set_string

   function yaml_set_integer(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      integer, intent(in) :: value
      logical :: success

      success = c_yaml_set_integer(node%ptr, trim(key)//c_null_char, int(value, c_int))
   end function yaml_set_integer

   function yaml_set_real(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(fp), intent(in) :: value
      logical :: success

      success = c_yaml_set_real(node%ptr, trim(key)//c_null_char, real(value, c_double))
   end function yaml_set_real

   function yaml_set_logical(node, key, value) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      logical, intent(in) :: value
      logical :: success

      success = c_yaml_set_logical(node%ptr, trim(key)//c_null_char, logical(value, c_bool))
   end function yaml_set_logical

   function yaml_save_file(node, filename) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: filename
      logical :: success

      success = c_yaml_save_file(node%ptr, trim(filename)//c_null_char)
   end function yaml_save_file

   !========================================================================
   ! High-level generic interface implementations
   !========================================================================

   subroutine yaml_get_string_generic(node, key, value, rc, default_value)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: value
      integer, intent(out), optional :: rc
      character(len=*), intent(in), optional :: default_value

      logical :: success
      integer :: local_rc

      success = yaml_get_string(node, key, value)

      if (success) then
         local_rc = 0
      else
         local_rc = -1
         if (present(default_value)) then
            value = default_value
            local_rc = 0
         endif
      endif

      if (present(rc)) rc = local_rc
   end subroutine yaml_get_string_generic

   subroutine yaml_get_integer_generic(node, key, value, rc, default_value)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      integer, intent(out) :: value
      integer, intent(out), optional :: rc
      integer, intent(in), optional :: default_value

      logical :: success
      integer :: local_rc

      success = yaml_get_integer(node, key, value)

      if (success) then
         local_rc = 0
      else
         local_rc = -1
         if (present(default_value)) then
            value = default_value
            local_rc = 0
         endif
      endif

      if (present(rc)) rc = local_rc
   end subroutine yaml_get_integer_generic

   subroutine yaml_get_real_sp_generic(node, key, value, rc, default_value)
      use iso_fortran_env, only: real32
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(kind=real32), intent(out) :: value
      integer, intent(out), optional :: rc
      real(kind=real32), intent(in), optional :: default_value

      logical :: success
      integer :: local_rc
      real(c_double) :: c_value

      ! Call C interface directly and convert to single precision
      success = c_yaml_get_real(node%ptr, trim(key)//c_null_char, c_value)
      if (success) then
         value = real(c_value, kind=real32)
         local_rc = 0
      else
         local_rc = -1
         if (present(default_value)) then
            value = default_value
            local_rc = 0
         endif
      endif

      if (present(rc)) rc = local_rc
   end subroutine yaml_get_real_sp_generic

   subroutine yaml_get_real_dp_generic(node, key, value, rc, default_value)
      use iso_fortran_env, only: real64
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(kind=real64), intent(out) :: value
      integer, intent(out), optional :: rc
      real(kind=real64), intent(in), optional :: default_value

      logical :: success
      integer :: local_rc
      real(c_double) :: c_value

      ! Call C interface directly and convert to double precision
      success = c_yaml_get_real(node%ptr, trim(key)//c_null_char, c_value)
      if (success) then
         value = real(c_value, kind=real64)
         local_rc = 0
      else
         local_rc = -1
         if (present(default_value)) then
            value = default_value
            local_rc = 0
         endif
      endif

      if (present(rc)) rc = local_rc
   end subroutine yaml_get_real_dp_generic

   subroutine yaml_get_logical_generic(node, key, value, rc, default_value)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      logical, intent(out) :: value
      integer, intent(out), optional :: rc
      logical, intent(in), optional :: default_value

      logical :: success
      integer :: local_rc

      success = yaml_get_logical(node, key, value)

      if (success) then
         local_rc = 0
      else
         local_rc = -1
         if (present(default_value)) then
            value = default_value
            local_rc = 0
         endif
      endif

      if (present(rc)) rc = local_rc
   end subroutine yaml_get_logical_generic

   subroutine yaml_set_string_generic(node, key, value, rc)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      integer, intent(out), optional :: rc

      logical :: success
      integer :: local_rc

      success = yaml_set_string(node, key, value)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
   end subroutine yaml_set_string_generic

   subroutine yaml_set_integer_generic(node, key, value, rc)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      integer, intent(in) :: value
      integer, intent(out), optional :: rc

      logical :: success
      integer :: local_rc

      success = yaml_set_integer(node, key, value)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
   end subroutine yaml_set_integer_generic

   subroutine yaml_set_real_dp_generic(node, key, value, rc)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(kind=fp), intent(in) :: value
      integer, intent(out), optional :: rc

      logical :: success
      integer :: local_rc

      success = yaml_set_real(node, key, value)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
   end subroutine yaml_set_real_dp_generic

   subroutine yaml_set_logical_generic(node, key, value, rc)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      logical, intent(in) :: value
      integer, intent(out), optional :: rc

      logical :: success
      integer :: local_rc

      success = yaml_set_logical(node, key, value)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
   end subroutine yaml_set_logical_generic

   subroutine yaml_get_string_array_generic(node, key, values, rc, actual_size)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: values(:)
      integer, intent(out), optional :: rc
      integer, intent(out), optional :: actual_size

      logical :: success
      integer :: local_rc, local_size

      success = yaml_get_string_array(node, key, values, local_size)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
      if (present(actual_size)) actual_size = local_size
   end subroutine yaml_get_string_array_generic

   subroutine yaml_get_integer_array_generic(node, key, values, rc, actual_size)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      integer, intent(out) :: values(:)
      integer, intent(out), optional :: rc
      integer, intent(out), optional :: actual_size

      logical :: success
      integer :: local_rc, local_size

      success = yaml_get_integer_array(node, key, values, local_size)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
      if (present(actual_size)) actual_size = local_size
   end subroutine yaml_get_integer_array_generic

   subroutine yaml_get_real_dp_array_generic(node, key, values, rc, actual_size)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(in) :: key
      real(kind=fp), intent(out) :: values(:)
      integer, intent(out), optional :: rc
      integer, intent(out), optional :: actual_size

      logical :: success
      integer :: local_rc, local_size

      success = yaml_get_real_array(node, key, values, local_size)
      local_rc = merge(0, -1, success)

      if (present(rc)) rc = local_rc
      if (present(actual_size)) actual_size = local_size
   end subroutine yaml_get_real_dp_array_generic

   function yaml_get_all_keys(node, keys, actual_count) result(success)
      type(yaml_node_t), intent(in) :: node
      character(len=*), intent(out) :: keys(:)
      integer, intent(out) :: actual_count
      logical :: success

      character(kind=c_char) :: c_keys(size(keys) * len(keys(1)))
      integer(c_int) :: c_actual_count
      integer :: i, j, key_len, start_pos

      key_len = len(keys(1))
      success = c_yaml_get_all_keys(node%ptr, c_keys, size(keys), key_len, c_actual_count)

      if (success) then
         actual_count = c_actual_count
         do i = 1, actual_count
            start_pos = (i - 1) * key_len + 1
            keys(i) = transfer(c_keys(start_pos:start_pos + key_len - 1), keys(i))
            ! Clean null characters from each key before returning to Fortran
            do j = 1, len(keys(i))
               if (ichar(keys(i)(j:j)) == 0) then
                  keys(i)(j:j) = ' '
               endif
            end do
         end do
      else
         actual_count = 0
      endif
   end function yaml_get_all_keys

   !========================================================================
   ! Safe conversion functions (no yaml-cpp error messages)
   !========================================================================

   subroutine safe_yaml_get_real(yaml_root, key, value, rc)
      implicit none
      type(yaml_node_t), intent(in) :: yaml_root
      character(len=*), intent(in) :: key
      real(fp), intent(out) :: value
      integer, intent(out) :: rc

      character(len=64) :: str_value
      integer :: iostat_val

      ! Try to read as string first to avoid yaml-cpp conversion errors
      call yaml_get(yaml_root, key, str_value, rc)
      if (rc /= 0) then
         return  ! Key not found
      endif

      ! Convert string to real
      read(str_value, *, iostat=iostat_val) value
      if (iostat_val /= 0) then
         rc = -1  ! Conversion failed
      else
         rc = 0   ! Success
      endif
   end subroutine safe_yaml_get_real

   subroutine safe_yaml_get_integer(yaml_root, key, value, rc)
      implicit none
      type(yaml_node_t), intent(in) :: yaml_root
      character(len=*), intent(in) :: key
      integer, intent(out) :: value
      integer, intent(out) :: rc

      character(len=64) :: str_value
      integer :: iostat_val

      ! Try to read as string first to avoid yaml-cpp conversion errors
      call yaml_get(yaml_root, key, str_value, rc)
      if (rc /= 0) then
         return  ! Key not found
      endif

      ! Convert string to integer
      read(str_value, *, iostat=iostat_val) value
      if (iostat_val /= 0) then
         rc = -1  ! Conversion failed
      else
         rc = 0   ! Success
      endif
   end subroutine safe_yaml_get_integer

   subroutine safe_yaml_get_logical(yaml_root, key, value, rc)
      implicit none
      type(yaml_node_t), intent(in) :: yaml_root
      character(len=*), intent(in) :: key
      logical, intent(out) :: value
      integer, intent(out) :: rc

      character(len=64) :: str_value
      character(len=64) :: lower_str

      ! Try to read as string first to avoid yaml-cpp conversion errors
      call yaml_get(yaml_root, key, str_value, rc)
      if (rc /= 0) then
         return  ! Key not found
      endif

      ! Convert to lowercase for comparison
      lower_str = trim(adjustl(str_value))
      call to_lowercase_internal(lower_str)

      ! Convert string to logical with various accepted formats
      select case (trim(lower_str))
       case ('true', 't', '1', 'yes', 'y', 'on')
         value = .true.
         rc = 0
       case ('false', 'f', '0', 'no', 'n', 'off')
         value = .false.
         rc = 0
       case default
         rc = -1  ! Conversion failed
      end select
   end subroutine safe_yaml_get_logical

   subroutine to_lowercase_internal(str)
      implicit none
      character(len=*), intent(inout) :: str
      integer :: i, ascii_val

      do i = 1, len_trim(str)
         ascii_val = ichar(str(i:i))
         if (ascii_val >= 65 .and. ascii_val <= 90) then  ! A-Z
            str(i:i) = char(ascii_val + 32)  ! Convert to lowercase
         endif
      end do
   end subroutine to_lowercase_internal

end module yaml_interface_mod
```


