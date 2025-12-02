!> \file CATChemNetCDF_Mod.F90
!! \brief High-level NetCDF I/O interface for CATChem with MPI support
!! \date 2025
!!
!! This module provides a simplified, high-level interface for reading and writing
!! NetCDF files in CATChem applications. It abstracts the low-level NetCDF operations
!! and provides convenient routines for common atmospheric data I/O patterns.
!!
!! \section mpi_compatibility MPI Compatibility
!!
!! This module provides both serial and MPI-parallel NetCDF I/O capabilities:
!!
!! **Serial Mode (Default):**
!! - Uses standard NetCDF Fortran interface (netcdf module)
!! - Each MPI rank operates independently on files
!! - Suitable for rank-local I/O operations
!! - Always available when NetCDF is enabled
!!
!! **Parallel Mode (Optional):**
!! - Uses parallel NetCDF-4 with MPI-IO (netcdf_par module)
!! - Supports collective I/O operations across MPI ranks
!! - Requires NetCDF-4 built with parallel HDF5 and MPI support
!! - Enable with NETCDF_PARALLEL and MPI_ENABLED preprocessor flags
!! - Access via open_mpi() method instead of open()
!!
!! \section build_requirements Build Requirements for MPI Support
!!
!! To enable MPI-parallel NetCDF support, ensure:
!! 1. NetCDF-4 library built with --enable-parallel and --enable-netcdf-4
!! 2. HDF5 library built with --enable-parallel
!! 3. MPI library (OpenMPI, MPICH, etc.)
!! 4. CMake flags: -DNETCDF_PARALLEL=ON -DMPI_ENABLED=ON
!!
!! \section usage_examples Usage Examples
!!
!! **Serial I/O (current default):**
!! \code{.f90}
!!   use CATChemNetCDF_Mod
!!   type(NetCDFFileType) :: nc_file
!!   real(fp), allocatable :: temperature(:,:,:)
!!
!!   call nc_file%open('input_file.nc', 'r')
!!   call nc_file%read_var('T', temperature)
!!   call nc_file%close()
!! \endcode
!!
!! **Parallel I/O with MPI:**
!! \code{.f90}
!!   use CATChemNetCDF_Mod
!!   use mpi
!!   type(NetCDFFileType) :: nc_file
!!   real(fp), allocatable :: temperature(:,:,:)
!!   integer :: mpi_comm, rc
!!
!!   if (mpi_netcdf_available()) then
!!      call nc_file%open_mpi('input_file.nc', 'r', MPI_COMM_WORLD, rc)
!!      call nc_file%read_var('T', temperature)  ! Collective I/O
!!      call nc_file%close()
!!   else
!!      ! Fallback to serial mode
!!      call nc_file%open('input_file.nc', 'r')
!!      call nc_file%read_var('T', temperature)
!!      call nc_file%close()
!!   endif
!! \endcode
!!
!! Features:
!! - Simple file opening/closing with automatic error handling
!! - Read/write operations for 1D, 2D, 3D, and 4D arrays
!! - Automatic dimension and variable discovery
!! - Metadata reading and writing
!! - Time series support for atmospheric data
!! - FV3 coordinate system support
!! - MPI-parallel collective I/O support (when available)
!! - Automatic fallback to serial mode when parallel NetCDF unavailable
!!
module CATChemNetCDF_Mod
#ifdef NETCDF_ENABLED
   use netcdf
#ifdef NETCDF_PARALLEL
   use netcdf_par
#endif
#endif
#ifdef MPI_ENABLED
   use mpi
#endif
   use precision_mod, only: fp
   use error_mod, only: error_type
   implicit none
   private

   public :: NetCDFFileType
   public :: NETCDF_READ, NETCDF_WRITE, NETCDF_APPEND
   public :: netcdf_available, mpi_netcdf_available

   ! File access modes
   integer, parameter :: NETCDF_READ = 1
   integer, parameter :: NETCDF_WRITE = 2
   integer, parameter :: NETCDF_APPEND = 3

   ! Maximum dimensions and variables
   integer, parameter :: MAX_DIMS = 10
   integer, parameter :: MAX_VARS = 1000
   integer, parameter :: MAX_NAME_LEN = 64
   integer, parameter :: MAX_STRING_LEN = 256

   !> NetCDF file interface type
   type :: NetCDFFileType
      private
      character(len=MAX_STRING_LEN) :: filename = ''
      integer :: ncid = -1
      integer :: mode = -1
      logical :: is_open = .false.

      ! MPI settings
      logical :: use_mpi = .false.
      integer :: mpi_comm = -1
      integer :: mpi_rank = -1
      integer :: mpi_size = -1

      ! File metadata
      integer :: ndims = 0
      integer :: nvars = 0
      integer :: ngatts = 0
      integer :: unlimdimid = -1

      ! Dimension information
      character(len=MAX_NAME_LEN) :: dim_names(MAX_DIMS)
      integer :: dim_sizes(MAX_DIMS)
      integer :: dim_ids(MAX_DIMS)

   contains
      ! File operations
      procedure :: open => netcdf_open
      procedure :: open_mpi => netcdf_open_mpi
      procedure :: close => netcdf_close
      procedure :: is_file_open => netcdf_is_open

      ! MPI utilities
      procedure :: set_mpi_comm => netcdf_set_mpi_comm
      procedure :: get_mpi_info => netcdf_get_mpi_info
      procedure :: is_mpi_enabled => netcdf_is_mpi_enabled

      ! Information queries
      procedure :: get_filename => netcdf_get_filename
      procedure :: get_dimensions => netcdf_get_dimensions
      procedure :: get_variables => netcdf_get_variables
      procedure :: has_variable => netcdf_has_variable
      procedure :: has_dimension => netcdf_has_dimension

      ! Variable reading (different data types and dimensions)
      procedure :: read_real_1d => netcdf_read_real_1d
      procedure :: read_real_2d => netcdf_read_real_2d
      procedure :: read_real_3d => netcdf_read_real_3d
      procedure :: read_real_4d => netcdf_read_real_4d
      procedure :: read_int_1d => netcdf_read_int_1d
      procedure :: read_string => netcdf_read_string
      generic :: read_var => read_real_1d, read_real_2d, read_real_3d, read_real_4d, &
         read_int_1d, read_string

      ! Variable writing
      procedure :: write_real_1d => netcdf_write_real_1d
      procedure :: write_real_2d => netcdf_write_real_2d
      procedure :: write_real_3d => netcdf_write_real_3d
      procedure :: write_real_4d => netcdf_write_real_4d
      procedure :: write_int_1d => netcdf_write_int_1d
      procedure :: write_string => netcdf_write_string
      generic :: write_var => write_real_1d, write_real_2d, write_real_3d, write_real_4d, &
         write_int_1d, write_string

      ! Attribute operations
      procedure :: read_attribute_real => netcdf_read_attr_real
      procedure :: read_attribute_string => netcdf_read_attr_string
      procedure :: write_attribute_real => netcdf_write_attr_real
      procedure :: write_attribute_string => netcdf_write_attr_string
      generic :: read_attr => read_attribute_real, read_attribute_string
      generic :: write_attr => write_attribute_real, write_attribute_string

      ! Specialized atmospheric data routines
      procedure :: read_fv3_coordinates => netcdf_read_fv3_coords
      procedure :: read_time_series => netcdf_read_time_series
      procedure :: read_meteorology => netcdf_read_meteorology
      procedure :: write_column_data => netcdf_write_column_data

      ! Utility functions
      procedure :: sync => netcdf_sync
      procedure :: get_error_message => netcdf_get_error_msg

   end type NetCDFFileType

contains

   !> Check if NetCDF support is available
   logical function netcdf_available()
#ifdef NETCDF_ENABLED
      netcdf_available = .true.
#else
      netcdf_available = .false.
#endif
   end function netcdf_available

   !> Check if parallel NetCDF with MPI support is available
   logical function mpi_netcdf_available()
#if defined(NETCDF_ENABLED) && defined(NETCDF_PARALLEL) && defined(MPI_ENABLED)
      mpi_netcdf_available = .true.
#else
      mpi_netcdf_available = .false.
#endif
   end function mpi_netcdf_available

   !> Set MPI communicator for parallel I/O
   subroutine netcdf_set_mpi_comm(this, mpi_comm, rc)
      class(NetCDFFileType), intent(inout) :: this
      integer, intent(in) :: mpi_comm
      integer, intent(out), optional :: rc

      integer :: local_rc, rank, size

      local_rc = 0

#if defined(MPI_ENABLED)
      this%mpi_comm = mpi_comm
      this%use_mpi = .true.

      ! Get MPI rank and size
      call MPI_Comm_rank(mpi_comm, rank, local_rc)
      if (local_rc == MPI_SUCCESS) then
         this%mpi_rank = rank
         call MPI_Comm_size(mpi_comm, size, local_rc)
         if (local_rc == MPI_SUCCESS) then
            this%mpi_size = size
         endif
      endif
#else
      local_rc = -1  ! MPI not available
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_set_mpi_comm

   !> Get MPI information
   subroutine netcdf_get_mpi_info(this, mpi_comm, mpi_rank, mpi_size, rc)
      class(NetCDFFileType), intent(in) :: this
      integer, intent(out), optional :: mpi_comm, mpi_rank, mpi_size
      integer, intent(out), optional :: rc

      if (present(mpi_comm)) mpi_comm = this%mpi_comm
      if (present(mpi_rank)) mpi_rank = this%mpi_rank
      if (present(mpi_size)) mpi_size = this%mpi_size
      if (present(rc)) rc = 0

   end subroutine netcdf_get_mpi_info

   !> Check if MPI is enabled for this file
   logical function netcdf_is_mpi_enabled(this)
      class(NetCDFFileType), intent(in) :: this
      netcdf_is_mpi_enabled = this%use_mpi
   end function netcdf_is_mpi_enabled

   !> Open NetCDF file with MPI support
   subroutine netcdf_open_mpi(this, filename, mode_str, mpi_comm, rc)
      class(NetCDFFileType), intent(inout) :: this
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: mode_str  ! 'r', 'w', 'a'
      integer, intent(in) :: mpi_comm
      integer, intent(out), optional :: rc

      integer :: local_rc, ncmode

      local_rc = 0

#if defined(NETCDF_ENABLED) && defined(NETCDF_PARALLEL) && defined(MPI_ENABLED)
      ! Set MPI communicator first
      call this%set_mpi_comm(mpi_comm, local_rc)
      if (local_rc /= 0) then
         if (present(rc)) rc = local_rc
         return
      endif

      ! Close file if already open
      if (this%is_open) then
         call this%close()
      endif

      this%filename = trim(filename)

      ! Determine NetCDF mode with parallel access
      select case (trim(mode_str))
       case ('r', 'read')
         this%mode = NETCDF_READ
         ncmode = IOR(NF90_NOWRITE, NF90_MPIIO)
       case ('w', 'write')
         this%mode = NETCDF_WRITE
         ncmode = IOR(IOR(NF90_CLOBBER, NF90_NETCDF4), NF90_MPIIO)
       case ('a', 'append')
         this%mode = NETCDF_APPEND
         ncmode = IOR(NF90_WRITE, NF90_MPIIO)
       case default
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      end select

      ! Open or create file with parallel access
      if (this%mode == NETCDF_WRITE) then
         local_rc = nf90_create_par(this%filename, ncmode, this%mpi_comm, MPI_INFO_NULL, this%ncid)
      else
         local_rc = nf90_open_par(this%filename, ncmode, this%mpi_comm, MPI_INFO_NULL, this%ncid)
      endif

      if (local_rc == NF90_NOERR) then
         this%is_open = .true.
         call this%read_file_info(local_rc)
      endif
#else
      ! Fall back to serial mode if parallel NetCDF not available
      call this%open(filename, mode_str, local_rc)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_open_mpi

   !> Close NetCDF file
   subroutine netcdf_close(this, rc)
      class(NetCDFFileType), intent(inout) :: this
      integer, intent(out), optional :: rc

      integer :: local_rc

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (this%is_open) then
         local_rc = nf90_close(this%ncid)
         this%is_open = .false.
         this%ncid = -1
         this%filename = ''
      endif
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_close

   !> Check if file is open
   logical function netcdf_is_open(this)
      class(NetCDFFileType), intent(in) :: this
      netcdf_is_open = this%is_open
   end function netcdf_is_open

   !> Get filename
   function netcdf_get_filename(this) result(filename)
      class(NetCDFFileType), intent(in) :: this
      character(len=MAX_STRING_LEN) :: filename
      filename = this%filename
   end function netcdf_get_filename

   !> Read file information (dimensions, variables, etc.)
   subroutine read_file_info(this, rc)
      class(NetCDFFileType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i

      rc = 0

#ifdef NETCDF_ENABLED
      ! Get file info
      rc = nf90_inquire(this%ncid, this%ndims, this%nvars, this%ngatts, this%unlimdimid)
      if (rc /= NF90_NOERR) return

      ! Read dimension information
      do i = 1, min(this%ndims, MAX_DIMS)
         this%dim_ids(i) = i - 1  ! NetCDF uses 0-based indexing
         rc = nf90_inquire_dimension(this%ncid, this%dim_ids(i), &
            this%dim_names(i), this%dim_sizes(i))
         if (rc /= NF90_NOERR) return
      enddo
#endif

   end subroutine read_file_info

   !> Get dimension information
   subroutine netcdf_get_dimensions(this, dim_names, dim_sizes, ndims, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(out), optional :: dim_names(:)
      integer, intent(out), optional :: dim_sizes(:)
      integer, intent(out), optional :: ndims
      integer, intent(out), optional :: rc

      integer :: local_rc, i, n

      local_rc = 0

      if (present(ndims)) ndims = this%ndims

      if (present(dim_names)) then
         n = min(size(dim_names), this%ndims)
         do i = 1, n
            dim_names(i) = this%dim_names(i)
         enddo
      endif

      if (present(dim_sizes)) then
         n = min(size(dim_sizes), this%ndims)
         do i = 1, n
            dim_sizes(i) = this%dim_sizes(i)
         enddo
      endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_get_dimensions

   !> Get variable names
   subroutine netcdf_get_variables(this, var_names, nvars, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(out), optional :: var_names(:)
      integer, intent(out), optional :: nvars
      integer, intent(out), optional :: rc

      integer :: local_rc, i, varid, n
      character(len=MAX_NAME_LEN) :: var_name

      local_rc = 0

      if (present(nvars)) nvars = this%nvars

#ifdef NETCDF_ENABLED
      if (present(var_names)) then
         n = min(size(var_names), this%nvars)
         do i = 1, n
            varid = i - 1  ! NetCDF uses 0-based indexing
            local_rc = nf90_inquire_variable(this%ncid, varid, name=var_name)
            if (local_rc == NF90_NOERR) then
               var_names(i) = var_name
            else
               exit
            endif
         enddo
      endif
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_get_variables

   !> Check if variable exists
   logical function netcdf_has_variable(this, var_name)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name

      integer :: varid, rc

      netcdf_has_variable = .false.

#ifdef NETCDF_ENABLED
      if (this%is_open) then
         rc = nf90_inq_varid(this%ncid, var_name, varid)
         netcdf_has_variable = (rc == NF90_NOERR)
      endif
#endif

   end function netcdf_has_variable

   !> Check if dimension exists
   logical function netcdf_has_dimension(this, dim_name)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: dim_name

      integer :: dimid, rc

      netcdf_has_dimension = .false.

#ifdef NETCDF_ENABLED
      if (this%is_open) then
         rc = nf90_inq_dimid(this%ncid, dim_name, dimid)
         netcdf_has_dimension = (rc == NF90_NOERR)
      endif
#endif

   end function netcdf_has_dimension

   !> Read 1D real variable
   subroutine netcdf_read_real_1d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(out), allocatable :: data(:)
      integer, intent(out), optional :: rc

      integer :: local_rc, varid, ndims_var, dimids(MAX_DIMS)
      integer :: dim_sizes(1)

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (.not. this%is_open) then
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get variable info
      local_rc = nf90_inq_varid(this%ncid, var_name, varid)
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inquire_variable(this%ncid, varid, ndims=ndims_var, dimids=dimids)
      if (local_rc /= NF90_NOERR .or. ndims_var /= 1) then
         local_rc = -2  ! Wrong number of dimensions
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get dimension size
      local_rc = nf90_inquire_dimension(this%ncid, dimids(1), len=dim_sizes(1))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      ! Allocate and read data
      if (allocated(data)) deallocate(data)
      allocate(data(dim_sizes(1)))

      local_rc = nf90_get_var(this%ncid, varid, data)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_real_1d

   !> Read 2D real variable
   subroutine netcdf_read_real_2d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(out), allocatable :: data(:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc, varid, ndims_var, dimids(MAX_DIMS)
      integer :: dim_sizes(2)

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (.not. this%is_open) then
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get variable info
      local_rc = nf90_inq_varid(this%ncid, var_name, varid)
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inquire_variable(this%ncid, varid, ndims=ndims_var, dimids=dimids)
      if (local_rc /= NF90_NOERR .or. ndims_var /= 2) then
         local_rc = -2
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get dimension sizes
      local_rc = nf90_inquire_dimension(this%ncid, dimids(1), len=dim_sizes(1))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif
      local_rc = nf90_inquire_dimension(this%ncid, dimids(2), len=dim_sizes(2))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      ! Allocate and read data
      if (allocated(data)) deallocate(data)
      allocate(data(dim_sizes(1), dim_sizes(2)))

      local_rc = nf90_get_var(this%ncid, varid, data)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_real_2d

   !> Read 3D real variable
   subroutine netcdf_read_real_3d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(out), allocatable :: data(:,:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc, varid, ndims_var, dimids(MAX_DIMS)
      integer :: dim_sizes(3)

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (.not. this%is_open) then
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get variable info
      local_rc = nf90_inq_varid(this%ncid, var_name, varid)
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inquire_variable(this%ncid, varid, ndims=ndims_var, dimids=dimids)
      if (local_rc /= NF90_NOERR .or. ndims_var /= 3) then
         local_rc = -2
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get dimension sizes
      local_rc = nf90_inquire_dimension(this%ncid, dimids(1), len=dim_sizes(1))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif
      local_rc = nf90_inquire_dimension(this%ncid, dimids(2), len=dim_sizes(2))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif
      local_rc = nf90_inquire_dimension(this%ncid, dimids(3), len=dim_sizes(3))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      ! Allocate and read data
      if (allocated(data)) deallocate(data)
      allocate(data(dim_sizes(1), dim_sizes(2), dim_sizes(3)))

      ! Enable collective I/O if using MPI
#if defined(NETCDF_ENABLED) && defined(NETCDF_PARALLEL)
      if (this%use_mpi) then
         local_rc = nf90_var_par_access(this%ncid, varid, NF90_COLLECTIVE)
         if (local_rc /= NF90_NOERR) then
            if (present(rc)) rc = local_rc
            return
         endif
      endif
#endif

      local_rc = nf90_get_var(this%ncid, varid, data)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_real_3d

   !> Read 4D real variable
   subroutine netcdf_read_real_4d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(out), allocatable :: data(:,:,:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc, varid, ndims_var, dimids(MAX_DIMS)
      integer :: dim_sizes(4)

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (.not. this%is_open) then
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get variable info and dimensions (similar to 3D case)
      local_rc = nf90_inq_varid(this%ncid, var_name, varid)
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inquire_variable(this%ncid, varid, ndims=ndims_var, dimids=dimids)
      if (local_rc /= NF90_NOERR .or. ndims_var /= 4) then
         local_rc = -2
         if (present(rc)) rc = local_rc
         return
      endif

      ! Get all dimension sizes
      local_rc = nf90_inquire_dimension(this%ncid, dimids(1), len=dim_sizes(1))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif
      local_rc = nf90_inquire_dimension(this%ncid, dimids(2), len=dim_sizes(2))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif
      local_rc = nf90_inquire_dimension(this%ncid, dimids(3), len=dim_sizes(3))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif
      local_rc = nf90_inquire_dimension(this%ncid, dimids(4), len=dim_sizes(4))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      ! Allocate and read data
      if (allocated(data)) deallocate(data)
      allocate(data(dim_sizes(1), dim_sizes(2), dim_sizes(3), dim_sizes(4)))

      local_rc = nf90_get_var(this%ncid, varid, data)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_real_4d

   !> Read integer 1D variable
   subroutine netcdf_read_int_1d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      integer, intent(out), allocatable :: data(:)
      integer, intent(out), optional :: rc

      integer :: local_rc, varid, ndims_var, dimids(MAX_DIMS)
      integer :: dim_sizes(1)

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (.not. this%is_open) then
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inq_varid(this%ncid, var_name, varid)
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inquire_variable(this%ncid, varid, ndims=ndims_var, dimids=dimids)
      if (local_rc /= NF90_NOERR .or. ndims_var /= 1) then
         local_rc = -2
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inquire_dimension(this%ncid, dimids(1), len=dim_sizes(1))
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      if (allocated(data)) deallocate(data)
      allocate(data(dim_sizes(1)))

      local_rc = nf90_get_var(this%ncid, varid, data)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_int_1d

   !> Read string variable
   subroutine netcdf_read_string(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      character(len=*), intent(out) :: data
      integer, intent(out), optional :: rc

      integer :: local_rc, varid

      local_rc = 0
      data = ''

#ifdef NETCDF_ENABLED
      if (.not. this%is_open) then
         local_rc = -1
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_inq_varid(this%ncid, var_name, varid)
      if (local_rc /= NF90_NOERR) then
         if (present(rc)) rc = local_rc
         return
      endif

      local_rc = nf90_get_var(this%ncid, varid, data)
#endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_string

   ! Write operations (placeholder implementations)
   subroutine netcdf_write_real_1d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(in) :: data(:)
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_real_1d

   subroutine netcdf_write_real_2d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(in) :: data(:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_real_2d

   subroutine netcdf_write_real_3d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(in) :: data(:,:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_real_3d

   subroutine netcdf_write_real_4d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(in) :: data(:,:,:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_real_4d

   subroutine netcdf_write_int_1d(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      integer, intent(in) :: data(:)
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_int_1d

   subroutine netcdf_write_string(this, var_name, data, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: data
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_string

   ! Attribute operations (simplified)
   subroutine netcdf_read_attr_real(this, var_name, attr_name, value, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name, attr_name
      real(fp), intent(out) :: value
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      value = 0.0_fp
      if (present(rc)) rc = local_rc
   end subroutine netcdf_read_attr_real

   subroutine netcdf_read_attr_string(this, var_name, attr_name, value, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name, attr_name
      character(len=*), intent(out) :: value
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      value = ''
      if (present(rc)) rc = local_rc
   end subroutine netcdf_read_attr_string

   subroutine netcdf_write_attr_real(this, var_name, attr_name, value, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name, attr_name
      real(fp), intent(in) :: value
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_attr_real

   subroutine netcdf_write_attr_string(this, var_name, attr_name, value, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name, attr_name
      character(len=*), intent(in) :: value
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_attr_string

   !> Read FV3 coordinate information
   subroutine netcdf_read_fv3_coords(this, ak, bk, rc)
      class(NetCDFFileType), intent(in) :: this
      real(fp), intent(out), allocatable :: ak(:), bk(:)
      integer, intent(out), optional :: rc

      integer :: local_rc

      local_rc = 0

      ! Try to read ak and bk arrays
      call this%read_var('ak', ak, local_rc)
      if (local_rc /= 0) then
         if (present(rc)) rc = local_rc
         return
      endif

      call this%read_var('bk', bk, local_rc)
      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_fv3_coords

   !> Read time series data
   subroutine netcdf_read_time_series(this, var_name, data, time_index, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(out), allocatable :: data(:,:,:)
      integer, intent(in), optional :: time_index
      integer, intent(out), optional :: rc

      integer :: local_rc

      local_rc = -999  ! Not fully implemented
      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_time_series

   !> Read meteorological data bundle
   subroutine netcdf_read_meteorology(this, temperature, pressure, humidity, &
      wind_u, wind_v, surface_pressure, rc)
      class(NetCDFFileType), intent(in) :: this
      real(fp), intent(out), allocatable, optional :: temperature(:,:,:)
      real(fp), intent(out), allocatable, optional :: pressure(:,:,:)
      real(fp), intent(out), allocatable, optional :: humidity(:,:,:)
      real(fp), intent(out), allocatable, optional :: wind_u(:,:,:)
      real(fp), intent(out), allocatable, optional :: wind_v(:,:,:)
      real(fp), intent(out), allocatable, optional :: surface_pressure(:,:)
      integer, intent(out), optional :: rc

      integer :: local_rc

      local_rc = 0

      if (present(temperature)) then
         call this%read_var('T', temperature, local_rc)
         if (local_rc /= 0 .and. present(rc)) then
            rc = local_rc
            return
         endif
      endif

      if (present(pressure)) then
         call this%read_var('P', pressure, local_rc)
         if (local_rc /= 0 .and. present(rc)) then
            rc = local_rc
            return
         endif
      endif

      if (present(humidity)) then
         call this%read_var('Q', humidity, local_rc)
         if (local_rc /= 0 .and. present(rc)) then
            rc = local_rc
            return
         endif
      endif

      if (present(wind_u)) then
         call this%read_var('U', wind_u, local_rc)
         if (local_rc /= 0 .and. present(rc)) then
            rc = local_rc
            return
         endif
      endif

      if (present(wind_v)) then
         call this%read_var('V', wind_v, local_rc)
         if (local_rc /= 0 .and. present(rc)) then
            rc = local_rc
            return
         endif
      endif

      if (present(surface_pressure)) then
         call this%read_var('PS', surface_pressure, local_rc)
      endif

      if (present(rc)) rc = local_rc

   end subroutine netcdf_read_meteorology

   !> Write column data
   subroutine netcdf_write_column_data(this, var_name, data, time_index, rc)
      class(NetCDFFileType), intent(in) :: this
      character(len=*), intent(in) :: var_name
      real(fp), intent(in) :: data(:)
      integer, intent(in), optional :: time_index
      integer, intent(out), optional :: rc

      integer :: local_rc
      local_rc = -999  ! Not implemented
      if (present(rc)) rc = local_rc
   end subroutine netcdf_write_column_data

   !> Sync data to disk
   subroutine netcdf_sync(this, rc)
      class(NetCDFFileType), intent(in) :: this
      integer, intent(out), optional :: rc

      integer :: local_rc

      local_rc = 0

#ifdef NETCDF_ENABLED
      if (this%is_open) then
         local_rc = nf90_sync(this%ncid)
      endif
#endif

      if (present(rc)) rc = local_rc
   end subroutine netcdf_sync

   !> Get error message
   function netcdf_get_error_msg(this, error_code) result(error_msg)
      class(NetCDFFileType), intent(in) :: this
      integer, intent(in) :: error_code
      character(len=MAX_STRING_LEN) :: error_msg

#ifdef NETCDF_ENABLED
      error_msg = nf90_strerror(error_code)
#else
      write(error_msg, '(A,I0)') 'NetCDF error (NetCDF not available): ', error_code
#endif

   end function netcdf_get_error_msg

end module CATChemNetCDF_Mod
