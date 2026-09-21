!> \file netcdf_mpi_example.F90
!! \brief Example demonstrating MPI-compatible NetCDF I/O
!! \date 2025
!!
!! This example shows how to use the CATChemNetCDF module in both
!! serial and MPI-parallel modes, with automatic fallback.
!!
program netcdf_mpi_example
   use CATChemNetCDF_Mod
   use catchem_bridge_precision, only: fp
#ifdef MPI_ENABLED
   use mpi
#endif
   implicit none

   type(NetCDFFileType) :: nc_file
   real(fp), allocatable :: temperature(:,:,:)
   integer :: rc, mpi_rank, mpi_size
   character(len=256) :: filename
   logical :: file_exists

#ifdef MPI_ENABLED
   ! Initialize MPI
   call MPI_Init(rc)
   call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, rc)
   call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, rc)
#else
   mpi_rank = 0
   mpi_size = 1
#endif

   ! Example filename
   filename = 'test_meteorology.nc'

   ! Check if file exists (rank 0 only for output)
   if (mpi_rank == 0) then
      inquire(file=trim(filename), exist=file_exists)

      write(*,*) '======================================='
      write(*,*) 'CATChem NetCDF MPI Compatibility Test'
      write(*,*) '======================================='
      write(*,*) 'NetCDF available: ', netcdf_available()
      write(*,*) 'MPI NetCDF available: ', mpi_netcdf_available()
      write(*,*) 'MPI ranks: ', mpi_size
      write(*,*) 'Test file: ', trim(filename)
      write(*,*) 'File exists: ', file_exists
      write(*,*) '======================================='
   endif

#ifdef MPI_ENABLED
   ! Synchronize all ranks
   call MPI_Barrier(MPI_COMM_WORLD, rc)
#endif

   ! Try to open file with MPI support if available
   if (mpi_netcdf_available() .and. mpi_size > 1) then
      if (mpi_rank == 0) then
         write(*,*) 'Opening file with MPI-parallel NetCDF...'
      endif

#ifdef MPI_ENABLED
      call nc_file%open_mpi(filename, 'r', MPI_COMM_WORLD, rc)
#endif

      if (rc == 0) then
         if (mpi_rank == 0) then
            write(*,*) 'Success: File opened with parallel NetCDF'
         endif
      else
         if (mpi_rank == 0) then
            write(*,*) 'Warning: Parallel NetCDF failed, trying serial mode...'
         endif
         call nc_file%open(filename, 'r', rc)
      endif
   else
      ! Use serial mode
      if (mpi_rank == 0) then
         write(*,*) 'Using serial NetCDF mode...'
      endif
      call nc_file%open(filename, 'r', rc)
   endif

   if (rc == 0 .and. nc_file%is_file_open()) then
      if (mpi_rank == 0) then
         write(*,*) 'File opened successfully'
         write(*,*) 'MPI mode enabled: ', nc_file%is_mpi_enabled()
      endif

      ! Try to read a variable (if it exists)
      if (nc_file%has_variable('T')) then
         call nc_file%read_var('T', temperature, rc)
         if (rc == 0) then
            if (mpi_rank == 0) then
               write(*,*) 'Temperature data read successfully'
               write(*,*) 'Array shape: ', shape(temperature)
            endif
         else
            if (mpi_rank == 0) then
               write(*,*) 'Failed to read temperature data'
            endif
         endif
      else
         if (mpi_rank == 0) then
            write(*,*) 'Variable T not found in file'
         endif
      endif

      ! Close file
      call nc_file%close()
      if (mpi_rank == 0) then
         write(*,*) 'File closed successfully'
      endif
   else
      if (mpi_rank == 0) then
         write(*,*) 'Failed to open file: ', trim(filename)
         write(*,*) 'Error code: ', rc
         write(*,*) 'This is expected if the test file does not exist'
      endif
   endif

   if (mpi_rank == 0) then
      write(*,*) '======================================='
      write(*,*) 'Test completed'
      write(*,*) '======================================='
   endif

#ifdef MPI_ENABLED
   call MPI_Finalize(rc)
#endif

end program netcdf_mpi_example
