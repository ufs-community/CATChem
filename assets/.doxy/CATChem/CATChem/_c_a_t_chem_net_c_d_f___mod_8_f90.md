

# File CATChemNetCDF\_Mod.F90



[**FileList**](files.md) **>** [**api**](dir_da61e3e9a357748887e3ca8d7c5a0c16.md) **>** [**CATChemNetCDF\_Mod.F90**](_c_a_t_chem_net_c_d_f___mod_8_f90.md)

[Go to the source code of this file](_c_a_t_chem_net_c_d_f___mod_8_f90_source.md)

_High-level NetCDF I/O interface for CATChem with MPI support._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**catchemnetcdf\_mod**](namespacecatchemnetcdf__mod.md) <br> |




















































## Detailed Description




**Date:**

2025


This module provides a simplified, high-level interface for reading and writing NetCDF files in CATChem applications. It abstracts the low-level NetCDF operations and provides convenient routines for common atmospheric data I/O patterns.


## MPI Compatibility



This module provides both serial and MPI-parallel NetCDF I/O capabilities:


**Serial Mode (Default):**
* Uses standard NetCDF Fortran interface (netcdf module)
* Each MPI rank operates independently on files
* Suitable for rank-local I/O operations
* Always available when NetCDF is enabled




**Parallel Mode (Optional):**
* Uses parallel NetCDF-4 with MPI-IO (netcdf\_par module)
* Supports collective I/O operations across MPI ranks
* Requires NetCDF-4 built with parallel HDF5 and MPI support
* Enable with NETCDF\_PARALLEL and MPI\_ENABLED preprocessor flags
* Access via open\_mpi() method instead of open()





## Build Requirements for MPI Support



To enable MPI-parallel NetCDF support, ensure:
* NetCDF-4 library built with enable-parallel and enable-netcdf-4
* HDF5 library built with enable-parallel
* MPI library (OpenMPI, MPICH, etc.)
* CMake flags: -DNETCDF\_PARALLEL=ON -DMPI\_ENABLED=ON





## Usage Examples



**Serial I/O (current default):** 
```Fortran
use catchemnetcdf_mod
type(NetCDFFileType) :: nc_file
real(fp), allocatable :: temperature(:,:,:)

call nc_file%open('input_file.nc', 'r')
call nc_file%read_var('T', temperature)
call nc_file%close()
```



**Parallel I/O with MPI:** 
```Fortran
use catchemnetcdf_mod
use mpi
type(NetCDFFileType) :: nc_file
real(fp), allocatable :: temperature(:,:,:)
integer :: mpi_comm, rc

if (mpi_netcdf_available()) then
   call nc_file%open_mpi('input_file.nc', 'r', mpi_comm_world, rc)
   call nc_file%read_var('T', temperature)  ! Collective I/O
   call nc_file%close()
else
   ! Fallback to serial mode
   call nc_file%open('input_file.nc', 'r')
   call nc_file%read_var('T', temperature)
   call nc_file%close()
endif
```



Features:
* Simple file opening/closing with automatic error handling
* Read/write operations for 1D, 2D, 3D, and 4D arrays
* Automatic dimension and variable discovery
* Metadata reading and writing
* Time series support for atmospheric data
* FV3 coordinate system support
* MPI-parallel collective I/O support (when available)
* Automatic fallback to serial mode when parallel NetCDF unavailable 





    

------------------------------
The documentation for this class was generated from the following file `src/api/CATChemNetCDF_Mod.F90`

