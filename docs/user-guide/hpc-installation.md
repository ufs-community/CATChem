# HPC Installation Guide

This guide covers installation of CATChem on High Performance Computing (HPC) systems.

## Overview

CATChem is designed to run efficiently on HPC systems with support for:

- **Message Passing Interface (MPI)** for distributed memory parallelism
- **OpenMP** for shared memory parallelism
- **Hybrid MPI+OpenMP** execution
- **GPU acceleration** with OpenACC/CUDA
- **Containerized deployment** with Singularity/Apptainer

## System Requirements

### Minimum Requirements

- **Compiler**: Modern Fortran compiler (Intel 2019+, GNU 9+, NVIDIA HPC SDK)
- **MPI**: OpenMPI 3.0+, Intel MPI 2019+, or Cray MPICH
- **Libraries**: NetCDF-Fortran, HDF5, ESMF (optional)
- **Memory**: 4GB per MPI process minimum
- **Storage**: 100MB for installation, additional for input/output data

### Recommended Requirements

- **Compiler**: Intel Fortran 2021+ or GNU 11+
- **MPI**: Latest stable version of chosen implementation
- **Memory**: 8-16GB per MPI process for operational runs
- **Network**: InfiniBand or high-speed Ethernet
- **Storage**: Parallel file system (Lustre, GPFS) for I/O intensive runs

## Common HPC Systems

### NOAA/NCEP Systems

**Hera (RDHPCS):**
```bash
# Load environment
module purge
module load intel/2021.3.0
module load impi/2021.3.0
module load netcdf/4.7.4
module load hdf5/1.10.6

# Configure build
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
```

**Orion (RDHPCS):**
```bash
# Load environment
module purge
module load intel/2020
module load impi/2020
module load netcdf/4.7.0
module load hdf5/1.10.5

# Configure build
export FC=mpiifort
export CC=mpiicc
```

### NCAR Systems

**Cheyenne:**
```bash
# Load environment
module purge
module load ncarenv/1.3
module load intel/19.1.1
module load mpt/2.22
module load netcdf/4.7.4

# Configure build
export FC=mpif90
export CC=mpicc
```

**Casper:**
```bash
# Load environment
module purge
module load ncarenv/1.3
module load gnu/9.1.0
module load openmpi/4.0.3
module load netcdf/4.7.4

# Configure build
export FC=mpif90
export CC=mpicc
```

### DOE Systems

**Summit (OLCF):**
```bash
# Load environment
module purge
module load gcc/9.3.0
module load spectrum-mpi/10.4.0.3-20210112
module load netcdf-fortran/4.5.3
module load hdf5/1.10.7

# Configure for GPU
export FC=mpif90
export CC=mpicc
export FCFLAGS="-acc -gpu=cc70 -Minfo=accel"
```

**Cori (NERSC):**
```bash
# Load environment
module purge
module load PrgEnv-intel/6.0.9
module load cray-netcdf/4.7.4.2
module load cray-hdf5/1.12.0.2

# Configure build
export FC=ftn
export CC=cc
```

## Installation Steps

### 1. Environment Setup

Create an environment script for your system:

```bash
# save as setup_catchem_env.sh
#!/bin/bash

# System-specific modules (example for generic cluster)
module purge
module load compiler/intel/2021.3.0
module load mpi/openmpi/4.1.1
module load netcdf-fortran/4.5.3
module load hdf5/1.12.1

# Build environment
export FC=mpif90
export CC=mpicc
export CXX=mpicxx

# NetCDF paths (adjust for your system)
export NETCDF_ROOT=$NETCDF_FORTRAN_ROOT
export HDF5_ROOT=$HDF5_ROOT

# Optimization flags
export FCFLAGS="-O3 -xHOST -ipo"
export CFLAGS="-O3 -xHOST -ipo"

# Parallel make
export MAKEFLAGS="-j 8"

echo "Environment configured for CATChem build"
```

### 2. Dependency Installation

**Option A: Use System Modules (Recommended)**
```bash
# Load system-provided libraries
source setup_catchem_env.sh
```

**Option B: Build Dependencies from Source**
```bash
# Build NetCDF-Fortran
cd $HOME/software
wget https://github.com/Unidata/netcdf-fortran/archive/v4.5.4.tar.gz
tar -xzf v4.5.4.tar.gz
cd netcdf-fortran-4.5.4

./configure --prefix=$HOME/software/netcdf-fortran \
           --enable-shared --enable-static \
           CPPFLAGS="-I$HDF5_ROOT/include -I$NETCDF_C_ROOT/include" \
           LDFLAGS="-L$HDF5_ROOT/lib -L$NETCDF_C_ROOT/lib"

make -j 8 && make install
```

### 3. Build CATChem

```bash
# Clone repository
git clone https://github.com/UFS-Community/CATChem.git
cd CATChem

# Setup environment
source setup_catchem_env.sh

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI=ON \
  -DOPENMP=ON \
  -DNETCDF_ROOT=$NETCDF_ROOT \
  -DHDF5_ROOT=$HDF5_ROOT

# Build
make -j 8

# Install (optional)
make install
```

### 4. Testing Installation

```bash
# Run basic tests
cd tests
./run_tests.sh

# Run parallel test
mpirun -np 4 ../build/bin/catchem_test_main

# Performance test
sbatch --ntasks=16 --time=10:00 test_performance.sh
```

## GPU Installation

### NVIDIA GPU Support

**Requirements:**
- CUDA Toolkit 11.0+
- NVIDIA HPC SDK 21.5+
- GPU-aware MPI (optional but recommended)

**Build Configuration:**
```bash
# Load GPU environment
module load nvhpc/22.2
module load cuda/11.7
module load openmpi/4.1.1+cuda

# Configure build
export FC=nvfortran
export CC=nvc
export FCFLAGS="-acc -gpu=cc70,cc80 -Minfo=accel"

# Build with GPU support
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DGPU=ON \
  -DACCELERATOR=OPENACC \
  -DMPI=ON
```

**GPU Testing:**
```bash
# Test GPU functionality
nvidia-smi
./test_gpu_kernels

# Run with GPU
export CUDA_VISIBLE_DEVICES=0,1,2,3
mpirun -np 4 --map-by ppr:1:socket:pe=1 \
  ./catchem_driver --gpu-enabled
```

## Optimized Configurations

### Intel Systems

**Build Script:**
```bash
#!/bin/bash
module load intel/2021.4.0
module load impi/2021.4.0
module load netcdf/4.7.4

export FC=mpiifort
export CC=mpiicc
export FCFLAGS="-O3 -xCORE-AVX512 -ipo -no-prec-div -fp-model fast=2"
export CFLAGS="-O3 -xCORE-AVX512 -ipo"

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI=ON \
  -DOPENMP=ON \
  -DOPTIMIZATION=AGGRESSIVE
```

### AMD Systems

**Build Script:**
```bash
#!/bin/bash
module load gcc/11.2.0
module load openmpi/4.1.1
module load netcdf-fortran/4.5.4

export FC=mpif90
export CC=mpicc
export FCFLAGS="-O3 -march=znver3 -mtune=znver3 -ffast-math"
export CFLAGS="-O3 -march=znver3 -mtune=znver3"

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI=ON \
  -DOPENMP=ON
```

### ARM Systems

**Build Script:**
```bash
#!/bin/bash
module load gcc/11.1.0
module load openmpi/4.1.0
module load netcdf-fortran/4.5.3

export FC=mpif90
export CC=mpicc
export FCFLAGS="-O3 -mcpu=neoverse-n1 -ffast-math"
export CFLAGS="-O3 -mcpu=neoverse-n1"

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI=ON \
  -DOPENMP=ON
```

## Container Installation

### Singularity/Apptainer

**Definition File:**
```singularity
Bootstrap: docker
From: ubuntu:22.04

%environment
    export PATH=/opt/catchem/bin:$PATH
    export LD_LIBRARY_PATH=/opt/catchem/lib:$LD_LIBRARY_PATH

%post
    # Install dependencies
    apt-get update && apt-get install -y \
        build-essential \
        gfortran \
        libopenmpi-dev \
        libnetcdff-dev \
        libhdf5-dev \
        cmake \
        git

    # Build CATChem
    cd /opt
    git clone https://github.com/UFS-Community/CATChem.git
    cd CATChem
    mkdir build && cd build

    cmake .. -DCMAKE_INSTALL_PREFIX=/opt/catchem \
             -DMPI=ON -DOPENMP=ON
    make -j 4 install

%runscript
    exec /opt/catchem/bin/catchem_driver "$@"
```

**Build and Run:**
```bash
# Build container
singularity build catchem.sif catchem.def

# Run container
singularity exec catchem.sif catchem_driver input.yml

# MPI run with container
mpirun singularity exec catchem.sif catchem_driver input.yml
```

## Job Submission Examples

### SLURM

**Basic Job Script:**
```bash
#!/bin/bash
#SBATCH --job-name=catchem
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=16
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --account=myproject

# Load environment
source setup_catchem_env.sh

# Set runtime environment
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Run CATChem
srun ./catchem_driver config.yml
```

**GPU Job Script:**
```bash
#!/bin/bash
#SBATCH --job-name=catchem-gpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=01:00:00
#SBATCH --partition=gpu

# Load GPU environment
module load nvhpc/22.2 cuda/11.7

# Set GPU environment
export CUDA_VISIBLE_DEVICES=0,1,2,3
export CUDA_MPS=1

# Run with GPU
srun --gpu-bind=closest ./catchem_driver --gpu config.yml
```

### PBS/Torque

**Job Script:**
```bash
#!/bin/bash
#PBS -N catchem
#PBS -l nodes=4:ppn=16
#PBS -l walltime=02:00:00
#PBS -q standard
#PBS -A myproject

cd $PBS_O_WORKDIR

# Load environment
source setup_catchem_env.sh

# Run CATChem
mpirun -hostfile $PBS_NODEFILE -np 64 ./catchem_driver config.yml
```

### LSF

**Job Script:**
```bash
#!/bin/bash
#BSUB -J catchem
#BSUB -n 64
#BSUB -R "span[ptile=16]"
#BSUB -W 02:00
#BSUB -q general
#BSUB -P myproject

# Load environment
source setup_catchem_env.sh

# Run CATChem
mpirun.lsf ./catchem_driver config.yml
```

## Performance Tuning

### MPI Configuration

**OpenMPI:**
```bash
export OMPI_MCA_btl=^openib
export OMPI_MCA_pml=ucx
export OMPI_MCA_btl_tcp_if_include=ib0
```

**Intel MPI:**
```bash
export I_MPI_FABRICS=shm:ofi
export I_MPI_OFI_PROVIDER=mlx
export I_MPI_PIN_DOMAIN=omp
```

### Memory Optimization

```bash
# Set memory limits
ulimit -s unlimited
export OMP_STACKSIZE=256M

# NUMA awareness
export OMP_PROC_BIND=true
export OMP_PLACES=cores
```

### I/O Optimization

```bash
# For Lustre file systems
lfs setstripe -c 8 -S 1M output_directory/

# For GPFS
mmchattr --set-pool-policy=data1 output_directory/
```

## Troubleshooting

### Common Build Issues

**Missing Dependencies:**
```bash
# Check for required libraries
ldd build/bin/catchem_driver

# Verify module loading
module list
echo $NETCDF_ROOT
```

**Compiler Issues:**
```bash
# Test compiler
$FC --version
$FC -show  # For wrapper compilers

# Test MPI
mpirun --version
mpirun -np 2 hostname
```

### Runtime Issues

**Memory Problems:**
```bash
# Check memory usage
/usr/bin/time -v ./catchem_driver config.yml

# Monitor with system tools
top -p $(pgrep catchem)
```

**Performance Issues:**
```bash
# Profile application
export TMPDIR=/tmp
mpirun -np 4 valgrind --tool=callgrind ./catchem_driver config.yml
```

## System-Specific Notes

### Cray Systems

- Use PrgEnv modules instead of direct compiler modules
- May need to set `CRAYPE_LINK_TYPE=dynamic`
- Use `ftn` and `cc` compiler wrappers

### IBM Power Systems

- Use IBM XL compilers for best performance
- Consider POWER9 specific optimizations
- GPU systems may require special MPI configuration

### Cloud Systems

**AWS ParallelCluster:**

- Use Intel or ARM optimized AMIs
- Configure EFA for high-performance networking
- Consider spot instances for cost optimization

**Google Cloud:**

- Use HPC-optimized machine types
- Configure with Cloud filestore for shared storage
- Consider preemptible instances

## Support

For HPC-specific installation issues:

1. **Check system documentation** for your specific HPC system
2. **Contact system administrators** for module and configuration help
3. **File issues** on the **[CATChem GitHub Issues](https://github.com/ufs-community/CATChem/issues)**
4. **Join discussions** on **[GitHub Discussions](https://github.com/ufs-community/CATChem/discussions)**

## References

- [Build System Documentation](build-system.md)
- [Performance Guide](advanced_topics/performance.md)
