## Prerequisites for MUMPS installation on Linux

We need the following tools before compiling MUMPS on Linux.

### Compilers 

The compilation of MUMPS requires both C and Fortran compilers. Although both C and Fortran compilers are included in GNU Compiler Collection (GCC), we recommend using Intel [C compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html#gs.gtmcma), [Fortran compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.gtma2f), and [MPI library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html#gs.gtmtr3).

On a Linux cluster, if Intel compilers and MPI library are not loaded, we may load them on a Lmod module system by typing
```shell
module load intel
module load intel intel-mpi
```
In non-Lmod module system, we may type (for Intel oneAPI),
```shell
source .../intel/oneapi/compiler/latest/env/vars.sh
source  .../intel/oneapi/mpi/latest/env/vars.sh
```
where `...` is the path to Intel folder.

On a local Linux machine, we can download the Intel installer [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html). Make sure to use custom installation and choose the following: 
	(a) Intel Inspector
	(b) Intel MPI Library
	(c) Intel oneAPI DPC++/C++ Compiler \& Intel C++ Compiler Classic
	(d) Intel Fortran Compiler \& Intel Fortran Compiler Classic        

After installing the Intel compiler, we need to set the environment variables specific to these compilers. Run the following commands:

```shell
source opt/intel/oneapi/compiler/latest/env/vars.sh
source opt/intel/oneapi/mpi/latest/env/vars.sh
```

These commands will configure the necessary environment variables for the Intel compilers, and Intel MPI.

### BLAS, LAPACK, and ScaLAPACK

MUMPS requires both BLAS, LAPACK, and ScaLAPACK libraries, which are standard libraries on a Linux cluster. These libraries are also included in many implementations, such as MKL and OpenBLAS. 

In the example below, we use MKL. Different clusters install MKL in different paths. We will need the correct path for the linker to find corresponding BLAS, LAPACK, and ScaLAPACK libraries. In the provided `Makefile.inc`, we assume the MKL path has been exported to an environment variable called `MKLROOT`. Here we show how to export the correct `MKLROOT`. The BLAS, LAPACK, and ScaLAPACK libraries can be found under `$(MKLROOT)/lib/intel64`. 

If we are using the Lmod module system and MKL is installed, we can use 

```shell
module load intel
module load intel-mkl
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out, if MKL is available. The `MKLROOT` should be under the path `.../mkl`, where `...` is the Intel directory path.

If Lmod module system is not used in our cluster, then we can also try to find the path like `.../intel/oneapi/mkl`, where `...` is the Intel directory path. Then we can assign the `MKLROOT` by 

```shell
source .../mkl/bin/mklvars.sh intel64
echo $MKLROOT
```

`MKLROOT` is assigned and should be printed out.

On a local machine, please visit [Intel oneAPI Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) to install MKL. During the installation process, ensure that we select custom installation and choose only the Intel oneAPI Math Kernel Library. This will ensure the necessary libraries are properly installed for MUMPS to function correctly.

After installing the MKL, we can set the environment variable `MKLROOT`. Run the following commands:

```shell
source /opt/intel/oneapi/mkl/latest/env/vars.sh
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out.

Note that MUMPS supports shared memory, multithreaded parallelism through the use of multithreaded
BLAS libraries. We provide `Makefile.inc` on Linux and they activate the OpenMP feature. We can use the environment variable `OMP_NUM_THREADS` to set the number of threads.