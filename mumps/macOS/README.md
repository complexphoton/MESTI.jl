## Prerequisites for MUMPS installation on macOS

We need the following tools before compiling MUMPS on macOS.

### Xcode Command Line Tools

Apple's Xcode Command Line Tools (CLT) include <code>make</code>, <code>ar</code>, <code>ranlib</code>, and Apple's implementation of BLAS and LAPACK, [vecLib](https://developer.apple.com/documentation/accelerate/veclib), within its Accelerate framework. We need those.

In the next step, we will install Homebrew, which will install CLT (if not already installed). So, nothing needs to be done here.

### Homebrew

Xcode and CLT do not include a Fortran compiler. Here, we use [Homebrew](https://brew.sh/) to install one; Homebrew can also be used for the optional installation of <code>openblas</code> and <code>cmake</code> below. If you already have Homebrew installed, run <code>brew update</code> to update it. If you don't have Homebrew installed, copy the following line and paste it in terminal.
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Follow instructions from the script to install Homebrew and then to add it to your PATH.

If CLT was not installed prior, Homebrew will install it as part of the script above. 

As described in the [Homebrew installation page](https://docs.brew.sh/Installation), this installs Homebrew to <code>/opt/homebrew</code> for an Apple Silicon Mac, <code>/usr/local</code> for an Intel Mac.

After installation, enter
```
brew doctor
```
to make sure there's no outstanding issues.

### GNU compiler collection

After installing Homebrew, enter
```
brew install gcc
```
in terminal. This will install the [GNU compiler collection (GCC)](https://gcc.gnu.org/), which includes the Fortran compiler <code>gfortran</code>. If you already installed GCC through Homebrew before, this will update GCC to the most current stable version.

### ScaLAPACK, MPI and LLVM

MUMPS requires the ScaLAPACK, MPI and OpenMP libraries, which can be installed via Homebrew
```
brew install scalapack open-mpi llvm
```

 The <code>clang</code> compiler from Apple does not support OpenMP by default and we need to use the <code>clang</code> compiler from LLVM instead. After installing LLVM, enter

```
export PATH=/opt/homebrew/opt/llvm/bin:$PATH
```
to override the <code>clang</code> compiler from Apple.

### METIS

In 3D system, because METIS ordering is more efficient than AMD ordering, we should install the METIS program for graph partitioning (not to be confused with MESTI).  We can use [METIS](https://github.com/scivision/METIS/tree/743ae96033f31907d89c80e3470c0325e9a97f7b) (version 5.1.0) program for graph partitioning. We can install them in the following steps:

(a) Downloading METIS (version 5.1.0)

```shell
wget https://github.com/scivision/METIS/blob/743ae96033f31907d89c80e3470c0325e9a97f7b/archive/metis-5.1.0.tar.gz
```

(b) Decompress metis-5.1.0.tar.gz

```shell
tar zxvf metis-5.1.0.tar.gz
```

(c) Setting METIS to double precision

```shell
sed -i "43s/32/64/" metis-5.1.0/include/metis.h
```

(d) Installing METIS

```shell
cd metis-5.1.0; make config; sudo make install;
```

Then, by default, the library file, header file, and binaries will be installed in `/usr/local/lib`, `/usr/local/include`, and `/usr/local/bin`.

In some rare cases, your machine cannot find METIS libraries by itself when you run Julia interface for MUMPS. You can append the METIS libraries to your `LD_LIBRARYP_PATH`

```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LMETISDIR
```

`LMETISDIR` is the path to the folder where the METIS library is.

### Running MUMPS in Julia

You may need to configure [MPI.jl](https://juliaparallel.org/MPI.jl/stable/configuration/) before running MUMPS in Julia. The steps are straightforward using MPIPreferences.jl. First, install MPIPreferences.jl by entering
```
julia --project -e 'using Pkg; Pkg.add("MPIPreferences")'
```
in terminal. Then run <code>MPIPreferences.use_system_binary()</code> in Julia or through the command line:
```
julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```

This should automatically find the open MPI installed above.