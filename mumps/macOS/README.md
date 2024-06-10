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

 The <code>clang</code> compiler from Apple does not support OpenMP and we need to use the <code>clang</code> compiler from LLVM instead. After installing LLVM, enter

```
export PATH=/opt/homebrew/opt/llvm/bin:$PATH
```
to override the <code>clang</code> compiler from Apple.

## Troubleshooting

### Issue 1: "clang: error: unsupported option '-fopenmp'"
This error occurs if you compile with the default <code>clang</code> compiler from Apple. Make sure you are using the <code>clang</code> compiler from LLVM by running the following command in terminal:
```
which clang
```
The output should be <code>/path/to/llvm/bin/clang</code>. If this is not the case, you need to override the <code>clang</code> compiler from Apple by running the following command in terminal:
```
export PATH=/opt/homebrew/opt/llvm/bin:$PATH
```

### Issue 2: The tests from MESTI.jl take a very long time to run
The tests in <code>MESTI.jl/test</code> and <code>MESTI.jl/examples</code> should complete in a few seconds. If the test keeps running without outputting any errors, ensure that MPI.jl is configured properly by following those [steps](../README.md#Running-MUMPS-in-Julia).

### Issue 3: "shmem: mmap: an error occurred while determining whether or not /var/folders/.../sm_segment.xxx could be created."
If you are using OpenMPI v5.0.3 or later, you may encounter this error while using MESTI.jl. The code can still complete but segmentation fault or deadlock may occur in some [cases](https://github.com/open-mpi/ompi/issues/12307). If those issues occur, you can downgrade OpenMPI to v5.0.2