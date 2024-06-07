# MESTI.jl

**MESTI** (<ins>M</ins>axwell's <ins>E</ins>quations <ins>S</ins>olver with <ins>T</ins>housands of <ins>I</ins>nputs) is an open-source software for full-wave electromagnetic simulations in frequency domain using finite-difference discretization on the [Yee lattice](https://meep.readthedocs.io/en/latest/Yee_Lattice).

MESTI implements the **augmented partial factorization (APF)** method described in [this paper](https://doi.org/10.1038/s43588-022-00370-6). While conventional methods solve Maxwell's equations on every element of the discretization basis set (which contains much more information than is typically needed), APF bypasses such intermediate solution step and directly computes the projected quantities of interest: a generalized scattering matrix given any list of input source profiles and any list of output projection profiles. It can jointly handle thousands of inputs without a loop over them, using fewer computing resources than what a conventional direct method uses to handle a single input. It is exact with no approximation beyonddiscretization.

The MESTI.jl package here provides all the features in the 2D MATLAB version [MESTI.m](https://github.com/complexphoton/MESTI.m) and additionally supports 3D vectorial systems, anisotropic *ε*, MPI parallelization, both single-precision and double-precision arithmetics, and can perform subpixel smoothing for the geometric shapes handled by [GeometryPrimitives.jl](https://github.com/stevengj/GeometryPrimitives.jl). It is written in Julia.

In the 3D case, MESTI.jl solves

$$
\left[ \nabla\times\left( \nabla\times \right) - \frac{\omega^2}{c^2}\bar{\bar{\varepsilon}}({\bf r}) \right] {\bf E}({\bf r}) = {\bf b}({\bf r}).
$$

For 2D systems in the transverse-magnetic (TM) polarization (*Ex*, *Hy*, *Hz*), MESTI.jl solves

$$
\left[ -\frac{\partial^2}{\partial y^2} -\frac{\partial^2}{\partial z^2} - \frac{\omega^2}{c^2}\varepsilon_{xx}(y,z) \right] E_x(y,z)= b(y,z),
$$

where **b**(**r**) or *b*(*y*,*z*) is the source profile.

MESTI.jl is a general-purpose solver with its interface written to provide maximal flexibility. It supports
 - Full 3D vectorial systems.
 - TM polarization for 2D systems.
 - Any tensor or scalar relative permittivity profile $\bar{\bar{\varepsilon}}({\bf r})$, which can be real-valued or complex-valued.
 - Open boundary modeled by a [perfectly matched layer (PML)](https://en.wikipedia.org/wiki/Perfectly_matched_layer) placed on any side(s), with both imaginary-coordinate and real-coordinate stretching (so the PML can accelerate the attenuation of evanescent waves in addition to attenuating the propagating waves).
- Periodic, Bloch periodic, perfect electrical conductor (PEC), and/or perfect magnetic conductor (PMC) boundary conditions.
- [Subpixel smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing) for the geometric shapes handled by [GeometryPrimitives.jl](https://github.com/stevengj/GeometryPrimitives.jl) through function [<code>mesti_subpixel_smoothing()</code>](https://github.com/complexphoton/MESTI.jl/blob/main/src/mesti_subpixel_smoothing.jl).
 - Any material dispersion $\bar{\bar{\varepsilon}}$(*ω*), since this is a frequency-domain method.
 - Any list of input source profiles (user-specified or automatically built).
 - Any list of output projection profiles (or no projection, in which case the complete field profiles are returned).
 - Real-valued or complex-valued frequency *ω*.
 - Automatic or manual choice between APF or a conventional direct solver (*e.g.*, to compute the full field profile) as the solution method.
 - Linear solver using MUMPS (requires installation) or the built-in routines in Julia (which uses UMFPACK).
 - Shared memory parallelism (with multithreaded BLAS and with OpenMP in MUMPS) and distributed memory parallelism (with MPI in MUMPS).
 - Single-precision or double-precision arithmetic.
## When to use MESTI.jl?

MESTI.jl can perform most linear-response computations for arbitrary structures, such as

- Scattering problems: transmission, reflection transport through complex media, waveguide bent, grating coupler, radar cross-section, controlled-source electromagnetic surveys, *etc*.
- Thermal emission.
- Local density of states.
- [Inverse design](https://github.com/complexphoton/metalens_inverse_design) based on the above quantities.

Since MESTI.jl can use the APF method to handle a large number of input states simultaneously, the computational advantage of MESTI.jl is the most pronounced in multi-input systems.

There are use cases that MESTI.jl can handle but is not necessarily the most efficient, such as
- Broadband response problems involving many frequencies but only a few input states. Time-domain methods like FDTD may be preferred as they can compute a broadband response without looping over the frequencies.
- Problems like plasmonics that require more than an order of magnitude difference in the discretization grid size at different regions of the structure. Finite-element methods may be preferred as they can handle varying spatial resolutions. (Finite-element methods can also adopt APF, but MESTI.jl uses finite difference with a fixed grid size.)
- Homogeneous structures with a small surface-to-volume ratio. Boundary element methods may be preferred as they only discretize the surface.

Problems that MESTI.jl does not handle:
- Nonlinear systems (*e.g.*, *χ*<sup>(2)</sup>, *χ*<sup>(3)</sup>, gain media).
- Magnetic systems (*e.g.*, spatially varying permeability *μ*).

For eigenmode computation, such as waveguide mode solver and photonic band structure computation, one can use [<code>mesti_build_fdfd_matrix.jl</code>](./src/mesti_build_fdfd_matrix.jl) to build the matrix and then compute its eigenmodes. However, we don't currently provide a dedicated function to do so.

## Installation
MESTI.jl is written and run in Julia programming language. Follow the standard process to download and install Julia. We can download Julia [here](https://julialang.org/downloads/). Suppose we installed the 1.9.3 version of Julia. After Julia is installed, we can add the path of your Julia to <code>PATH</code> through the terminal by 

```shell
export PATH=".../julia-1.9.3/bin/"
```

where  <code>...</code> is the path to your Julia.

Before installing MESTI.jl, the user first need to install the parallel version of the sparse linear solver [MUMPS](https://mumps-solver.org/index.php). Without MUMPS, MESTI.jl can still run but cannot use the APF method and will only use a conventional method with the built-in linear solver, which can be orders of magnitude slower and uses much more memory (especially in 3D and for large 2D systems). See this [MUMPS installation](./mumps) page for steps to install MUMPS. For this MESTI (v0.4.3) version, it is compatible with MUMPS versions between v5.3.3 and v5.6.2. It is not compatible with MUMPS versions v5.7.0 or v5.7.1. 

After the MUMPS installation, if you have a clean Julia environment (*i.e.* have not installed any Julia package before or have not installed new version of [Makie.jl](https://github.com/MakieOrg/Makie.jl)(v0.20 and v0.21) and [GeometryPrimitives.jl](https://github.com/stevengj/GeometryPrimitives.jl)(v0.5.0), you can install MESTI.jl (v0.4.3) by opening the command-line interface of Julia and typing:  

```julia
import Pkg; Pkg.add("MESTI")
```

On the other hand, if you have installed these two Julia packages: the new version of [Makie.jl](https://github.com/MakieOrg/Makie.jl)(v0.20 and v0.21) and [GeometryPrimitives.jl](https://github.com/stevengj/GeometryPrimitives.jl)(v0.5.0). You will encounter some compatibility problems and the installation will fail. In this case, please run 

```julia
import Pkg; Pkg.add(Pkg.PackageSpec(;name="Makie", version="0.19.12")); Pkg.add("GeometryPrimitives"); Pkg.add("MESTI")
```

to downgrade them to the old version compatible with our MESTI.jl and then install MESTI.jl(v0.4.3).

After installing MESTI.jl, if the user happens to run the command to upgrade every Julia package to the latest version, such as *Pkg.update()*, it would install the latest version of Makie.jl and make MESTI.jl not compatible with it. So please rerun the line above to downgrade to the compatible version and reinstall MESTI.jl.

To see whether you install the right version, you can always check the version of MESTI by opening the command-line interface of Julia and typing: 

```julia
import Pkg; Pkg.status("MESTI")
```

If you want our latest MESTI.jl, it should show

```julia
[8d7f31fa] MESTI v0.4.3
```

## Tests

After compiling MUMPS and installation MESTI.jl, run <code>[install_packages.jl](./test/install_packages.jl)</code> to install other Julia packages used in the tests and examples.

Now we are ready to run the following test scripts

- <code>[basic_solve.jl](./mumps/basic_solve.jl)</code>
- <code>[schur_complement.jl](./mumps/schur_complement.jl)</code>

If any of them does not run successfully, please look back at the compilation of MUMPS or the Julia interface to see if there were serious warning messages.

Then, run <code>[runtests.jl](./test/runtests.jl)</code> in the [test](./test) folder. This script runs three tests:

- <code>[matrix_solver_test.jl](./test/matrix_solver_test.jl)</code>
- <code>[interface_t_r_test.jl](./test/interface_t_r_test.jl)</code>
- <code>[unitary_test.jl](./test/unitary_test.jl)</code>

If all pass, congratulations! You are done and be able to run MESTI.jl with MUMPS solver.

## Usage Summary 

The function <code>mesti(syst, B, C, D)</code> provides the most flexibility. Structure <code>syst</code> specifies the permittivity profile, boundary conditions in *x*, *y*, and *z*, which side(s) to put PML with what parameters, the wavelength, and the discretization grid size. Any list of input source profiles can be specified with matrix <code>B</code>, each column of which specifies one source profile **b**(**r**). Any list of output projection profiles can be specified with matrix <code>C</code>. Matrix <code>D</code> is optional (treated as zero when not specified) and subtracts the baseline contribution; see the [APF paper](https://doi.org/10.1038/s43588-022-00370-6) for details.

The function <code>mesti2s(syst, in, out)</code> deals specifically with scattering problems in two-sided or one-sided geometries where $\bar{\bar{\varepsilon}}(\bf r)$  consists of an inhomogeneous scattering region with homogeneous spaces on the low (*-z*) and high (*+z*) sides, light is incident from those sides, the boundary condition in *z* is outgoing, and the boundary condition in *y* and *x* is closed (*e.g.*, periodic or PEC). The user only needs to specify the input and output sides or channel indices or wavefronts through <code>in</code> and <code>out</code>. The function <code>mesti2s()</code> automatically builds the source matrix <code>B</code>, projection matrix <code>C</code>, baseline matrix <code>D</code>, and calls <code>mesti()</code> for the computation. Flux normalization in *z* is applied automatically so, when the PML is thick enough, the full scattering matrix is unitary when $\bar{\bar{\varepsilon}}$ is Hermitian everywhere [*i.e.*, $\bar{\bar{\varepsilon}}({\bf r})=\bar{\bar{\varepsilon}}^\dagger({\bf r})$ for all ${\bf r}$].

To compute the complete field profiles, simply omit the argument <code>C</code> or <code>out</code>.

The solution method, the linear solver to use, and other options can be specified with a structure <code>opts</code> as an optional input argument to <code>mesti()</code> or <code>mesti2s()</code>; see documentation for details. They are chosen automatically when not explicitly specified.

The function [<code>mesti_build_channels()</code>](./src/mesti_build_channels.jl) can be used to build the input and/or output matrices when using <code>mesti()</code> or to determine which channels are of interest when using <code>mesti2s()</code>.

The function [<code>mesti_subpixel_smoothing()</code>](./src/mesti_subpixel_smoothing.jl) can be used to build the permittivity profile with subpixel smoothing.

## Documentation

Detailed documentation is given in comments at the beginning of the function files:
 - [<code>mesti_main.jl</code>](./src/mesti_main.jl) for <code>mesti()</code> 
 - [<code>mesti2s.jl</code>](./src/mesti2s.jl) for <code>mesti2s()</code> 
 - [<code>mesti_build_channels.jl</code>](./src/mesti_build_channels.jl) for <code>mesti_build_channels()</code> 
 - [<code>mesti_subpixel_smoothing</code>](./src/mesti_subpixel_smoothing.jl) for <code>mesti_subpixel_smoothing()</code> 

For example, typing <code>? mesti2s</code> in Julia brings up the documentation for <code>mesti2s()</code>.

## Multithreading and MPI

MESTI.jl can use both distributed memory parallelization across nodes/sockets through MPI and shared memory parallelization within one node/socket through multithreading (if MUMPS was compiled with multithreading enabled). The multithreading speed-up comes from using a multithreaded BLAS library and L0-threads feature in MUMPS. Parts of the MUMPS utilize multithreading with L0-threads through OpenMP directives. With APF, most of the computing time is spent on factorization within MUMPS (*e.g.*, see Fig 2d of the [APF paper](https://doi.org/10.1038/s43588-022-00370-6)). The factorization and solving stages within MUMPS are parallelized. The building and analyzing stages are not performance-critical and are not parallelized.

L0-threads (see the MUMPS Users' guide for details) in multithreading enhances the time performance, but marginally increases the memory usage. It is enabled by default when it is in 1D, 2D, or 2.5D ((width in *x*)*(width in *y*)/(thickness in *z*) $\geqslant 100$). In full 3D ((width in *x*)*(width in *y*)/(thickness in *z*) $< 100$), memory usage is critical and we do not use L0-threads by default. We can change from the default by setting the field `opts.use_L0_threads` in the input argument `opts`.

In MUMPS, multithreading is more efficient than MPI, both in speed and in memory usage. So, we should maximize multhreading before using MPI. For example, if we use one node with a single socket having 8 cores (where the 8 cores sharing the same memory), we should use one MPI process (*i.e.*, no MPI) with 8 threads, instead of 8 MPI processes with one thread each. As another example, if we use 3 nodes, each node has 2 sockets, and each socket has 4 cores sharing the same memory of that socket (so, 24 cores in total), we should use 6 MPI processes (one per socket) with 4 threads per MPI process, instead of 24 MPI processes with one thread each.

The default number of threads is the number of cores available on the machine (either the number of physical cores, or the number of cores requested when running the job with a scheduler like Slurm on a cluster). Therefore, we only need to launch MESTI with the number of MPI processes equaling the total number of sockets.

We can set the number of threads to be different from the default by setting the environment variable <code>OMP_NUM_THREADS</code> or the field <code>opts.nthreads_OMP</code> of the optional input argument <code>opts</code>.

To use MPI, we should prepare the script to construct the system on the main processor and call worker processors when needed. An example script and its corresponding submission script on a cluster are provided in the [MPI](./MPI) folder. 

To check the actual number of MPI processes and threads used in MUMPS, set <code>opts.verbal_solver = true</code> in the input argument and look at the standard output from MUMPS. For example, the following output

```text:Output
      executing #MPI =      2 and #OMP =      4
```

shows 2 MPI processes with 4 threads each.

## Examples

Examples in the [examples](./examples) folder illustrate the usage and the main functionalities of MESTI. Each example has its own folder, with its <code>.jl</code> and <code>.ipynb</code> script, auxiliary files specific to that example, and a <code>README.md</code> page that shows the example script with its outputs:

- [Open channel in a disordered system](./examples/2d_open_channel_through_disorder): 2D, using <code>mesti2s()</code>, transmission matrix & field profile with customized wavefronts.
- [Phase-conjugated focusing in disordered system](./examples/2d_focusing_inside_disorder_with_phase_conjugation): 2D, using <code>mesti()</code> and <code>mesti2s()</code>, customized source & field profile with customized wavefronts.
- [Reflection matrix in Gaussian-beam basis](./examples/2d_reflection_matrix_Gaussian_beams): 2D, using <code>mesti()</code>, reflection matrix in customized basis for a fully open system.
- [Metalens focusing via angular spectrum propagation](./examples/2d_metalens_focusing_via_angular_spectrum_propagation): 2D, using <code>mesti()</code>, focusing field calculated using angular spectrum propagation.

Also see the following repository:
- [metalens_inverse_design](https://github.com/complexphoton/metalens_inverse_design): Using MESTI.jl to perform multi-angle inverse design of a wide-field-of-view metalens.

## Gallery

Here are some animations from the examples above:

1. Open channel propagating through disorder
   <img src="https://github.com/complexphoton/MESTI.jl/assets/44913081/68421516-db18-4793-9ef5-304079671113" width="1080" height="720">
2. Focusing phase-conjugated light through disorder
   <img src="https://github.com/complexphoton/MESTI.jl/assets/44913081/a1ddcb6a-f43a-44c1-a0b0-6a611e57bd98" width="1080" height="720">
3. Reflection matrix of a scatterer in Gaussian-beam basis
   <img src="https://github.com/complexphoton/MESTI.jl/assets/44913081/27a99177-8ee1-4106-8bce-4e4166a1214c" width="864" height="576">
4. [Inverse designed wide-field-of-view metalens](https://github.com/complexphoton/metalens_inverse_design)
   <img src="https://github.com/complexphoton/MESTI.jl/assets/109620064/cfe707e5-445b-40b8-ab9a-6275041b1a7a" width="540" height="795">
## Acknowledgment

We thank [William Sweeney](https://github.com/wrs28) for granting us permission to integrate his MUMPS-julia interface, [MUMPS3.jl](https://github.com/wrs28/MUMPS3.jl/tree/5.3.3-update), into this package. The files bearing the mumps3 prefix in the [src](./src) directory have been adopted from MUMPS3.jl.

## Reference & Credit

For more information on the theory, capability, and benchmarks (*e.g.*, scaling of computing time, memory usage, and accuracy), please see:

- Ho-Chun Lin, Zeyu Wang, and Chia Wei Hsu. [Fast multi-source nanophotonic simulations using augmented partial factorization](https://doi.org/10.1038/s43588-022-00370-6). *Nature Computational Science* **2**, 815–822 (2022).

```bibtex
@article{2022_Lin_NCS,
  title = {Fast multi-source nanophotonic simulations using augmented partial factorization},
  author = {Lin, Ho-Chun and Wang, Zeyu and Hsu, Chia Wei},
  journal = {Nat. Comput. Sci.},
  volume = {2},
  issue = {12},
  pages = {815--822},
  year = {2022},
  month = {Dec},
  doi = {10.1038/s43588-022-00370-6}
}
```

Please cite this paper when you use MESTI.
