# MESTI.jl

**MESTI** (<ins>M</ins>axwell's <ins>E</ins>quations <ins>S</ins>olver with <ins>T</ins>housands of <ins>I</ins>nputs) is an open-source software for full-wave electromagnetic simulations in frequency domain using finite-difference discretization on the [Yee lattice](https://meep.readthedocs.io/en/latest/Yee_Lattice).

MESTI implements the **augmented partial factorization (APF)** method described in [this paper](https://doi.org/10.1038/s43588-022-00370-6). While conventional methods solve Maxwell's equations on every element of the discretization basis set (which contains much more information than is typically needed), APF bypasses such intermediate solution step and directly computes the projected quantities of interest: a generalized scattering matrix given any list of input source profiles and any list of output projection profiles. It can jointly handle thousands of inputs without a loop over them, using fewer computing resources than what a conventional direct method uses to handle a single input. It is exact with no approximation beyond discretization.

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
- [Subpixel smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing) for the geometric shapes handled by [GeometryPrimitives.jl](https://github.com/stevengj/GeometryPrimitives.jl) through function [<code>mesti_subpixel_smoothing()</code>](./src/mesti_subpixel_smoothing.jl).
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

Before installing MESTI.jl, the user should first install the parallel version of the sparse linear solver [MUMPS](https://mumps-solver.org/index.php). Without MUMPS, MESTI.jl can still run but will only use the built-in linear solver, which is significantly slower and uses much more memory (especially in 3D and for large 2D systems). See this [MUMPS installation](./mumps) page for steps to install MUMPS.

To install MESTI.jl, simply open the command-line interface of Julia and type:  

```julia
import Pkg; Pkg.add("MESTI")
```

After installing MESTI.jl, you may also install other packages used in the tests and examples by running <code>[install_packages.jl](./mumps/install_packages.jl)</code>.

## Usage Summary 

The function <code>mesti(syst, B, C, D)</code> provides the most flexibility. Structure <code>syst</code> specifies the permittivity profile, boundary conditions in *x*, *y*, and *z*, which side(s) to put PML with what parameters, the wavelength, and the discretization grid size. Any list of input source profiles can be specified with matrix <code>B</code>, each column of which specifies one source profile **b**(**r**). Any list of output projection profiles can be specified with matrix <code>C</code>. Matrix <code>D</code> is optional (treated as zero when not specified) and subtracts the baseline contribution; see the [APF paper](https://doi.org/10.1038/s43588-022-00370-6) for details.

The function <code>mesti2s(syst, in, out)</code> deals specifically with scattering problems in two-sided or one-sided geometries where $\bar{\bar{\varepsilon}}(\bf r)$  consists of an inhomogeneous scattering region with homogeneous spaces on the low (*-z*) and high (*+z*) sides, light is incident from those sides, the boundary condition in *z* is outgoing, and the boundary condition in *y* and *x* is closed (*e.g.*, periodic or PEC). The user only needs to specify the input and output sides or channel indices or wavefronts through <code>in</code> and <code>out</code>. The function <code>mesti2s()</code> automatically builds the source matrix <code>B</code>, projection matrix <code>C</code>, baseline matrix <code>D</code>, and calls <code>mesti()</code> for the computation. 

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

## Examples

Examples in the [examples](./examples) folder illustrate the usage and the main functionalities of MESTI. Each example has its own folder, with its <code>.jl</code> and <code>.ipynb</code> script, auxiliary files specific to that example, and a <code>README.md</code> page that shows the example script with its outputs:

- [Open channel in a disordered system](./examples/2d_open_channel_through_disorder): 2D, using <code>mesti2s()</code>, transmission matrix & field profile with customized wavefronts.
- [Phase-conjugated focusing in disordered system](./examples/2d_focusing_phase_conjugated_light_through_disorder): 2D, using <code>mesti()</code> and <code>mesti2s()</code>, customized source & field profile with customized wavefronts.
- [Reflection matrix in Gaussian-beam basis](./examples/2d_reflection_matrix_Gaussian_beams): 2D, using <code>mesti()</code>, reflection matrix in customized basis for a fully open system.

Also see the following repository:
- [metalens_inverse_design](https://github.com/complexphoton/metalens_inverse_design): Using MESTI.jl to perform multi-angle inverse design of a wide-field-of-view metalens.

## Gallery

Here are some animations from the examples above:

1. Open channel propagating through disorder
   <img src="./examples/2d_open_channel_through_disorder/disorder_open_channel.gif" width="540" height="360"> 
2. Focusing phase-conjugated light through disorder
   <img src="./examples/2d_focusing_inside_disorder_with_phase_conjugation/phase_conjugated_focusing.gif" width="540" height="360"> 
3. Reflection matrix of a scatterer in Gaussian-beam basis:
   <img src="./examples/2d_reflection_matrix_Gaussian_beams/reflection_matrix_Gaussian_beams.gif" width="432" height="288">

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
