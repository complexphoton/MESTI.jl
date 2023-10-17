# MESTI.jl

**MESTI** (Maxwell's Equations Solver with Thousands of Inputs) is an open-source software for full-wave electromagnetic simulations in frequency domain using finite-difference discretization on the [Yee lattice](https://meep.readthedocs.io/en/latest/Yee_Lattice).

MESTI implements the **augmented partial factorization (APF)** method described in [this paper](https://doi.org/10.1038/s43588-022-00370-6). While conventional methods solve Maxwell's equations on every element of the discretization basis set (which contains much more information than is typically needed), APF bypasses such intermediate solution step and directly computes the information of interest: a generalized scattering matrix given any list of input source profiles and any list of output projection profiles. It can jointly handle thousands of inputs without a loop over them, using fewer computing resources than what a conventional direct method uses to handle a single input. It is exact with no approximation beyond discretization.

Compared to [MESTI.m](https://github.com/complexphoton/MESTI.m) which is a 2D MATLAB version of MESTI, MESTI.jl here uses Julia with single-precision or double-precision arithmetic and considers either full 3D systems 

$$
\left[ \nabla\times\left( \nabla\times \right) - \frac{\omega^2}{c^2}\bar{\bar{\varepsilon}}(\bf r) \right] {\bf E}(\bf r) = {\bf b}(\bf r),
$$

or 2D systems in transverse-magnetic (TM) polarization (*Ex*, *Hy*, *Hz*) with

$$
\left[ -\frac{\partial^2}{\partial y^2} -\frac{\partial^2}{\partial z^2} - \frac{\omega^2}{c^2}\varepsilon_{xx}(y,z) \right] E_x(y,z)= b(y,z),
$$

where (**b**(**r**) or *b*(*y*,*z*)) is the source profile.

Note that the coordinate notation switches between MESTI.m and MESTI.jl: (x, y, z) in MESTI.jl corresponding to (z, y, x) in MESTI.m.

MESTI.jl is a general-purpose solver with its interface written to provide maximal flexibility. It supports
 - Full 3D system supporting both *s* and *p* polarizations.
 - TM polarization for 2D system.
 - Any relative permittivity profile $\bar{\bar{\varepsilon}}(\bf r)$ (or *ε*<sub>*xx*</sub>(*y*,*z*)), real-valued or complex-valued. Users can optionally average the interface pixels for [subpixel smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing) before calling MESTI.jl. It can be done by another code SBPSM.jl, which would be released soon.
 - Infinite open spaces can be described with a [perfectly matched layer (PML)](https://en.wikipedia.org/wiki/Perfectly_matched_layer) placed on any side(s), which also allows for infinite substrates, waveguides, photonic crystals, *etc*. The PML implemented in MESTI includes both imaginary-coordinate and real-coordinate stretching, so it can accelerate the attenuation of evanescent waves in addition to attenuating the propagating waves.
 - Any material dispersion $\bar{\bar{\varepsilon}}$(*ω*) can be used since this is in frequency domain.
 - Any list of input source profiles (user-specified or automatically built).
 - Any list of output projection profiles (or no projection, in which case the complete field profiles are returned).
 - Periodic, Bloch periodic, perfect electrical conductor (PEC), and/or perfect magnetic conductor (PMC) boundary conditions.
 - Real-valued or complex-valued frequency *ω*.
 - Automatic or manual choice between APF, conventional direct solver (*e.g.*, to compute the full field profile) as the solution method.
 - Linear solver using MUMPS (requires installation) or the built-in routines in Julia (which uses UMFPACK).
 - Shared memory parallelism (with multithreaded BLAS and with OpenMP in MUMPS) and distributed memory parallelism (with MPI in MUMPS).
 - Single-precision and double-precision arithmetic.

## When to use MESTI?

MESTI.jl can perform most linear-response computations for arbitrary structures, such as

- Scattering problems: transmission, reflection transport through complex media, waveguide bent, grating coupler, radar cross-section, controlled-source electromagnetic surveys, *etc*.
- Thermal emission.
- Local density of states.
- [Inverse design](https://github.com/complexphoton/APF_inverse_design) based on the above quantities.

Since MESTI can use the APF method to handle a large number of input states simultaneously, the computational advantage of MESTI is the most pronounced in multi-input systems.

There are use cases that MESTI.jl can handle but is not necessarily the most efficient, such as
- Broadband response problems involving many frequencies but only a few input states. Time-domain methods like FDTD may be preferred as they can compute a broadband response without looping over frequencies.
- Problems like plasmonics that require more than an order of magnitude difference in the discretization grid size at different regions of the structure. Finite-element methods may be preferred as they can handle varying spatial resolutions. (Finite-element methods can also adopt APF, but MESTI uses finite difference with a fixed grid size.)
- Homogeneous structures with a small surface-to-volume ratio. Boundary element methods may be preferred as they only discretize the surface.

Problems that MESTI.jl currently does not handle:
- Nonlinear systems (*e.g.*, *χ*<sup>(2)</sup>, *χ*<sup>(3)</sup>, gain media).
- Magnetic systems (*e.g.*, spatially varying permeability *μ*)

For eigenmode computation, such as waveguide mode solver and photonic band structure computation, one can use [<code>mesti_build_fdfd_matrix.jl</code>](./src/mesti_build_fdfd_matrix.jl) to build the matrix and then compute its eigenmodes. However, we don't currently provide a dedicated function to do so.

## Installation

To install MESTI.jl, simply open Julia REPL and type: 

<code>import Pkg; Pkg.add(rul="https://github.com/complexphoton/MESTI.jl")</code> 

or 

<code>import Pkg; Pkg.add("MESTI")</code> 

You may also install other necessary packages by running <code>[install_packages.jl](./mumps/install_packages.jl)</code>

However, to use the APF method, the user needs to install the parallel version of [MUMPS](https://graal.ens-lyon.fr/MUMPS/index.php) and its Julia interface [MUMPS3](https://github.com/wrs28/MUMPS3.jl/tree/5.3.3-update). Without MUMPS, MESTI will still run but will only use other methods, which generally take longer and use more memory. So, MUMPS installation is strongly recommended for large-scale multi-input simulations or whenever efficiency is important. See this [MUMPS installation](./mumps) page for steps to install MUMPS.

## Usage Summary 

The function [<code>mesti(syst, B, C, D)</code>](./src/mesti_main.m) provides the most flexibility. Structure <code>syst</code> specifies the polarization to use, permittivity profile, boundary conditions in *x* and *y*, which side(s) to put PML with what parameters, the wavelength, and the discretization grid size. Any list of input source profiles can be specified with matrix <code>B</code>, each column of which specifies one source profile **b**(**r**). Any list of output projection profiles can be specified with matrix <code>C</code>. Matrix <code>D</code> is optional (treated as zero when not specified) and subtracts the baseline contribution; see [this paper](https://doi.org/10.1038/s43588-022-00370-6) for details.

The function [<code>mesti2s(syst, in, out)</code>](./src/mesti2s.m) deals specifically with scattering problems in two-sided or one-sided geometries where $\bar{\bar{\varepsilon}}(\bf r)$  consists of an inhomogeneous scattering region with homogeneous spaces on the low (*-z*) and high (*+z*) sides, light is incident from the low side and/or high side, the boundary condition in *z* is outgoing, and the boundary condition in *y* and *x* is closed (*e.g.*, periodic or PEC). The user only needs to specify the input and output sides or channel indices or wavefronts through <code>in</code> and <code>out</code>. The function <code>mesti2s()</code> automatically builds the source matrix <code>B</code>, projection matrix <code>C</code>, baseline matrix <code>D</code>, and calls <code>mesti()</code> for the computation. Flux normalization in *z* is applied automatically and exactly, so the full scattering matrix is always unitary when  $\bar{\bar{\varepsilon}}(\bf r)=(\bar{\bar{\varepsilon}}(\bf r))^\dagger$ . 

To compute the complete field profiles, simply omit the argument <code>C</code> or <code>out</code>.

The solution method, the linear solver to use, and other options can be specified with a structure <code>opts</code> as an optional input argument to <code>mesti()</code> or <code>mesti2s()</code>; see documentation for details. They are chosen automatically when not explicitly specified.

The function [<code>mesti_build_channels()</code>](./src/mesti_build_channels.jl) can be used to build the input and/or output matrices when using <code>mesti()</code>, or to determine which channels are of interest when using <code>mesti2s()</code>.

## Documentation

Detailed documentation is given in comments at the beginning of the function files:
 - [<code>mesti_main.jl</code>](./src/mesti_main.jl) for <code>mesti()</code> 
 - [<code>mesti2s.jl</code>](./src/mesti2s.jl) for <code>mesti2s()</code> 
 - [<code>mesti_build_channels.jl</code>](./src/mesti_build_channels.jl) for <code>mesti_build_channels()</code> 

For example, typing <code>? mesti2s</code> in Julia brings up the documentation for <code>mesti2s()</code>.


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
