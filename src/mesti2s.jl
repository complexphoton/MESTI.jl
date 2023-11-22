# Export composite data types
export channel_type
export channel_index
export wavefront

# Export a function mesti2s()
export mesti2s

mutable struct channel_type
    # A composite data type to specify channel type
    # See also: mesti and mesti2s   
    side::String
    polarization::String   
    channel_type() = new()    
end

mutable struct channel_index
    # A composite data type to specify channel index
    # See also: mesti and mesti2s   

    # Used in 3D systems
    ind_low_s::Vector{Int64}
    ind_low_p::Vector{Int64}
    ind_high_s::Vector{Int64}
    ind_high_p::Vector{Int64}
    
    # Used in 2D systems
    ind_low::Vector{Int64}
    ind_high::Vector{Int64}
    
    channel_index() = new()
end

mutable struct wavefront
    # A composite data type to specify wavefront, which is linear combinations of channels
    # See also: mesti and mesti2s   

    # Used in 3D systems
    v_low_s::Union{Array{Int64,2}, Array{Float64,2}, Array{ComplexF64,2}}
    v_low_p::Union{Array{Int64,2}, Array{Float64,2}, Array{ComplexF64,2}}
    v_high_s::Union{Array{Int64,2}, Array{Float64,2}, Array{ComplexF64,2}}
    v_high_p::Union{Array{Int64,2}, Array{Float64,2}, Array{ComplexF64,2}} 
    
    # Used in 2D systems
    v_low::Union{Array{Int64,2}, Array{Float64,2}, Array{ComplexF64,2}}
    v_high::Union{Array{Int64,2}, Array{Float64,2}, Array{ComplexF64,2}}
    
    wavefront() = new()
end   


"""
    MESTI2S Solves frequency-domain scattering problems in a two-sided geometry.
        ---3D field profile---
        (Ex, Ey, Ez, channels, info) = MESTI2S(syst, input) returns the spatial field profiles
        of Ex(x,y,z), Ey(x,y,z), and Ez(x,y,z) satisfying
            [(∇ × ∇ ×) - (omega/c)^2*epsilon(x,y,z)]*[Ex(x,y,z); Ey(x,y,z); Ez(x,y,z)] = 
            i*omega*mu_0*[Jx(x,y,z); Jy(x,y,z); Jz(x,y,z)], where 
        The relative permittivity profile epsilon(x,y,z), frequency omega, and boundary conditions 
        are specified by structure "syst". epsilon(x,y,z) must have homogeneous spaces on the low (-z) 
        and high (+z) sides, with an outgoing boundary in z for the scattered waves and a closed
        boundary in x and y.  
            Note that relative permittivity profile epsilon(x,y,z) is a rank 2 tenor: 
                [epsilon_xx, epsilon_xy, epsilon_xz; 
                epsilon_yx, epsilon_yy, epsilon_yz; 
                epsilon_zx, epsilon_zy, epsilon_zz] in general. 
        Users can specify the diagonal terms only (epsilon_xx(x,y,z), epsilon_yy(x,y,z), and epsilon_zz(x,y,z)) 
        or all of them.
        The incident wavefronts from low and/or high are specified by variable "input".
            The returned "Ex", "Ey", "Ez" is a 4D array, such as Ex(:,:,:,i), 
        being the field profile Ex given the i-th input source profile. Same data structure for Ey and Ez. 
        The information of the computation is returned in structure "info".

        ---2D TM field profile---
        (Ex, channels, info) = MESTI2S(syst, input) returns the spatial
        field profiles of Ex(y,z) for scattering problems of 2D transverse-magnetic (TM) fields satisfying
        [- (d/dy)^2 - (d/dz)^2 - (omega/c)^2*epsilon_xx(y,z)] Ex(y,z) = 0,
        where Ex is each a superposition of incident and scattered waves.
        The returned "Ex" is a 3D array, with Ex(:,:,i) being the total (incident + scattered) field profile
        of Ex given the i-th incident wavefront.
        
        ---Scattering matrix S---
        (S, channels, info) = MESTI2S(syst, input, output) returns the scattering matrix
        S, where "input" and "output" specify either the list of input/output channels or
        the input/output wavefronts. When the MUMPS solver is available,
        this is typically done by computing the Schur complement of an augmented
        matrix K through a partial factorization.
        
        (Ex, Ey, Ez, channels, info) = MESTI2S(syst, input, opts), 
        (Ex, channels, info) = MESTI2S(syst, input, opts), and
        (S, channels, info) = MESTI2S(syst, input, output, opts) allow detailed options to
        be specified with structure "opts".
    
        In mesti2s(), for 3D cases, the boundary condition in x and y must be closed 
        e.g., periodic or PEC. For 2D cases, the boundary condition in y must be closed. 
        Given the closed boundary, the set of transverse modes forms a
        complete and orthonormal basis of propagating and evanescent channels.
        The inputs and outputs are specified in the basis of these propagating 
        channels, with coefficients normalized with respect to the flux in the
        longitudinal (z) direction. Properties of those channels are given by
        mesti_build_channels().

        When in 3D cases an open boundary in x or y is of interest (in 2D cases an open boundary 
        in y is of interest), the appropriate input/output basis depends on the problem, 
        so the user should use the more general function mesti() and will need to specify 
        the input and output matrices B and C.
        
        MESTI only considers nonmagnetic materials.
        
        This file builds the input and output channels using mesti_build_channels(),
        builds the matrices B and C, and then calls mesti() to solve the scattering
        problems.
        
        === Input Arguments ===
        syst (Syst struct; required):
            A structure that specifies the system, used to build the FDFD matrix A.
            It contains the following fields:
            syst.epsilon_xx (numeric array or matrix, real or complex):
                For 3D systems, an nx_Ex-by-ny_Ex-by-nz_Ex array discretizing the relative permittivity
                profile epsilon_xx(x,y,z). Specifically, syst.epsilon_xx(n,m,l) is the scalar
                epsilon_xx(n,m,l) averaged over a cube with volume (syst.dx)^3 centered at
                the point (x_n, y_m, z_l) where Ex(x,y,z) is located on the Yee lattice.
                Outside the scattering region (with z < 0 and z > H), epsilon_xx(x,y,z)
                is given by scalars syst.epsilon_low and syst.epsilon_high for the
                homogeneous spaces on the two sides. Note that nx_Ez = 0 corresponds
                to H = 0 (ie, no scattering region) and is allowed.
                For 2D TM fields, an ny_Ex-by-nz_Ex matrix discretizing the relative permittivity
                profile epsilon_xx(y,z).
            syst.epsilon_xy (numeric array or matrix, real or complex, optional):
                For 3D, an nx_Ez-by-ny_Ez-by-nz_Ex array discretizing the relative permittivity
                profile epsilon_xy(x,y,z). Specifically, syst.epsilon_xy(n,m,l) is the scalar
                epsilon_xy(n,m,l) averaged over a cube with volume (syst.dx)^3 centered at
                the low corner on the Yee lattice (n,m,l).
            syst.epsilon_xz (numeric array or nothing, real or complex, optional):    
                Discretizing the relative permittivity profile epsilon_xz(x,y,z), 
                analogous to syst.epsilon_xy.
            syst.epsilon_yx (numeric array or nothing, real or complex, optional):    
                Discretizing the relative permittivity profile epsilon_yx(x,y,z), 
                analogous to syst.epsilon_xy.
            syst.epsilon_yy (numeric array or nothing, real or complex, required for 3D):    
                Discretizing the relative permittivity profile epsilon_yy(x,y,z), 
                analogous to syst.epsilon_xx.
            syst.epsilon_yz (numeric array or nothing, real or complex, optional):    
                Discretizing the relative permittivity profile epsilon_yz(x,y,z), 
                analogous to syst.epsilon_xy.
            syst.epsilon_zx (numeric array or nothing, real or complex, optional):    
                Discretizing the relative permittivity profile epsilon_zx(x,y,z), 
                analogous to syst.epsilon_xy.
            syst.epsilon_zy (numeric array or nothing, real or complex, optional):    
                Discretizing the relative permittivity profile epsilon_zy(x,y,z), 
                analogous to syst.epsilon_xy.
            syst.epsilon_zz (numeric array or nothing, real or complex, required for 3D):    
                Discretizing the relative permittivity profile epsilon_zz(x,y,z), 
                analogous to syst.epsilon_xx.
            syst.epsilon_low (real scalar; required):
                Relative permittivity of the homogeneous space on the low.
            syst.epsilon_high (real scalar or nothing; optional):
                Relative permittivity of the homogeneous space on the high. If
                syst.epsilon_high is not given or is nothing, for 3D cases the system
                will be one-sided, terminated on the high with a PEC boundary with 
                Ez(n,m,l) = 0 at l = nz_Ez + 1, and for 2D TM caes the system
                will be one-sided, terminated on the high with a PEC boundary with 
                Ex(m,l) = 0 at l = nx_Ex + 1.
            syst.length_unit (string; optional):
                Length unit, such as micron, nm, or some reference wavelength. This
                code only uses dimensionless quantities, so syst.length_unit is never
                used. This syst.length_unit is meant to help the user interpret the
                units of (x,y,z), dx, wavelength, kx_B, ky_B, etc.
            syst.wavelength (numeric scalar, real or complex; required):
                Vacuum wavelength 2*pi*c/omega, in units of syst.length_unit.
            syst.dx (positive scalar; required):
                Discretization grid size, in units of syst.length_unit.    
            syst.xBC (string or nothing; required unless syst.kx_B is specified):
                Boundary condition (BC) at the two ends in x direction, effectively
                specifying Ex(n,m,l) at n=0 and n=nx_Ex+1,
                           Ey(n,m,l) at n=0 and n=nx_Ey+1,
                           Ez(n,m,l) at n=0 and n=nx_Ez+1.    
                one pixel beyond the computation domain. Available choices are:
                "Bloch"    - Ex(n+nx_Ex,m,l) = Ex(n,m,l)*exp(1i*syst.kx_B*nx_Ex*syst.dx),
                             Ey(n+nx_Ey,m,l) = Ey(n,m,l)*exp(1i*syst.kx_B*nx_Ey*syst.dx),
                             Ez(n+nx_Ez,m,l) = Ez(n,m,l)*exp(1i*syst.kx_B*nx_Ez*syst.dx).   
                "periodic" - equivalent to "Bloch" with syst.kx_B = 0
                "PEC"      - Ex(0,m,l) = Ex(1,m,l); Ex(nx_Ex+1,m,l) = Ez(nx_Ex,m,l),
                             Ey(0,m,l) = Ey(nx_Ey+1,m,l) = 0,
                             Ez(0,m,l) = Ez(nx_Ez+1,m,l) = 0.   
                "PMC"      - Ex(0,m,l) = Ex(nx_Ex+1,m,l) = 0,
                             Ey(0,m,l) = Ey(1,m,l); Ey(nx_Ey+1,m,l) = Ey(nx_Ey,m,l),
                             Ez(0,m,l) = Ez(1,m,l); Ez(nx_Ez+1,m,l) = Ez(nx_Ez,m,l).    
                "PECPMC"   - Ex(0,m,l) = Ex(1,m,l); Ex(nx_Ex+1,m,l) = 0,
                             Ey(0,m,l) = 0; Ey(nx_Ey+1,m,l) = Ey(nx_Ey,m,l),
                             Ez(0,m,l) = 0; Ez(nx_Ez+1,m,l) = Ez(nx_Ez,m,l),    
                "PMCPEC"   - Ex(0,m,l) = 0; Ex(nx_Ex+1,m,l) = Ex(nx_Ex,m,l),
                             Ey(0,m,l) = Ey(1,m,l); Ey(nx_Ey+1,m,l) = 0,
                             Ez(0,m,l) = Ez(1,m,l); Ez(nx_Ez+1,m,l) = 0.
                where PEC stands for perfect electric conductor and PMC stands for perfect
                magnetic conductor.
                    Note that this xBC also defines a complete and orthonormal set of
                transverse modes, upon which the input and output channels in 
                arguments "input" and "output" are defined; mesti2s() does not support PML
                in x direction because a closed boundary is necessary for defining
                such a transverse basis.
                    Here, syst.xBC is required, with no default choice (except when
                syst.kx_B is given, in which case syst.xBC = "Bloch" is automatically
                used).
                    For 2D TM fields, xBC = nothing.
            syst.yBC (string; optional):
                Boundary condition in y direction, analogous to syst.xBC.    
            syst.kx_B (real scalar; optional):
                Bloch wave number in x direction, in units of 1/syst.length_unit.
                syst.kx_B is only used when syst.xBC = "Bloch". It is allowed to
                specify a complex-valued syst.kx_B, but a warning will be displayed.
            syst.ky_B (real scalar; optional):
                Bloch wave number in y direction, analogous to syst.kx_B.
            syst.zPML (a vector of PML structure; required):
                Parameters of the perfectly matched layer (PML) used to simulate an
                open boundary, which attenuates outgoing waves with minimal reflection. 
                    Note that in mesti2s(), the PML is placed in the homogeneous spaces
                specified by syst.epsilon_low and syst.epsilon_high, outside of the
                scattering region specified by syst.epsilon_xx/Ey/Ez. (This is different 
                from the more general function mesti(), where syst.epsilon_xx/Ey/Ez
                specifies the entire simulation domain, so PML is placed inside 
                syst.epsilon_xx/Ey/Ez.)
                    To use the same PML on both sides or to use PML on the low of a
                one-sided geometry, set syst.zPML to be a scalar structure with the
                following fields:    
                    npixels (positive integer scalar; required): Number of PML pixels.
                        This number of pixels is added in addition to the
                        scattering region.
                    power_sigma (non-negative scalar; optional): 
                        Power of the polynomial grading for the conductivity sigma; 
                        defaults to 3.
                    sigma_max_over_omega (non-negative scalar; optional):
                        Conductivity at the end of the PML; defaults to
                        0.8*(power_sigma+1)/((2*pi/wavelength)*dx*sqrt(epsilon_bg)).
                        where epsilon_bg is the average relative permittivity along the
                        last slice of the PML. This is used to attenuate propagating
                        waves.
                    power_kappa (non-negative scalar; optional): 
                        Power of the polynomial grading for the real-coordinate-stretching 
                        factor kappa; defaults to 3.
                    kappa_max (real scalar no smaller than 1; optional):
                        Real-coordinate-stretching factor at the end of the PML;
                        defaults to 15. This is used to accelerate the attenuation of
                        evanescent waves. kappa_max = 1 means no real-coordinate
                        stretching.
                    power_alpha (non-negative scalar; optional): 
                        Power of the polynomial grading for the CFS alpha factor; 
                        defaults to 1.
                    alpha_max_over_omega (non-negative scalar; optional): 
                        Complex-frequency-shifting (CFS) factor at the beginning 
                        of the PML. This is typically used in time-domain simulations 
                        to suppress late-time (low-frequency) reflections. 
                        We don't use it by default (alpha_max_over_omega = 0) 
                        since we are in frequency domain.
                We use the following PML coordinate-stretching factor:
                    s(p) = kappa(p) + sigma(p)./(alpha(p) - i*omega)
                with
                    sigma(p)/omega = sigma_max_over_omega*(p.^power_sigma),
                    kappa(p) = 1 + (kappa_max-1)*(p.^power_kappa),
                    alpha(p)/omega = alpha_max_over_omega*((1-p).^power_alpha),
                where omega is frequency, and p goes linearly from 0 at the beginning
                of the PML to 1 at the end of the PML. 
                    When there is just one PML element in the vector, 
                    syst.zPML = [PML],
                then the same PML parameters are used in the both sides of two-sided geometry, 
                or the PML parameters are used in low side of one-sided geometry.
                    When multiple sets of PML parameters are used in the two-sided geometry 
                (e.g., a thinner PML on one side, a thicker PML on another side), these
                parameters can be specified with a vector of PML strcture.
                    syst.zPML = [PML_low, PML_high],
                with PML_low and PML_high each being a structure containing the above
                fields; they can specify different PML parameters on low and high sides.
                    With real-coordinate stretching, PML can attenuate evanescent waves
                more efficiently than free space, so there is no need to place free
                space in front of PML.
                    The PML thickness should be chosen based on the acceptable level of
                reflectivity given the discretization resolution and the range of wave
                numbers (i.e., angles) involved; more PML pixels gives lower
                reflectivity. Typically 10-40 pixels are sufficient.
        input (channel_type, channel_index, or wavefront structure; required):
            The set of input channels or input wavefronts
                To specify all propagating channels on one side or on both sides for 
                polarizations, use ''input = channel_type()'', then it contains the following fields:
                    input.side (string, required): specify all propagating channels on 
                        sides. Available choices are:
                            "low"    - specify all (input.polarization)-polarization 
                                    propagating channels on the low side
                            "high"   - specify all (input.polarization)-polarization 
                                    propagating channels on the high side
                            "both"   - specify all (input.polarization)-polarization  
                                    propagating channels on the both sides
                    input.polarization (string, optional): specify polarizations on 
                        sides for 3D systems. Available choices are:
                            "s"      - specify s-polarization channels on the input.side
                            "p"      - specify p-polarization channels on the input.side
                            "both"   - specify both polarizations (s and p) channels on
                                    on the input.side
                        By default, input.polarization = "both";                  
                    To specify a subset of the propagating channels use ''input = channel_index()'',
                    then it contains the following fields:
                        For 3D systems:
                        input.ind_low_s (integer vector): Vector containing the indices of
                            propagating channels incident on the low side with s-polarization
                        input.ind_low_p (integer vector): Vector containing the indices of
                            propagating channels incident on the low side with p-polarization
                        input.ind_high_s (integer vector): Vector containing the indices of
                            propagating channels incident on the high side with s-polarization
                        input.ind_high_p (integer vector): Vector containing the indices of
                            propagating channels incident on the high side with p-polarization
                        One can provide only one or more of them.
                            The user can first use mesti_build_channels() to get the indices, wave
                        numbers, and transverse profiles of the propagating channels; based on
                        that, the user can specify the list of channels of interest through
                        input.ind_low_s, input.ind_low_p, input.ind_high_s, and input.ind_high_p 
                        above, or a list of customized wavefronts through input.v_low_s, 
                        input.v_low_p, input.v_high_s, or input.v_high_p below.
                        For 2D TM case:
                        input.ind_low (integer vector): Vector containing the indices of
                            propagating channels incident on the low side.
                        input.ind_high (integer vector): Vector containing the indices of
                            propagating channels incident on the high side.
                    To specify a custom input wavefronts, a superposition of multiple 
                    propagating channels, use ''input = wavefront()'', then it contains the following fields:
                        For 3D systems:
                        input.v_low_s (numeric matrix): Matrix where each column specifies the
                            coefficients of s-polarization propagating channels on the low 
                            for one input wavefront from the low; the wavefront is a superposition 
                            of all propagating channels with the superposition coefficients 
                            given by that column of input.v_low_s. size(input.v_low_s, 1) must equal 
                            N_prop_low, the total number of propagating channels on the low; 
                            size(input.v_low_s, 2) is the number of input wavefronts.
                        input.v_low_p (numeric matrix): Analogous to to input.v_low_s, but 
                            specifying input p-polarization wavefronts from the low instead.
                        input.v_high_s (numeric matrix): Analogous to to input.v_low_s, but 
                            specifying input s-polarization wavefronts from the high instead.
                        input.v_high_p (numeric matrix): Analogous to to input.v_low_s, but 
                            specifying input p-polarization wavefronts from the high instead.    
                    Note that the input wavefronts from the low and the input wavefronts
                    from the high and two polarizations are treated as separate inputs. In other 
                    words, each input either comes from the low or comes from the high on one
                    polarization. If an input with incidence from both sides or both polarization 
                    is of interest, the user can manually superimpose results from the separate
                    computations.
                        For 2D TM case:
                        input.v_low (numeric matrix): Matrix where each column specifies the
                            coefficients of propagating channels on the low for one input
                            wavefront from the low, with the superposition coefficients given by that 
                            column of input.v_low. size(input.v_low, 1) must equal N_prop_low, the total 
                            number of propagating channels on the low; size(input.v_L, 2) is the number
                            of input wavefronts.
                        input.v_high (numeric matrix): Analogous to to input.v_low, but specifying
                            input wavefronts from the high instead.
        output (channel_type, channel_index, wavefront structure, or nothing; optional):
            The set of output channels or output wavefronts.
                When out = nothing or when out is omitted, no output projection is used, and 
            the spatial field profiles Ex(x,y,z), Ey(x,y,z), and Ez(x,y,z) corresponding to 
            the set of inputs are returned.
                When out is given, the scattering matrix is returned, with the output
            basis of the scattering matrix specified by "output". In this case, out
            follows the same syntax as the argument "input".    
                To specify all propagating channels on one side or on both sides for 
                polarizations, use ''output = channel_type()'', then it contains the following fields:
                    output.side (string, required): specify all propagating channels on sides. 
                        Available choices are:
                            "low"    - specify all propagating channels on the low side for 
                                    polarization specified by output.polarization
                            "high"   - specify all propagating channels on the high side for 
                                    polarization specified by output.polarization
                            "both"   - specify all propagating channels on the both sides
                                    (low and high) for polarization specified by 
                                    output.polarization
                    output.polarization (string, optional): specify polarizations on sides for 3D systems. 
                        Available choices are:
                            "s"      - specify s-polarization channels on the output.side
                            "p"      - specify p-polarization channels on the output.side
                            "both"   - specify both polarizations (s and p) channels on
                                    on the output.side
                        By default, output.polarization = "both";                  
                    To specify a subset of the propagating channels use ''output = channel_index()'',
                    then it contains the following fields:
                        For 3D systems:
                        output.ind_low_s (integer vector): Vector containing the indices of
                            propagating channels incident on the low side with s-polarization
                        output.ind_low_p (integer vector): Vector containing the indices of
                            propagating channels incident on the low side with p-polarization
                        output.ind_high_s (integer vector): Vector containing the indices of
                            propagating channels incident on the high side with s-polarization
                        output.ind_high_p (integer vector): Vector containing the indices of
                            propagating channels incident on the high side with p-polarization
                        One can provide only one or more of them.
                    The user can first use mesti_build_channels() to get the indices, wave
                    numbers, and transverse profiles of the propagating channels; based on
                    that, the user can specify the list of channels of interest through
                    output.ind_low_s, output.ind_low_p, output.ind_high_s, and output.ind_high_p 
                    above, or a list of customized wavefronts through output.v_low_s, output.v_low_p, 
                    output.v_high_s, or output.v_high_p below.
                        For 2D TM case:
                        output.ind_low (integer vector): Vector containing the indices of
                            propagating channels incident on the low side.
                        output.ind_high (integer vector): Vector containing the indices of
                            propagating channels incident on the high side.
                    To specify a custom output wavefronts, a superposition of multiple 
                propagating channels, use ''output = wavefront()'', then it contains the following fields:
                    For 3D systems:
                    output.v_low_s (numeric matrix): Matrix where each column specifies the
                        coefficients of s-polarization propagating channels on the low 
                        for one output wavefront from the low; the wavefront is a superposition 
                        of all propagating channels with the superposition coefficients 
                        given by that column of output.v_low_s. size(output.v_low_s, 1) must equal 
                        N_prop_low, the total number of propagating channels on the low; 
                        size(output.v_low_s, 2) is the number of output wavefronts.
                    output.v_low_p (numeric matrix): Analogous to to output.v_low_s, but 
                        specifying output p-polarization wavefronts from the low instead.
                    output.v_high_s (numeric matrix): Analogous to to output.v_low_s, but 
                        specifying output s-polarization wavefronts from the high instead.
                    output.v_high_p (numeric matrix): Analogous to to output.v_low_s, but 
                        specifying output p-polarization wavefronts from the high instead.        
                Here, each row of the scattering matrix corresponds to projection onto an
            output wavefront specified by a column of output.v_low_s, output.v_low_p, 
            output.v_high_s, or output.v_high_p.
                    For 2D TM case:
                    input.v_low (numeric matrix): Matrix where each column specifies the
                    coefficients of propagating channels on the low for one input
                    wavefront from the low, with the superposition coefficients given by that column of input.v_low.
                    size(input.v_low, 1) must equal N_prop_low, the total number of
                    propagating channels on the low; size(input.v_L, 2) is the number
                    of input wavefronts.
                    input.v_high (numeric matrix): Analogous to to input.v_low, but specifying
                    input wavefronts from the high instead.
        opts (Opts structure; optional):
            A structure that specifies the options of computation. 
            It can contain the following fields (all optional):
            opts.verbal (boolean scalar; optional, defaults to true):
                Whether to print system information and timing to the standard output.
            opts.nz_low (non-negative integer scalar; optional, defaults to 0):
                Number of pixels of homogeneous space on the low (syst.epsilon_low) to
                include when returning the spatial field profile; not used for
                scattering matrix computations.
            opts.nz_high (non-negative integer scalar; optional, defaults to 0):
                Number of pixels of homogeneous space on the high (syst.epsilon_high) to
                include when returning the spatial field profile; not used for
                scattering matrix computations. Note that opts.nx_high can still be used
                in one-sided geometries where syst.epsilon_high is not given; the field
                profile on the high is simply zero in such case.
            opts.solver (string; optional):
                The solver used for sparse matrix factorization. Available choices
                are (case-insensitive):
                    "MUMPS"  - (default when MUMPS is available) To use MUMPS, MUMPS must
	                    be installed and the Julia environment variable "MUMPS_PREFIX"
	                    should be specified.
                    "JULIA" -  (default when MUMPS is not available) Uses the built-in 
                            lu() function in JULIA, which uses UMFPACK. 
                MUMPS is faster and uses less memory than lu(), and is required for
                the APF method.
            opts.method (string; optional):
                The solution method. Available choices are (case-insensitive):
                    "APF" - Augmented partial factorization. C*inv(A)*B is obtained
                            through the Schur complement of an augmented matrix
                            K = [A,B;C,0] using a partial factorization. Must have
                            opts.solver = "MUMPS". This is the most efficient method,
                            but it cannot be used for computing the full field profile
                            X=inv(A)*B or with iterative refinement.
                    "FG"  - Factorize and group. Factorize A=L*U, and obtain C*inv(A)*B
                            through C*inv(U)*inv(L)*B with optimized grouping. Must
                            have opts.solver = "JULIA". This is slightly better than
                            "FS" when MUMPS is not available, but it cannot be used for
                            computing the full field profile X=inv(A)*B.
                    "FS"  - Factorize and solve. Factorize A=L*U, solve for X=inv(A)*B
                            with forward and backward substitutions, and project with
                            C as C*inv(A)*B = C*X. Here, opts.solver can be either
                            "MUMPS" or "JULIA", and it can be used for computing
                            the full field profile X=inv(A)*B or with iterative
                            refinement.
                    "C*inv(U)*inv(L)*B"   - Same as "FG".    
                    "factorize_and_solve" - Same as "FS".
                By default, if C is given and opts.iterative_refinement = false, then
                "APF" is used when opts.solver = "MUMPS", and "C*inv(U)*inv(L)*B" is
                used when opts.solver = "JULIA". Otherwise, "factorize_and_solve" is
                used.
            opts.symmetrize_K (boolean scalar; optional):
                Whether or not to pad input and/or output channels and perform
                permutations to make matrix K = [A,B;C,0] symmetric when computing its
                Schur complement, which lowers computing time and memory usage. Such
                channel padding and permutation is reversed afterwards and does not
                affect what mesti2s() returns.
            opts.clear_syst (boolean scalar; optional, defaults to false):
                When opts.clear_syst = true, variable "syst" will be cleared in the
                caller's workspace to reduce peak memory usage. This can be used when
                syst.epsilon_xx, syst.epsilon_yy, and syst.epsilon_zz take up
                significant memory and are not needed after calling mesti().
            opts.clear_memory (boolean scalar; optional, defaults to true):
                Whether or not to clear variables inside mesti() to reduce peak memory
                usage.
            opts.verbal_solver (boolean scalar; optional, defaults to false):
                Whether to have the solver print detailed information to the standard
                output. Note the behavior of output from MUMPS depends on compiler.
            opts.use_single_precision_MUMPS (boolean scalar; optional, defaults to true):
                Whether to use single precision version of MUMPS; used only when 
                opts.solver = "MUMPS". Using single precision version of MUMPS can 
                reduce memory usage and computing time.
            opts.use_METIS (boolean scalar; optional, defaults to false in 2D and to true in 3D):
                Whether to use METIS (instead of the default AMD) to compute the
                ordering in MUMPS. Using METIS can sometimes reduce memory usage
                and/or factorization and solve time, but it typically takes longer at
                the analysis (i.e., ordering) stage in 2D. In 3D METIS is general better 
                than AMD.
            opts.nrhs (positive integer scalar; optional):
                The number of right-hand sides (number of columns of the input matrix
                B) to consider simultaneously, used only when opts.method =
                "factorize_and_solve" and C is given. Defaults to 1 if
                opts.iterative_refinement = true, 10 if opts.solver = "MUMPS" with
                opts.iterative_refinement = false, 4 otherwise.
            opts.store_ordering (boolean scalar; optional, defaults to false):
                Whether to store the ordering sequence (permutation) for matrix A or
                matrix K; only possible when opts.solver = "MUMPS". If
                opts.store_ordering = true, the ordering will be returned in
                info.ordering.
            opts.ordering (positive integer vector; optional):
                A user-specified ordering sequence for matrix A or matrix K, used only
                when opts.solver = "MUMPS". Using the ordering from a previous
                computation can speed up (but does not eliminate) the analysis stage.
                The matrix size must be the same, and the sparsity structure should be
                similar among the previous and the current computation.
            opts.analysis_only (boolean scalar; optional, defaults to false):
                When opts.analysis_only = true, the factorization and solution steps
                will be skipped, and S = nothing will be returned. The user can use
                opts.analysis_only = true with opts.store_ordering = true to return
                the ordering for A or K; only possible when opts.solver = "MUMPS".
            opts.nthreads_OMP (positive integer scalar; optional):
                Number of OpenMP threads used in MUMPS; overwrites the OMP_NUM_THREADS
                environment variable.
            opts.parallel_dependency_graph (logical scalar; optional):
                If MUMPS is multithread, whether to use parallel dependency graph in MUMPS.
                This typically improve the time performance, but marginally increase 
                the memory usage.
            opts.iterative_refinement (boolean scalar; optional, defaults to false):
                Whether to use iterative refinement in MUMPS to lower round-off
                errors. Iterative refinement can only be used when opts.solver =
                "MUMPS" and opts.method = "factorize_and_solve" and C is given, in
                case opts.nrhs must equal 1. When iterative refinement is used, the
                relevant information will be returned in info.itr_ref_nsteps,
                info.itr_ref_omega_1, and info.itr_ref_omega_2.
            opts.use_continuous_dispersion (boolean scalar; optional):
                Whether to use the dispersion equation of the continuous wave equation
                when building the input/output channels. Defaults to false, in which
                case the finite-difference dispersion is used.
            opts.n0 (real numeric scalar; optional, defaults to 0):
                Center of the 1D transverse mode profile with periodic or Bloch periodic
                boundary condition, u(m,a) = exp(i*kx(a)*syst.dx*(m-m0))/sqrt(nx),
                where kx(a) = syst.kx_B + a*(2*pi/nx*syst.dx) and nx = nx_Ex = nx_Ey = nx_Ez.
            opts.m0 (real numeric scalar; optional, defaults to 0):
                Center of the 1D transverse mode profile with periodic or Bloch periodic
                boundary conditions, analogous to opts.n0.
            opts.use_BLR (logical scalar; optional, defaults to false):
                Whether to use block low-rank approximation in MUMPS to possibly lower computational
                cost (but in most case it does not). It can only be used when opts.solver = "MUMPS".
            opts.threshold_BLR (positive real scalar; optional):
                The dropping parameter controls the accuracy of the block low-rank approximations. 
                It can only be used when opts.solver = "MUMPS" and opts.use_BLR = true.
                Please refer to the section of BLR API in MUMPS userguide.
            opts.icntl_36 (positive integer scalar; optional):
                It controls the choice of the BLR factorization variant. 
                It can only be used when opts.solver = "MUMPS" and opts.use_BLR = true.
                Please refer to the section of BLR API in MUMPS userguide.
            opts.icntl_38 (positive integer scalar; optional):
                It estimated compression rate of LU factors.
                It can only be used when opts.solver = "MUMPS" and opts.use_BLR = true.
                Please refer to the section of BLR API in MUMPS userguide.

        === Output Arguments ===
        ---When users specify out---
        field_profiles (4D array):
            For field-profile computations (i.e., when "out" is not given), the
            returned field_profiles is a 4D array containing the field profiles, such
            that field_profiles(:,:,:,a) is the total (incident + scattered) field of
            Ex/Ey/Ez corresponding to the a-th input wavefront, with
                Recall that nz_Ez = nz_Ex + 1 = nz_Ey + 1 = H/syst.dx + 1 for a 
            two-sided geometry; nz_Ez = nz_Ex = nz_Ey = H/syst.dx for a one-sided geometry.
                By default, opts.nz_low = opts.nz_high = 0, and field_profiles only contain
            field_profiles in the scattering region. When opts.nz_low and opts.nz_high are
            specified, field_profiles will also contain opts.nz_low pixels in the homogeneous 
            space on the low, opts.nz_high pixel in the homogeneous space on the high.
                To be specific, the field_profiles contains:
                For 3D systems:
                    Ex (4D array):
                    electrical field profile for Ex component
                    (nx_Ex, ny_Ex, nz_Ex, number of inputs) = size(Ex)
                    Ey (4D array):
                    electrical field profile for Ex component
                    (nx_Ey, ny_Ey, nz_Ey, number of inputs) = size(Ey)
                    Ez (4D array):
                    electrical field profile for Ex component
                    (nx_Ez, ny_Ez, nz_Ez, number of inputs) = size(Ez)
                For 2D TM case:
                    Ex (3D array):
                    electrical field profile for Ex component
                    (ny_Ex, nz_Ex, number of inputs) = size(Ex)
        ---or when users  do not specify out---
        S (numeric matrix):
            For scattering-matrix computations (i.e., when ''output'' is given), S is the
            scattering matrix, such that S(b,a) is the flux-normalized field
            coefficient in the b-th propagating output channel (or the b-th output
            wavefront) given incident wave in the a-th propagating input channel (or
            the a-th input wavefront).
                When all propagating channels on one side or on both sides are
            requested, e.g. with ''input = channel_type()'', ''input.side = high'' or 
            ''output = channel_type()'', ''output.side = both'',  matrix S includes 
            channels.low.N_prop propagating channels on the low and/or
            2*channels.high.N_prop on the high, with longitudinal wave
            numbers given by vectors channels.low.kzdx_prop and channels.high.kzdx_prop.
                The phases of the elements of the scattering matrix depend on the
            reference plane. For channels on the low, the reference plane is at
            z = 0. For channels on the high, the reference plane is at z = H.
                When a subset of the propagating channels are requested, in 3D
            ''input = channel_index()'' and ''output = channel_index()''
            with in.ind_low_s, in.ind_low_p, in.ind_high_s, in.ind_high_p, 
            out.ind_low_s, out.ind_low_p, out.ind_low_p and/or out.ind_high_p, 
            matrix S includes such subset of the propagating channels. 
            In 2D, it used with in.ind_low, in.ind_high, out.ind_low, and out.ind_high.
                When the input wavefronts and/or output basis, i.e. ''input = wavefront()''
            and/or ''output = wavefront()'' is specified in 3D by input.v_low_s, input.v_low_p, 
            input.v_high_s, input.v_high_p, output.v_low_s, output.v_low_p, output.v_high_s, 
            and/or output.v_high_p. In 2D, it used with input.v_low, input.v_high, output.v_low, 
            and output.v_high. Matrix S is the full scattering matrix with the input and/or output 
            basis changed.
        channels (Channel structure):
            A structure returned by function mesti_build_channels() that contains
            properties of the propagating and evanescent channels of the homogeneous
            spaces on the low and high. Type "? mesti_build_channels" for more
            information.   
        info (Info structure):
            A structure that contains the following fields:
            info.opts (Opts structure):
                The final "opts" used, excluding the user-specified matrix ordering.
            info.timing_init (non-negative scalar):
            info.timing_build (non-negative scalar):
            info.timing_analyze (non-negative scalar):
            info.timing_factorize (non-negative scalar):
            info.timing_total (non-negative scalar):    
            info.zPML (two-element cell array);
                PML parameters on the sides of z direction.
            info.ordering_method (string; optional):
                Ordering method used in MUMPS.
            info.ordering (positive integer vector; optional):
                Ordering sequence returned by MUMPS when opts.store_ordering = true.
            info.itr_ref_nsteps (integer vector; optional):
                Number of steps of iterative refinement for each input, if
                opts.iterative_refinement = true; 0 means no iterative refinement.
            info.itr_ref_omega_1 (real vector; optional):
                Scaled residual omega_1 at the end of iterative refinement for each
                input; see MUMPS user guide section 3.3.2 for definition.
            info.itr_ref_omega_2 (real vector; optional):
                Scaled residual omega_2 at the end of iterative refinement for each
                input; see MUMPS user guide section 3.3.2 for definition.
            
        See also: mesti_build_channels, mesti2s
"""
function mesti2s(syst::Syst, input::Union{channel_type, channel_index, wavefront}, output::Union{channel_type, channel_index, wavefront, Nothing}, opts::Union{Opts, Nothing})
    # Make deepcopy of them to avoid mutating input argument 
    syst = deepcopy(syst); input = deepcopy(input); output = deepcopy(output); opts = deepcopy(opts) 
    
    ## Part 1.1: Check validity of syst, assign default values to its fields, and parse BC and PML specifications
    t0 = time() 

    if ~(isdefined(syst, :epsilon_low) && isa(syst.epsilon_low, Real)); throw(ArgumentError("Input argument syst must have field \"epsilon_low\" and be a real scalar.")); end    
    if ~isdefined(syst, :wavelength); throw(ArgumentError("Input argument syst must have field \"wavelength\".")); end
    if ~isdefined(syst, :dx); throw(ArgumentError("Input argument syst must have field \"dx\".")); end
    if ~(syst.dx > 0); throw(ArgumentError("syst.dx must be a positive scalar.")); end
    
    # Take care of the 2D TM cases
    if ndims(syst.epsilon_xx) == 2
        use_2D_TM = true
        if (isdefined(syst, :epsilon_yy) && ~isa(syst.epsilon_yy, Nothing)) || (isdefined(syst, :epsilon_zz) && ~isa(syst.epsilon_zz, Nothing)) || (isdefined(syst, :epsilon_xy) && ~isa(syst.epsilon_xy, Nothing)) || (isdefined(syst, :epsilon_xz) && ~isa(syst.epsilon_xz, Nothing)) || (isdefined(syst, :epsilon_yx) && ~isa(syst.epsilon_yx, Nothing)) || (isdefined(syst, :epsilon_yz) && ~isa(syst.epsilon_yz, Nothing)) || (isdefined(syst, :epsilon_zx) && ~isa(syst.epsilon_zx, Nothing)) || (isdefined(syst, :epsilon_zy) && ~isa(syst.epsilon_zy, Nothing))
            throw(ArgumentError("Only syst.epsilon_xx is required for 2D TM fields Ex(y,z), but other components should not be given or should be nothing"))
        end
    else
        use_2D_TM = false
    end
    
    if ~use_2D_TM && ~(isdefined(syst, :epsilon_xx) && isdefined(syst, :epsilon_yy) && isdefined(syst, :epsilon_zz))
        throw(ArgumentError("Input argument syst must have field \"epsilon_xx\", \"epsilon_yy\", and \"epsilon_zz\" for 3D systems."))
    end    
    
    # Check that the user did not accidentally use options only in mesti()
    if isdefined(syst, :PML)
        throw(ArgumentError("syst.PML is not used in mesti2s(); use syst.zPML instead."))
    elseif isdefined(syst, :zBC)
        throw(ArgumentError("syst.zBC is not supported in mesti2s(); use mesti() if the boundary condition behind PML is not PEC."))                
    elseif isdefined(syst, :kz_B)
        throw(ArgumentError("syst.kz_B is not supported in mesti2s(); use mesti() if Bloch periodic BC in z is needed."))
    elseif isdefined(syst, :PML_type) && (lowercase(syst.PML_type) != "upml")
        throw(ArgumentError("syst.PML_type = \"UPML\" is the only optiona in mesti2s(); use mesti() if SC-PML is needed."))
    end
   
    # syst.epsilon_high is an optional argument; determines whether the system is one-sided
    if ~isdefined(syst, :epsilon_high) || isa(syst.epsilon_high, Nothing)
        two_sided = false
        syst.epsilon_high = nothing       
    elseif ~isa(syst.epsilon_high, Real) 
        throw(ArgumentError("syst.epsilon_high must be a real scalar, if given."))    
    else
        two_sided = true        
    end    
    
    if ~use_2D_TM
        # Number of grid points in x, y, and z for Ex, Ey, and Ez
        (nx_Ex, ny_Ex, nz_Ex) = size(syst.epsilon_xx)
        (nx_Ey, ny_Ey, nz_Ey) = size(syst.epsilon_yy)
        (nx_Ez, ny_Ez, nz_Ez) = size(syst.epsilon_zz)

        if maximum([nx_Ex*ny_Ex, nx_Ey*ny_Ey, nx_Ez*ny_Ez]) == 0 
            throw(ArgumentError("Total number of pixel cannot be 0 for all Ex, Ey, and Ez components in the tranverse xy plane."))
        end        
    else
        # Number of grid points in y and z for Ex
        (ny_Ex, nz_Ex) = size(syst.epsilon_xx)
    end
    
    dn = 0.5;  # the source/detection is half a pixel away from z=0 and z=d
    
    if ~use_2D_TM
        # Check boundary condition in x
        if isdefined(syst, :kx_B) 
            if isdefined(syst, :xBC) && lowercase(syst.xBC) != "bloch"
                throw(ArgumentError("When syst.kx_B is given, syst.xBC must be \"Bloch\" if specified."))
            end
            syst.xBC = "Bloch"
            # mesti_build_channels() uses kx_B*periodicity as the input arguments xBC for Bloch BC
            xBC = (syst.kx_B)*(nx_Ex*syst.dx) # dimensionless
        else
            if ~isdefined(syst, :xBC)
                throw(ArgumentError("Input argument syst must have non-empty field \"xBC\" when syst.kx_B is not given."))      
            elseif ~(lowercase(syst.xBC) in ["bloch", "periodic", "pec", "pmc", "pecpmc", "pmcpec"])
                throw(ArgumentError("syst.xBC = \"$(syst.xBC)\" is not a supported option; type ''? mesti2s'' for supported options."))
            elseif lowercase(syst.xBC) == "bloch"
                throw(ArgumentError("syst.xBC = \"Bloch\" but syst.kx_B is not given."))
            end
            xBC = syst.xBC
        end
    end
    
    # Check boundary condition in y
    if isdefined(syst, :ky_B) 
        if isdefined(syst, :yBC) && lowercase(syst.yBC) != "bloch"
            throw(ArgumentError("When syst.ky_B is given, syst.yBC must be \"Bloch\" if specified."))
        end
        syst.yBC = "Bloch"
        # mesti_build_channels() uses ky_B*periodicity as the input arguments yBC for Bloch BC
        if ~use_2D_TM
            yBC = (syst.ky_B)*(ny_Ey*syst.dx) # dimensionless
        else
            yBC = (syst.ky_B)*(ny_Ex*syst.dx) # dimensionless
        end
    else
        if ~isdefined(syst, :yBC)
            throw(ArgumentError("Input argument syst must have non-empty field \"yBC\" when syst.ky_B is not given."))       
        elseif ~(lowercase(syst.yBC) in ["bloch", "periodic", "pec", "pmc", "pecpmc", "pmcpec"])
            throw(ArgumentError("syst.yBC = \"$(syst.yBC)\" is not a supported option; type ''? mesti2s'' for supported options."))
        elseif lowercase(syst.yBC) == "bloch"
            throw(ArgumentError("syst.yBC = \"Bloch\" but syst.ky_B is not given."))
        end
        yBC = syst.yBC
    end

    # Check zPML
    if ~isdefined(syst, :zPML) 
        throw(ArgumentError("syst.zPML must be given."))
    elseif isa(syst.zPML, Nothing)
        throw(ArgumentError("syst.zPML must be a PML structure or a vector of PML structure.")) 
    elseif isa(syst.zPML, PML)
        syst.zPML = [syst.zPML] 
    end

    if length(syst.zPML) > 2 || length(syst.zPML) == 0
        throw(ArgumentError("The length of syst.zPML must be 1 or 2."))          
    end
    
    if ~(length(syst.zPML)==2) && two_sided
        if isdefined(syst.zPML, :side) && syst.zPML.side != "both"
            throw(ArgumentError("syst.zPML.direction should be assigned as \"both\", if given, when the user just provided only one PML object is provided for a two-sided geometry."))
        end
        # Apply the same zPML on both sides    
        syst.zPML = [deepcopy(syst.zPML[1]), deepcopy(syst.zPML[1])]
        syst.zPML[1].side = "-"; syst.zPML[2].side = "+" 
    elseif length(syst.zPML)==2 && ~two_sided
        throw(ArgumentError("For a one-sided geometry, the boundary condition on the high side is PEC; syst.zPML only specifies PML on the low side and must be a PML structure or a vector containing only one PML structure."))
    end

    # Convert BC to take care of lowercase or uppercase
    if ~use_2D_TM
        xBC = convert_BC(xBC, "x")
    end
    yBC = convert_BC(yBC, "y")
    
    if ~use_2D_TM
        # Check number of grid points with the boundary conditions
        if nx_Ey != nx_Ez; throw(ArgumentError("Number of grids along x provided by syst.epsilon_yy and syst.epsilon_zz should be same.")); end
        if ny_Ex != ny_Ez; throw(ArgumentError("Number of grids along y provided by syst.epsilon_xx and syst.epsilon_zz should be same.")); end
        if nz_Ex != nz_Ey; throw(ArgumentError("Number of grids along z provided by syst.epsilon_xx and syst.epsilon_yy should be same.")); end
        check_BC_and_grid(xBC, nx_Ex, nx_Ey, nx_Ez, "x")
        check_BC_and_grid(yBC, ny_Ex, ny_Ey, ny_Ez, "y")
        check_BC_and_grid("PEC", nz_Ex, nz_Ey, nz_Ez, "z")

        if (isdefined(syst, :epsilon_xy) && ~isa(syst.epsilon_xy, Nothing) && ~(size(syst.epsilon_xy) == (nx_Ez, ny_Ez, nz_Ex)))
            throw(ArgumentError("The size of syst.epsilon_xy should be should be (size(syst.epsilon_zz, 1), size(syst.epsilon_zz, 2), size(syst.epsilon_xx, 3)) = ($(size(syst.epsilon_zz, 1)), $(size(syst.epsilon_zz, 2)), $(size(syst.epsilon_xx, 3)))."))
        end
        if (isdefined(syst, :epsilon_xz) && ~isa(syst.epsilon_xz, Nothing) && ~(size(syst.epsilon_xz) == (nx_Ey, ny_Ex, nz_Ey)))
            throw(ArgumentError("The size of syst.epsilon_xz should be should be (size(syst.epsilon_yy, 1), size(syst.epsilon_xx, 2), size(syst.epsilon_yy, 3)) = ($(size(syst.epsilon_yy, 1)), $(size(syst.epsilon_xx, 2)), $(size(syst.epsilon_yy, 3)))."))
        end
        if (isdefined(syst, :epsilon_yx) && ~isa(syst.epsilon_yx, Nothing) && ~(size(syst.epsilon_yx) == (nx_Ez, ny_Ez, nz_Ey)))
            throw(ArgumentError("The size of syst.epsilon_yx should be should be (size(syst.epsilon_zz, 1), size(syst.epsilon_zz, 2), size(syst.epsilon_yy, 3)) = ($(size(syst.epsilon_zz, 1)), $(size(syst.epsilon_zz, 2)), $(size(syst.epsilon_yy, 3)))."))
        end
        if (isdefined(syst, :epsilon_yz) && ~isa(syst.epsilon_yz, Nothing) && ~(size(syst.epsilon_yz) == (nx_Ey, ny_Ex, nz_Ex)))
            throw(ArgumentError("The size of syst.epsilon_yz should be should be (size(syst.epsilon_yy, 1), size(syst.epsilon_xx, 2), size(syst.epsilon_xx, 3)) = ($(size(syst.epsilon_yy, 1)), $(size(syst.epsilon_xx, 2)), $(size(syst.epsilon_xx, 3)))."))
        end
        if (isdefined(syst, :epsilon_zx) && ~isa(syst.epsilon_zx, Nothing) && ~(size(syst.epsilon_zx) == (nx_Ey, ny_Ez, nz_Ey)))
            throw(ArgumentError("The size of syst.epsilon_zx should be should be (size(syst.epsilon_yy, 1), size(syst.epsilon_zz, 2), size(syst.epsilon_yy, 3)) = ($(size(syst.epsilon_yy, 1)), $(size(syst.epsilon_zz, 2)), $(size(syst.epsilon_yy, 3)))."))
        end
        if (isdefined(syst, :epsilon_zy) && ~isa(syst.epsilon_zy, Nothing) && ~(size(syst.epsilon_zy) == (nx_Ez, ny_Ex, nz_Ex)))
            throw(ArgumentError("The size of syst.epsilon_zy should be should be (size(syst.epsilon_zz, 1), size(syst.epsilon_xx, 2), size(syst.epsilon_xx, 3)) = ($(size(syst.epsilon_zz, 1)), $(size(syst.epsilon_xx, 2)), $(size(syst.epsilon_xx, 3)))."))
        end
    end

    # nz_extra is the number of homogeneous-space pixels to be added to syst.epsilon_xx, epsilon_yx, epsilon_zx, epsilon_xy, syst.epsilon_yy, epsilon_zy, epsilon_xz, epsilon_yz  and syst.epsilon_zz in z direction.
    # The two elements of nz_extra corresponds to the low and high sides.
    if two_sided
        n_sides = 2
        nz_extra = Vector{Int64}(undef, 2)
        str_zBC = ["PML", "PML"] # to be used for printing later
    else
        n_sides = 1
        nz_extra = Vector{Int64}(undef, 2)
        nz_extra[2] = 0
        str_zBC = ["PML", "PEC"]
    end
        
    syst.PML = Vector{PML}() # to be used in mesti()
        
    # Loop over syst.zPML and handle PML parameters, if specified
    str_sides = ["-", "+"]
    for ii = 1:n_sides
        PML_ii = syst.zPML[ii]
        if isdefined(PML_ii, :direction) && PML_ii.direction != "z"
            throw(ArgumentError("syst.zPML[$(ii)].direction must be \"z\", if given."))
        else
            PML_ii.direction = "z" # to be used in mesti()
        end
        if isdefined(PML_ii, :side) && PML_ii.side != str_sides[ii]
            throw(ArgumentError("syst.zPML[$(ii)].side must be $(str_sides[ii]), if given."))
        else
            PML_ii.side = str_sides[ii] # to be used in mesti()
        end
        
        # Number of PML pixels must be given
        # Other fields are optional and will be checked in mesti_build_fdfd_matrix()
        if ~isdefined(PML_ii, :npixels)
            throw(ArgumentError("syst.zPML[$(ii)] must contain field \"npixels\"."))
        end
        if ~(PML_ii.npixels >= 0)
            throw(ArgumentError("syst.zPML[$(ii)].npixels must be a non-negative integer scalar."))
        end

        # Defaults to no spacer
        if ~isdefined(PML_ii, :npixels_spacer)
            PML_ii.npixels_spacer = 0
        end
        if ~(PML_ii.npixels_spacer >= 0)
            throw(ArgumentError("syst.zPML[$(ii)].npixels_spacer must be a non-negative integer scalar, if given."))
        end
        
        # number of pixels in z to be added outside of syst.epsilon_xx, syst.epsilon_yy, and syst.epsilon_zz 
        nz_extra[ii] = 1 + PML_ii.npixels + PML_ii.npixels_spacer

        PML_ii.npixels_spacer = nothing # these will not be used in mesti()
        push!(syst.PML, PML_ii) # add this PML structure to the vector of PML structure
    end
    syst.zPML = nothing # mesti() will throw warning if syst.zPML is given
    
    if ~use_2D_TM
        # Total number of pixels in z 
        nz_tot_Ex = nz_Ex + sum(nz_extra)
        nz_tot_Ey = nz_Ey + sum(nz_extra)
        nz_tot_Ez = nz_Ez + sum(nz_extra)

        # Total number of grid points for Ex, Ey, and Ez    
        nt_tot_Ex = nx_Ex*ny_Ex*nz_tot_Ex
        nt_tot_Ey = nx_Ey*ny_Ey*nz_tot_Ey
        nt_tot_Ez = nx_Ez*ny_Ez*nz_tot_Ez
    else
        # Total number of pixels in z 
        nz_tot_Ex = nz_Ex + sum(nz_extra)
        
        # Total number of grid points for Ex    
        nt_tot_Ex = ny_Ex*nz_tot_Ex
    end
    
    ## Part 1.2: Check the arguments opts and assign default values 
    
    ## We continue from here next time
    if opts == nothing
        opts = Opts()
    end
    
    # Check that the user did not accidentally use options only in mesti()
    if isdefined(opts, :prefactor) && ~isa(opts.prefactor, Nothing) 
        throw(ArgumentError("opts.prefactor is not used in mesti2s(); the -2i prefactor is automatically included."))
    end    
    
    # We return the scattering matrix S=C*inv(A)*B if output is given from input argument; else we return the spatial field profiles
    # This opts.return_field_profile is not specified by the user; it will be returned as info.opts.return_field_profile to help debugging
    opts.return_field_profile = isa(output, Nothing)
    
    # By default, for field-profile computation, we only return result within the scattering region, setting nz_low = nz_high = 0.
    if opts.return_field_profile
        if ~isdefined(opts, :nz_low) || isa(opts.nz_low, Nothing)
            opts.nz_low = 0
        elseif ~(opts.nz_low >= 0)
            throw(ArgumentError("opts.nz_low must be a non-negative integer scalar, if given."))
        end
        if ~isdefined(opts, :nz_high) || isa(opts.nz_high, Nothing)
            opts.nz_high = 0;
        elseif ~(opts.nz_high >= 0 )
            throw(ArgumentError("opts.nz_high must be a non-negative integer scalar, if given."))
        end
    else
        if isdefined(opts, :nz_low) && ~isa(opts.nz_low, Nothing)
            @warn "opts.nz_low is not used for scattering matrix computation; will be ignored."
            opts.nz_low = nothing                
        end
        if isdefined(opts, :nz_high) && ~isa(opts.nz_high, Nothing)
            @warn ("opts.nz_high is not used for scattering matrix computation; will be ignored.")
            opts.nz_high = nothing
        end
    end

    # Turn on verbal output by default
    if ~isdefined(opts, :verbal)
        opts.verbal = true
    elseif ~isa(opts.verbal, Bool)
        throw(ArgumentError("opts.verbal must be a boolean, if given."))
    end

    # By default, we don't clear syst in the caller's workspace
    if ~isdefined(opts, :clear_syst) 
        opts.clear_syst = false
    elseif ~isa(opts.clear_syst, Bool)
        throw(ArgumentError("opts.clear_syst must be a boolean, if given."))
    end
    
    # By default, we will clear internal variables    
    if ~isdefined(opts, :clear_memory) 
        opts.clear_memory = true
    elseif ~isa(opts.clear_memory, Bool)
        throw(ArgumentError("opts.clear_memory must be a boolean, if given."))
    end

    # Use MUMPS for opts.solver when it is available
    MUMPS_available = haskey(ENV,"MUMPS_PREFIX")
    solver_specified = true   
    if ~isdefined(opts, :solver) || isa(opts.solver, Nothing)
        solver_specified = false
        if MUMPS_available
            opts.solver = "MUMPS"
        else
            opts.solver = "JULIA"
        end
    else
        opts.solver = uppercase(opts.solver)
        if ~(opts.solver in ["MUMPS", "JULIA"])
            throw(ArgumentError("opts.solver = \"$(opts.solver)\" is not a supported option; use \"MUMPS\" or \"JULIA\"."))
        elseif opts.solver == "MUMPS" && ~MUMPS_available
            throw(ArgumentError("opts.solver = \"$(opts.solver)\" but the Julia environment variable \"MUMPS_PREFIX\" is not specified."))
        end
    end

    method_specified = true;
    if ~isdefined(opts, :method) || isa(opts.method, Nothing)
        method_specified = false;        
        # By default, if the argument ''output'' is not given or if
        # opts.iterative_refinement = true, then "factorize_and_solve" is used.
        # Otherwise, then "APF" is used when opts.solver = "MUMPS", and
        # "C*inv(U)*inv(L)*B" is used when opts.solver = "JULIA".
        if opts.return_field_profile || (isdefined(opts, :iterative_refinement) && opts.iterative_refinement == true)
            opts.method = "factorize_and_solve"
        elseif isdefined(opts, :store_ordering) && opts.store_ordering == true
            opts.method = "APF"
        elseif isdefined(opts, :analysis_only) && opts.analysis_only == true
            opts.method = "APF"
        else
            if opts.solver == "MUMPS"
                opts.method = "APF"
            else # opts.solver == "JULIA"
                opts.method = "C*inv(U)*inv(L)*B"
            end
        end
    elseif ~((lowercase(opts.method) in ["apf", "fs", "factorize_and_solve", "fg", "c*inv(u)*inv(l)*b"]))
        throw(ArgumentError("opts.method = \"$(opts.method)\" is not a supported option; use \"APF\", \"factorize_and_solve\", or \"C*inv(U)*inv(L)*B\"."))
    elseif opts.return_field_profile && ~((lowercase(opts.method) in ["fs", "factorize_and_solve"]))
        throw(ArgumentError("opts.method = \"$(opts.method)\" cannot be used for field-profile computations; use opts.method = \"factorize_and_solve\" instead."))
    elseif isdefined(opts, :iterative_refinement) && opts.iterative_refinement == true && ~((lowercase(opts.method) in ["fs", "factorize_and_solve"]))
        throw(ArgumentError("opts.method = \"$(opts.method)\" cannot be used when opts.iterative_refinement = true; use opts.method = \"factorize_and_solve\" instead."))   
    elseif lowercase(opts.method) == "fs" || lowercase(opts.method) == "factorize_and_solve"
        opts.method = "factorize_and_solve"  # opts.method = "FS" is short for opts.method = "factorize_and_solve"
    elseif lowercase(opts.method) == "apf"
        opts.method = "APF"
    elseif lowercase(opts.method) == "fg" || lowercase(opts.method) == "c*inv(u)*inv(l)*b"
        opts.method = "C*inv(U)*inv(L)*B"        
    end

    if opts.method == "APF" && opts.solver == "JULIA"
        throw(ArgumentError("opts.method = \"APF\" requires opts.solver = \"MUMPS\"."))
    end
    if opts.method == "C*inv(U)*inv(L)*B" && opts.solver == "MUMPS"
        throw(ArgumentError("opts.method = \"C*inv(U)*inv(L)*B\" requires opts.solver = \"JULIA\"."))
    end
     
    # opts.use_continuous_dispersion, opts.n0, and opts.m0 will be initialized/checked in mesti_build_channels()
    if ~isdefined(opts, :use_continuous_dispersion)
        opts.use_continuous_dispersion = false # Use finite-difference dispersion by default
    end
    if ~isdefined(opts, :n0); opts.n0 = 0; end
    if ~isdefined(opts, :m0); opts.m0 = 0; end
    
    # Use METIS in 3D and use AMD in 2D by default
    if ~isdefined(opts, :use_METIS) && ~use_2D_TM
        opts.use_METIS = true # If this is 3D system, we use METIS ordering by default
    elseif ~isdefined(opts, :use_METIS) && use_2D_TM
        opts.use_METIS = false # If this is 2D system, we use AMD ordering by default        
    elseif ~isa(opts.use_METIS, Bool)
        throw(ArgumentError("opts.use_METIS must be a boolean, if given."))   
    end

    # opts.symmetrize_K will be checked/initialized later


    # The following fields of opts will be checked/initialized in mesti_matrix_solver():
    #    opts.solver
    #    opts.method
    #    opts.verbal_solver
    #    opts.use_single_precision_MUMPS
    #    opts.nrhs
    #    opts.store_ordering
    #    opts.ordering
    #    opts.analysis_only
    #    opts.nthreads_OMP
    #    opts.parallel_dependency_graph    
    #    opts.iterative_refinement
    #    
    #    opts.use_BLR
    #    opts.threshold_BLR
    #    opts.icntl_36
    #    opts.icntl_38

    # Set up the homogeneous-space channels on the two sides
    k0dx = (2*pi/syst.wavelength)*(syst.dx)
    if ~two_sided
        # For convenience, we set syst.epsilon_high so it is not nothing, and so we can still use channels.low below
        syst.epsilon_high = NaN
    end 
    if ~use_2D_TM
        channels = mesti_build_channels(nx_Ex, nx_Ey, xBC, ny_Ex, ny_Ey, yBC, k0dx, syst.epsilon_low, syst.epsilon_high, opts.use_continuous_dispersion, opts.n0, opts.m0)
    else
       channels = mesti_build_channels(ny_Ex, yBC, k0dx, syst.epsilon_low, syst.epsilon_high, opts.use_continuous_dispersion, opts.m0) 
    end

    if opts.clear_memory
        GC.gc()
    end
    
    N_prop_low = channels.low.N_prop
    if two_sided
        N_prop_high = channels.high.N_prop
    end        

    opts.n0 = nothing # this will not be used in mesti()
    opts.m0 = nothing # this will not be used in mesti()
    opts.use_continuous_dispersion = nothing # this will not be used in mesti()
    
    if opts.verbal
        # print basic system info
        @printf("===System size=== \n")
        if ~use_2D_TM
            @printf("nx_Ex = %d, ny_Ex = %d; nz_Ex = %d => %d \n", nx_Ex, ny_Ex, nz_Ex, nz_tot_Ex)       
            @printf("nx_Ey = %d, ny_Ey = %d; nz_Ey = %d => %d \n", nx_Ey, ny_Ey, nz_Ey, nz_tot_Ey)
            @printf("nx_Ez = %d, ny_Ez = %d; nz_Ez = %d => %d \n", nx_Ez, ny_Ez, nz_Ez, nz_tot_Ez)
        else
            @printf("ny_Ex = %d; nz_Ex = %d => %d for Ex(y,z) \n", ny_Ex, nz_Ex, nz_tot_Ex)
        end
        if two_sided
            @printf("[N_prop_low, N_prop_high] = [%d, %d] per polarization\n", N_prop_low, N_prop_high)    
        else
            @printf("one-sided; N_prop_low = %d per polarization\n", N_prop_low)    
        end
        if ~use_2D_TM
            @printf("xBC = %s; yBC = %s; zBC = [%s, %s]", syst.xBC, syst.yBC, str_zBC[1], str_zBC[2])
            if syst.xBC == "Bloch"; @printf(", (kx_B = %.4f)", syst.kx_B); end
        else
            @printf("yBC = %s; zBC = [%s, %s]", syst.yBC, str_zBC[1], str_zBC[2])
        end        
        if syst.yBC == "Bloch"; @printf(", (ky_B = %.4f)", syst.ky_B); end
        @printf("\n")
    end
        
    t1 = time(); timing_init = t1-t0; # Initialization time
        
    ## Part 2.1: Parse the argument ''input''
    if opts.verbal; @printf("Building B,C... "); end
        
    use_ind_in = ~isa(input, wavefront) # whether to use indices of input channels (if not, wavefront coefficients will be used)
    
    if ~use_2D_TM
        if use_ind_in
            # Indices of input channels, in vectors
            ind_in_low_s = Vector{Integer}()
            ind_in_low_p = Vector{Integer}()   
            if two_sided
                ind_in_high_s = Vector{Integer}()
                ind_in_high_p = Vector{Integer}()
            end
        else
            # Coefficients of wavefronts, in matrices     
            v_in_low_s = spzeros(ComplexF64, N_prop_low, 0)
            v_in_low_p = spzeros(ComplexF64, N_prop_low, 0)
            if two_sided
                v_in_high_s = spzeros(ComplexF64, N_prop_high, 0)
                v_in_high_p = spzeros(ComplexF64,N_prop_high, 0)
            end        
        end
    else
        if use_ind_in
            # Indices of input channels, in vectors
            ind_in_low = Vector{Integer}()
            if two_sided
                ind_in_high = Vector{Integer}()
            end
        else
            # Coefficients of wavefronts, in matrices     
            v_in_low = spzeros(ComplexF64, N_prop_low, 0)
            if two_sided
                v_in_high = spzeros(ComplexF64, N_prop_high, 0)
            end        
        end
    end
       
    # Flags of using s/p polarization on low/high side.
    use_low_s = false;  use_low_p = false
    use_high_s = false; use_high_p = false    
    
    if isa(input, channel_type)        
        if ~isdefined(input, :side) 
            throw(ArgumentError("When ''input'' is a channel_type strucutre, it must provide field \"side\"."))
        else
            if ~(input.side in ["low", "high", "both"])
                throw(ArgumentError("input.side = \"$(input.side)\" is not a supported option; use \"low\", \"high\", or \"both\"."))
            end
        end
        if ~use_2D_TM
            if isdefined(input, :polarization) && ~(input.polarization in ["s", "p", "both"])
                throw(ArgumentError("input.polarization, if given, must be \"s\", \"p\", or \"both\"."))
            elseif ~isdefined(input, :polarization)
                # Pick the polarization to use for users
                input.polarization = "both"
            end
        end
        
        # Take all input channels on the low and/or high
        if input.side == "low" || input.side == "both"
            if ~use_2D_TM
                if input.polarization == "s" || input.polarization == "both"
                    use_low_s = true
                    ind_in_low_s = 1:N_prop_low 
                end
                if input.polarization == "p" || input.polarization == "both"
                    use_low_p = true                
                    ind_in_low_p = 1:N_prop_low  
                end
            else
                ind_in_low = 1:N_prop_low
            end
        end
        if input.side == "high" || input.side == "both"
            if two_sided
                if ~use_2D_TM
                    if input.polarization == "s" || input.polarization == "both"
                        use_high_s = true                    
                        ind_in_high_s = 1:N_prop_high
                    end
                    if input.polarization == "p" || input.polarization == "both"
                        use_high_p = true                    
                        ind_in_high_p = 1:N_prop_high  
                    end
                else
                    ind_in_high = 1:N_prop_high
                end
            else
                throw(ArgumentError("input.side = \"$(input.side)\" cannot be used in a one-sided geometry."))
            end
        end
    elseif isa(input, channel_index)
        # Use the user-specified set of input channels on the low and/or high        
        if ~use_2D_TM && ~(isdefined(input, :ind_low_s) || isdefined(input, :ind_low_p) || isdefined(input, :ind_high_s) || isdefined(input, :ind_high_p))
            throw(ArgumentError("When ''input'' is a channel_index strucutre, it must provide field at least one field of \"ind_low_s\", \"ind_low_p\", \"ind_high_s\", and \"ind_high_p\" for 3D systems."))
        elseif use_2D_TM && ~(isdefined(input, :ind_low) || isdefined(input, :ind_high))
            throw(ArgumentError("When ''input'' is a channel_index strucutre, it must provide field at least one field of \"ind_low\", and \"ind_high\" for 2D systems."))
        end
        # Used in 3D systems
        if isdefined(input, :ind_low_s)             
            ind_in_low_s = input.ind_low_s
            if ~(minimum(ind_in_low_s) > 0 && maximum(ind_in_low_s) <= N_prop_low)
                throw(ArgumentError("input.ind_low_s, when specified, must be an array of positive integers not exceeding N_prop_low = $(N_prop_low)."))
            end
            use_low_s = true       
        end
        if isdefined(input, :ind_low_p) 
            ind_in_low_p = input.ind_low_p
            if ~(minimum(ind_in_low_p) > 0 && maximum(ind_in_low_p) <= N_prop_low)
                throw(ArgumentError("input.ind_low_p, when specified, must be an array of positive integers not exceeding N_prop_low = $(N_prop_low)."))
            end
            use_low_p = true                   
        end
        if isdefined(input, :ind_high_s)
            if two_sided
                ind_in_high_s = input.ind_high_s  
                if ~(minimum(ind_in_high_s) > 0 && maximum(ind_in_high_s) <= N_prop_high)
                    throw(ArgumentError("input.ind_high_s, when specified, must be an array of positive integers not exceeding N_prop_high = $(N_prop_high)."))
                end
                use_high_s = true                   
            else
                throw(ArgumentError("input.ind_high_s cannot be used in a one-sided geometry."))
            end
        end
        if isdefined(input, :ind_high_p)
            if two_sided
                ind_in_high_p = input.ind_high_p        
                if ~(minimum(ind_in_high_p) > 0 && maximum(ind_in_high_p) <= N_prop_high)
                    throw(ArgumentError("input.ind_high_p, when specified, must be an array of positive integers not exceeding N_prop_high = $(N_prop_high)."))
                end
                use_high_p = true                                   
            else
                throw(ArgumentError("input.ind_high_p cannot be used in a one-sided geometry."))
            end
        end
        # Used in 2D systems
        if isdefined(input, :ind_low)             
            ind_in_low = input.ind_low
            if ~(minimum(ind_in_low) > 0 && maximum(ind_in_low) <= N_prop_low)
                throw(ArgumentError("input.ind_low, when specified, must be an array of positive integers not exceeding N_prop_low = $(N_prop_low)."))
            end   
        end
        if isdefined(input, :ind_high)
            if two_sided
                ind_in_high = input.ind_high  
                if ~(minimum(ind_in_high) > 0 && maximum(ind_in_high) <= N_prop_high)
                    throw(ArgumentError("input.ind_high, when specified, must be an array of positive integers not exceeding N_prop_high = $(N_prop_high)."))
                end               
            else
                throw(ArgumentError("input.ind_high cannot be used in a one-sided geometry."))
            end
        end
    else # Use user-specified input wavefronts
        if ~use_2D_TM && ~(isdefined(input, :v_low_s) || isdefined(input, :v_low_p) || isdefined(input, :v_high_s) || isdefined(input, :v_high_p))
            throw(ArgumentError("When ''input'' is a wavefront strucutre, it must provide field at least one field of \"v_low_s\", \"v_low_p\", \"v_high_s\", and \"v_high_p\"."))
        elseif use_2D_TM && ~(isdefined(input, :v_low) || isdefined(input, :v_high))
            throw(ArgumentError("When ''input'' is a wavefront strucutre, it must provide field at least one field of \"v_low\", and \"v_high\"."))
        end
        # Used in 3D systems
        if isdefined(input, :v_low_s)
            use_low_s = true            
            v_in_low_s = input.v_low_s
            if ~(size(v_in_low_s, 1) == N_prop_low)
                throw(ArgumentError("input.v_low_s, when specified, must be a numeric array with size (N_prop_low, M_in_low_s) where N_prop_low = $(N_prop_low)."))
            end
        end
        if isdefined(input, :v_low_p)
            use_low_p = true                        
            v_in_low_p = input.v_low_p
            if ~(size(v_in_low_p, 1) == N_prop_low)
                throw(ArgumentError("input.v_low_p, when specified, must be a numeric array with size (N_prop_low, M_in_low_p) where N_prop_low = $(N_prop_low)."))
            end
        end           
        if isdefined(input, :v_high_s)
            if two_sided
                use_high_s = true                        
                v_in_high_s = input.v_high_s
                if ~(size(v_in_high_s, 1) == N_prop_high)
                    throw(ArgumentError("input.v_high_s, when specified, must be a numeric matrix with size (N_prop_high, M_in_high_s) where N_prop_high = $(N_prop_high)."))
                end
            else
                throw(ArgumentError("input.v_high_s cannot be used in a one-sided geometry."))
            end
        end
        if isdefined(input, :v_high_p)
            if two_sided
                use_high_p = true                                        
                v_in_high_p = input.v_high_p
                if ~(size(v_in_high_p, 1) == N_prop_high)
                    throw(ArgumentError("input.v_high_p, when specified, must be a numeric matrix with size (N_prop_high, M_in_high_p) where N_prop_high = $(N_prop_high)."))
                end
            else
                throw(ArgumentError("input.v_high_p cannot be used in a one-sided geometry."))
            end
        end
        # Used in 2D systems
        if isdefined(input, :v_low)          
            v_in_low = input.v_low
            if ~(size(v_in_low, 1) == N_prop_low)
                throw(ArgumentError("input.v_low, when specified, must be a numeric array with size (N_prop_low, M_in_low) where N_prop_low = $(N_prop_low)."))
            end
        end
        if isdefined(input, :v_high)
            if two_sided                                     
                v_in_high = input.v_high
                if ~(size(v_in_high, 1) == N_prop_high)
                    throw(ArgumentError("input.v_high, when specified, must be a numeric matrix with size (N_prop_high, M_in_high) where N_prop_high = $(N_prop_high)."))
                end
            else
                throw(ArgumentError("input.v_high cannot be used in a one-sided geometry."))
            end
        end
    end

    # Number of inputs
    if ~use_2D_TM
        if use_ind_in
            M_in_low_s = length(ind_in_low_s)
            M_in_low_p = length(ind_in_low_p)        
            M_in_low = M_in_low_s + M_in_low_p
            if two_sided
                M_in_high_s = length(ind_in_high_s)
                M_in_high_p = length(ind_in_high_p)   
                M_in_high = M_in_high_s + M_in_high_p
            end
        else
            M_in_low_s = size(v_in_low_s, 2)
            M_in_low_p = size(v_in_low_p, 2)
            M_in_low = M_in_low_s + M_in_low_p                
            if two_sided
                M_in_high_s = size(v_in_high_s, 2)
                M_in_high_p = size(v_in_high_p, 2)
                M_in_high = M_in_high_s + M_in_high_p
            end              
        end
    else
        if use_ind_in
            M_in_low = length(ind_in_low)     
            if two_sided
                M_in_high = length(ind_in_high)
            end
        else
            M_in_low = size(v_in_low, 2)            
            if two_sided
                M_in_high = size(v_in_high, 2)
            end              
        end
    end
        
    ## Part 2.2: Parse the argument ''output''
    if output != nothing
        use_ind_out = ~isa(output, wavefront) # whether to use indices of output channels (if not, wavefront coefficients will be used)
        
        if ~use_2D_TM
            if use_ind_out
                # Indices of output channels, in row vectors
                ind_out_low_s = Vector{Integer}()
                ind_out_low_p = Vector{Integer}()  
                if two_sided
                    ind_out_high_s = Vector{Integer}()
                    ind_out_high_p = Vector{Integer}()
                end
            else
                # Coefficients of wavefronts, in matrices     
                v_out_low_s = spzeros(ComplexF64, N_prop_low, 0)
                v_out_low_p = spzeros(ComplexF64, N_prop_low, 0)
                if two_sided
                    v_out_high_s = spzeros(ComplexF64, N_prop_high, 0)
                    v_out_high_p = spzeros(ComplexF64, N_prop_high, 0)
                end        
            end
        else
            if use_ind_out
                # Indices of output channels, in row vectors
                ind_out_low = Vector{Integer}()
                if two_sided
                    ind_out_high = Vector{Integer}()
                end
            else
                # Coefficients of wavefronts, in matrices     
                v_out_low = spzeros(ComplexF64, N_prop_low, 0)
                if two_sided
                    v_out_high = spzeros(ComplexF64, N_prop_high, 0)
                end        
            end
        end
        
        if isa(output, channel_type)        
            if ~isdefined(output, :side) 
                throw(ArgumentError("When ''output'' is a channel_type strucutre, it must provide field \"side\"."))
            else
                if ~(output.side in ["low", "high", "both"])
                    throw(ArgumentError("output.side = \"$(output.side)\" is not a supported option; use \"low\", \"high\", or \"both\"."))
                end
            end
            if ~use_2D_TM
                if isdefined(output, :polarization) && ~(output.polarization in ["s", "p", "both"])
                    throw(ArgumentError("output.polarization, if given, must be \"s\", \"p\", or \"both\"."))
                elseif ~isdefined(output, :polarization)
                    # Pick the polarization to use when uers specify the side;
                    output.polarization = "both"
                end
            end

            # Take all output channels on the low and/or high
            if output.side == "low" || output.side == "both"
                if ~use_2D_TM
                    if output.polarization == "s" || output.polarization == "both"
                        use_low_s = true
                        ind_out_low_s = 1:N_prop_low
                    end
                    if output.polarization == "p" || output.polarization == "both"
                        use_low_p = true                
                        ind_out_low_p = 1:N_prop_low
                    end
                else
                    ind_out_low = 1:N_prop_low
                end
            end
            if output.side == "high" || output.side == "both"
                if two_sided
                    if ~use_2D_TM
                        if output.polarization == "s" || output.polarization == "both"
                            use_high_s = true                    
                            ind_out_high_s = 1:N_prop_high
                        end
                        if output.polarization == "p" || output.polarization == "both"
                            use_high_p = true                    
                            ind_out_high_p = 1:N_prop_high 
                        end
                    else
                        ind_out_high = 1:N_prop_high
                    end
                else
                    throw(ArgumentError("output.side = \"$(output.side)\" cannot be used in a one-sided geometry."))
                end
            end
        elseif isa(output, channel_index)
            # Use the user-specified set of output channels on the low and/or high        
            if ~use_2D_TM && ~(isdefined(output, :ind_low_s) || isdefined(output, :ind_low_p) || isdefined(output, :ind_high_s) || isdefined(output, :ind_high_p))
                throw(ArgumentError("When ''output'' is a channel_index strucutre, it must provide field at least one field of \"ind_low_s\", \"ind_low_p\", \"ind_high_s\", and \"ind_high_p\" for 3D systems."))
            elseif use_2D_TM && ~(isdefined(output, :ind_low) || isdefined(output, :ind_high))
                throw(ArgumentError("When ''output'' is a channel_index strucutre, it must provide field at least one field of \"ind_low\", and \"ind_high\" for 2D systems."))
            end
            # Used in 3D systems
            if isdefined(output, :ind_low_s)             
                ind_out_low_s = output.ind_low_s
                if ~(minimum(ind_out_low_s) > 0 && maximum(ind_out_low_s) <= N_prop_low)
                    throw(ArgumentError("output.ind_low_s, when specified, must be an array of positive integers not exceeding N_prop_low = $(N_prop_low)."))
                end
                use_low_s = true       
            end
            if isdefined(output, :ind_low_p) 
                ind_out_low_p = output.ind_low_p
                if ~(minimum(ind_out_low_p) > 0 && maximum(ind_out_low_p) <= N_prop_low)
                    throw(ArgumentError("output.ind_low_p, when specified, must be an array of positive integers not exceeding N_prop_low = $(N_prop_low)."))
                end
                use_low_p = true                   
            end
            if isdefined(output, :ind_high_s)
                if two_sided
                    ind_out_high_s = output.ind_high_s  
                    if ~(minimum(ind_out_high_s) > 0 && maximum(ind_out_high_s) <= N_prop_high)
                        throw(ArgumentError("output.ind_high_s, when specified, must be an array of positive integers not exceeding N_prop_high = $(N_prop_high)."))
                    end
                    use_high_s = true                   
                else
                    throw(ArgumentError("output.ind_high_s cannot be used in a one-sided geometry."))
                end
            end
            if isdefined(output, :ind_high_p)
                if two_sided
                    ind_out_high_p = output.ind_high_p        
                    if ~(minimum(ind_out_high_p) > 0 && maximum(ind_out_high_p) <= N_prop_high)
                        throw(ArgumentError("output.ind_high_p, when specified, must be an array of positive integers not exceeding N_prop_high = $(N_prop_high)."))
                    end
                    use_high_p = true                                   
                else
                    throw(ArgumentError("output.ind_high_p cannot be used in a one-sided geometry."))
                end
            end
            # Used in 2D systems
            if isdefined(output, :ind_low)             
                ind_out_low = output.ind_low
                if ~(minimum(ind_out_low) > 0 && maximum(ind_out_low) <= N_prop_low)
                    throw(ArgumentError("output.ind_low, when specified, must be an array of positive integers not exceeding N_prop_low = $(N_prop_low)."))
                end 
            end
            if isdefined(output, :ind_high)
                if two_sided
                    ind_out_high = output.ind_high 
                    if ~(minimum(ind_out_high) > 0 && maximum(ind_out_high) <= N_prop_high)
                        throw(ArgumentError("output.ind_high, when specified, must be an array of positive integers not exceeding N_prop_high = $(N_prop_high)."))
                    end                 
                else
                    throw(ArgumentError("output.ind_high cannot be used in a one-sided geometry."))
                end
            end
        else # Use user-specified output wavefronts
            if ~use_2D_TM && ~(isdefined(output, :v_low_s) || isdefined(output, :v_low_p) || isdefined(output, :v_high_s) || isdefined(output, :v_high_p))
                throw(ArgumentError("When ''output'' is a wavefront strucutre, it must provide field at least one field of \"v_low_s\", \"v_low_p\", \"v_high_s\", and \"v_high_p\"."))
            elseif use_2D_TM && ~(isdefined(output, :v_low) || isdefined(output, :v_high))
                throw(ArgumentError("When ''output'' is a wavefront strucutre, it must provide field at least one field of \"v_low\", and \"v_high\"."))
            end              
            
            # Used in 3D systems
            if isdefined(output, :v_low_s)
                v_out_low_s = output.v_low_s
                if ~(size(v_out_low_s, 1) == N_prop_low)
                    throw(ArgumentError("output.v_low_s, when specified, must be a numeric array with size (N_prop_low, M_out_low_s) where N_prop_low = $(N_prop_low)."))
                end
            end
            if isdefined(output, :v_low_p)
                v_out_low_p = output.v_low_p
                if ~(size(v_out_low_p, 1) == N_prop_low)
                    throw(ArgumentError("output.v_low_p, when specified, must be a numeric array with size (N_prop_low, M_out_low_p) where N_prop_low = $(N_prop_low)."))
                end
            end           
            if isdefined(output, :v_high_s)
                if two_sided
                    v_out_high_s = output.v_high_s
                    if ~(size(v_out_high_s, 1) == N_prop_high)
                        throw(ArgumentError("output.v_high_s, when specified, must be a numeric matrix with size (N_prop_high, M_out_high_s) where N_prop_high = $(N_prop_high)."))
                    end
                else
                    throw(ArgumentError("output.v_high_s cannot be used in a one-sided geometry."))
                end
            end
            if isdefined(output, :v_high_p)
                if two_sided
                    v_out_high_p = output.v_high_p
                    if ~(size(v_out_high_p, 1) == N_prop_high)
                        throw(ArgumentError("output.v_high_p, when specified, must be a numeric matrix with size (N_prop_high, M_out_high_p) where N_prop_high = $(N_prop_high)."))
                    end
                else
                    throw(ArgumentError("output.v_high_p cannot be used in a one-sided geometry."))
                end
            end
            # Used in 2D systems
            if isdefined(output, :v_low)
                v_out_low = output.v_low
                if ~(size(v_out_low, 1) == N_prop_low)
                    throw(ArgumentError("output.v_low, when specified, must be a numeric array with size (N_prop_low, M_out_low) where N_prop_low = $(N_prop_low)."))
                end
            end
            if isdefined(output, :v_high)
                if two_sided
                    v_out_high = output.v_high
                    if ~(size(v_out_high, 1) == N_prop_high)
                        throw(ArgumentError("output.v_high, when specified, must be a numeric matrix with size (N_prop_high, M_out_high) where N_prop_high = $(N_prop_high)."))
                    end
                else
                    throw(ArgumentError("output.v_high cannot be used in a one-sided geometry."))
                end
            end
        end

        # Number of outputs
        if ~use_2D_TM
            if use_ind_out
                M_out_low_s = length(ind_out_low_s)
                M_out_low_p = length(ind_out_low_p)        
                M_out_low = M_out_low_s + M_out_low_p
                if two_sided
                    M_out_high_s = length(ind_out_high_s)
                    M_out_high_p = length(ind_out_high_p)   
                    M_out_high = M_out_high_s + M_out_high_p
                end
            else
                M_out_low_s = size(v_out_low_s, 2)
                M_out_low_p = size(v_out_low_p, 2)
                M_out_low = M_out_low_s + M_out_low_p                
                if two_sided
                    M_out_high_s = size(v_out_high_s, 2)
                    M_out_high_p = size(v_out_high_p, 2)
                    M_out_high = M_out_high_s + M_out_high_p
                end              
            end
        else
            if use_ind_out
                M_out_low = length(ind_out_low)
                if two_sided
                    M_out_high = length(ind_out_high)
                end
            else
                M_out_low = size(v_out_low, 2)          
                if two_sided
                    M_out_high = size(v_out_high, 2)
                end              
            end
        end
    end

    ## Part 2.3: Build the source B_Ej_low/B_Ej_high and the projection C_Ej_B/C_Ej_T on the two surfaces
        
    # Incorporate line 1383 to line 1463 in mesti2s.m

    # For 2D TM, we symmetrize the set of input channels and output channels, if possible
    # Bloch periodic boundary condition with ky_B*periodicity != 0 or pi breaks the symmetry of A
    # mesti_build_channels() currently does not return channels.low.ind_prop_conj when yBC == pi (even though such permutation does exist), so we cannot symmetrize K when yBC == pi.
    # When ny == 1, even though matrix A is symmetric, there is still no permutation that flips the sign of ky, so we still set is_symmetric_A to false.
    if !use_2D_TM || (isa(yBC, Number) && yBC != 0) # && yBC ~= pi && ny > 1
        is_symmetric_A = false
    else
        is_symmetric_A = true
    end
    # opts.symmetrize_K can only be used when all of the following are met: (1) input argument ''output'' is given, (2) ''input'' and ''output'' are not specified as wavefronts, (3) opts.method = "APF", and (4) the boundary condition in y is not Bloch periodic.
    if !isa(output, Nothing) && use_ind_in && use_ind_out && opts.method == "APF" && is_symmetric_A                
        # By default, opts.symmetrize_K = true if it's possible            
        if !isdefined(opts, :symmetrize_K) || isa(opts.symmetrize_K, Nothing)
            opts.symmetrize_K = true
        elseif !(isa(opts.clear_syst, Bool))
            throw(ArgumentError("opts.symmetrize_K must be a logical scalar, if given."))
        end
        use_transpose_B = opts.symmetrize_K
        opts.symmetrize_K = nothing
    else
        # symmetrize K is not possible
        if isdefined(opts, :symmetrize_K)
            if isequal(opts.symmetrize_K, true)
                @warn("opts.symmetrize_K = true is only available when use_2D_TM = true, isa(output, Nothing) = false, use_ind_in = true, use_ind_out = true, opts.solver = \"MUMPS\", opts.method = \"APF\", and syst.yBC is not \"Bloch\". Here use_2D_TM = $(use_2D_TM), isa(output, Nothing) = $(isa(output, Nothing)), use_ind_in = $(use_ind_in), use_ind_out = $(use_ind_out), opts.solver = \"$(opts.solver)\", opts.method = \"$(opts.method)\", syst.yBC = \"$(syst.yBC)\"; opts.symmetrize_K will be ignored")
            end
            opts.symmetrize_K = nothing
        end
        use_transpose_B = false
    end

    # Remove opts.solver and opts.method so that mesti_matrix_solver() will know they were not specified by the user.
    if !solver_specified
        opts.solver = nothing
    end
    if !method_specified
        opts.method = nothing
    end
    
    # No need to build C if we symmetrize K
    if opts.return_field_profile || use_transpose_B
        build_C = false;
    else
        build_C = true;
    end    
    
    # We only need the transverse functions for Ex, Ey, and derivative of the transverse functions for Ez; see Ref XXX ......
    if ~use_2D_TM
        @cast u_prop_low_Ex[(n,m),N_prop_low] := channels.u_x_m(channels.low.kydx_prop)[m,N_prop_low] * channels.u_x_n(channels.low.kxdx_prop)[n,N_prop_low]
        @cast u_prop_low_Ey[(n,m),N_prop_low] := channels.u_y_m(channels.low.kydx_prop)[m,N_prop_low] * channels.u_y_n(channels.low.kxdx_prop)[n,N_prop_low]
        @cast u_prop_low_dEz_over_dx[(n,m),N_prop_low] := channels.u_z_m(channels.low.kydx_prop)[m,N_prop_low] * channels.du_z_n(channels.low.kxdx_prop)[n,N_prop_low]    
        @cast u_prop_low_dEz_over_dy[(n,m),N_prop_low] := channels.du_z_m(channels.low.kydx_prop)[m,N_prop_low] * channels.u_z_n(channels.low.kxdx_prop)[n,N_prop_low]
        if two_sided
            if (syst.epsilon_high == syst.epsilon_low)
                u_prop_high_Ex = u_prop_low_Ex
                u_prop_high_Ey = u_prop_low_Ey
                u_prop_high_dEz_over_dx = u_prop_low_dEz_over_dx
                u_prop_high_dEz_over_dy = u_prop_low_dEz_over_dy
            else
                @cast u_prop_high_Ex[(n,m),N_prop_high] := channels.u_x_m(channels.high.kydx_prop)[m,N_prop_high] * channels.u_x_n(channels.high.kxdx_prop)[n,N_prop_high]
                @cast u_prop_high_Ey[(n,m),N_prop_high] := channels.u_y_m(channels.high.kydx_prop)[m,N_prop_high] * channels.u_y_n(channels.high.kxdx_prop)[n,N_prop_high]
                @cast u_prop_high_dEz_over_dx[(n,m),N_prop_high] := channels.u_z_m(channels.high.kydx_prop)[m,N_prop_high] * channels.du_z_n(channels.high.kxdx_prop)[n,N_prop_high]
                @cast u_prop_high_dEz_over_dy[(n,m),N_prop_high] := channels.du_z_m(channels.high.kydx_prop)[m,N_prop_high] * channels.u_z_n(channels.high.kxdx_prop)[n,N_prop_high]
            end
        end
    else
        @cast u_prop_low_Ex[m,N_prop_low] := channels.u_x_m(channels.low.kydx_prop)[m,N_prop_low]
        if two_sided
            if (syst.epsilon_high == syst.epsilon_low)
                u_prop_high_Ex = u_prop_low_Ex
            else
                @cast u_prop_high_Ex[m,N_prop_high] := channels.u_x_m(channels.high.kydx_prop)[m,N_prop_high]
            end
        end
    end
    # Here we build:
    # (1) the blocks of input matrix B on the low and high surfaces: B_Ej_low and B_Ej_high
    # (2) the blocks of output matrix C on the low and high surfaces: C_Ej_low and C_Ej_high
    # B_Ej_low, C_Ej_low, B_Ej_high, C_Ej_high are all dense matrices with size(..., 1) = nx_Ej*ny_Ej. They are inputs/outputs on a slice surface in xy-plane.
    # Inputs/outputs are placed at one pixel outside syst.epsilon_xx/epsilon_yy/epsilon_zz.
    # A line source including back propagation of B_Ex_low = -2im*sqrt(nu)*alpha_x(a,b,sigma,+)*u_x(n,m,a,b)*exp(1im*kz(a,b)*dx*(-1/2))+2/sqrt(nu)*alpha_z(a,b,sigma,+)*(u_z(n+1,m,a,b)-u_z(n,m,a,b))*cos(kz(a,b)*dx/2)*exp(1im*kz(a,b)*dx*(-1/2)) and B_Ey_low = -2im*sqrt(nu)*alpha_y(a,b,sigma,+)*u_y(n,m,a,b)*exp(1im*kz(a,b)*dx*(-1/2))+2/sqrt(nu)*alpha_z(a,b,sigma,+)*(u_z(n,m+1,a,b)-u_z(n,m,a,b))*cos(kz(a,b)*dx/2)*exp(1im*kz(a,b)*dx*(-1/2)) at l=0 will generate an z-flux-normalized incident field of Ex, Ey, and Ez. Ex = alpha_x(a,b,sigma,+/-)/sqrt(nu)*alpha_x(a,b,sigma,+/-)*u_x(n,m,a,b)*exp(1im*kz(a,b)*dx*(|l|-1/2)), where + for l>= 0 and - for l< 0. Ey = alpha_y(a,b,sigma,+/-)/sqrt(nu)*alpha_y(a,b,sigma,+/-)*u_x(n,m,a,b)*exp(1im*kz(a,b)*dx*(|l|-1/2)), where + for l>= 0 and - for l< 0. Ez = alpha_y(a,b,sigma,+/-)/sqrt(nu)*alpha_y(a,b,sigma,+/-)*u_x(n,m,a,b)*exp(1im*kz(a,b)*dx*|l|), where + for l>= 0 and - for l< 0. Note nu = sin(kzdx).
    # We will multiple the -2i prefactor at the end.
    # The flux-normalized output projection is \sum\limits_{j = x,y}{sqrt(nu)*conj(alpha_j(c,d,sigma,+/-)*u_x(n,m,c,d)*exp((+/-)1im*kz(c,d)*dx*(1/2)))+1/sqrt(nu)*conj(1im*alpha_z(c,d,sigma,+/-)*(u_z(n+delta_ix,m+delta_iy,c,d)-u_z(n,m,c,d))*cos(kz(c,d)*dx/2)))*exp(1im*kz(c,d)*dx*(-1/2)}; it will be transposed in mesti(). We also need to shift the phase by back propagating.
    # Therefore, we want C_Ex_low = sqrt(nu)*conj(alpha_x(c,d,sigma,-)*u_x(n,m,c,d)*exp(-1im*kz(c,d)*dx*(1/2)))+1/sqrt(nu)*conj(1im*alpha_z(c,d,sigma,-)*(u_z(n+1,m,c,d)-u_z(n,m,c,d))*cos(kz(c,d)*dx/2)))*exp(1im*kz(c,d)*dx*(-1/2) and C_Ey_low = sqrt(nu)*conj(alpha_y(c,d,sigma,-)*u_x(n,m,c,d)*exp(-1im*kz(c,d)*dx*(1/2)))+1/sqrt(nu)*conj(1im*alpha_z(c,d,sigma,-)*(u_z(n+1,m,c,d)-u_z(n,m,c,d))*cos(kz(c,d)*dx/2)))*exp(1im*kz(c,d)*dx*(-1/2).
    # Note that the complex conjugation only applies to u; the sqrt(nu) prefactor is not conjugated. (At real-valued frequency, nu is real-valued, so this doesn't matter. But at complex-valued frequency, nu is complex-valued, and we should not conjugate it. Note that the transverse basis is complete and orthonormal even when the frequency is complex, so the output projection doesn't need to be modified when the frequency is complex.)
    # When the input/output channels are specified by channel indices, we will multiply the exp(-i*kzdx(a,b)*dn) and exp(-i*kzdx(c,d)*dn) prefactors at the end, after C*inv(A)*B is computed.
    # When the input/output wavefronts are specified, we take superpositions of the channels using the v coefficients, with the sqrt(nu)*exp(-i*kzdx*dn) prefactors included.
    if ~use_2D_TM
        # Construct coefficient of two polarization
        # alpha_x_low_s stands for x-component coefficients for s-polarization waves.
        # z-component coefficients for s-polarization waves are always zeros, so we don't initialize it.
        alpha_x_low_s = zeros(N_prop_low,0);  alpha_y_low_s = zeros(N_prop_low,0)
        alpha_x_low_p = zeros(N_prop_low,0);  alpha_y_low_p = zeros(N_prop_low,0);  alpha_z_low_p = zeros(N_prop_low,0)
        alpha_x_high_s = zeros(N_prop_high,0); alpha_y_high_s = zeros(N_prop_high,0)
        alpha_x_high_p = zeros(N_prop_high,0); alpha_y_high_p = zeros(N_prop_high,0); alpha_z_high_p = zeros(N_prop_high,0)

        # Here we define kappa and use them to construct alpha, which is the x, y, and z coefficents for polarization basis.
        # Note thes kappa_x, kappa_y, kappa_z are different from the real-coordinate-stretching factor kappa in PML parameters.
        # Please refer to our upcoming 3D vectorial wave paper or email authors if you are interested in the theory.  
        kappa_x_low = sin.((channels.low.kxdx_prop)/2)
        kappa_y_low = sin.((channels.low.kydx_prop)/2)
        kappa_z_low = sin.((channels.low.kzdx_prop)/2)

        if use_low_s
            denominator = sqrt.(kappa_x_low.^2+kappa_y_low.^2)
            alpha_x_low_s = -kappa_y_low./denominator
            alpha_y_low_s = kappa_x_low./denominator
            # For the case, kappa_x_low = 0 and kappa_y_low = 0. We choose s-polarization along y-direction.
            alpha_x_low_s[isnan.(alpha_x_low_s)] .= 0 
            alpha_y_low_s[isnan.(alpha_y_low_s)] .= 1
        end

        # Here we use the convention that the alpha_low is alpha_+ (+ means wave propagation along positive z-direction) on low side and alpha_high is also alpha_+ (+ means wave propagation along positive z-direction) on high side.
        # We should make the sign of alpha_z following by the convention.
        if use_low_p
            denominator = sqrt.((abs.(kappa_x_low.*kappa_z_low)).^2+(abs.(kappa_y_low.*kappa_z_low)).^2+(abs.(kappa_x_low.^2+kappa_y_low.^2)).^2)
            alpha_x_low_p = kappa_x_low.*kappa_z_low./denominator
            alpha_y_low_p = kappa_y_low.*kappa_z_low./denominator
            alpha_z_low_p = (-1)*(kappa_x_low.^2+kappa_y_low.^2)./denominator
            # For the case, kappa_x_low = 0 and kappa_y_low = 0. We choose s-polarization along x-direction.
            alpha_x_low_p[isnan.(alpha_x_low_p)] .= 1
            alpha_y_low_p[isnan.(alpha_y_low_p)] .= 0
            alpha_z_low_p[isnan.(alpha_z_low_p)] .= 0
        end                

        if two_sided
            if syst.epsilon_high == syst.epsilon_low
                kappa_x_high = kappa_x_low
                kappa_y_high = kappa_y_low
                kappa_z_high = kappa_z_low                
            else
                kappa_x_high = sin.((channels.high.kxdx_prop)/2)
                kappa_y_high = sin.((channels.high.kydx_prop)/2)
                kappa_z_high = sin.((channels.high.kzdx_prop)/2)
            end
            if use_high_s
                denominator = sqrt.(kappa_x_high.^2+kappa_y_high.^2)
                alpha_x_high_s = -kappa_y_high./denominator
                alpha_y_high_s = kappa_x_high./denominator
                alpha_x_high_s[isnan.(alpha_x_high_s)] .= 0
                alpha_y_high_s[isnan.(alpha_y_high_s)] .= 1
            end
            if use_high_p
                denominator = sqrt.((abs.(kappa_x_high.*kappa_z_high)).^2+(abs.(kappa_y_high.*kappa_z_high)).^2+(abs.(kappa_x_high.^2+kappa_y_high.^2)).^2)                
                alpha_x_high_p = kappa_x_high.*kappa_z_high./denominator
                alpha_y_high_p = kappa_y_high.*kappa_z_high./denominator
                alpha_z_high_p = (-1)*(kappa_x_high.^2+kappa_y_high.^2)./denominator
                alpha_x_high_p[isnan.(alpha_x_high_p)] .= 1
                alpha_y_high_p[isnan.(alpha_y_high_p)] .= 0
                alpha_z_high_p[isnan.(alpha_z_high_p)] .= 0
            end
        end
    end

    if use_transpose_B # when use_2D_TM && opts.symmetrize_K = true 
        # Here, we pad channels and/or permutate them such that C = transpose(B); this makes matrix K = [A,B;C,0] symmetric.
        # To have C=transpose(B), the complex conjugate of the transverse field profiles of the list of output channels must equal the transverse field profiles of the list of input channels, and the list of input channels and the list of output channels must have the same prefactor nu.
        # So, we expand the list of input channels (ind_in_low) to include the conjugate pairs of the output channels (channels.low.ind_prop_conj(ind_out_low)). The conjugate pairs correspond to flipping the sign of ky, and they share the same nu (which only depends on ky^2).
        # We only build B_low and B_high here, from which matrix B will be built in mesti(). C_low and C_high are not needed since matrix C will not be used until in mesti_matrix_solver().
        # We only need to keep the unique channel indices, since ind_in_low and channels.low.ind_prop_conj(ind_out_low) are likely to contain the same indices.
        # ind_in_out_low satisfies ind_L(ind_in_out_low) = [ind_in_low, channels.L.ind_prop_conj(ind_out_low)]. The computations will be done in ind_low, so later we can use ind_in_out_low to retrieve the original lists of input channels (with the first half of ind_in_out_low) and the original list of output channels (with the second half of ind_in_out_low).

        # The following four lines are equivalent to [ind_low, ~, ind_in_out_low] = unique([ind_in_low, channels.low.ind_prop_conj(ind_out_low)]) in MATLAB
        ind_low = sort(unique(vcat(ind_in_low, channels.low.ind_prop_conj[ind_out_low])))
        ind_in_out_low_dict = Dict(ind_low[ii] => ii for ii in 1:length(ind_low))
        ind_in_out_low = [ind_in_out_low_dict[ii] for ii in vcat(ind_in_low, channels.low.ind_prop_conj[ind_out_low])]
        ind_in_out_low_dict = nothing
        
        B_Ex_low = u_prop_low_Ex[:,ind_low]
        
        if two_sided
            ind_high = sort(unique(vcat(ind_in_high, channels.high.ind_prop_conj[ind_out_high])))
            ind_in_out_high_dict = Dict(ind_high[ii] => ii for ii in 1:length(ind_high))
            ind_in_out_high = [ind_in_out_high_dict[ii] for ii in vcat(ind_in_high, channels.high.ind_prop_conj[ind_out_high])]
            ind_in_out_high_dict = nothing
            
            B_Ex_high = u_prop_high_Ex[:,ind_high]
        end
    else # without opts.symmetrize_K                
        # Build up matrices B_low_sigma_Ej and B_high_sigma_Ej
        # Then combine them into B_Ej_low and B_Ej_high
        if use_ind_in # input channels specified by ind_in_low_s, ind_in_low_p, ind_in_high_s and ind_in_high_p
                      # or input channels specified by ind_in_low and ind_in_high
            if ~use_2D_TM
                B_low_s_Ex = u_prop_low_Ex[:,ind_in_low_s].*reshape(channels.low.sqrt_nu_prop[ind_in_low_s],1,:)
                B_low_s_Ey = u_prop_low_Ey[:,ind_in_low_s].*reshape(channels.low.sqrt_nu_prop[ind_in_low_s],1,:)
                B_low_p_Ex = u_prop_low_Ex[:,ind_in_low_p].*reshape(channels.low.sqrt_nu_prop[ind_in_low_p],1,:)
                B_low_p_Ey = u_prop_low_Ey[:,ind_in_low_p].*reshape(channels.low.sqrt_nu_prop[ind_in_low_p],1,:)
                B_low_p_dEz_over_dx = u_prop_low_dEz_over_dx[:,ind_in_low_p].*reshape(1 ./channels.low.sqrt_nu_prop[ind_in_low_p].*cos.(channels.low.kzdx_prop[ind_in_low_p]/2),1,:)
                B_low_p_dEz_over_dy = u_prop_low_dEz_over_dy[:,ind_in_low_p].*reshape(1 ./channels.low.sqrt_nu_prop[ind_in_low_p].*cos.(channels.low.kzdx_prop[ind_in_low_p]/2),1,:)            
                B_Ex_low = reshape([B_low_s_Ex.*reshape(alpha_x_low_s[ind_in_low_s],1,:) B_low_p_Ex.*reshape(alpha_x_low_p[ind_in_low_p],1,:)+1im*B_low_p_dEz_over_dx.*reshape(alpha_z_low_p[ind_in_low_p],1,:)], nx_Ex, ny_Ex, 1, :)
                B_Ey_low = reshape([B_low_s_Ey.*reshape(alpha_y_low_s[ind_in_low_s],1,:) B_low_p_Ey.*reshape(alpha_y_low_p[ind_in_low_p],1,:)+1im*B_low_p_dEz_over_dy.*reshape(alpha_z_low_p[ind_in_low_p],1,:)], nx_Ey, ny_Ey, 1, :)
                if two_sided
                    B_high_s_Ex = u_prop_high_Ex[:,ind_in_high_s].*reshape(channels.high.sqrt_nu_prop[ind_in_high_s],1,:)
                    B_high_s_Ey = u_prop_high_Ey[:,ind_in_high_s].*reshape(channels.high.sqrt_nu_prop[ind_in_high_s],1,:)
                    B_high_p_Ex = u_prop_high_Ex[:,ind_in_high_p].*reshape(channels.high.sqrt_nu_prop[ind_in_high_p],1,:)
                    B_high_p_Ey = u_prop_high_Ey[:,ind_in_high_p].*reshape(channels.high.sqrt_nu_prop[ind_in_high_p],1,:)
                    B_high_p_dEz_over_dx = u_prop_high_dEz_over_dx[:,ind_in_high_p].*reshape(1 ./channels.high.sqrt_nu_prop[ind_in_high_p].*cos.(channels.high.kzdx_prop[ind_in_high_p]/2),1,:)
                    B_high_p_dEz_over_dy = u_prop_high_dEz_over_dy[:,ind_in_high_p].*reshape(1 ./channels.high.sqrt_nu_prop[ind_in_high_p].*cos.(channels.high.kzdx_prop[ind_in_high_p]/2),1,:)
                    B_Ex_high = reshape([B_high_s_Ex.*reshape(alpha_x_high_s[ind_in_high_s],1,:) B_high_p_Ex.*reshape(alpha_x_high_p[ind_in_high_p],1,:)+1im*B_high_p_dEz_over_dx.*reshape(alpha_z_high_p[ind_in_high_p],1,:)], nx_Ex, ny_Ex, 1, :)
                    B_Ey_high = reshape([B_high_s_Ey.*reshape(alpha_y_high_s[ind_in_high_s],1,:) B_high_p_Ey.*reshape(alpha_y_high_p[ind_in_high_p],1,:)+1im*B_high_p_dEz_over_dy.*reshape(alpha_z_high_p[ind_in_high_p],1,:)], nx_Ey, ny_Ey, 1, :)      
                end
            else
                B_Ex_low = u_prop_low_Ex[:,ind_in_low]
                if two_sided
                   B_Ex_high = u_prop_high_Ex[:,ind_in_high]
                end
            end
        else # input wavefronts specified by v_in_low_s, v_in_low_p, v_in_high_s, and v_in_high_p.
             # or input wavefronts specified by v_in_low and v_in_high
            if ~use_2D_TM
                B_low_s_Ex = u_prop_low_Ex*((channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_in_low_s.*alpha_x_low_s)) # use implicit expansion
                B_low_s_Ey = u_prop_low_Ey*((channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_in_low_s.*alpha_y_low_s))
                B_low_p_Ex = u_prop_low_Ex*((channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_in_low_p.*alpha_x_low_p))
                B_low_p_Ey = u_prop_low_Ey*((channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_in_low_p.*alpha_y_low_p))
                B_low_p_dEz_over_dx = u_prop_low_dEz_over_dx*((1 ./channels.low.sqrt_nu_prop.*cos.(channels.low.kzdx_prop/2).*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_in_low_p.*alpha_z_low_p))
                B_low_p_dEz_over_dy = u_prop_low_dEz_over_dy*((1 ./channels.low.sqrt_nu_prop.*cos.(channels.low.kzdx_prop/2).*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_in_low_p.*alpha_z_low_p))
                B_Ex_low = reshape([B_low_s_Ex B_low_p_Ex+1im*B_low_p_dEz_over_dx], nx_Ex, ny_Ex, 1, :)
                B_Ey_low = reshape([B_low_s_Ey B_low_p_Ey+1im*B_low_p_dEz_over_dy], nx_Ey, ny_Ey, 1, :) 
                if two_sided
                    B_high_s_Ex = u_prop_high_Ex*((channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_in_high_s.*alpha_x_high_s)) # use implicit expansion
                    B_high_s_Ey = u_prop_high_Ey*((channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_in_high_s.*alpha_y_high_s))
                    B_high_p_Ex = u_prop_high_Ex*((channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_in_high_p.*alpha_x_high_p))
                    B_high_p_Ey = u_prop_high_Ey*((channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_in_high_p.*alpha_y_high_p))
                    B_high_p_dEz_over_dx = u_prop_high_dEz_over_dx*((1 ./channels.high.sqrt_nu_prop.*cos.(channels.high.kzdx_prop/2).*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_in_high_p.*alpha_z_high_p))
                    B_high_p_dEz_over_dy = u_prop_high_dEz_over_dy*((1 ./channels.high.sqrt_nu_prop.*cos.(channels.high.kzdx_prop/2).*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_in_high_p.*alpha_z_high_p))
                    B_Ex_high = reshape([B_high_s_Ex B_high_p_Ex+1im*B_high_p_dEz_over_dx], nx_Ex, ny_Ex, 1, :)
                    B_Ey_high = reshape([B_high_s_Ey B_high_p_Ey+1im*B_high_p_dEz_over_dy], nx_Ey, ny_Ey, 1, :)
                end
            else
                B_Ex_low = u_prop_low_Ex*(channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop).*v_in_low)
                if two_sided
                    B_Ex_high = u_prop_high_Ex*(channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop).*v_in_high)
                end
            end
        end        
        # Build up component of matrices C_low_sigma_Ej and C_high_sigma_Ej
        # Then combine them into C_Ej_low and C_Ej_high        
        if build_C
            if use_ind_out # output channels specified by ind_out_low_s, ind_out_low_p, ind_out_high_s and ind_out_high_p
                           # or input channels specified by ind_out_low and ind_out_high
                if ~use_2D_TM
                    C_low_s_Ex = conj(u_prop_low_Ex[:,ind_out_low_s]).*reshape(channels.low.sqrt_nu_prop[ind_out_low_s],1,:)
                    C_low_s_Ey = conj(u_prop_low_Ey[:,ind_out_low_s]).*reshape(channels.low.sqrt_nu_prop[ind_out_low_s],1,:)
                    C_low_p_Ex = conj(u_prop_low_Ex[:,ind_out_low_p]).*reshape(channels.low.sqrt_nu_prop[ind_out_low_p],1,:)
                    C_low_p_Ey = conj(u_prop_low_Ey[:,ind_out_low_p]).*reshape(channels.low.sqrt_nu_prop[ind_out_low_p],1,:)
                    C_low_p_dEz_over_dx = conj(u_prop_low_dEz_over_dx[:,ind_out_low_p]).*reshape(1 ./channels.low.sqrt_nu_prop[ind_out_low_p].*cos.(channels.low.kzdx_prop[ind_out_low_p]/2),1,:)
                    C_low_p_dEz_over_dy = conj(u_prop_low_dEz_over_dy[:,ind_out_low_p]).*reshape(1 ./channels.low.sqrt_nu_prop[ind_out_low_p].*cos.(channels.low.kzdx_prop[ind_out_low_p]/2),1,:)
                    C_Ex_low = reshape([C_low_s_Ex.*conj(reshape(alpha_x_low_s[ind_out_low_s],1,:)) C_low_p_Ex.*conj(reshape(alpha_x_low_p[ind_out_low_p],1,:))-1im*C_low_p_dEz_over_dx.*conj(reshape(alpha_z_low_p[ind_out_low_p],1,:))], nx_Ex, ny_Ex, 1, :)
                    C_Ey_low = reshape([C_low_s_Ey.*conj(reshape(alpha_y_low_s[ind_out_low_s],1,:)) C_low_p_Ey.*conj(reshape(alpha_y_low_p[ind_out_low_p],1,:))-1im*C_low_p_dEz_over_dy.*conj(reshape(alpha_z_low_p[ind_out_low_p],1,:))], nx_Ey, ny_Ey, 1, :)
                    if two_sided
                        C_high_s_Ex = conj(u_prop_high_Ex[:,ind_out_high_s]).*reshape(channels.high.sqrt_nu_prop[ind_out_high_s],1,:)
                        C_high_s_Ey = conj(u_prop_high_Ey[:,ind_out_high_s]).*reshape(channels.high.sqrt_nu_prop[ind_out_high_s],1,:)
                        C_high_p_Ex = conj(u_prop_high_Ex[:,ind_out_high_p]).*reshape(channels.high.sqrt_nu_prop[ind_out_high_p],1,:)
                        C_high_p_Ey = conj(u_prop_high_Ey[:,ind_out_high_p]).*reshape(channels.high.sqrt_nu_prop[ind_out_high_p],1,:)
                        C_high_p_dEz_over_dx = conj(u_prop_high_dEz_over_dx[:,ind_out_high_p]).*reshape(1 ./channels.high.sqrt_nu_prop[ind_out_high_p].*cos.(channels.high.kzdx_prop[ind_out_high_p]/2),1,:)
                        C_high_p_dEz_over_dy = conj(u_prop_high_dEz_over_dy[:,ind_out_high_p]).*reshape(1 ./channels.high.sqrt_nu_prop[ind_out_high_p].*cos.(channels.high.kzdx_prop[ind_out_high_p]/2),1,:)
                        C_Ex_high = reshape([C_high_s_Ex.*conj(reshape(alpha_x_high_s[ind_out_high_s],1,:)) C_high_p_Ex.*conj(reshape(alpha_x_high_p[ind_out_high_p],1,:))-1im*C_high_p_dEz_over_dx.*conj(reshape(alpha_z_high_p[ind_out_high_p],1,:))], nx_Ex, ny_Ex, 1, :) 
                        C_Ey_high = reshape([C_high_s_Ey.*conj(reshape(alpha_y_high_s[ind_out_high_s],1,:)) C_high_p_Ey.*conj(reshape(alpha_y_high_p[ind_out_high_p],1,:))-1im*C_high_p_dEz_over_dy.*conj(reshape(alpha_z_high_p[ind_out_high_p],1,:))], nx_Ey, ny_Ey, 1, :)
                    end
              else
                    C_Ex_low = conj(u_prop_low_Ex[:,ind_out_low])
                    if two_sided
                       C_Ex_high = conj(u_prop_high_Ex[:,ind_out_high])
                    end
              end
            else # output wavefronts specified by v_out_low_s, v_out_low_p, v_out_high_s, and v_out_high_p
                 # u_prop_low_Ex is a large nx_Ex*ny_Ex-by-N_prop_L matrix; instead of conjugating it, we conjugate the smaller nx_Ex*ny_Ey-by-M_out_L matrix.
                if ~use_2D_TM
                    C_low_s_Ex = conj(u_prop_low_Ex*(conj(channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_out_low_s.*alpha_x_low_s))) # use implicit expansion
                    C_low_s_Ey = conj(u_prop_low_Ey*(conj(channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_out_low_s.*alpha_y_low_s)))
                    C_low_p_Ex = conj(u_prop_low_Ex*(conj(channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_out_low_p.*alpha_x_low_p)))                
                    C_low_p_Ey = conj(u_prop_low_Ey*(conj(channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_out_low_p.*alpha_y_low_p)))
                    C_low_p_dEz_over_dx = conj(u_prop_low_dEz_over_dx*(conj(1 ./channels.low.sqrt_nu_prop.*cos.(channels.low.kzdx_prop/2).*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_out_low_p.*alpha_z_low_p)))
                    C_low_p_dEz_over_dy = conj(u_prop_low_dEz_over_dy*(conj(1 ./channels.low.sqrt_nu_prop.*cos.(channels.low.kzdx_prop/2).*exp.((-1im*dn)*channels.low.kzdx_prop)).*(v_out_low_p.*alpha_z_low_p)))
                    C_Ex_low = reshape([C_low_s_Ex C_low_p_Ex-1im*C_low_p_dEz_over_dx], nx_Ex, ny_Ex, 1, :)
                    C_Ey_low = reshape([C_low_s_Ey C_low_p_Ey-1im*C_low_p_dEz_over_dy], nx_Ey, ny_Ey, 1, :)
                    if two_sided
                        C_high_s_Ex = conj(u_prop_high_Ex*(conj(channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_out_high_s.*alpha_x_high_s))) # use implicit expansion
                        C_high_s_Ey = conj(u_prop_high_Ey*(conj(channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_out_high_s.*alpha_y_high_s)))
                        C_high_p_Ex = conj(u_prop_high_Ex*(conj(channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_out_high_p.*alpha_x_high_p)))
                        C_high_p_Ey = conj(u_prop_high_Ey*(conj(channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_out_high_p.*alpha_y_high_p)))
                        C_high_p_dEz_over_dx = conj(u_prop_high_dEz_over_dx*(conj(1 ./channels.high.sqrt_nu_prop.*cos.(channels.high.kzdx_prop/2).*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_out_high_p.*alpha_z_high_p)))
                        C_high_p_dEz_over_dy = conj(u_prop_high_dEz_over_dy*(conj(1 ./channels.high.sqrt_nu_prop.*cos.(channels.high.kzdx_prop/2).*exp.((-1im*dn)*channels.high.kzdx_prop)).*(v_out_high_p.*alpha_z_high_p)))
                        C_Ex_high = reshape([C_high_s_Ex C_high_p_Ex-1im*C_high_p_dEz_over_dx], nx_Ex, ny_Ex, 1, :) 
                        C_Ey_high = reshape([C_high_s_Ey C_high_p_Ey-1im*C_high_p_dEz_over_dy], nx_Ey, ny_Ey, 1, :)
                    end
               else # or output wavefronts specified by v_out_low and v_out_high
                    # u_prop_low_Ex is a large ny_Ex-by-N_prop_L matrix; instead of conjugating it, we conjugate the smaller ny_Ex-by-M_out_L matrix
                    C_Ex_low = conj(u_prop_low_Ex*(conj(channels.low.sqrt_nu_prop.*exp.((-1im*dn)*channels.low.kzdx_prop)).*v_out_low))
                    if two_sided
                        C_Ex_high = conj(u_prop_high_Ex*(conj(channels.high.sqrt_nu_prop.*exp.((-1im*dn)*channels.high.kzdx_prop)).*v_out_high))                
                    end
               end
            end
        end
    end
    
    if opts.clear_memory
        if ~use_2D_TM
            u_prop_low_Ex = nothing; u_prop_low_Ey = nothing
            u_prop_low_dEz_over_dx = nothing; u_prop_low_dEz_over_dy = nothing
            u_prop_high_Ex = nothing; u_prop_high_Ey = nothing
            u_prop_high_dEz_over_dx = nothing; u_prop_high_dEz_over_dy = nothing

            alpha_x_low_s = nothing; alpha_y_low_s = nothing; alpha_x_low_p = nothing; alpha_y_low_p = nothing
            alpha_x_high_s = nothing; alpha_y_high_s = nothing; alpha_x_high_p = nothing; alpha_y_high_p = nothing
            B_low_s_Ex = nothing; B_low_s_Ey = nothing; B_low_p_Ex = nothing; B_low_p_Ey = nothing
            B_low_p_dEz_over_dx = nothing; B_low_p_dEz_over_dy = nothing
            B_high_s_Ex = nothing; B_high_s_Ey = nothing; B_high_p_Ex = nothing; B_high_p_Ey = nothing
            B_high_p_dEz_over_dx = nothing; B_high_p_dEz_over_dy = nothing
            C_low_s_Ex = nothing; C_low_s_Ey = nothing; C_low_p_Ex = nothing; C_low_p_Ey = nothing
            C_low_p_dEz_over_dx = nothing; C_low_p_dEz_over_dy = nothing
            C_high_s_Ex = nothing; C_high_s_Ey = nothing; C_high_p_Ex = nothing; C_high_p_Ey = nothing
            C_high_p_dEz_over_dx = nothing; C_high_p_dEz_over_dy = nothing
        else
            u_prop_low_Ex = nothing; u_prop_high_Ex = nothing;
        end
        GC.gc()
    end
    
    t2 = time(); timing_build_BC = t2 - t1
    if opts.verbal; @printf("elapsed time: %7.3f secs\n", timing_build_BC); end
            
    ## Part 3: call mesti() to do the computation
    # Add syst.epsilon_low and syst.epsilon_high, to be used in mesti()
    if two_sided
        if ~use_2D_TM
            syst.epsilon_xx = cat(syst.epsilon_low*ones(nx_Ex,ny_Ex,nz_extra[1]), syst.epsilon_xx, syst.epsilon_high*ones(nx_Ex,ny_Ex,nz_extra[2]), dims=3)
            syst.epsilon_yy = cat(syst.epsilon_low*ones(nx_Ey,ny_Ey,nz_extra[1]), syst.epsilon_yy, syst.epsilon_high*ones(nx_Ey,ny_Ey,nz_extra[2]), dims=3)
            syst.epsilon_zz = cat(syst.epsilon_low*ones(nx_Ez,ny_Ez,nz_extra[1]), syst.epsilon_zz, syst.epsilon_high*ones(nx_Ez,ny_Ez,nz_extra[2]), dims=3)
            if (isdefined(syst, :epsilon_xy) && ~isa(syst.epsilon_xy, Nothing)) 
                syst.epsilon_xy = cat(zeros(nx_Ez,ny_Ez,nz_extra[1]), syst.epsilon_xy, zeros(nx_Ez,ny_Ez,nz_extra[2]), dims=3)                
            end
            if (isdefined(syst, :epsilon_xz) && ~isa(syst.epsilon_xz, Nothing))             
                syst.epsilon_xz = cat(zeros(nx_Ey,ny_Ex,nz_extra[1]), syst.epsilon_xz, zeros(nx_Ey,ny_Ex,nz_extra[2]), dims=3)                
            end
            if (isdefined(syst, :epsilon_yx) && ~isa(syst.epsilon_yx, Nothing)) 
                syst.epsilon_yx = cat(zeros(nx_Ez,ny_Ez,nz_extra[1]), syst.epsilon_yx, zeros(nx_Ez,ny_Ez,nz_extra[2]), dims=3)                
            end
            if (isdefined(syst, :epsilon_yz) && ~isa(syst.epsilon_yz, Nothing)) 
                syst.epsilon_yz = cat(zeros(nx_Ey,ny_Ex,nz_extra[1]), syst.epsilon_yz, zeros(nx_Ey,ny_Ex,nz_extra[2]), dims=3)                
            end
            if (isdefined(syst, :epsilon_zx) && ~isa(syst.epsilon_zx, Nothing)) 
                syst.epsilon_zx = cat(zeros(nx_Ey,ny_Ez,nz_extra[1]), syst.epsilon_zx, zeros(nx_Ey,ny_Ez,nz_extra[2]), dims=3)                
            end
            if (isdefined(syst, :epsilon_zy) && ~isa(syst.epsilon_zy, Nothing)) 
                syst.epsilon_zy = cat(zeros(nx_Ez,ny_Ex,nz_extra[1]), syst.epsilon_zy, zeros(nx_Ez,ny_Ex,nz_extra[2]), dims=3)                
            end
        else
            syst.epsilon_xx = cat(syst.epsilon_low*ones(ny_Ex,nz_extra[1]), syst.epsilon_xx, syst.epsilon_high*ones(ny_Ex,nz_extra[2]), dims=2)
        end
    else
        if ~use_2D_TM
            syst.epsilon_xx = cat(syst.epsilon_low*ones(nx_Ex,ny_Ex,nz_extra[1]), syst.epsilon_xx, dims=3)
            syst.epsilon_yy = cat(syst.epsilon_low*ones(nx_Ey,ny_Ey,nz_extra[1]), syst.epsilon_yy, dims=3)
            syst.epsilon_zz = cat(syst.epsilon_low*ones(nx_Ez,ny_Ez,nz_extra[1]), syst.epsilon_zz, dims=3)
            if (isdefined(syst, :epsilon_xy) && ~isa(syst.epsilon_xy, Nothing)) 
                syst.epsilon_xy = cat(syst.epsilon_low*ones(nx_Ez,ny_Ez,nz_extra[1]), syst.epsilon_xy, dims=3)
            end
            if (isdefined(syst, :epsilon_xz) && ~isa(syst.epsilon_xz, Nothing))             
                syst.epsilon_xz = cat(syst.epsilon_low*ones(nx_Ey,ny_Ex,nz_extra[1]), syst.epsilon_xz, dims=3)
            end
            if (isdefined(syst, :epsilon_yx) && ~isa(syst.epsilon_yx, Nothing)) 
                syst.epsilon_yx = cat(syst.epsilon_low*ones(nx_Ez,ny_Ez,nz_extra[1]), syst.epsilon_yx, dims=3)
            end
            if (isdefined(syst, :epsilon_yz) && ~isa(syst.epsilon_yz, Nothing)) 
                syst.epsilon_yz = cat(syst.epsilon_low*ones(nx_Ey,ny_Ex,nz_extra[1]), syst.epsilon_yz, dims=3)
            end
            if (isdefined(syst, :epsilon_zx) && ~isa(syst.epsilon_zx, Nothing)) 
                syst.epsilon_zx = cat(syst.epsilon_low*ones(nx_Ey,ny_Ez,nz_extra[1]), syst.epsilon_zx, dims=3)
            end
            if (isdefined(syst, :epsilon_zy) && ~isa(syst.epsilon_zy, Nothing)) 
                syst.epsilon_zy = cat(syst.epsilon_low*ones(nx_Ez,ny_Ex,nz_extra[1]), syst.epsilon_zy, dims=3)
            end
        else
            syst.epsilon_xx = cat(syst.epsilon_low*ones(ny_Ex,nz_extra[1]), syst.epsilon_xx, dims=2)
        end
    end
    syst.epsilon_low = nothing # mesti() will throw warning if syst.epsilon_low is given
    syst.epsilon_high = nothing # mesti() will throw warning if syst.epsilon_high is given

    # Whether we use PML or have a one-sided geometry, the BC in z direction will be Dirichlet when building matrix A in mesti()
    syst.zBC = "PEC"

    # location of inputs/outputs on the low surface (l = 0)
    l_low = nz_extra[1]
    
    # location of inputs/outputs on the high surface (l = nz_Ex+1 = nz_Ey+1)
    l_high = nz_extra[1] + nz_Ex + 1
        
    # Specify inputs
    if ~use_2D_TM                        
        B_Ex = Source_struct()
        B_Ey = Source_struct()
        B_Ez = Source_struct()
        B_Ez.isempty = true # The source only has x and y components; see Ref XXX ......
        if ~two_sided
            B_Ex.pos = [[1, 1, l_low, nx_Ex, ny_Ex, 1]]
            B_Ey.pos = [[1, 1, l_low, nx_Ey, ny_Ey, 1]]
            B_Ex.data = [B_Ex_low]
            B_Ey.data = [B_Ey_low]
        else
            B_Ex.pos = [[1, 1, l_low, nx_Ex, ny_Ex, 1], [1, 1, l_high, nx_Ex, ny_Ex, 1]]
            B_Ey.pos = [[1, 1, l_low, nx_Ey, ny_Ey, 1], [1, 1, l_high, nx_Ey, ny_Ey, 1]]
            B_Ex.data = [B_Ex_low, B_Ex_high]
            B_Ey.data = [B_Ey_low, B_Ey_high]
        end
        B = [B_Ex, B_Ey, B_Ez]
    else
        B_Ex = Source_struct()
        if ~two_sided
            B_Ex.pos = [[1, l_low, ny_Ex, 1]]
            B_Ex.data = [B_Ex_low]
        else
            B_Ex.pos = [[1, l_low, ny_Ex, 1], [1, l_high, ny_Ex, 1]]
            B_Ex.data = [B_Ex_low, B_Ex_high]
        end
        B = [B_Ex]
    end
                            
    # Specify outputs
    if build_C
        if ~use_2D_TM
            C_Ex = Source_struct()
            C_Ey = Source_struct()
            C_Ez = Source_struct()
            C_Ez.isempty = true # The output projection only has x and y components; see Ref XXX ......        
            if ~two_sided
                C_Ex.pos = [[1, 1, l_low, nx_Ex, ny_Ex, 1]]
                C_Ey.pos = [[1, 1, l_low, nx_Ey, ny_Ey, 1]]
                C_Ex.data = [C_Ex_low]
                C_Ey.data = [C_Ey_low]
            else
                C_Ex.pos = [[1, 1, l_low, nx_Ex, ny_Ex, 1], [1, 1, l_high, nx_Ex, ny_Ex, 1]]
                C_Ey.pos = [[1, 1, l_low, nx_Ey, ny_Ey, 1], [1, 1, l_high, nx_Ey, ny_Ey, 1]]
                C_Ex.data = [C_Ex_low, C_Ex_high]
                C_Ey.data = [C_Ey_low, C_Ey_high]
            end
            C = [C_Ex, C_Ey, C_Ez]
        else
            C_Ex = Source_struct()
            if ~two_sided
                C_Ex.pos = [[1, l_low, ny_Ex, 1]]
                C_Ex.data = [C_Ex_low]
            else
                C_Ex.pos = [[1, l_low, ny_Ex, 1], [1, l_high, ny_Ex, 1]]
                C_Ex.data = [C_Ex_low, C_Ex_high]
            end
            C = [C_Ex]
        end
    elseif use_transpose_B
        C = "transpose(B)"
    else
        C = nothing
    end
    
    if ~use_2D_TM
        # Find out the index of trivial and nontrivial channels for input and output
        ind_in_trivial_ch = Vector{Int64}()
        ind_in_nontrivial_ch = Vector{Int64}()
        for ii = 1:M_in_low
            if iszero(B_Ex_low[:,:,:,ii]) && iszero(B_Ey_low[:,:,:,ii])
                push!(ind_in_trivial_ch,ii)
            else
                push!(ind_in_nontrivial_ch,ii)            
            end
        end
        if two_sided
            for ii = 1:M_in_high
                if iszero(B_Ex_high[:,:,:,ii]) && iszero(B_Ey_high[:,:,:,ii])
                    push!(ind_in_trivial_ch,M_in_low+ii)
                else
                    push!(ind_in_nontrivial_ch,M_in_low+ii)                            
                end
            end
        end

        if build_C
            ind_out_trivial_ch = Vector{Int64}()
            ind_out_nontrivial_ch = Vector{Int64}()
            for ii = 1:M_out_low
                if iszero(C_Ex_low[:,:,:,ii]) && iszero(C_Ey_low[:,:,:,ii])
                    push!(ind_out_trivial_ch,ii)
                else
                    push!(ind_out_nontrivial_ch,ii)                            
                end
            end
            if two_sided
                for ii = 1:M_out_high
                    if iszero(C_Ex_high[:,:,:,ii]) && iszero(C_Ey_high[:,:,:,ii])
                        push!(ind_out_trivial_ch,M_out_low+ii)                    
                    else
                        push!(ind_out_nontrivial_ch,M_out_low+ii)                            
                    end
                end
            end
        end
    end
    
    if opts.clear_memory
        if ~use_2D_TM
            B_Ex_low = nothing; B_Ey_low = nothing; B_Ex_high = nothing; B_Ey_high = nothing
            B_Ex = nothing; B_Ey = nothing; B_Ez = nothing
            C_Ex_low = nothing; C_Ey_low = nothing; C_Ex_high = nothing; C_Ey_high = nothing
            C_Ex = nothing; C_Ey = nothing; C_Ez = nothing
        else
            B_Ex_low = nothing; B_Ex_high = nothing; B_Ex = nothing
            C_Ex_low = nothing; C_Ex_high = nothing; C_Ex = nothing
        end
        GC.gc()
    end

    # D will be subtracted later
    D = nothing

    # variables syst, B, and C can be cleared from mesti() since we don't need them anymore
    opts.clear_syst = opts.clear_memory
    opts.clear_BC = opts.clear_memory

    if ~isdefined(opts, :analysis_only) || isa(opts.analysis_only, Nothing)
        opts.analysis_only = false       
    end 

    # Main computation happens here
    # Note that we should no longer use syst, B, and C beyond this point since they may be cleared inside mesti()
    (S, info) = mesti(syst, B, C, D, opts)
    if opts.clear_memory
        GC.gc()
    end
    
    # Include the -2i prefactor that should have been in the input matrix B
    S = (-2im)*S
    
    if opts.return_field_profile && ~opts.analysis_only
        M = size(S,2)
        if ~use_2D_TM
            Ex = reshape(S[1:nt_tot_Ex, :], nx_Ex, ny_Ex, nz_tot_Ex, M)
            Ey = reshape(S[nt_tot_Ex+1:nt_tot_Ex+nt_tot_Ey, :], nx_Ey, ny_Ey, nz_tot_Ey, M)
            Ez = reshape(S[nt_tot_Ex+nt_tot_Ey+1:nt_tot_Ex+nt_tot_Ey+nt_tot_Ez, :], nx_Ez, ny_Ez, nz_tot_Ez, M)
        else
            Ex = reshape(S[1:nt_tot_Ex, :], ny_Ex, nz_tot_Ex, M)
        end
    end
    
    if ~use_2D_TM
        info.ind_in_trivial_ch  = ind_in_trivial_ch
        info.ind_in_nontrivial_ch  = ind_in_nontrivial_ch

        if build_C
            info.ind_out_trivial_ch = ind_out_trivial_ch
            info.ind_out_nontrivial_ch = ind_out_nontrivial_ch
        end
    end
    
    info.opts.prefactor = nothing
    info.opts.symmetrize_K = use_transpose_B    
            
    info.timing_init  = info.timing_init  + timing_init
    info.timing_build = info.timing_build + timing_build_BC # combine with build time for A and K

    if info.opts.analysis_only
        t3 = time()
        info.timing_total = t3 - t0
        if opts.verbal; @printf("          Total elapsed time: %7.3f secs\n", info.timing.total); end
        return S, channels, info
    end

    t1 = time()
                    
    ## Part 4: wrap up
                        
    # Recover the original list of input and output channels if we symmetrized K = [A,B;C,0]
    if use_transpose_B # when opts.symmetrize_K = true
        # Indices for the original list of input channels on the low
        ind_in = ind_in_out_low[1:M_in_low]

        # Indices for the original list of output channels on the low
        # There is no need to use ind_prop_conj again since the use of C = transpose(B) already compensates the previous use of ind_prop_conj.
        ind_out = ind_in_out_low[M_in_low .+ (1:M_out_low)]

        if two_sided
            # Include channels on the high; note ind_in_out_B and ind_in_out_T are column vectors
            ind_in  = [ind_in; length(ind_low) .+ ind_in_out_high[1:M_in_high]]
            ind_out = [ind_out; length(ind_low) .+ ind_in_out_high[M_in_high .+ (1:M_out_high)]]
        end
        
        S = S[ind_out, ind_in]
    end        

    # Include the exp(-i*kzdx*dn) prefactors that should have been in the input matrix B
    if use_ind_in # input channels specified by ind_in_low_s, ind_in_low_p, ind_in_high_s, or ind_in_high_p
                # or input channels specified by ind_in_low and ind_in_high
        if ~use_2D_TM
            prefactor = [exp.((-1im*dn)*channels.low.kzdx_prop[ind_in_low_s]); exp.((-1im*dn)*channels.low.kzdx_prop[ind_in_low_p])]
            if two_sided
                prefactor = [prefactor; exp.((-1im*dn)*channels.high.kzdx_prop[ind_in_high_s]); exp.((-1im*dn)*channels.high.kzdx_prop[ind_in_high_p])]
            end
            if opts.return_field_profile
                Ex = Ex.*reshape(prefactor, 1, 1, 1, :) # use implicit expansion
                Ey = Ey.*reshape(prefactor, 1, 1, 1, :)
                Ez = Ez.*reshape(prefactor, 1, 1, 1, :)
            else
                S = S.*reshape(prefactor,1,:) # use implicit expansion
            end
        else
            prefactor = channels.low.sqrt_nu_prop[ind_in_low].*exp.((-1im*dn)*channels.low.kzdx_prop[ind_in_low])
            if two_sided
                prefactor = [prefactor; channels.high.sqrt_nu_prop[ind_in_high].*exp.((-1im*dn)*channels.high.kzdx_prop[ind_in_high])]
            end
            if opts.return_field_profile
                Ex = Ex.*reshape(prefactor, 1, 1, :) # use implicit expansion
            else
                S = S.*reshape(prefactor,1,:) # use implicit expansion
            end
        end
    end            
    
    if ~opts.return_field_profile
        # Include the exp(-i*kzdx*dn) prefactors that should have been in the output matrix C
        if use_ind_out # output channels specified by ind_out_low_s, ind_out_low_p, ind_out_high_s, or ind_out_high_p 
                       # or output channels specified by ind_out_low and ind_out_high
            if ~use_2D_TM
                prefactor = [exp.((-1im*dn)*channels.low.kzdx_prop[ind_out_low_s]); exp.((-1im*dn)*channels.low.kzdx_prop[ind_out_low_p])]
                if two_sided
                    prefactor = [prefactor; exp.((-1im*dn)*channels.high.kzdx_prop[ind_out_high_s]); exp.((-1im*dn)*channels.high.kzdx_prop[ind_out_high_p])]
                end
            else
                prefactor = channels.low.sqrt_nu_prop[ind_out_low].*exp.((-1im*dn)*channels.low.kzdx_prop[ind_out_low])
                if two_sided
                    prefactor = [prefactor; channels.high.sqrt_nu_prop[ind_out_high].*exp.((-1im*dn)*channels.high.kzdx_prop[ind_out_high])]
                end
            end
            S = prefactor.*S; # use implicit expansion
        end
                                
        # Subtract D = C*inv(A_0)*B - S_0 where A_0 is a reference system and S_0 is its scattering matrix
        # The form of D is summation phase shift by dn over x, y, and z components; see Ref XXX ......    
        # When user-specified input and output wavefronts are used, we have D_low = (v_out_low')*diag(D)*v_in_low 
        if ~use_2D_TM
            phase_factor = repeat(exp.((-1im*2*dn)*channels.low.kzdx_prop), 2)
            if use_ind_in
                D_low = spdiagm(N_prop_low*2, N_prop_low*2, 0 => vec(phase_factor))
                D_low = D_low[:, [ind_in_low_s; N_prop_low.+ind_in_low_p]] 
            else
                D_low = phase_factor .* [[v_in_low_s; spzeros(N_prop_low, M_in_low_s)] [spzeros(N_prop_low, M_in_low_p); v_in_low_p]]
            end
            if use_ind_out
                D_low = D_low[[ind_out_low_s; N_prop_low.+ind_out_low_p],:]
            else
                D_low = [[v_out_low_s; spzeros(N_prop_low, M_out_low_s)] [spzeros(N_prop_low, M_out_low_p); v_out_low_p]]' * D_low
            end        
            if two_sided && (M_in_high != 0 || M_out_high != 0)
                phase_factor = repeat(exp.((-1im*2*dn)*channels.high.kzdx_prop), 2)
                if use_ind_in
                    D_high = spdiagm(N_prop_high*2, N_prop_high*2, 0 => vec(phase_factor))
                    D_high = D_high[:, cat(ind_in_high_s, N_prop_high.+ind_in_high_p, dims=1)]
                else
                    D_high = phase_factor .* [[v_in_high_s; spzeros(N_prop_high, M_in_high_s)] [spzeros(N_prop_high, M_in_high_p); v_in_high_p]]
                end
                if use_ind_out
                    D_high = D_high[[ind_out_high_s; N_prop_high.+ind_out_high_p], :]
                else
                    D_high = [[v_out_high_s; spzeros(N_prop_high, M_out_high_s)] [spzeros(N_prop_high, M_out_high_p); v_out_high_p]]'*D_high  
                end
                S = S - [[D_low; spzeros(M_out_high, M_in_low)] [spzeros(M_out_low, M_in_high); D_high]]
            else
                S = S - D_low  
            end
       else
            phase_factor = exp.((-1im*2*dn)*channels.low.kzdx_prop)
            if use_ind_in
                D_low = spdiagm(N_prop_low, N_prop_low, 0 => vec(phase_factor))
                D_low = D_low[:, ind_in_low]
            else
                D_low = phase_factor .* v_in_low
            end
            if use_ind_out
                D_low = D_low[ind_out_low, :]
            else
                D_low = v_out_low' * D_low
            end
            if two_sided && (M_in_high != 0 || M_out_high != 0)
                phase_factor = exp.((-1im*2*dn)*channels.high.kzdx_prop)
                if use_ind_in
                    D_high = spdiagm(N_prop_high, N_prop_high, 0 => vec(phase_factor))
                    D_high = D_high[:, ind_in_high]
                else
                    D_high = phase_factor .* v_in_high
                end
                if use_ind_out
                    D_high = D_high[ind_out_high, :]
                else
                    D_high = v_out_high' * D_high
                end
                S = S - [[D_low; spzeros(M_out_high, M_in_low)] [spzeros(M_out_low, M_in_high); D_high]]
            else
                S = S - D_low
            end
       end
    else # when opts.return_field_profile = true
        # Ex returned by mesti() has size [nx_Ex, ny_Ex, nz_tot_Ex, M_in_low+M_in_high] where nz_Ex_tot = nz_Ex + sum(nz_extra). Ey and Ez follow the same.
        # Here, we remove the PML.npixels and npixels_spacer pixels
        nz_remove_low = nz_extra[1] - 1 # we keep the surface pixel
        if two_sided
            nz_remove_high = nz_extra[2] - 1
        else
            nz_remove_high = 0
        end
        if ~use_2D_TM
            Ex = Ex[:, :, (1+nz_remove_low):(nz_tot_Ex-nz_remove_high), :]
            Ey = Ey[:, :, (1+nz_remove_low):(nz_tot_Ey-nz_remove_high), :]
            Ez = Ez[:, :, (1+nz_remove_low):(nz_tot_Ez-nz_remove_high), :]
        else
            Ex = Ex[:, (1+nz_remove_low):(nz_tot_Ex-nz_remove_high), :]
        end
                
        # At this point, Ex, Ey, and Ez have this number of pixels in z:
        # Ex, two-sided: 1+nz_Ex+1
        # Ex, one-sided: 1+nz_Ex
        # Ey, two-sided: 1+nz_Ey+1
        # Ey, one-sided: 1+nz_Ey
        # Ez, two-sided: 1+nz_Ez+1
        # Ez, one-sided: 1+nz_Ez
        
        if opts.nz_low <= 1 && (~two_sided || opts.nz_high <= 1) # nothing to add (except possibly zeros)
            if ~(opts.nz_low == 1 && ((two_sided && opts.nz_high == 1) || (~two_sided && opts.nz_high == 0)))
                # Remove pixels such that we have nz_low pixels of syst.epsilon_low and nz_high pixels of syst.epsilon_high
                if two_sided
                    if ~use_2D_TM
                        Ex = Ex[:, :, (2-opts.nz_low):(1+nz_Ex+opts.nz_high), :]
                        Ey = Ey[:, :, (2-opts.nz_low):(1+nz_Ey+opts.nz_high), :]
                        Ez = Ez[:, :, (2-opts.nz_low):(1+nz_Ez+opts.nz_high), :]
                    else
                        Ex = Ex[:, (2-opts.nz_low):(1+nz_Ex+opts.nz_high), :]
                    end
                else
                    if ~use_2D_TM
                        Ex = Ex[:, :, (2-opts.nz_low):(1+nz_Ex), :]
                        Ey = Ey[:, :, (2-opts.nz_low):(1+nz_Ey), :]
                        Ez = Ez[:, :, (2-opts.nz_low):(1+nz_Ez), :]
                    else
                        Ex = Ex[:, (2-opts.nz_low):(1+nz_Ex), :]
                    end
                    if opts.nz_high > 0
                        # Add nz_high slices of zero on the high
                        if ~use_2D_TM
                            Ex = cat(Ex, zeros(nx_Ex, ny_Ex, opts.nz_high, M_in_low), dims=3)
                            Ey = cat(Ey, zeros(nx_Ey, ny_Ey, opts.nz_high, M_in_low), dims=3)
                            Ez = cat(Ez, zeros(nx_Ez, ny_Ez, opts.nz_high, M_in_low), dims=3)
                        else
                            Ex = cat(Ex, zeros(ny_Ex, opts.nz_high, M_in_low), dims=2)
                        end
                    end
                end
            end
        else
            if opts.verbal; @printf("            ... "); end
            
            if ~use_2D_TM
                u_Ex = kron(channels.u_x_m(channels.kydx), channels.u_x_n(channels.kxdx))
                u_Ey = kron(channels.u_y_m(channels.kydx), channels.u_y_n(channels.kxdx))
                u_Ez = kron(channels.u_z_m(channels.kydx), channels.u_z_n(channels.kxdx))
                #u_dEz_over_dx = kron(channels.u_z_m(channels.kydx), channels.du_z_n(channels.kxdx))
                #u_dEz_over_dy = kron(channels.du_z_m(channels.kydx), channels.u_z_n(channels.kxdx))

                kappa_x_all_low = sin.(reshape(reapt(channels.kxdx_all,1,size(channels.kydx_all)),:)/2)  
                kappa_y_all_low = sin.(reshape(reapt(transpose(channels.kydx_all),size(channels.kxdx_all),1),:)/2)    
                kappa_z_all_low = sin.((channels.low.kzdx_all)/2)     

                denominator = sqrt.(kappa_x_all_low.^2+kappa_y_all_low.^2)
                alpha_x_all_low_s = -kappa_y_all_low./denominator
                alpha_y_all_low_s = kappa_x_all_low./denominator
                alpha_z_all_low_s = zeros(size(channels.low.kzdx_all,1),:)
                
                alpha_x_all_low_s[isnan.(alpha_x_all_low_s)] .= 0 
                alpha_y_all_low_s[isnan.(alpha_y_all_low_s)] .= 1

                denominator = sqrt.((abs.(kappa_x_all_low.*kappa_z_all_low)).^2+(abs.(kappa_y_all_low.*kappa_z_all_low)).^2+(abs.(kappa_x_all_low.^2+kappa_y_all_low.^2)).^2)
                alpha_x_all_low_p = kappa_x_all_low.*kappa_z_all_low./denominator
                alpha_y_all_low_p = kappa_y_all_low.*kappa_z_all_low./denominator
                alpha_z_all_low_p = (-1)*(-1)*(kappa_x_all_low.^2+kappa_y_all_low.^2)./denominator # for propagation along -z direction

                alpha_x_all_low_p[isnan.(alpha_x_all_low_p)] .= 1   
                alpha_y_all_low_p[isnan.(alpha_y_all_low_p)] .= 0  
                alpha_z_all_low_p[isnan.(alpha_z_all_low_p)] .= 0                

                # u where the a-th column is the a-th channels; it includes all the propagating and evanescent channels.
                # u_low = [[u_Ex.*reshape(alpha_x_all_low_s,1,:); u_Ey.*reshape(alpha_y_all_low_s,1,:)]
                #          [u_Ex.*reshape(alpha_x_all_low_p,1,:)+1im*u_dEz_over_dx*reshape(cos.(channels.low.kzdx_all/2).*alpha_z_all_low_p./sin.(channels.low.kzdx_all),1,:); 
                #           u_Ey.*reshape(alpha_y_all_low_p,1,:)+1im*u_dEz_over_dy*reshape(cos.(channels.low.kzdx_all/2).*alpha_z_all_low_p./sin.(channels.low.kzdx_all),1,:)]]
                # u_low = [[u_Ex.*reshape(alpha_x_all_low_s,1,:); u_Ey.*reshape(alpha_y_all_low_s,1,:); u_Ez.*reshape(alpha_z_all_low_s,1,:)]
                #          [u_Ex.*reshape(alpha_x_all_low_p,1,:); u_Ey.*reshape(alpha_y_all_low_p,1,:); u_Ez.*reshape(alpha_z_all_low_p,1,:)]]
                u_low_s = [u_Ex.*reshape(alpha_x_all_low_s,1,:); u_Ey.*reshape(alpha_y_all_low_s,1,:); u_Ez.*reshape(alpha_z_all_low_s,1,:)]
                u_low_p = [u_Ex.*reshape(alpha_x_all_low_p,1,:); u_Ey.*reshape(alpha_y_all_low_p,1,:); u_Ez.*reshape(alpha_z_all_low_p,1,:)]

                # We use u' to project the field on the low or high surface onto the complete and orthonormal set of transverse modes.            
                u_low_s_prime = u_low_s'
                u_low_p_prime = u_low_p'

                if two_sided
                    kappa_x_all_high = kappa_x_all_low  
                    kappa_y_all_high = kappa_y_all_low    
                    kappa_z_all_high = sin.((channels.high.kzdx_all)/2)     

                    alpha_x_all_high_s = alpha_x_all_low_s
                    alpha_y_all_high_s = alpha_x_all_low_s
                    alpha_z_all_high_s = zeros(size(channels.high.kzdx_all,1),:)
                    
                    denominator = sqrt.((abs.(kappa_x_all_high.*kappa_z_all_high)).^2+(abs.(kappa_y_all_high.*kappa_z_all_high)).^2+(abs.(kappa_x_all_high.^2+kappa_y_all_high.^2)).^2)
                    alpha_x_all_high_p = kappa_x_all_high.*kappa_z_all_high./denominator
                    alpha_y_all_high_p = kappa_y_all_high.*kappa_z_all_high./denominator
                    alpha_z_all_high_p = (-1)*(kappa_x_all_high.^2+kappa_y_all_high.^2)./denominator # for propagation along +z direction
                    alpha_x_all_high_p[isnan.(alpha_x_all_high_p)] .= 1
                    alpha_y_all_high_p[isnan.(alpha_y_all_high_p)] .= 0
                    alpha_z_all_high_p[isnan.(alpha_z_all_high_p)] .= 0 

                    # u_high = [[u_Ex.*reshape(alpha_x_all_high_s,1,:); u_Ey.*reshape(alpha_y_all_high_s,1,:)]
                    #           [u_Ex.*reshape(alpha_x_all_high_p,1,:)+1im*u_dEz_over_dx*reshape(cos.(channels.high.kzdx_all/2).*alpha_z_all_high_p./sin.(channels.high.kzdx_all),1,:);
                    #            u_Ey.*reshape(alpha_y_all_high_p,1,:)+1im*u_dEz_over_dy*reshape(cos.(channels.high.kzdx_all/2).*alpha_z_all_high_p./sin.(channels.high.kzdx_all),1,:)]] 
                    # u_high = [[u_Ex.*reshape(alpha_x_all_high_s,1,:); u_Ey.*reshape(alpha_y_all_high_s,1,:); u_Ez.*reshape(alpha_z_all_high_s,1,:)]
                    #           [u_Ex.*reshape(alpha_x_all_high_p,1,:); u_Ey.*reshape(alpha_y_all_high_p,1,:); u_Ez.*reshape(alpha_z_all_high_p,1,:)]]
                    u_high_s = [u_Ex.*reshape(alpha_x_all_high_s,1,:); u_Ey.*reshape(alpha_y_all_high_s,1,:); u_Ez.*reshape(alpha_z_all_high_s,1,:)]
                    u_high_p = [u_Ex.*reshape(alpha_x_all_high_p,1,:); u_Ey.*reshape(alpha_y_all_high_p,1,:); u_Ez.*reshape(alpha_z_all_high_p,1,:)]

                    u_high_s_prime = u_high_s'
                    u_high_p_prime = u_high_p'    
                end
            else
                u = channels.u_x_m(channels.kydx_all)
                u_prime = u'
            end
            
            nz_low_extra = opts.nz_low - 1

            if nz_low_extra == -1
                # Remove the single pixel of syst.epsilon_low on the low
                if ~use_2D_TM
                    Ex = Ex[:, :, 2:end, :]
                    Ey = Ey[:, :, 2:end, :]
                    Ez = Ez[:, :, 2:end, :]
                else
                    Ex = Ex[:, 2:end, :]
                end
            elseif nz_low_extra > 0
                if ~use_2D_TM
                    # Add pixels such that we have opts.nz_low pixels of syst.epsilon_low on the low
                    Ex_low = zeros(ComplexF64, nx_Ex, ny_Ex, nz_low_extra, size(Ex,4))
                    Ey_low = zeros(ComplexF64, nx_Ey, ny_Ey, nz_low_extra, size(Ey,4))
                    Ez_low = zeros(ComplexF64, nx_Ez, ny_Ez, nz_low_extra, size(Ez,4))                
                    # The l below is l = l - l_low, where l_low is the location of the inputs/outputs on the low surface.
                    l = (-nz_low_extra):1:-1

                    # What is the case nx_Ex != nx_Ey != nx_Ez?
                    # What is the case ny_Ex != ny_Ey != ny_Ez?
                    
                    nx = maximum([nx_Ex,nx_Ey])
                    ny = maximum([ny_Ex,ny_Ey])
                    #kz_z = [repeat(reshape(channels.low.kzdx_all, nx*ny, 1), 2, 1)].*reshape(l, 1, :) # kz*z; 2*nx*ny-by-nz_low_extra matrix through implicit expansion
                    kz_z = reshape(channels.low.kzdx_all, nx*ny, 1).*reshape(l, 1, :) # kz*z; nx*ny-by-nz_low_extra matrix through implicit expansion
                    exp_pikz = exp.( 1im*kz_z) # exp(+i*kz*z)
                    exp_mikz = exp.(-1im*kz_z) # exp(-i*kz*z)
                    c_in_s = zeros(ComplexF64, nx*ny, 1)
                    c_in_p = zeros(ComplexF64, nx*ny, 1)

                    l_low = 1 # index for the inputs/outputs on the low surface

                    if M_in_low > 0
                        prefactor = reshape(exp.((-1im*dn)*channels.low.kzdx_prop)./channels.low.sqrt_nu_prop, N_prop_low, 1)
                    end
                    ###
                    for ii = 1:M_in_low_s # s-polarized input from low
                        #c = u_low_prime*[reshape(Ex[:, :, l_low, ii], :); reshape(Ey[:, :, l_low, ii], :)] # c is a 2*nx*ny-by-1 column vector of transverse mode coefficients
                        c_s = u_low_s_prime*[reshape(Ex[:, :, l_low, ii], :); reshape(Ey[:, :, l_low, ii], :); reshape(Ez[:, :, l_low, ii], :)] # c_s is a nx*ny-by-1 column vector of transverse mode coefficients for s-polarization
                        c_p = u_low_p_prime*[reshape(Ex[:, :, l_low, ii], :); reshape(Ey[:, :, l_low, ii], :); reshape(Ez[:, :, l_low, ii], :)] # c_p is a nx*ny-by-1 column vector of transverse mode coefficients for p-polarization
        
                        # c_in_s is the s-polarized incident wavefront at n_low; note we need to back propagate dn pixel from z=0
                        c_in_s[:] .= 0
                        if use_ind_in
                            c_in_s[channels.low.ind_prop[ind_in_low_s[ii]]] = prefactor[ind_in_low_s[ii]]
                        else
                            c_in_s[channels.low.ind_prop] = prefactor.*v_in_low_s[:,ii]                                
                        end
                        c_out_s = c_s - c_in_s
                        c_out_p = c_p - c_in_p

                        Ex_low[:,:,:,ii] = reshape(u_Ex.*reshape(alpha_x_all_low_s,1,:)*(c_in_s.*exp_pikz + c_out_s.*exp_mikz) + u_Ex.*reshape(alpha_x_all_low_p,1,:)*(c_in_p.*exp_pikz + c_out_p.*exp_mikz), nx_Ex,ny_Ex,nz_low_extra)
                        Ey_low[:,:,:,ii] = reshape(u_Ey.*reshape(alpha_y_all_low_s,1,:)*(c_in_s.*exp_pikz + c_out_s.*exp_mikz) + u_Ey.*reshape(alpha_y_all_low_p,1,:)*(c_in_p.*exp_pikz + c_out_p.*exp_mikz), nx_Ey,ny_Ey,nz_low_extra)
                        Ez_low[:,:,:,ii] = reshape(u_Ez.*reshape(alpha_z_all_low_s,1,:)*(c_in_s.*exp_pikz + c_out_s.*exp_mikz) + u_Ez.*reshape(alpha_z_all_low_p,1,:)*(c_in_p.*exp_pikz + c_out_p.*exp_mikz), nx_Ez,ny_Ez,nz_low_extra)                       
                    end
                    for ii = 1:M_in_low_p # p-polarized input from low
                        #c = u_low_prime*[reshape(Ex[:, :, l_low, M_in_low_s+ii], :); reshape(Ey[:, :, l_low, M_in_low_s+ii], :)] # c is a 2*nx*ny-by-1 column vector of transverse mode coefficients
                        c_s = u_low_s_prime*[reshape(Ex[:, :, l_low, M_in_low_s+ii], :); reshape(Ey[:, :, l_low, M_in_low_s+ii], :); reshape(Ez[:, :, l_low, M_in_low_s+ii], :)] # c_s is a nx*ny-by-1 column vector of transverse mode coefficients for s-polarization
                        c_p = u_low_p_prime*[reshape(Ex[:, :, l_low, M_in_low_s+ii], :); reshape(Ey[:, :, l_low, M_in_low_s+ii], :); reshape(Ez[:, :, l_low, M_in_low_s+ii], :)] # c_p is a nx*ny-by-1 column vector of transverse mode coefficients for p-polarization

                        c_in_p[:] .= 0
                        if use_ind_in
                            c_in_p[channels.low.ind_prop[ind_in_low_p[ii]]] = prefactor[ind_in_low_p[ii]]
                        else
                            c_in_p[channels.low.ind_prop] = prefactor.*v_in_low_p[:,ii]
                        end
                        c_out_s = c_s - c_in_s
                        c_out_p = c_p - c_in_p

                        Ex_low[:,:,:,M_in_low_s+ii] = reshape(u_Ex.*reshape(alpha_x_all_low_s,1,:)*(c_in_s.*exp_pikz + c_out_s.*exp_mikz) + u_Ex.*reshape(alpha_x_all_low_p,1,:)*(c_in_p.*exp_pikz + c_out_p.*exp_mikz), nx_Ex,ny_Ex,nz_low_extra)
                        Ey_low[:,:,:,M_in_low_s+ii] = reshape(u_Ey.*reshape(alpha_y_all_low_s,1,:)*(c_in_s.*exp_pikz + c_out_s.*exp_mikz) + u_Ey.*reshape(alpha_y_all_low_p,1,:)*(c_in_p.*exp_pikz + c_out_p.*exp_mikz), nx_Ey,ny_Ey,nz_low_extra)
                        Ez_low[:,:,:,M_in_low_s+ii] = reshape(u_Ez.*reshape(alpha_z_all_low_s,1,:)*(c_in_s.*exp_pikz + c_out_s.*exp_mikz) + u_Ez.*reshape(alpha_z_all_low_p,1,:)*(c_in_p.*exp_pikz + c_out_p.*exp_mikz), nx_Ez,ny_Ez,nz_low_extra)                       
                    end

                    if two_sided
                        for ii = 1:M_in_high_s # s-polarized input from high
                            # c_out_s = c_s because there is no input on the high side
                            c_out_s = u_low_s_prime*[reshape(Ex[:, :, l_low, M_in_low+ii], :); reshape(Ey[:, :, l_low, M_in_low+ii], :); reshape(Ez[:, :, l_low, M_in_low+ii], :)]
                            c_out_p = u_low_p_prime*[reshape(Ex[:, :, l_low, M_in_low+ii], :); reshape(Ey[:, :, l_low, M_in_low+ii], :); reshape(Ez[:, :, l_low, M_in_low+ii], :)]
                            
                            Ex_low[:,:,:,M_in_low+ii] = reshape(u_Ex.*reshape(alpha_x_all_low_s,1,:)*(c_out_s.*exp_mikz) + u_Ex.*reshape(alpha_x_all_low_p,1,:)*(c_out_p.*exp_mikz), nx_Ex,ny_Ex,nz_low_extra)
                            Ey_low[:,:,:,M_in_low+ii] = reshape(u_Ey.*reshape(alpha_y_all_low_s,1,:)*(c_out_s.*exp_mikz) + u_Ey.*reshape(alpha_y_all_low_p,1,:)*(c_out_p.*exp_mikz), nx_Ey,ny_Ey,nz_low_extra)
                            Ez_low[:,:,:,M_in_low+ii] = reshape(u_Ez.*reshape(alpha_z_all_low_s,1,:)*(c_out_s.*exp_mikz) + u_Ez.*reshape(alpha_z_all_low_p,1,:)*(c_out_p.*exp_mikz), nx_Ez,ny_Ez,nz_low_extra)
                        end
                        for ii = 1:M_in_high_p # p-polarized input from high
                            # c_out_s = c_s because there is no input on the high side
                            c_out_s = u_low_s_prime*[reshape(Ex[:, :, l_low, M_in_low+M_in_high_s+ii], :); reshape(Ey[:, :, l_low, M_in_low+M_in_high_s+ii], :); reshape(Ez[:, :, l_low, M_in_low+M_in_high_s+ii], :)]
                            c_out_p = u_low_p_prime*[reshape(Ex[:, :, l_low, M_in_low+M_in_high_s+ii], :); reshape(Ey[:, :, l_low, M_in_low+M_in_high_s+ii], :); reshape(Ez[:, :, l_low, M_in_low+M_in_high_s+ii], :)]

                            Ex_low[:,:,:,M_in_low+M_in_high_s+ii] = reshape(u_Ex.*reshape(alpha_x_all_low_s,1,:)*(c_out_s.*exp_mikz) + u_Ex.*reshape(alpha_x_all_low_p,1,:)*(c_out_p.*exp_mikz), nx_Ex,ny_Ex,nz_low_extra)
                            Ey_low[:,:,:,M_in_low+M_in_high_s+ii] = reshape(u_Ey.*reshape(alpha_y_all_low_s,1,:)*(c_out_s.*exp_mikz) + u_Ey.*reshape(alpha_y_all_low_p,1,:)*(c_out_p.*exp_mikz), nx_Ey,ny_Ey,nz_low_extra)
                            Ez_low[:,:,:,M_in_low+M_in_high_s+ii] = reshape(u_Ez.*reshape(alpha_z_all_low_s,1,:)*(c_out_s.*exp_mikz) + u_Ez.*reshape(alpha_z_all_low_p,1,:)*(c_out_p.*exp_mikz), nx_Ez,ny_Ez,nz_low_extra)
                        end
                    end
                    Ex = cat(Ex_low, Ex, dims=3)
                    Ey = cat(Ey_low, Ey, dims=3)
                    Ez = cat(Ez_low, Ez, dims=3)  
                else
                    # Add pixels such that we have opts.nz_low pixels of syst.epsilon_low on the low
                    Ex_low = zeros(ComplexF64, ny_Ex, nz_low_extra, size(Ex,3))
                    # The l below is l = l - l_low, where l_low is the location of the inputs/outputs on the low surface.
                    l = (-nz_low_extra):1:-1
                    # For the back propagation, incident wavefront would just use propagating channels and output wavefront would use all channels
                    kz_z = reshape(channels.low.kzdx_all, ny_Ex, 1).*reshape(l, 1, :) # kz*z; ny_Ex-by-nz_low_extra matrix through implicit expansion
                    exp_mikz = exp.(-1im*kz_z) # exp(-i*kz*z)
                    kz_z_prop = reshape(channels.low.kzdx_prop, channels.low.N_prop, 1).*reshape(l, 1, :) # kz*z; channels.low.N_prop-by-nz_low_extra matrix through implicit expansion
                    exp_pikz_prop = exp.( 1im*kz_z_prop); # exp(+i*kz*z)
                    c_in = zeros(ComplexF64, ny_Ex, 1)
                    c_in_prop = zeros(ComplexF64, channels.low.N_prop, 1)
                    l_low = 1 # index for the inputs/outputs on the low surface
                    if M_in_low > 0
                        prefactor = reshape(exp.((-1im*dn)*channels.low.kzdx_prop)./channels.low.sqrt_nu_prop, N_prop_low, 1)
                    end
                    for ii = 1:M_in_low # input from low
                        c = u_prime*Ex[:, l_low, ii] # c is a ny_Ex-by-1 column vector of transverse mode coefficients
                        # c_in is the incident wavefront at n_low; note we need to back propagate dn pixel from z=0
                        c_in[:] .= 0
                        if use_ind_in
                            c_in_prop[ind_in_low[ii]] = prefactor[ind_in_low[ii]]
                        else
                            c_in_prop = prefactor.*v_in_low[:,ii]                                
                        end
                        c_in[channels.low.ind_prop] = c_in_prop
                        c_out = c - c_in
                        Ex_low[:,:,ii] = u[:,channels.low.ind_prop]*(c_in_prop.*exp_pikz_prop) + u*(c_out.*exp_mikz)               
                    end
                    if two_sided
                        for ii = 1:M_in_high # input from high
                            c_out = u_prime*Ex[:, l_low, M_in_low+ii] # c_out = c because there is no input on the high side
                            Ex_low[:,:,M_in_low+ii] = u*(c_out.*exp_mikz)
                        end
                    end
                    Ex = cat(Ex_low, Ex, dims=2)
                end
            end
                    
            if two_sided
                nz_high_extra = opts.nz_high - 1
                if nz_high_extra == -1
                    # Remove the single pixel of syst.epsilon_high on the high
                    if ~use_2D_TM
                        Ex = Ex[:, :, 1:(end-1), :]
                        Ey = Ey[:, :, 1:(end-1), :]
                        Ez = Ez[:, :, 1:(end-1), :]
                    else
                        Ex = Ex[:, 1:(end-1), :]
                    end
                elseif nz_high_extra > 0
                    if ~use_2D_TM
                        # Add pixels such that we have opts.nz_high pixels of syst.epsilon_high on the high
                        Ex_high = zeros(ComplexF64, nx_Ex, ny_Ex, nz_high_extra, size(Ex,4))
                        Ey_high = zeros(ComplexF64, nx_Ey, ny_Ey, nz_high_extra, size(Ey,4))
                        Ez_high = zeros(ComplexF64, nx_Ez, ny_Ez, nz_high_extra, size(Ez,4))                
                        # The n below is l = l - l_high, where l_high is the location of the inputs/outputs on the high surface.
                        l = 1:nz_high_extra

                        # What is the case nx_Ex != nx_Ey != nx_Ez?
                        # What is the case ny_Ex != ny_Ey != ny_Ez?
                        
                        nx = maximum([nx_Ex,nx_Ey])
                        ny = maximum([ny_Ex,ny_Ey])
                        kz_z = reshape(channels.high.kzdx_all, nx*ny, 1).*reshape(l, 1, :) # kz*z; nx*ny-by-nz_low_extra matrix through implicit expansion
                        exp_pikz = exp.( 1im*kz_z) # exp(+i*kz*z)
                        exp_mikz = exp.(-1im*kz_z) # exp(-i*kz*z)
                        c_in_s = zeros(ComplexF64, nx*ny, 1)
                        c_in_p = zeros(ComplexF64, nx*ny, 1)
                        l_high = size(Ex, 3) # index for the inputs/outputs on the high surface

                        for ii = 1:M_in_low_s # s-polarized input from low
                            # c_out_s = c_s because there is no input on the low side
                            c_out_s = u_low_s_prime*[reshape(Ex[:, :, l_high, ii], :); reshape(Ey[:, :, l_high, ii], :); reshape(Ez[:, :, l_high+1, ii], :)] 
                            c_out_p = u_low_p_prime*[reshape(Ex[:, :, l_high, ii], :); reshape(Ey[:, :, l_high, ii], :); reshape(Ez[:, :, l_high+1, ii], :)] 

                            Ex_high[:,:,:,ii] = reshape(u_Ex.*reshape(alpha_x_all_high_s,1,:)*(c_out_s.*exp_pikx) + u_Ex.*reshape(alpha_x_all_high_p,1,:)*(c_out_p.*exp_pikx), nx_Ex,ny_Ex,nz_high_extra)
                            Ey_high[:,:,:,ii] = reshape(u_Ey.*reshape(alpha_y_all_high_s,1,:)*(c_out_s.*exp_pikx) + u_Ey.*reshape(alpha_y_all_high_p,1,:)*(c_out_p.*exp_pikx), nx_Ey,ny_Ey,nz_high_extra)
                            Ez_high[:,:,:,ii] = reshape(u_Ez.*reshape(alpha_z_all_high_s,1,:)*(c_out_s.*exp_pikx) + u_Ez.*reshape(alpha_z_all_high_p,1,:)*(c_out_p.*exp_pikx), nx_Ez,ny_Ez,nz_high_extra)
                        end                            

                        for ii = 1:M_in_low_p # p-polarized input from low
                            # c_out_s = c_s because there is no input on the low side
                            c_out_s = u_low_s_prime*[reshape(Ex[:, :, l_high, M_in_low_s+ii], :); reshape(Ey[:, :, l_high, M_in_low_s+ii], :); reshape(Ez[:, :, l_high+1, M_in_low_s+ii], :)] 
                            c_out_p = u_low_p_prime*[reshape(Ex[:, :, l_high, M_in_low_s+ii], :); reshape(Ey[:, :, l_high, M_in_low_s+ii], :); reshape(Ez[:, :, l_high+1, M_in_low_s+ii], :)] 

                            Ex_high[:,:,:,M_in_low_s+ii] = reshape(u_Ex.*reshape(alpha_x_all_high_s,1,:)*(c_out_s.*exp_pikx) + u_Ex.*reshape(alpha_x_all_high_p,1,:)*(c_out_p.*exp_pikx), nx_Ex,ny_Ex,nz_high_extra)
                            Ey_high[:,:,:,M_in_low_s+ii] = reshape(u_Ey.*reshape(alpha_y_all_high_s,1,:)*(c_out_s.*exp_pikx) + u_Ey.*reshape(alpha_y_all_high_p,1,:)*(c_out_p.*exp_pikx), nx_Ey,ny_Ey,nz_high_extra)
                            Ez_high[:,:,:,M_in_low_s+ii] = reshape(u_Ez.*reshape(alpha_z_all_high_s,1,:)*(c_out_s.*exp_pikx) + u_Ez.*reshape(alpha_z_all_high_p,1,:)*(c_out_p.*exp_pikx), nx_Ez,ny_Ez,nz_high_extra)
                        end

                        if M_in_high_s > 0 || M_in_high_p > 0
                            prefactor = reshape(exp.((-1im*dn)*channels.high.kzdx_prop)./channels.high.sqrt_nu_prop, N_prop_high, 1)                                
                        end
                        for ii = 1:M_in_high_s # s-polarized input from high
                            c_s = u_high_s_prime*[reshape(Ex[:, :, l_high, M_in_low+ii], :); reshape(Ey[:, :, l_high, M_in_low+ii], :); reshape(Ez[:, :, l_high+1, M_in_low+ii], :)] # c_s is a nx*ny-by-1 column vector of transverse mode coefficients
                            c_p = u_high_p_prime*[reshape(Ex[:, :, l_high, M_in_low+ii], :); reshape(Ey[:, :, l_high, M_in_low+ii], :); reshape(Ez[:, :, l_high+1, M_in_low+ii], :)] # c_p is a nx*ny-by-1 column vector of transverse mode coefficients

                            c_in_s[:] .= 0
                            if use_ind_in
                                c_in_s[channels.high.ind_prop[ind_in_high_s[ii]]] = prefactor[ind_in_high_s[ii]]
                            else
                                c_in_s[channels.high.ind_prop] = prefactor.*v_in_high_s[:,ii]
                            end
                            c_out_s = c_s - c_in_s
                            c_out_p = c_p - c_in_p

                            Ex_high[:,:,:,M_in_low+ii] = reshape(u_Ex.*reshape(alpha_x_all_high_s,1,:)*(c_in_s.*exp_mikz + c_out_s.*exp_pikz) + u_Ex.*reshape(alpha_x_all_high_p,1,:)*(c_in_p.*exp_mikz + c_out_p.*exp_pikz), nx_Ex,ny_Ex,nz_high_extra)
                            Ey_high[:,:,:,M_in_low+ii] = reshape(u_Ey.*reshape(alpha_y_all_high_s,1,:)*(c_in_s.*exp_mikz + c_out_s.*exp_pikz) + u_Ey.*reshape(alpha_y_all_high_p,1,:)*(c_in_p.*exp_mikz + c_out_p.*exp_pikz), nx_Ey,ny_Ey,nz_high_extra)
                            Ez_high[:,:,:,M_in_low+ii] = reshape(u_Ez.*reshape(alpha_z_all_high_s,1,:)*(c_in_s.*exp_mikz + c_out_s.*exp_pikz) + u_Ez.*reshape(alpha_z_all_high_p,1,:)*(c_in_p.*exp_mikz + c_out_p.*exp_pikz), nx_Ez,ny_Ez,nz_high_extra)
                        end

                        for ii = 1:M_in_high_p # p-polarized input from high
                            c_s = u_high_s_prime*[reshape(Ex[:, :, l_high, M_in_low+M_in_high_s+ii], :); reshape(Ey[:, :, l_high, M_in_low+M_in_high_s+ii], :); reshape(Ez[:, :, l_high+1, M_in_low+M_in_high_s+ii], :)] # c_s is a nx*ny-by-1 column vector of transverse mode coefficients
                            c_p = u_high_p_prime*[reshape(Ex[:, :, l_high, M_in_low+M_in_high_s+ii], :); reshape(Ey[:, :, l_high, M_in_low+M_in_high_s+ii], :); reshape(Ez[:, :, l_high+1, M_in_low+M_in_high_s+ii], :)] # c_p is a nx*ny-by-1 column vector of transverse mode coefficients
                            
                            c_in_p[:] .= 0
                            if use_ind_in
                                c_in_p[channels.high.ind_prop[ind_in_high_s[ii]]] = prefactor[ind_in_high_p[ii]]
                            else
                                c_in_p[channels.high.ind_prop] = prefactor.*v_in_high_p[:,ii]
                            end
                            c_out_s = c_s - c_in_s
                            c_out_p = c_p - c_in_p

                            Ex_high[:,:,:,M_in_low+M_in_high_s+ii] = reshape(u_Ex.*reshape(alpha_x_all_high_s,1,:)*(c_in_s.*exp_mikz + c_out_s.*exp_pikz) + u_Ex.*reshape(alpha_x_all_high_p,1,:)*(c_in_p.*exp_mikz + c_out_p.*exp_pikz), nx_Ex,ny_Ex,nz_high_extra)
                            Ey_high[:,:,:,M_in_low+M_in_high_s+ii] = reshape(u_Ey.*reshape(alpha_y_all_high_s,1,:)*(c_in_s.*exp_mikz + c_out_s.*exp_pikz) + u_Ey.*reshape(alpha_y_all_high_p,1,:)*(c_in_p.*exp_mikz + c_out_p.*exp_pikz), nx_Ey,ny_Ey,nz_high_extra)
                            Ez_high[:,:,:,M_in_low+M_in_high_s+ii] = reshape(u_Ez.*reshape(alpha_z_all_high_s,1,:)*(c_in_s.*exp_mikz + c_out_s.*exp_pikz) + u_Ez.*reshape(alpha_z_all_high_p,1,:)*(c_in_p.*exp_mikz + c_out_p.*exp_pikz), nx_Ez,ny_Ez,nz_high_extra)
                        end
                        Ex = cat(Ex, Ex_high, dims=3)
                        Ey = cat(Ey, Ey_high, dims=3)
                        Ez = cat(Ez, Ez_high, dims=3)                
                    else
                        # Add pixels such that we have opts.nz_high pixels of syst.epsilon_high on the high
                        Ex_high = zeros(ComplexF64, ny_Ex, nz_high_extra, size(Ex,3))
                        # The n below is l = l - l_high, where l_high is the location of the inputs/outputs on the high surface.
                        l = 1:nz_high_extra
                        # For the back propagation, incident wavefront would just use propagating channels and output wavefront would use all channels                
                        kz_z = reshape(channels.high.kzdx_all, ny_Ex, 1).*reshape(l, 1, :) # kz*z; ny_Ex-by-nz_low_extra matrix through implicit expansion
                        exp_pikz = exp.( 1im*kz_z) # exp(+i*kz*z)
                        kz_z_prop = reshape(channels.high.kzdx_prop, channels.high.N_prop, 1).*reshape(l, 1, :) # kz*z; channels.high.N_prop-by-nz_high_extra matrix through implicit expansion
                        exp_mikz_prop = exp.(-1im*kz_z_prop) # exp(-i*kz*z)
                        c_in = zeros(ComplexF64, ny_Ex, 1)
                        c_in_prop = zeros(ComplexF64, channels.high.N_prop, 1)
                        l_high = size(Ex, 2) # index for the inputs/outputs on the high surface
                        for ii = 1:M_in_low # input from low
                            c_out = u_prime*Ex[:, l_high, ii] # c_out = c because there is no input on the high side
                            Ex_high[:,:,ii] = u*(c_out.*exp_pikz)
                        end                
                        if M_in_high > 0
                                prefactor = reshape(exp.((-1im*dn)*channels.high.kzdx_prop)./channels.high.sqrt_nu_prop, N_prop_high, 1)       
                        end                     
                        for ii = 1:M_in_high # input from high
                            c = u_prime*Ex[:, l_high, M_in_low+ii] # c is a ny_Ex-by-1 column vector of transverse mode coefficients
                            # c_in is the incident wavefront at n_R; note we need to back propagate dn pixel from z=d
                            c_in[:] .= 0
                            if use_ind_in
                                c_in_prop[ind_in_high[ii]] = prefactor[ind_in_high[ii]]
                            else
                                c_in_prop = prefactor.*v_in_high[:,ii]
                            end
                            c_in[channels.high.ind_prop] = c_in_prop
                            c_out = c - c_in
                            Ex_high[:,:,M_in_low+ii] = u[:,channels.high.ind_prop]*(c_in_prop.*exp_mikz_prop) + u*(c_out.*exp_pikz)
                        end
                        Ex = cat(Ex, Ex_high, dims=2)
                    end
                end
            elseif opts.nz_high > 0 # one-sided
               # Add nz_high slices of zero on the high
               if ~use_2D_TM
                   Ex = cat(Ex, zeros(nx_Ex, ny_Ex, opts.nz_high, M_in_low), dims=3)
                   Ey = cat(Ey, zeros(nx_Ey, ny_Ey, opts.nz_high, M_in_low), dims=3)
                   Ez = cat(Ez, zeros(nx_Ez, ny_Ez, opts.nz_high, M_in_low), dims=3)
               else
                   Ex = cat(Ex, zeros(ny_Ex, opts.nz_high, M_in_low), dims=2)
               end
            end
                                    
            t2 = time()
            if opts.verbal; @printf("elapsed time: %7.3f secs\n", t2-t1); end
        end        
    end

    t3 = time()
    info.timing_solve = info.timing_solve + t3 - t1 # Add the post-processing time
    info.timing_total = t3 - t0

    if opts.verbal; @printf("          Total elapsed time: %7.3f secs\n", info.timing_total); end
    
    if opts.return_field_profile
        if ~use_2D_TM
            return Ex, Ey, Ez, channels, info
        else
            return Ex, channels, info
        end
    else
        return S, channels, info                    
    end
end

# The following are mesti functions to take different number of input arguments, but all of them will
# call the mesti main function.


# When syst, and input are specified; return the components of the field profiles.
function mesti2s(syst::Syst, input::Union{channel_type, channel_index, wavefront})
    # Check if 2D TM fields are required
    if ndims(syst.epsilon_xx) == 2
        use_2D_TM = true
    else
        use_2D_TM = false
    end
    
    if ~use_2D_TM
        (Ex, Ey, Ez, channels, info) = mesti2s(syst, input, nothing, nothing)
        return Ex, Ey, Ez, channels, info
    else
        (Ex, channels, info) = mesti2s(syst, input, nothing, nothing)
        return Ex, channels, info
    end
end


# When syst, input, and opts are specified; return the components of the field profiles.
function mesti2s(syst::Syst, input::Union{channel_type, channel_index, wavefront}, opts::Union{Opts,Nothing})
    # Check if 2D TM fields are required
    if ndims(syst.epsilon_xx) == 2
        use_2D_TM = true
    else
        use_2D_TM = false
    end
    
    if ~use_2D_TM
        (Ex, Ey, Ez, channels, info) = mesti2s(syst, input, nothing, opts)
        return Ex, Ey, Ez, channels, info
    else
        (Ex, channels, info) = mesti2s(syst, input, nothing, opts)
        return Ex, channels, info
    end
end

# When syst, input, and output are specified; return the scattering matrix.
function mesti2s(syst::Syst, input::Union{channel_type, channel_index, wavefront}, output::Union{channel_type, channel_index, wavefront,Nothing})
    (S, channels, info) = mesti2s(syst, input, output, nothing)
    return S, channels, info
end
