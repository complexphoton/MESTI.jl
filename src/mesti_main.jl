# Export composite data types
export Source_struct
export Syst

# Export a function mesti()
export mesti

mutable struct Source_struct
    # A composite data type to specfiy source
    # See also: mesti and mesti2s   
   pos::Vector{Vector{Int64}}
   data::Union{Vector{Array{Int64,2}},Vector{Array{Float64,2}},Vector{Array{ComplexF64,2}},
               Vector{Array{Int64,3}},Vector{Array{Float64,3}},Vector{Array{ComplexF64,3}},
               Vector{Array{Int64,4}},Vector{Array{Float64,4}},Vector{Array{ComplexF64,4}}} 
   ind::Vector{Vector{Int64}}
   isempty::Integer
    
   Source_struct() = new()
end

mutable struct Syst
    # A composite data type to specfiy system
    # See also: mesti and mesti2s   

    # Below are used in both mesti() and mesti2s()
    epsilon_xx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Matrix{Int64},Matrix{Float64},Matrix{ComplexF64},Nothing}
    epsilon_xy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}
    epsilon_xz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}

    epsilon_yx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}
    epsilon_yy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}
    epsilon_yz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}
 
    epsilon_zx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}
    epsilon_zy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}    
    epsilon_zz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}

    length_unit::String
    wavelength::Number
    dx::Real
    xBC::Union{String,Nothing}
    yBC::String
    kx_B::Union{Number,Nothing}
    ky_B::Number

    # Below are used in mesti() only
    zBC::String    
    kz_B::Number
    PML::Union{PML,Vector{PML}}
    PML_type::String    
    
    # Below are used in mesti2s() only
    epsilon_low::Union{Real,Nothing}
    epsilon_high::Union{Real,Nothing}    
    zPML::Union{PML,Vector{PML},Nothing}

    Syst() = new()
end 

"""
    MESTI Multi-source frequency-domain electromagnetic simulations.
        ---3D field profile---
        (Ex, Ey, Ez, info) = MESTI(syst, B) returns the spatial field profiles
        of Ex(x,y,z), Ey(x,y,z), and Ez(x,y,z) satisfying
            [(∇ × ∇ ×)  - (omega/c)^2*epsilon(x,y,z)]*[Ex(x,y,z); Ey(x,y,z); Ez(x,y,z)] = 
            i*omega*mu_0*[Jx(x,y,z); Jy(x,y,z); Jz(x,y,z)], where 
        The relative permittivity profile epsilon(x,y,z), frequency omega, and boundary conditions 
        are specified by structure "syst".
            Note that relative permittivity profile epsilon(x,y,z) is a rank 2 tenor: 
                [epsilon_xx, epsilon_xy, epsilon_xz; 
                epsilon_yx, epsilon_yy, epsilon_yz; 
                epsilon_zx, epsilon_zy, epsilon_zz] in general. 
        Users can specify the diagonal terms only (epsilon_xx(x,y,z), epsilon_yy(x,y,z), and epsilon_zz(x,y,z)) 
        or all of them.
        Each column of matrix "B" specifies a distinct input source profile.
        The returned "Ex", "Ey", "Ez" is a 4D array, such as Ex(:,:,:,i), 
        being the field profile Ex given the i-th input source profile. Same data structure for Ey and Ez. 
        The information of the computation is returned in structure "info".
        
        MESTI uses finite-difference discretization on the Yee lattice, after which
        the differential operator becomes an (nt_Ex*nt_Ey*nt_Ez)-by-(nt_Ex*nt_Ey*nt_Ez) sparse matrix A
        where nt_Ex is the number of total grid sites for Ex and nt_Ex = nx_Ex*ny_Ex*nz_Ex
        (nx_Ex, ny_Ex, nz_Ex) is the number of sites Ex is discretized onto. Same notation for Ey and Ez.
        Ex = reshape((inv(A)*B)[nt_Ex,:], nx_Ex, ny_Ex, nz_Ex, :).
        Ey = reshape((inv(A)*B)[nt_Ey,:], nx_Ey, ny_Ey, nz_Ey, :).
        Ez = reshape((inv(A)*B)[nt_Ez,:], nx_Ez, ny_Ez, nz_Ez, :).
            
        ---2D TM field profile---
        (Ex, info) = MESTI(syst, B) returns the spatial field profiles
        of Ex(y,z) for 2D transverse-magnetic (TM) fields satisfying
        [- (d/dy)^2 - (d/dz)^2 - (omega/c)^2*epsilon(y,z)] Ex(y,z) = i*omega*mu_0*Jx(y,z).
            The returned 'Ex' is a 3D array, with Ex(:,:,i) being the field profile
        of Ex given the i-th input source profile. The information of the computation is returned in structure 'info'.
            MESTI uses finite-difference discretization on the Yee lattice, after which
        the differential operator becomes an (ny_Ex*nz_Ex)-by-(ny_Ex*nz_Ex) sparse matrix A
        where (ny_Ex, nz_Ex) is the number of sites Ex is discretized onto, and
        Ex = reshape(inv(A)*B, ny_Ex, nz_Ex, []).

        ---Generalized scattering matrix S---
        (S, info) = MESTI(syst, B, C) returns S = C*inv(A)*B where the solution
        inv(A)*B is projected onto the output channels or locations of interest
        through matrix C; each row of matrix "C" is a distinct output projection
        profile, discretized into a 1-by-(nt_Ex*nt_Ey*nt_Ez) vector for 3D or
        1-by-(ny_Ex*nz_Ex) vector for 2D TM fields in the same order as matrix A.
        When the MUMPS is available, this is done by computing the Schur complement of an 
        augmented matrix K = [A,B;C,0] through a partial factorization.
        
        (S, info) = MESTI(syst, B, C, D) returns S = C*inv(A)*B - D. This can be
        used for the computation of scattering matrices, where S is the scattering
        matrix, and matrix D can be derived analytically or computed as D =
        C*inv(A0)*B - S0 from a reference system A0 for which the scattering matrix
        S0 is known.
        
        (Ex, Ey, Ez, info) = MESTI(syst, B, opts),
        (Ex, info) = MESTI(syst, B, opts),
        (S, info) = MESTI(syst, B, C, opts), and
        (S, info) = MESTI(syst, B, C, D, opts) allow detailed options to be
        specified with structure "opts".
        
        MESTI only considers nonmagnetic materials.
        
        This file checks and parses the parameters, and it can build matrices B and
        C from its nonzero elements specified by the user (see details below). It
        calls function mesti_build_fdfd_matrix() to build matrix A and function
        mesti_matrix_solver!() to compute C*inv(A)*B or inv(A)*B, where most of the
        computation is done.
        
        === Input Arguments ===
        syst (Syst struct; required):
            A structure that specifies the system, used to build the FDFD matrix A.
            It contains the following fields:
            syst.epsilon_xx (numeric array or matrix, real or complex, required):
                For 3D systems, an nx_Ex-by-ny_Ex-by-nz_Ex array discretizing the relative permittivity
                profile epsilon_xx(x,y,z). Specifically, syst.epsilon_xx(n,m,l) is the scalar
                epsilon_xx(n,m,l) averaged over a cube with volume (syst.dx)^3 centered at
                the point (x_n, y_m, z_l) where Ex(x,y,z) is located on the Yee lattice.
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
            syst.length_unit (string; optional):
                Length unit, such as micron, nm, or some reference wavelength. This
                code only uses dimensionless quantities, so syst.length_unit is never
                used. This syst.length_unit is meant to help the user interpret the
                units of (x,y,z), dx, wavelength, kx_B, ky_B, kz_B, etc.
            syst.wavelength (numeric scalar, real or complex; required):
                Vacuum wavelength 2*pi*c/omega, in units of syst.length_unit.
            syst.dx (positive scalar; required):
                Discretization grid size, in units of syst.length_unit.
            syst.PML (PML structure or a vector of PML structure; optional):
                Parameters of the perfectly matched layer (PML) used to simulate an
                open boundary. Note that PML is not a boundary condition; it is a
                layer placed within the simulation domain (just before the boundary)
                that attenuates outgoing waves with minimal reflection.
                    In mesti(), the PML starts from the interior of the system
                specified by syst.epsilon_ij (i = x,y,z and j = x,y,z), 
                and ends at the first or last pixel inside syst.epsilon_ij 
                (i = x,y,z and j = x,y,z). (Note: this is different from the function 
                mesti2s() that handles two-sided geometries, where the homogeneous spaces
                on the low and high are specified separately through syst.epsilon_low 
                and syst.epsilon_high, and where PML is placed in such homogeneous space, 
                outside of the syst.epsilon_ij (i = x,y,z and j = x,y,z).
                    When only one set of PML parameters is used in the system (as is
                the most common), such parameters can be specified with a scalar PML
                structure syst.PML that contains the following fields:
                    npixels (positive integer scalar; required): 
                        Number of PML pixels.
                        Note this is within syst.epsilon_ij (i = x,y,z and j = x,y,z),
                        not in addition to.
                    direction (string; optional): 
                        Direction(s) where PML is placed. Available choices are (case-insensitive):
                            "all" - (default) PML in x, y, and z directions for 3D and PML in y and z directions for 2D TM
                            "x"   - PML in x direction
                            "y"   - PML in y direction
                            "z"   - PML in z direction    
                    side (string; optional): 
                        Side(s) where PML is placed.Available choices are (case-insensitive):
                            "both" - (default) PML on both sides
                            "-"    - one-sided PML; end at the first pixel
                            "+"    - one-sided PML; end at the last pixel
                    power_sigma (non-negative scalar; optional): 
                        Power of the polynomial grading for the conductivity sigma; defaults to 3.
                    sigma_max_over_omega (non-negative scalar; optional):
                        Conductivity at the end of the PML; defaults to
                        0.8*(power_sigma+1)/((2*pi/wavelength)*dx*sqrt(epsilon_bg)).
                        where epsilon_bg is the average relative permittivity along the
                        last slice of the PML. This is used to attenuate propagating
                        waves.
                    power_kappa (non-negative scalar; optional): 
                        Power of the polynomial grading for the real-coordinate-stretching factor
                        kappa; defaults to 3.
                    kappa_max (real scalar no smaller than 1; optional):
                        Real-coordinate-stretching factor at the end of the PML;
                        defaults to 15. This is used to accelerate the attenuation of
                        evanescent waves. kappa_max = 1 means no real-coordinate
                        stretching.
                    power_alpha (non-negative scalar; optional): 
                        Power of the polynomial grading for the CFS alpha factor; defaults to 1.
                    alpha_max_over_omega (non-negative scalar; optional): 
                        Complex-frequency-shifting (CFS) factor at the beginning of the PML.
                        This is typically used in time-domain simulations to suppress
                        late-time (low-frequency) reflections. We don't use it by
                        default (alpha_max_over_omega = 0) since we are in frequency
                        domain.
                We use the following PML coordinate-stretching factor:
                    s(p) = kappa(p) + sigma(p)./(alpha(p) - i*omega)
                with
                    sigma(p)/omega = sigma_max_over_omega*(p.^power_sigma),
                    kappa(p) = 1 + (kappa_max-1)*(p.^power_kappa),
                    alpha(p)/omega = alpha_max_over_omega*((1-p).^power_alpha),
                where omega is frequency, and p goes linearly from 0 at the beginning
                of the PML to 1 at the end of the PML. 
                    If isdefined(syst, :PML) = false, which means no PML on any side. PML is
                only placed on the side(s) specified by syst.PML.
                    When multiple sets of PML parameters are used in the system (e.g.,
                a thinner PML on one side, a thicker PML on another side), these
                parameters can be specified with a vector of PML strcture.
                    syst.PML = [PML_1, PML_2, ...],
                with PML_1 and PML_2 each being a structure containing the above
                fields; they can specify different PML parameters on different sides.
                Each side cannot be specified more than once.
                    With real-coordinate stretching, PML can attenuate evanescent waves
                more efficiently than free space, so there is no need to place free
                space in front of PML.
                    The PML thickness should be chosen based on the acceptable level of
                reflectivity given the discretization resolution and the range of wave
                numbers (i.e., angles) involved; more PML pixels gives lower
                reflectivity. Typically 10-40 pixels are sufficient.
            syst.PML_type (string; optional):
                Type of PML. Available choices are (case-insensitive):
                    "UPML"   - (default) uniaxial PML
                    "SC-PML" - stretched-coordinate PML
                The two are mathematically equivalent, but using SC-PML has lower 
                condition number.
            syst.xBC (string or nothing; optional):
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
                By default, syst.xBC = "Bloch" if syst.kx_B is given; otherwise syst.xBC = "PEC"
                The choice of syst.xBC has little effect on the numerical accuracy
                when PML is used.
                For 2D TM fields, xBC = nothing.
            syst.yBC (string; optional):
                Boundary condition in y direction, analogous to syst.xBC.
            syst.zBC (string; optional):
                Boundary condition in z direction, analogous to syst.xBC.
            syst.kx_B (real scalar or nothing; optional):
                Bloch wave number in x direction, in units of 1/syst.length_unit.
                syst.kx_B is only used when syst.xBC = "Bloch". It is allowed to
                specify a complex-valued syst.kx_B, but a warning will be displayed.
                For 2D TM fields, kx_B = nothing.
            syst.ky_B (real scalar; optional):
                Bloch wave number in y direction, analogous to syst.kx_B.
            syst.kz_B (real scalar; optional):
                Bloch wave number in z direction, analogous to syst.kx_B.    
        B (numeric matrix or vector of source structure; required):
            Matrix specifying the input source profiles B in the C*inv(A)*B - D or
            C*inv(A)*B or inv(A)*B returned. When the input argument B is a matrix,
            it is directly used, and size(B,1) must equal nt_Ex*nt_Ey*nt_Ez for 3D system and nt_Ex for 2D TM case; 
            each column of B specifies a source profile, placed on the grid points of E.
                Note that matrix A is (syst.dx)^2 times the differential operator and
            is unitless, so each column of B is (syst.dx)^2 times the i*omega*mu_0*[Jx;Jy;Jz] on
            the right-hand side of the differential equation and has the same unit as E.
                Instead of specifying matrix B directly, one can specify only its
            nonzero parts, from which mesti() will build the sparse matrix B. To do
            so, B in the input argument should be set as a source structure vector; 
            (Bx_struct, By_struct, Bz_struct) here we refer to such source structure 
            as Bx_struct as source for x component to distinguish it from the
            resulting matrix B. If for every column of matrix B, all of its nonzero
            elements are spatially located within a cuboid (e.g., block sources),
            one can use the following fields:
                Bx_struct.pos (six-element integer vector or four-element integer vector): Bx_struct.pos =
                    [n1, m1, l1, d, w, h] specifies the location and the size of the
                    cuboid on Ex-grid. Here, (n1, m1, l1) is the index of the (x, y, z) coordinate of
                    the smaller-x, smaller-y, and smaller-z corner of the cuboid, 
                    at the location of f(n1, m1, l1) where f = Ex; (d, w, h) is the depth, width, and height 
                    of the cuboid, such that (n2, m2, l2) = (n1+d-1, m1+w-1, l1+h-1) is the index
                    of the higher-index corner of the cuboid.
                    For 2D TM case, users only need to specify [m1, l1, w, h].
                Bx_struct.data (2D, 3D or 4D numeric array): nonzero elements of matrix B for x component
                    within the cuboid specified by Bx_struct.pos.
                    When it is a 4D array (general 3D system), Bx_struct.data[n',m',l',a] is the a-th input
                    source at the location of f(n=n1+n'-1, m=m1+m'-1, l=l1+l'-1), which becomes
                    B[n+(m-1)*nx_Ex+(l-1)*nx_Ex*ny_Ex, a]. In other words, Bx_struct.data[:,:,:,a] gives the
                    sources at the cuboid f(n1+(0:(d-1)), m1+(0:(w-1)), l1+(0:(h-1))). So,
                    size(Bx_struct.data) must equal (d, w, h, number of inputs).
                        When it is a 3D array (2D TM system), Bx_struct.data(m',l',a) is the a-th input
                    source at the location of f(m=m1+m'-1, l=l1+l'-1), which becomes
                    Bx(m+(l-1)*ny, a). In other words, Bx_struct.data(:,:,a) gives the
                    sources at the rectangle f(m1+(0:(w-1)), l1+(0:(h-1))). So,
                    size(Bx_struct.data, [1,2]) must equal [w, h], and
                    size(Bx_struct.data, 3) is the number of inputs.
                    Alternatively, Bx_struct.data can be a 2D array that is
                    equivalent to reshape(data_in_4D_array, d*w*h, :) or reshape(data_in_3D_array, w*h, :), in which case
                    size(Bx_struct.data, 2) is the number of inputs; in this case,
                    Bx_struct.data can be a sparse matrix, and its sparsity will be
                    preserved when building matrix B.
                By_struct.pos (six-element integer vector): By_struct.pos =
                    [n1, m1, l1, d, w, h] specifies the location and the size of the
                    cuboid on Ey-grid, analogous to Bx_struct.pos.
                By_struct.data (2D or 4D numeric array): nonzero elements of matrix B for y component
                    within the cuboid specified by By_struct.pos, analogous to Bx_struct.data.
                Bz_struct.pos (six-element integer vector): Bz_struct.pos =
                    [n1, m1, l1, d, w, h] specifies the location and the size of the
                    cuboid on Ez-grid, analogous to Bx_struct.pos.
                Bz_struct.data (2D or 4D numeric array): nonzero elements of matrix B for z component
                    within the cuboid specified by Bz_struct.pos, analogous to Bx_struct.data.
                If different inputs are located within different cuboids (e.g.,
            inputs from plane sources on the low and separate inputs from plane
            sources on the high), Bx_struct can be a structure array with multiple
            elements (e.g., Bx_struct.pos[1] and Bx_struct.data[1] specify plane sources
            on the low; Bx_struct.pos[2] and Bx_struct.data[2] specify plane sources on
            the high); these inputs are treated separately, and the total number of
            inputs is size(Bx_struct.data[1], 4) + size(Bx_struct.data[2], 4) + ... +
            size(Bx_struct.data[end], 4).
                If the nonzero elements of matrix B do not have cuboid shapes in
            space [e.g., for total-field/scattered-field (TF/SF) simulations], one
            can use a source structure with the following fields:
                Bx_struct.ind (integer vector): linear indices of the spatial
                    locations of the nonzero elements of matrix B for x component, such that
                    f(Bx_struct.ind) are the points where the source is placed. Such linear 
                    indices can be constructed from Base._sub2ind().
                Bx_struct.data (2D numeric matrix): nonzero elements of matrix B for x component 
                    at the locations specified by Bx_struct.ind. Specifically,
                    Bx_struct.data[i,a] is the a-th input source at the location of
                    f(Bx_struct.ind[i]), which becomes B[Bx_struct.ind[i], a]. So,
                    size(Bx_struct.data, 1) must equal length(Bx_struct.ind), and
                    size(Bx_struct.data, 2) is the number of inputs.
                By_struct.ind (integer vector): linear indices of the spatial
                    locations of the nonzero elements of matrix B for y component, analogous to Bx_struct.ind.
                By_struct.data (2D numeric matrix): nonzero elements of matrix B for y component 
                    at the locations specified by By_struct.ind, analogous to Bx_struct.data.
                Bz_struct.ind (integer vector): linear indices of the spatial
                    locations of the nonzero elements of matrix B for z component, analogous to Bx_struct.ind.
                Bz_struct.data (2D numeric matrix): nonzero elements of matrix B for z component 
                    at the locations specified by Bz_struct.ind, analogous to Bx_struct.data.
                Similarly, one can use Bx_struct.ind[1], Bx_struct.ind[2] etc together
            with Bx_struct.data[1], Bx_struct.data[2] etc to specify inputs at
            different sets of locations. Every element of the structure array must
            have the same fields (e.g., one cannot specify Bx_struct.pos[1] and
            Bx_struct.ind[2]), so the more general Bx_struct.ind syntax should be used
            when some of the inputs are rectangular and some are not.
        C (numeric matrix, vector of source structure, or "transpose(B)"; optional):
            Matrix specifying the output projections in the C*inv(A)*B - D or
            C*inv(A)*B returned. When the input argument C is a matrix, it is
            directly used, and size(C,2) must equal nt_Ex*nt_Ey*nt_Ez for 3D system and nt_Ex for 2D TM case;
            each row of C specifies a projection profile, placed on the grid points of E.
                Scattering matrix computations often have C = transpose(B); if that
            is the case, the user can set C = "transpose(B)" as a string,
            and it will be replaced by transpose(B) in the code.
                For field-profile computations, the user can simply omit C from the
            input arguments, as in mesti(syst, B). If opts is needed, the user can use
            mesti(syst, B, opts).
                Similar to B, here one can specify only the nonzero parts of the
            output matrix C, from which mesti() will build the sparse matrix C. The
            syntax of vector of source structure (Cx_struct, Cy_struct, Cz_struct) is 
            the same as for B, summarized below. If for every row of matrix
            C, all of its nonzero elements are spatially located withing a cuboid
            (e.g., projection of fields on a line), one can set the input argument C
            to be a structure array (referred to as Cx_struct below) with the
            following fields:
                Cx_struct.pos (six-element integer vector or four-element integer vector): Cx_struct.pos =
                    [n1, m1, l1, d, w, h] specifies the location and the size of the
                    cuboid on Ex-grid. Here, (n1, m1, l1) is the index of the (x, y, z) coordinate of
                    the smaller-x, smaller-y, and smaller-z corner of the cuboid, 
                    at the location of f(n1, m1, l1) where f = Ex; (d, w, h) is the depth, width, and height 
                    of the cuboid, such that (n2, m2, l2) = (n1+d-1, m1+w-1, l1+h-1) is the index
                    of the higher-index corner of the cuboid.
                    For 2D TM case, users only need to specify [m1, l1, w, h].
                Cx_struct.data (2D, 3D or 4D numeric array): nonzero elements of matrix C for x component
                    within the cuboid specified by Cx_struct.pos.
                    When it is a 4D array (general 3D system), Cx_struct.data[n',m',l',a] is the a-th input
                    source at the location of f(n=n1+n'-1, m=m1+m'-1, l=l1+l'-1), which becomes
                    C[n+(m-1)*nx_Ex+(l-1)*nx_Ex*ny_Ex, a]. In other words, Cx_struct.data[:,:,:,a] gives the
                    sources at the cuboid f(n1+(0:(d-1)), m1+(0:(w-1)), l1+(0:(h-1))). So,
                    size(Cx_struct.data) must equal (d, w, h, number of inputs).
                        When it is a 3D array (2D system), Cx_struct.data(m',l',a) is the a-th input
                    source at the location of f(m=m1+m'-1, l=l1+l'-1), which becomes
                    Cx(m+(l-1)*ny, a). In other words, Cx_struct.data(:,:,a) gives the
                    sources at the rectangle f(m1+(0:(w-1)), l1+(0:(h-1))). So,
                    size(Cx_struct.data, [1,2]) must equal [w, h], and
                    size(Cx_struct.data, 3) is the number of inputs.
                    Alternatively, Cx_struct.data can be a 2D array that is
                    equivalent to reshape(data_in_4D_array, d*w*h, :) or reshape(data_in_3D_array, w*h, :), in which case
                    size(Cx_struct.data, 2) is the number of inputs; in this case,
                    Cx_struct.data can be a sparse matrix, and its sparsity will be
                    preserved when building matrix C.
                Cy_struct.pos (six-element integer vector): Cy_struct.pos =
                    [n1, m1, l1, d, w, h] specifies the location and the size of the
                    cuboid on Ey-grid, analogous to Cx_struct.pos.
                Cy_struct.data (2D or 4D numeric array): nonzero elements of matrix C for y component
                    within the cuboid specified by Cy_struct.pos, analogous to Cx_struct.data.
                Cz_struct.pos (six-element integer vector): Cz_struct.pos =
                    [n1, m1, l1, d, w, h] specifies the location and the size of the
                    cuboid on Ez-grid, analogous to Cx_struct.pos.
                Cz_struct.data (2D or 4D numeric array): nonzero elements of matrix C for z component
                    within the cuboid specified by Cz_struct.pos, analogous to Cx_struct.data.
                If the nonzero elements of matrix C do not have rectangular shapes in
            space [e.g., for near-field-to-far-field transformations], one can set C
            to a structure array with the following fields:
                Cx_struct.ind (integer vector): linear indices of the spatial
                    locations of the nonzero elements of matrix C for x component, such that
                    f(Cx_struct.ind) are the points where the source is placed.
                Cx_struct.data (2D numeric matrix): nonzero elements of matrix C for x component 
                    at the locations specified by Cx_struct.ind. Specifically,
                    Cx_struct.data[i,a] is the a-th input source at the location of
                    f(Cx_struct.ind[i]), which becomes C[Cx_struct.ind[i], a]. So,
                    size(Cx_struct.data, 1) must equal length(Cx_struct.ind), and
                    size(Cx_struct.data, 2) is the number of inputs.
                Cy_struct.ind (integer vector): linear indices of the spatial
                    locations of the nonzero elements of matrix C for y component, analogous to Cx_struct.ind.
                Cy_struct.data (2D numeric matrix): nonzero elements of matrix C for y component 
                    at the locations specified by Cy_struct.ind, analogous to Cx_struct.data.
                Cz_struct.ind (integer vector): linear indices of the spatial
                    locations of the nonzero elements of matrix C for z component, analogous to Cx_struct.ind.
                Cz_struct.data (2D numeric matrix): nonzero elements of matrix C for z component 
                    at the locations specified by Cz_struct.ind, analogous to Cx_struct.data.
                Like in Bx_struct, one can use structure arrays with multiple elements
            to specify outputs at different spatial locations.
        D (numeric matrix; optional):
            Matrix D in the C*inv(A)*B - D returned, which specifies the baseline
            contribution; size(D,1) must equal size(C,1), and size(D,2) must equal
            size(B,2). When D is not an input argument, it will not be subtracted from C*inv(A)*B. 
            For field-profile computations where C is not an input argument, D should not be an input
            argument, either.
        opts (Opts structure; optional):
            A structure that specifies the options of computation. 
            It can contain the following fields (all optional):
            opts.verbal (boolean scalar; optional, defaults to true):
                Whether to print system information and timing to the standard output.
            opts.prefactor (numeric scalar, real or complex; optional):
                When opts.prefactor is given, mesti() will return
                opts.prefactor*C*inv(A)*B - D or opts.prefactor*C*inv(A)*B or
                opts.prefactor*inv(A)*B. Such prefactor makes it easier to use C =
                transpose(B) to take advantage of reciprocity. Defaults to 1.
            opts.exclude_PML_in_field_profiles (boolean scalar; optional, defaults to false):
                When opts.exclude_PML_in_field_profiles = true, the PML pixels
                (specified by syst.PML.npixels) are excluded from the returned
                field_profiles on each side where PML is used; otherwise the full
                field profiles are returned. Only used for field-profile computations
                (i.e., when the output projection matrix C is not given).    
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
            opts.clear_BC (boolean scalar; optional, defaults to false):
                When opts.clear_BC = true, variables "B" and "C" will be cleared in
                the caller's workspace to reduce peak memory usage. Can be used when B
                and/or C take up significant memory and are not needed after calling
                mesti(). However, it is not implemented and does not work for current 
                julia version.
            opts.clear_syst (boolean scalar; optional, defaults to false):
                When opts.clear_syst = true, variable "syst" will be cleared in the
                caller's workspace to reduce peak memory usage. This can be used when
                syst.epsilon_xx, syst.epsilon_yy, and syst.epsilon_zz take up
                significant memory and are not needed after calling mesti(). However, 
                it is not implemented and does not work for current julia version.
            opts.clear_memory (boolean scalar; optional, defaults to true):
                Whether or not to clear variables inside mesti() to reduce peak memory
                usage. If opts.clear_memory == true, we call GC.gc() inside mesti().
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
        ---When users specify C---
        S (full numeric matrix):
            The generalized scattering matrix S = C*inv(A)*B or S = C*inv(A)*B - D.
        ---or when users do not specify C---
            When opts.exclude_PML_in_field_profiles = false,
            For 3D system,
        Ex (4D array):
            Electrical field profile for Ex component
            (nx_Ex, ny_Ex, nz_Ex, number of inputs) = size(Ex)
        Ey (4D array):
            Electrical field profile for Ey component
            (nx_Ey, ny_Ey, nz_Ey, number of inputs) = size(Ey)
        Ez (4D array):
            Electrical field profile for Ez component
            (nx_Ez, ny_Ez, nz_Ez, number of inputs) = size(Ez)
            ---For 2D TM case---
            Ex (3D array):
                Electrical field profile for Ex component
                (ny_Ex, nz_Ex, number of inputs) = size(Ex)
        When opts.exclude_PML_in_field_profiles = true, the PML pixels
            (specified by syst.PML.npixels) are excluded from Ex, Ey, and Ez on each
            side where PML is used.    
        info (Info structure):
            A structure that contains the following fields:
            info.opts (Opts structure):
                The final "opts" used, excluding the user-specified matrix ordering.
            info.timing_init (non-negative scalar):
                Timing in the "init" stages
            info.timing_build (non-negative scalar):
                Timing in the "building" stages
            info.timing_analyze (non-negative scalar):
                Timing in the "analysis" stages
            info.timing_factorize (non-negative scalar):
                Timing in the "factorizing" stages
            info.timing_total (non-negative scalar):    
                Timing in total
            info.xPML (two-element cell array; optional);
                PML parameters on the sides of x direction, if used.
            info.yPML (two-element cell array; optional);
                PML parameters on the sides of y direction, if used.
            info.zPML (two-element cell array; optional);
                PML parameters on the sides of z direction, if used.    
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
    
       See also: mesti_build_fdfd_matrix, mesti_matrix_solver!, mesti2s
"""
function mesti(syst::Syst, B::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2},Vector{Source_struct}}, C::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2}, Vector{Source_struct},String,Nothing}, D::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2},Array{Int32,2},Array{Float32,2},Array{ComplexF32,2},Nothing}, opts::Union{Opts,Nothing})
    
    if ~(stacktrace()[2].func == :mesti2s)
        # Make deepcopy of them to avoid mutating input argument 
        syst = deepcopy(syst); opts = deepcopy(opts)
    end
    
    ## Part 1.1: Check validity of syst, assign default values to its fields, and parse BC and PML specifications
    t0 = time() 
    
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
    
    if ~use_2D_TM
        if (isdefined(syst, :epsilon_yy) && ~isa(syst.epsilon_yy, Nothing)) || (isdefined(syst, :epsilon_zz) && ~isa(syst.epsilon_zz, Nothing)) || (isdefined(syst, :epsilon_xy) && ~isa(syst.epsilon_xy, Nothing)) || (isdefined(syst, :epsilon_xz) && ~isa(syst.epsilon_xz, Nothing)) || (isdefined(syst, :epsilon_yx) && ~isa(syst.epsilon_yx, Nothing)) || (isdefined(syst, :epsilon_yz) && ~isa(syst.epsilon_yz, Nothing)) || (isdefined(syst, :epsilon_zx) && ~isa(syst.epsilon_zx, Nothing)) || (isdefined(syst, :epsilon_zy) && ~isa(syst.epsilon_zy, Nothing))
            include_off_diagonal_epsilon = true
            if ~isdefined(syst, :epsilon_xy); syst.epsilon_xy = nothing; end
            if ~isdefined(syst, :epsilon_xz); syst.epsilon_xz = nothing; end
            if ~isdefined(syst, :epsilon_yx); syst.epsilon_yx = nothing; end
            if ~isdefined(syst, :epsilon_yz); syst.epsilon_yz = nothing; end
            if ~isdefined(syst, :epsilon_zx); syst.epsilon_zx = nothing; end
            if ~isdefined(syst, :epsilon_zy); syst.epsilon_zy = nothing; end
        else
            include_off_diagonal_epsilon = false
        end   
        
        # Number of grid points in x, y, and z for Ex, Ey, and Ez
        (nx_Ex, ny_Ex, nz_Ex) = size(syst.epsilon_xx)
        (nx_Ey, ny_Ey, nz_Ey) = size(syst.epsilon_yy)
        (nx_Ez, ny_Ez, nz_Ez) = size(syst.epsilon_zz)

        # Total number of grid points for Ex, Ey, and Ez    
        nt_Ex = nx_Ex*ny_Ex*nz_Ex
        nt_Ey = nx_Ey*ny_Ey*nz_Ey
        nt_Ez = nx_Ez*ny_Ez*nz_Ez
        
    else
        # Consider 2D TM field: Ex(y,z)
        # Number of grid points in y and z for Ex
        (ny_Ex, nz_Ex) = size(syst.epsilon_xx)
        
        # Total number of grid points for Ex
        nt_Ex = ny_Ex*nz_Ex
    end

    # Check that the user did not accidentally use options only in mesti2s()
    if isdefined(syst, :epsilon_low) && ~isa(syst.epsilon_low, Nothing)
        @warn "syst.epsilon_low is not used in mesti(); will be ignored."
        syst.epsilon_low = nothing
    end
    if isdefined(syst, :epsilon_high) && ~isa(syst.epsilon_high, Nothing)
        @warn "syst.epsilon_high is not used in mesti(); will be ignored."
        syst.epsilon_high = nothing        
    end
    if isdefined(syst, :zPML) && ~isa(syst.zPML, Nothing)
        @warn "syst.zPML is not used in mesti(); will be ignored."
        syst.zPML = nothing        
    end
    
    if ~use_2D_TM
        # Check boundary condition in x
        if isdefined(syst, :kx_B) 
            if isdefined(syst, :xBC) && lowercase(syst.xBC) != "bloch"
                throw(ArgumentError("When syst.kx_B is given, syst.xBC must be \"Bloch\" if specified."))
            end
            syst.xBC = "Bloch"
            # mesti_build_fdfd_matrix() uses (kx_B,ky_B,kz_B)*periodicity as the input arguments xBC, yBC, and zBC for Bloch BC
            xBC = (syst.kx_B)*(nx_Ex*syst.dx) # dimensionless
        else
            # Defaults to Dirichlet boundary condition unless syst.kx_B is given
            if ~isdefined(syst, :xBC)
                syst.xBC = "PEC"
            elseif ~(lowercase(syst.xBC) in ["bloch", "periodic", "pec", "pmc", "pecpmc", "pmcpec"])
                throw(ArgumentError("syst.xBC = \"$(syst.xBC)\" is not a supported option; type ''? mesti'' for supported options."))
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
        # mesti_build_fdfd_matrix() uses (kx_B,ky_B,kz_B)*periodicity as the input arguments xBC, yBC, and zBC for Bloch BC
        if ~use_2D_TM
            yBC = (syst.ky_B)*(ny_Ey*syst.dx) # dimensionless
        else
            yBC = (syst.ky_B)*(ny_Ex*syst.dx) # dimensionless
        end
    else
        # Defaults to Dirichlet boundary condition unless syst.ky_B is given
        if ~isdefined(syst, :yBC)
            syst.yBC = "PEC"
        elseif ~(lowercase(syst.yBC) in ["bloch", "periodic", "pec", "pmc", "pecpmc", "pmcpec"])
            throw(ArgumentError("syst.yBC = \"$(syst.yBC)\" is not a supported option; type ''? mesti'' for supported options."))
        elseif lowercase(syst.yBC) == "bloch"
            throw(ArgumentError("syst.yBC = \"Bloch\" but syst.ky_B is not given."))
        end
        yBC = syst.yBC
    end
        
    # Check boundary condition in z
    if isdefined(syst, :kz_B) 
        if isdefined(syst, :zBC) && lowercase(syst.zBC) != "bloch"
            throw(ArgumentError("When syst.kz_B is given, syst.zBC must be \"Bloch\" if specified."))
        end
        syst.zBC = "Bloch"
        # mesti_build_fdfd_matrix() uses (kx_B,ky_B,kz_B)*periodicity as the input arguments xBC, yBC, and zBC for Bloch BC
        if ~use_2D_TM
            zBC = (syst.kz_B)*(nz_Ez*syst.dx) # dimensionless
        else
            zBC = (syst.kz_B)*(nz_Ex*syst.dx) # dimensionless
        end
    else
        # Defaults to Dirichlet boundary condition unless syst.kz_B is given
        if ~isdefined(syst, :zBC)
            syst.zBC = "PEC"
        elseif ~(lowercase(syst.zBC) in ["bloch", "periodic", "pec", "pmc", "pecpmc", "pmcpec"])
            throw(ArgumentError("syst.zBC = \"$(syst.zBC)\" is not a supported option; type ''? mesti'' for supported options."))
        elseif lowercase(syst.zBC) == "bloch"
            throw(ArgumentError("syst.zBC = \"Bloch\" but syst.kz_B is not given."))
        end
        zBC = syst.zBC
    end

    # Convert BC to take care of lowercase or uppercase
    yBC = convert_BC(yBC, "y")
    zBC = convert_BC(zBC, "z")        
    
    if ~use_2D_TM
        xBC = convert_BC(xBC, "x")
        
        # Check number of grid points with the boundary conditions
        if nx_Ey != nx_Ez; throw(ArgumentError("Number of grids along x provided by syst.epsilon_yy and syst.epsilon_zz should be same.")); end
        if ny_Ex != ny_Ez; throw(ArgumentError("Number of grids along y provided by syst.epsilon_xx and syst.epsilon_zz should be same.")); end
        if nz_Ex != nz_Ey; throw(ArgumentError("Number of grids along z provided by syst.epsilon_xx and syst.epsilon_yy should be same.")); end
        check_BC_and_grid(xBC, nx_Ex, nx_Ey, nx_Ez, "x")
        check_BC_and_grid(yBC, ny_Ex, ny_Ey, ny_Ez, "y")
        check_BC_and_grid(zBC, nz_Ex, nz_Ey, nz_Ez, "z")

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
    
    # Defaults to no PML anywhere
    if ~isdefined(syst, :PML) 
        syst.PML = [PML(0)]
    elseif isa(syst.PML, PML) 
        # If user specifies PML structure instead of vector of PML strcture, transform it to vector of PML strcture. 
        syst.PML = [syst.PML]
    end
               
    # Parse the user-specified PML parameters to PML on the six sides
    # PML_list = [xPML_low, xPML_high, yPML_low, yPML_high, zPML_low, zPML_high]
    PML_list = Vector{PML}(undef, 6)
    str_sides = ["-x", "+x", "-y", "+y", "-z", "+z"]
    
    use_PML = false
    
    for ii = 1:length(syst.PML)
        PML_ii = syst.PML[ii]
              
        # Check that the user did not accidentally use options only in mesti2s()
        if isdefined(PML_ii, :npixels_spacer) && ~isa(PML_ii.npixels_spacer, Nothing)
            @warn "syst.PML[$(ii)] field \"npixels_spacer\" is not used in mesti() and will be ignored."
            PML_ii.npixels_spacer = nothing
        end
        
        # Number of PML pixels must be given
        # Other fields are optional and will be checked in mesti_build_fdfd_matrix()
        if ~isdefined(PML_ii, :npixels)
            throw(ArgumentError("syst.PML[$(ii)] must contain field \"npixels\"."))
        end

        if PML_ii.npixels != 0
            use_PML = true
        end
                
        # If PML is specified, we put it on both x, y, and z directions by default for 3D, and both y and z directions by default for 2D
        if ~isdefined(PML_ii, :direction)
            PML_ii.direction = "all"
        elseif ~(lowercase(PML_ii.direction) in ["all", "x", "y", "z"])
            throw(ArgumentError("syst.PML[$(ii)].direction = \"$(PML_ii.direction)\" is not a supported option; use \"all\", \"x\", \"y\", or \"z\".")) 
        end

        # If PML is specified, we put it on both sides by default
        if ~isdefined(PML_ii, :side)
            PML_ii.side = "both"
        elseif ~(lowercase(PML_ii.side) in ["both", "-", "+"])
            throw(ArgumentError("syst.PML[$(ii)].side = \"$(PML_ii.side)\" is not a supported option; use \"both\", \"-\", or \"+\"."))
        end
                    
        # Convert [PML_ii.direction and PML_ii.side] to a list of the PML locations
        # 1=xPML_low, 2=xPML_high, 3=yPML_low, 4=yPML_high, 5=zPML_low, 6=zPML_high 
        if lowercase(PML_ii.direction) == "all" # x & y & z
            if lowercase(PML_ii.side) == "both" # low & high
                if ~use_2D_TM
                    ind_ii = [1,2,3,4,5,6]
                else
                    # For 2D TM fields Ex(y,z), only put PML on y and z directions
                    ind_ii = [3,4,5,6]
                end
            elseif PML_ii.side == "-"
                if ~use_2D_TM
                    ind_ii = [1,3,5]
                else
                    # For 2D TM fields Ex(y,z), only put PML on y and z directions
                    ind_ii = [3,5]
                end
            else # PML_ii.side == "+"
                if ~use_2D_TM
                    ind_ii = [2,4,6]
                else
                    # For 2D TM fields Ex(y,z), only put PML on y and z directions
                    ind_ii = [4,6]
                end
            end
        elseif lowercase(PML_ii.direction) == "x"
            if lowercase(PML_ii.side) == "both" # low & high
                ind_ii = [1,2]
            elseif PML_ii.side == "-"
                ind_ii = 1
            else # PML_ii.side = "+"
                ind_ii = 2
            end                  
        elseif lowercase(PML_ii.direction) == "y"
            if lowercase(PML_ii.side) == "both" # low & high
                ind_ii = [3,4]
            elseif lowercase(PML_ii.side) == "-"
                ind_ii = 3
            else # PML_ii.side = "+"
                ind_ii = 4
            end                                    
        else # PML_ii.direction == "z"
            if lowercase(PML_ii.side) == "both" # low & high
                ind_ii = [5,6]
            elseif lowercase(PML_ii.side) == "-"
                ind_ii = 5
            else # PML_ii.side = "+"
                ind_ii = 6
            end                                                            
        end
                                                          
        # Specify PML at those locations
        for jj = 1:length(ind_ii)
            ind_side = ind_ii[jj]
            # Check that PML has not been specified at that location yet
            if isassigned(PML_list,ind_side)
                throw(ArgumentError("PML on $(str_sides[ind_side]) side is specified more than once in syst.PML."))
            end
            PML_list[ind_side] = deepcopy(PML_ii)
        end
    end
    for ii = 1:length(PML_list)
        if ~isassigned(PML_list, ii)  
            PML_list[ii] = PML(0)            
        end
        PML_list[ii].side = string(str_sides[ii][1])
        PML_list[ii].direction = string(str_sides[ii][2])        
    end
    
    # Convert to three separate PML vector for mesti_build_fdfd_matrix()
    if ~use_2D_TM
        xPML = PML_list[1:2] # [xPML_low, xPML_high]
    end
    yPML = PML_list[3:4] # [yPML_low, yPML_high]
    zPML = PML_list[5:6] # [zPML_low, zPML_high]
    
    # Use UPML by default as it produces a symmetric matrix A (unless Bloch periodic boundary is used)
    if ~isdefined(syst, :PML_type)
        syst.PML_type = "UPML"
    elseif ~(lowercase(syst.PML_type) in ["upml", "sc-pml", "scpml"])
        throw(ArgumentError("syst.PML_type = \"$(syst.PML_type)\" is not a supported option; use \"UPML\" or \"SC-PML\"."))
    end
    if lowercase(syst.PML_type) == "upml"
        use_UPML = true
    else
        use_UPML = false
    end
                
    ## Part 1.2: Check validity of the other input arguments and assign default values
    if opts == nothing
        opts = Opts()
    end
        
    # opts.return_field_profile is only used internally (but will be returned within info.opts)    
    if C == nothing
        opts.return_field_profile = true
    elseif isa(C, SparseMatrixCSC) || isa(C, Array)
        opts.return_field_profile = false
        use_transpose_B = false
    elseif C == "transpose(B)"
        opts.return_field_profile = false
        use_transpose_B = true
    else
        throw(ArgumentError("Input argument C must be a numeric matrix or a vector of Source_struct object or \"transpose(B)\", if given."))
    end

    # Check that the user did not accidentally use options only in mesti2s()
    if isdefined(opts, :is_symmetric_A) && ~isa(opts.is_symmetric_A, Nothing)
        @warn "opts.is_symmetric_A is not used in mesti(); will be ignored."
        opts.is_symmetric_A = nothing
    end                
    
    if isdefined(opts, :use_continuous_dispersion) && ~isa(opts.use_continuous_dispersion, Nothing)
        @warn "opts.use_continuous_dispersion is not used in mesti(); will be ignored."
        opts.use_continuous_dispersion = nothing        
    end    
    if isdefined(opts, :n0) && ~isa(opts.n0, Nothing) 
        @warn "opts.n0 is not used in mesti(); will be ignored."
        opts.n0 = nothing        
    end
    if isdefined(opts, :m0) && ~isa(opts.m0, Nothing) 
        @warn "opts.m0 is not used in mesti(); will be ignored."
        opts.m0 = nothing        
    end
    
    # We do not check opts.nz_low and opts.nz_high anymore, because they would be used before and after mesti2s() call mesti()
    #if isdefined(opts, :nz_low) && ~isa(opts.nz_low, Nothing)
    #    @warn "opts.nz_low is not used in mesti(); will be ignored."
    #    opts.nz_low = nothing
    #end                
    #if isdefined(opts, :nz_high) && ~isa(opts.nz_high, Nothing)
    #    @warn "opts.nz_high is not used in mesti(); will be ignored."
    #    opts.nz_high = nothing
    #end                    
    
    # Turn on verbal output by default
    if ~isdefined(opts, :verbal)
        opts.verbal = true
    elseif ~isa(opts.verbal, Bool)
        throw(ArgumentError("opts.verbal must be a boolean scalar, if given."))
    end

    # Defaults the prefactor to 1
    if ~isdefined(opts, :prefactor) || isa(opts.prefactor, Nothing)
        opts.prefactor = 1
    end

    # By default, we don't exclude the PML pixels from the returned field_profiles.
    if opts.return_field_profile
        if ~isdefined(opts, :exclude_PML_in_field_profiles) || isa(opts.exclude_PML_in_field_profiles, Nothing)
            opts.exclude_PML_in_field_profiles = false
        elseif ~isa(opts.exclude_PML_in_field_profiles, Bool)
            throw(ArgumentError("opts.exclude_PML_in_field_profiles must be a logical scalar, if given."))
        end
    else
        if isdefined(opts, :exclude_PML_in_field_profiles) && ~isa(opts.exclude_PML_in_field_profiles, Nothing)
            @warn "opts.exclude_PML_in_field_profiles is not used when output projection C is given; will be ignored."
            opts.exclude_PML_in_field_profiles = nothing
        end
    end    
    
    # By default, we don't clear syst in the caller's workspace
    if ~isdefined(opts, :clear_syst) 
        opts.clear_syst = false
    elseif ~isa(opts.clear_syst, Bool)
        throw(ArgumentError("opts.clear_syst must be a boolean scalar, if given."))
    end
       
    # By default, we don't clear B and C in the caller's workspace
    if ~isdefined(opts, :clear_BC) 
        opts.clear_BC = false
    elseif ~isa(opts.clear_BC, Bool)
        throw(ArgumentError("opts.clear_BC must be a boolean scalar, if given."))
    end
                
    # By default, we will clear internal variables to save memory; this is only used in mesti_matrix_solver!()
    # Note that mesti_matrix_solver!() defaults opts.clear_memory to false because some users that use mesti_matrix_solver!() directly may want to keep the input arguments A,B,C after calling it. But the opts.clear_memory in mesti() here only deals with the variables internal to mesti() and mesti_matrix_solver!(); it doesn't deal with the input arguments provided by the user (which are specified by opts.clear_syst and opts.clear_BC), so it is safe to default opts.clear_memory to true here.
    if ~isdefined(opts, :clear_memory) 
        opts.clear_memory = true
    elseif ~isa(opts.clear_memory, Bool)
        throw(ArgumentError("opts.clear_memory must be a boolean scalar, if given."))
    end

    # Use METIS in 3D and use AMD in 2D by default
    if ~isdefined(opts, :use_METIS) && ~use_2D_TM
        opts.use_METIS = true # If this is 3D system, we use METIS ordering by default
    elseif ~isdefined(opts, :use_METIS) && use_2D_TM
        opts.use_METIS = false # If this is 2D system, we use AMD ordering by default        
    elseif ~isa(opts.use_METIS, Bool)
        throw(ArgumentError("opts.use_METIS must be a boolean, if given."))   
    end                                
    
    # The following fields of opts will be checked/initialized in mesti_matrix_solver!():
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
    #    opts.use_BLR
    #    opts.threshold_BLR
    #    opts.icntl_36
    #    opts.icntl_38
    
    if opts.verbal
        # print basic system info if the calling function is not mesti2s()
        if stacktrace()[2].func == :mesti2s
            called_from_mesti2s = true;
            @printf("            ... ")
        else            
            called_from_mesti2s = false
            @printf("===System size=== \n")
            if ~use_2D_TM
                @printf("nx_Ex = %d, ny_Ex = %d; nz_Ex = %d \n", nx_Ex, ny_Ex, nz_Ex)
                @printf("nx_Ey = %d, ny_Ey = %d; nz_Ey = %d \n", nx_Ey, ny_Ey, nz_Ey)
                @printf("nx_Ez = %d, ny_Ez = %d; nz_Ez = %d \n", nx_Ez, ny_Ez, nz_Ez)
            else
                @printf("ny_Ex = %d; nz_Ex = %d for Ex(y,z) \n", ny_Ex, nz_Ex)
            end
            if use_PML
                @printf("%s on ", syst.PML_type)
                for ind_side = 1:6
                    if PML_list[ind_side].npixels != 0
                        @printf("%s ", str_sides[ind_side])
                    end
                end
                @printf("sides; ")
            else
                @printf("no PML; ")
            end
            if ~use_2D_TM
                @printf("xBC = %s", syst.xBC)
                if lowercase(syst.xBC) == "bloch"; @printf(" (kx_B = %.4f)", syst.kx_B); end
            end
            @printf("; yBC = %s", syst.yBC)
            if lowercase(syst.yBC) == "bloch"; @printf(" (ky_B = %.4f)", syst.ky_B); end
            @printf("; zBC = %s", syst.zBC)
            if lowercase(syst.zBC) == "bloch"; @printf(" (kz_B = %.4f)", syst.kz_B); end       
            @printf("\nBuilding B,C... ")
        end   
    end             
    t1 = time(); timing_init = t1-t0 # Initialization time
                
    ## Part 2.1: Build matrices B and C
    # matrices has fields A, B, and C to store matrices A, B, and C.
    matrices = Matrices()
    
    # Build the input matrix B from its nonzero elements specified by user
    component = ["x", "y", "z"] 
    if ~use_2D_TM 
        nx_list = [nx_Ex, nx_Ey, nx_Ez]
        ny_list = [ny_Ex, ny_Ey, ny_Ez]
        nz_list = [nz_Ex, nz_Ey, nz_Ez]
        nt_list = [nt_Ex, nt_Ey, nt_Ez]
    else
        nx_list = 1     # For computational use only
        ny_list = ny_Ex
        nz_list = nz_Ex
        nt_list = nt_Ex
    end
    
    if isa(B,Vector{Source_struct})
        if ~use_2D_TM
            if length(B) != 3
                throw(ArgumentError("The length of B must equal 3, when B is a vector of Source_struct."))
            end
            B_ii = Vector(undef, 3)
        else
            if length(B) != 1
                throw(ArgumentError("The length of B must equal 1 for 2D case."))
            end
            B_ii = Vector(undef, 1)
        end
        
        # Loop over Bx, By, and Bz for 3D
        # Only have Bx for 2D TM
        for ii = 1:length(B)
            B_struct = B[ii]   
            global M
            if isdefined(B_struct, :isempty) && ~isa(B_struct.isempty, Bool)
                throw(ArgumentError("If input argument B[$(ii)] have field \"isempty\" and it must be a boolean scalar."))
            elseif isdefined(B_struct, :isempty) && B_struct.isempty && (isdefined(B_struct, :data) || isdefined(B_struct, :ind) || isdefined(B_struct, :pos))
                throw(ArgumentError("If input argument B[$(ii)].isempty = true and B[$(ii)] should not have the field \"pos\", \"ind\", or \"data\"."))
            elseif ~isdefined(B_struct, :isempty)
                 B_struct.isempty = false
            end
            if ~B_struct.isempty
                if ~isdefined(B_struct, :data)
                    throw(ArgumentError("Input argument B[$(ii)] must have field \"data\" when B[$(ii)] is a Source_struct."))
                elseif ~isdefined(B_struct, :pos) && ~isdefined(B_struct, :ind)
                    throw(ArgumentError("Input argument B[$(ii)] must have field \"pos\" or \"ind\" when B[$(ii)] is a Source_struct."))
                elseif isdefined(B_struct, :pos) && isdefined(B_struct, :ind)
                    throw(ArgumentError("Input argument B[$(ii)] cannot have both field \"pos\" and field \"ind\" when B[$(ii)] is a Source_struct."))
                end     
                if isdefined(B_struct, :pos)
                    if length(B_struct.pos) != length(B_struct.data)
                        throw(ArgumentError("The length of B[$(ii)].pos must equal the length of B[$(ii)].data."))   
                    end
                    # Loop over different source locations
                    for jj = 1:length(B_struct.pos)
                        pos = B_struct.pos[jj]
                        if use_2D_TM 
                            if use_2D_TM && ~(length(pos) == 4 && minimum(pos)>0)
                                throw(ArgumentError("B[$(ii)].pos[$(jj)] must be a positive integer vector with 4 elements for 2D case."))
                            elseif ~use_2D_TM && ~(length(pos)==6 && minimum(pos)>0)
                                throw(ArgumentError("B[$(ii)].pos[$(jj)] must be a positive integer vector with 6 elements for 3D case."))
                            end
                        end
                    end
                else
                    if length(B_struct.ind) != length(B_struct.data)
                        throw(ArgumentError("The length of B[$(ii)].ind must equal the length of B[$(ii)].data."))   
                    end                
                end

                # We first pick the most efficient way to build matrix B.
                # If all of the following are satisfied: (1) length(B_struct) is small, (2) B_struct[ii].pos is used, and (3) the cuboid specified by B_struct[ii].pos is a single vertical slice in x-y plane or if it spans the full cross area of nx_list[ii]*ny_list[ii], then we will stack reshaped B_struct(ii).data with zeros. This avoids the overhead of building B with index-value pairs.
                # If any of the above is not satisfied, we will build B with index-value pairs.
                use_iv_pairs = false
                if length(B_struct.data) > 10
                    use_iv_pairs = true
                elseif ~isdefined(B_struct, :pos)
                    use_iv_pairs = true
                else
                    for jj = 1:length(B_struct.pos)
                        pos = B_struct.pos[jj]
                        if (~use_2D_TM && ~((pos[4] == nx_list[ii] && pos[5] == ny_list[ii]) || pos[6] == 1)) || (use_2D_TM && ~(pos[3] == ny_list[ii] || pos[4] == 1))
                            use_iv_pairs = true
                        end
                    end
                end
                
                if use_iv_pairs
                    # Construct matrix B_ii from the complete set of index-value pairs
                    N_tot = 0; # total number of nonzero elements in B_ii
                    for jj = 1:length(B_struct.data)
                        N_tot = N_tot + length(B_struct.data[jj])
                    end
                    ind_list = zeros(N_tot)
                    a_list = zeros(N_tot)
                    v_list = zeros(N_tot)
                    N = 0
                    M = 0
                else
                    # Build matrix B incrementally
                    B_ii[ii] = spzeros(nt_list[ii], 0)
                end
                
                # Loop over different positions of the input source
                for jj = 1:length(B_struct.data)
                    data = B_struct.data[jj]
                    if isdefined(B_struct, :pos)
                        # B_struct(ii).pos specifies a cuboid inside (nx, ny, nz); (n1, m1, l1) and (n2, m2, l2) are its two diagonal corners
                        # For 2D TM field Ex(y,z), we only consider y, z directions; so we set n1 = n2 = 1 (pos[4] = 1), and the two diagonal corners of this cuboid becomes (1, m1, l1) and (1, m2, l2)
                        pos = B_struct.pos[jj]
                        if use_2D_TM && length(pos) == 4
                            pos = [1, pos[1], pos[2], 1, pos[3], pos[4]]
                        end
                        n1 = pos[1] # first index in x
                        m1 = pos[2] # first index in y
                        l1 = pos[3] # first index in z           
                        n2 = n1 + pos[4] - 1 # last index in x
                        m2 = m1 + pos[5] - 1 # last index in y
                        l2 = l1 + pos[6] - 1 # last index in z
                        nxyz_data = pos[4]*pos[5]*pos[6]; # number of elements in this cuboid
                        if n1 > nx_list[ii]
                            if ii == 1; throw(ArgumentError("B[$ii].pos[$jj][1] = $(n1) exceeds nx_E$(component[ii]) = $(nx_Ex).")); end
                            if ii == 2; throw(ArgumentError("B[$ii].pos[$jj][1] = $(n1) exceeds nx_E$(component[ii]) = $(nx_Ey).")); end
                            if ii == 3; throw(ArgumentError("B[$ii].pos[$jj][1] = $(n1) exceeds nx_E$(component[ii]) = $(nx_Ez).")); end
                        elseif m1 > ny_list[ii]
                            temp = 2 - use_2D_TM
                            if ii == 1; throw(ArgumentError("B[$ii].pos[$jj][$temp] = $(m1) exceeds ny_E$(component[ii]) = $(ny_Ex).")); end
                            if ii == 2; throw(ArgumentError("B[$ii].pos[$jj][$temp] = $(m1) exceeds ny_E$(component[ii]) = $(ny_Ey).")); end
                            if ii == 3; throw(ArgumentError("B[$ii].pos[$jj][$temp] = $(m1) exceeds ny_E$(component[ii]) = $(ny_Ez).")); end
                        elseif l1 > nz_list[ii]
                            temp = 3 - use_2D_TM
                            if ii == 1; throw(ArgumentError("B[$ii].pos[$jj][$temp] = $(l1) exceeds nz_E$(component[ii]) = $(nz_Ex).")); end
                            if ii == 2; throw(ArgumentError("B[$ii].pos[$jj][$temp] = $(l1) exceeds nz_E$(component[ii]) = $(nz_Ey).")); end
                            if ii == 3; throw(ArgumentError("B[$ii].pos[$jj][$temp] = $(l1) exceeds nz_E$(component[ii]) = $(nz_Ez).")); end
                        elseif n2 > nx_list[ii]
                            if ii == 1; throw(ArgumentError("B[$ii].pos[$jj][1] + B[$ii].pos[$jj][4] - 1 = $(n2) exceeds nx_E$(component[ii]) = $(nx_Ex).")); end
                            if ii == 2; throw(ArgumentError("B[$ii].pos[$jj][1] + B[$ii].pos[$jj][4] - 1 = $(n2) exceeds nx_E$(component[ii]) = $(nx_Ey).")); end
                            if ii == 3; throw(ArgumentError("B[$ii].pos[$jj][1] + B[$ii].pos[$jj][4] - 1 = $(n2) exceeds nx_E$(component[ii]) = $(nx_Ez).")); end
                        elseif m2 > ny_list[ii]
                            temp = 2 - use_2D_TM
                            if ii == 1; throw(ArgumentError("B[$ii].pos[$jj][$temp] + B[$ii].pos[$jj][$(temp*2+1)] - 1 = $(m2) exceeds ny_E$(component[ii]) = $(ny_Ex).")); end
                            if ii == 2; throw(ArgumentError("B[$ii].pos[$jj][$temp] + B[$ii].pos[$jj][$(temp*2+1)] - 1 = $(m2) exceeds ny_E$(component[ii]) = $(ny_Ey).")); end
                            if ii == 3; throw(ArgumentError("B[$ii].pos[$jj][$temp] + B[$ii].pos[$jj][$(temp*2+1)] - 1 = $(m2) exceeds ny_E$(component[ii]) = $(ny_Ez).")); end
                        elseif l2 > nz_list[ii]
                            temp = 3 - use_2D_TM                            
                            if ii == 1; throw(ArgumentError("B[$ii].pos[$jj][$temp] + B[$ii].pos[$jj][$(temp*2)] - 1 = $(l2) exceeds nz_E$(component[ii]) = $(nz_Ex).")); end
                            if ii == 2; throw(ArgumentError("B[$ii].pos[$jj][$temp] + B[$ii].pos[$jj][$(temp*2)] - 1 = $(l2) exceeds nz_E$(component[ii]) = $(nz_Ey).")); end
                            if ii == 3; throw(ArgumentError("B[$ii].pos[$jj][$temp] + B[$ii].pos[$jj][$(temp*2)] - 1 = $(l2) exceeds nz_E$(component[ii]) = $(nz_Ez).")); end
                        elseif (~use_2D_TM && ~(ndims(data) == 2 || ndims(data) == 4)) || (use_2D_TM && ~(ndims(data) == 2 || ndims(data) == 3))
                            throw(ArgumentError("B[$ii].data[$jj] must be a 2D or 4D numeric array for 3D systems, or a 2D or 3D numeric array for 2D systems, when B[$ii].pos[$jj] is given."))
                        end

                        if ndims(data) == 2
                            M_ii = size(data, 2) # number of inputs
                        else
                            if ~use_2D_TM
                                M_ii = size(data, 4) # number of inputs
                            else
                                M_ii = size(data, 3) # number of inputs
                            end
                        end
                        if use_iv_pairs
                            # convert to linear indices
                            n_list = repeat((n1:n2), 1, pos(5), pos(6))
                            m_list = repeat(transpose(m1:m2), pos(4), 1, pos(6))
                            l_list = repeat(reshape((l1:l2),1,1,:), pos(4), pos(5), 1)
                            #ind = LinearIndices((nx_list[ii], ny_list[ii], nz_list[ii]))[CartesianIndex.(n_list, m_list, l_list)]
                            ind = Base._sub2ind((nx_list[ii], ny_list[ii], nz_list[ii]), n_list, m_list, l_list)
                        end
                    else
                        # B_struct.ind[jj] specifies linear indices of [n, m, l]
                        ind = B_struct.ind[jj]
                        if ~(minimum(ind)>0 && maximum(ind)<= nt_list[ii] && length(ind)==length(unique(ind)))
                            throw(ArgumentError("B[$ii][$jj].ind must be a vector with no repeated elements, where every element is an integer between 1 and nx_E$(component[ii])*ny_E$(component[ii])*nz_E$(component[ii]) = $(nt_list[ii])."))
                        elseif size(data, 1) != length(ind)
                            throw(ArgumentError("size(B[$ii][$jj].data, 1) = %d does not match length(B[$ii][$jj].ind) = $(length(ind))"))
                        end
                        nxyz_data = length(ind) # number of linear indices
                        M_ii = size(data, 2) # number of inputs
                    end
                    if use_iv_pairs
                        # Build index-value pairs: (ind_list, a_list, v_list)
                        N_ii = nxyz_data*M_ii; # number of nonzero elements in the jj-th part of matrix B_ii
                        ind_temp = N .+ (1:N_ii)
                        ind_list[ind_temp] = repeat(ind, M_ii) # spatial index
                        a_list[ind_temp] = reshape(repeat(M.+(1:M_ii), nxyz_data), N_ii) # input index
                        v_list[ind_temp] = reshape(data, N_ii)
                        N = N + N_ii # number of nonzero elements in matrix B up to the ii-th part
                        M = M + M_ii # number of columns in matrix B up to the ii-th part
                    else                
                        # If the cuboid of B_struct.pos is a single vertical slice or if it spans the full height of ny, we can stack reshaped B_struct[ii].data with zeros to avoid the overhead of building B with index-value pairs. But then B must be built incrementally, which is slow when length(B_struct) is large.
                        #nxyz_before = LinearIndices((nx_list[ii], ny_list[ii], nz_list[ii]))[CartesianIndex(n1, m1, l1)] - 1
                        nxyz_before = Base._sub2ind((nx_list[ii], ny_list[ii], nz_list[ii]), n1, m1, l1) - 1
                        #nxyz_after = nt_list[ii] - LinearIndices((nx_list[ii], ny_list[ii], nz_list[ii]))[CartesianIndex(n2, m2, l2)]
                        nxyz_after = nt_list[ii] - Base._sub2ind((nx_list[ii], ny_list[ii], nz_list[ii]), n2, m2, l2)
                        B_ii[ii] = [B_ii[ii] [spzeros(nxyz_before, M_ii); reshape(data, nxyz_data, M_ii); spzeros(nxyz_after, M_ii)]]
                    end                
                end
                if use_iv_pairs
                    B_ii[ii] = sparse(ind_list, a_list, v_list, nt_list[ii], M)
                    if opts.clear_BC
                        data = nothing; B_struct = nothing; m_list = nothing; n_list = nothing; l_list = nothing
                        ind = nothing; ind_list = nothing; a_list = nothing; v_list = nothing; ind_temp = nothing
                        GC.gc()
                    end
                else
                    M = size(B_ii[ii],2)
                end
                if opts.clear_BC
                    data = nothing; B_struct = nothing
                    GC.gc()
                end
            end
        end
        if ~isassigned(B_ii)
            throw(ArgumentError("B cannot be no source for all x, y, and z components."))            
        end
        for ii = 1:length(B_ii)
            if ~isassigned(B_ii, ii)  
                B_ii[ii] = spzeros(nt_list[ii], M)
            end
        end
        if ~use_2D_TM
            if ~(size(B_ii[1],2) == size(B_ii[2],2) && size(B_ii[1],2) == size(B_ii[3],2))
                throw(ArgumentError("B cannot have different number of inputs for x, y, or z components."))                     
            end
            matrices.B = [B_ii[1]; B_ii[2]; B_ii[3]] # B=[Bx; By; Bz]
        else
            matrices.B = B_ii[1] # B = Bx for 2D TM field
        end
        if opts.clear_BC
            B_ii = nothing
            GC.gc()
        end
    end
    
    if isa(C,Vector{Source_struct})
        if ~use_2D_TM
            if length(C) != 3
                throw(ArgumentError("The length of C must equal 3, when C is a vector of Source_struct."))
            end
            C_ii = Vector(undef, 3)
        else
            if length(C) != 1
                throw(ArgumentError("The length of C must equal 1 for 2D case."))
            end
            C_ii = Vector(undef, 1)
        end
       
        # Loop over Cx, Cy, and Cz
        # Only have Cx for 2D TM
        for ii = 1:length(C)
            C_struct = C[ii]   
            global M
            if isdefined(C_struct, :isempty) && ~isa(C_struct.isempty, Bool)
                throw(ArgumentError("If input argument C[$(ii)] have field \"isempty\" and it must be a boolean scalar."))
            elseif isdefined(C_struct, :isempty) && C_struct.isempty && (isdefined(C_struct, :data) || isdefined(C_struct, :ind) || isdefined(C_struct, :pos))
                throw(ArgumentError("If input argument C[$(ii)].isempty = true and C[$(ii)] should not have the field \"pos\", \"ind\", or \"data\"."))
            elseif ~isdefined(C_struct, :isempty)
                 C_struct.isempty = false
            end
            if ~C_struct.isempty            
                if ~isdefined(C_struct, :data)
                    throw(ArgumentError("Input argument C[$(ii)] must have field \"data\" when C[$(ii)] is a Source_struct."))
                elseif ~isdefined(C_struct, :pos) && ~isdefined(C_struct, :ind)
                    throw(ArgumentError("Input argument C[$(ii)] must have field \"pos\" or \"ind\" when C[$(ii)] is a Source_struct."))
                elseif isdefined(C_struct, :pos) && isdefined(C_struct, :ind)
                    throw(ArgumentError("Input argument C[$(ii)] cannot have both field \"pos\" and field \"ind\" when C[$(ii)] is a Source_struct."))
                end     
                if isdefined(C_struct, :pos)
                    if length(C_struct.pos) != length(C_struct.data)
                        throw(ArgumentError("The length of C[$(ii)].pos must equal the length of C[$(ii)].data."))   
                    end
                    # Loop over different source locations                
                    for jj = 1:length(C_struct.pos)
                        pos = C_struct.pos[jj]
                        if use_2D_TM 
                            if use_2D_TM && ~(length(pos) == 4 && minimum(pos)>0)
                                throw(ArgumentError("C[$(ii)].pos[$(jj)] must be a positive integer vector with 4 elements for 2D case."))
                            elseif ~use_2D_TM && ~(length(pos)==6 && minimum(pos)>0)
                                throw(ArgumentError("C[$(ii)].pos[$(jj)] must be a positive integer vector with 6 elements for 3D case."))
                            end
                        end
                    end
                else
                    if length(C_struct.ind) != length(C_struct.data)
                        throw(ArgumentError("The length of C[$(ii)].ind must equal the length of C[$(ii)].data."))   
                    end                
                end

                # We first pick the most efficient way to build matrix C.
                # If all of the following are satisfied: (1) length(C_struct) is small, (2) C_struct[ii].pos is used, and (3) the cuboid specified by C_struct[ii].pos is a single vertical slice in x-y plane or if it spans the full cross area of nx_list[ii]*ny_list[ii], then we will stack reshaped C_struct(ii).data with zeros. This avoids the overhead of building C with index-value pairs.
                # If any of the above is not satisfied, we will build C with index-value pairs.
                use_iv_pairs = false
                if length(C_struct.data) > 10
                    use_iv_pairs = true
                elseif ~isdefined(C_struct, :pos)
                    use_iv_pairs = true
                else
                    for jj = 1:length(C_struct.pos)
                        pos = C_struct.pos[jj]
                        if (~use_2D_TM && ~((pos[4] == nx_list[ii] && pos[5] == ny_list[ii]) || pos[6] == 1)) || (use_2D_TM && ~(pos[3] == ny_list[ii] || pos[4] == 1))
                            use_iv_pairs = true
                        end
                    end
                end
                if use_iv_pairs
                    # Construct matrix C_ii from the complete set of index-value pairs
                    N_tot = 0; # total number of nonzero elements in C_ii
                    for jj = 1:length(C_struct.data)
                        N_tot = N_tot + length(C_struct.data[jj])
                    end
                    ind_list = zeros(N_tot)
                    a_list = zeros(N_tot)
                    v_list = zeros(N_tot)
                    N = 0
                    M = 0
                else
                    # Build matrix C incrementally
                    C_ii[ii] = spzeros(0, nt_list[ii])                    
                end
                # Loop over different positions of the output projection
                for jj = 1:length(C_struct.data)
                    data = C_struct.data[jj]
                    if isdefined(C_struct, :pos)
                        # C_struct(ii).pos specifies a cuboid inside (nx, ny, nz); (n1, m1, l1) and (n2, m2, l2) are its two diagonal corners
                        # For 2D TM field Ex(y,z), we only consider y, z directions; so we set n1 = n2 = 1 (pos[4] = 1), and the two diagonal corners of this cuboid becomes (1, m1, l1) and (1, m2, l2)
                        pos = C_struct.pos[jj]
                        if use_2D_TM && length(pos) == 4
                            pos = [1, pos[1], pos[2], 1, pos[3], pos[4]]
                        end
                        n1 = pos[1] # first index in x
                        m1 = pos[2] # first index in y
                        l1 = pos[3] # first index in z           
                        n2 = n1 + pos[4] - 1 # last index in x
                        m2 = m1 + pos[5] - 1 # last index in y
                        l2 = l1 + pos[6] - 1 # last index in z
                        nxyz_data = pos[4]*pos[5]*pos[6]; # number of elements in this cuboid
                        if n1 > nx_list[ii]
                            if ii == 1; throw(ArgumentError("C[$ii].pos[$jj][1] = $(n1) exceeds nx_E$(component[ii]) = $(nx_Ex).")); end
                            if ii == 2; throw(ArgumentError("C[$ii].pos[$jj][1] = $(n1) exceeds nx_E$(component[ii]) = $(nx_Ey).")); end
                            if ii == 3; throw(ArgumentError("C[$ii].pos[$jj][1] = $(n1) exceeds nx_E$(component[ii]) = $(nx_Ez).")); end
                        elseif m1 > ny_list[ii]
                            temp = 2 - use_2D_TM
                            if ii == 1; throw(ArgumentError("C[$ii].pos[$jj][$temp] = $(m1) exceeds ny_E$(component[ii]) = $(ny_Ex).")); end
                            if ii == 2; throw(ArgumentError("C[$ii].pos[$jj][$temp] = $(m1) exceeds ny_E$(component[ii]) = $(ny_Ey).")); end
                            if ii == 3; throw(ArgumentError("C[$ii].pos[$jj][$temp] = $(m1) exceeds ny_E$(component[ii]) = $(ny_Ez).")); end
                        elseif l1 > nz_list[ii]
                            temp = 3 - use_2D_TM
                            if ii == 1; throw(ArgumentError("C[$ii].pos[$jj][$temp] = $(l1) exceeds nz_E$(component[ii]) = $(nz_Ex).")); end
                            if ii == 2; throw(ArgumentError("C[$ii].pos[$jj][$temp] = $(l1) exceeds nz_E$(component[ii]) = $(nz_Ey).")); end
                            if ii == 3; throw(ArgumentError("C[$ii].pos[$jj][$temp] = $(l1) exceeds nz_E$(component[ii]) = $(nz_Ez).")); end
                        elseif n2 > nx_list[ii]
                            if ii == 1; throw(ArgumentError("C[$ii].pos[$jj][1] + C[$ii].pos[$jj][4] - 1 = $(n2) exceeds nx_E$(component[ii]) = $(nx_Ex).")); end
                            if ii == 2; throw(ArgumentError("C[$ii].pos[$jj][1] + C[$ii].pos[$jj][4] - 1 = $(n2) exceeds nx_E$(component[ii]) = $(nx_Ey).")); end
                            if ii == 3; throw(ArgumentError("C[$ii].pos[$jj][1] + C[$ii].pos[$jj][4] - 1 = $(n2) exceeds nx_E$(component[ii]) = $(nx_Ez).")); end
                        elseif m2 > ny_list[ii]
                            temp = 2 - use_2D_TM
                            if ii == 1; throw(ArgumentError("C[$ii].pos[$jj][$temp] + C[$ii].pos[$jj][$(temp*2+1)] - 1 = $(m2) exceeds ny_E$(component[ii]) = $(ny_Ex).")); end
                            if ii == 2; throw(ArgumentError("C[$ii].pos[$jj][$temp] + C[$ii].pos[$jj][$(temp*2+1)] - 1 = $(m2) exceeds ny_E$(component[ii]) = $(ny_Ey).")); end
                            if ii == 3; throw(ArgumentError("C[$ii].pos[$jj][$temp] + C[$ii].pos[$jj][$(temp*2+1)] - 1 = $(m2) exceeds ny_E$(component[ii]) = $(ny_Ez).")); end
                        elseif l2 > nz_list[ii]
                            temp = 3 - use_2D_TM                            
                            if ii == 1; throw(ArgumentError("C[$ii].pos[$jj][$temp] + C[$ii].pos[$jj][$(temp*2)] - 1 = $(l2) exceeds nz_E$(component[ii]) = $(nz_Ex).")); end
                            if ii == 2; throw(ArgumentError("C[$ii].pos[$jj][$temp] + C[$ii].pos[$jj][$(temp*2)] - 1 = $(l2) exceeds nz_E$(component[ii]) = $(nz_Ey).")); end
                            if ii == 3; throw(ArgumentError("C[$ii].pos[$jj][$temp] + C[$ii].pos[$jj][$(temp*2)] - 1 = $(l2) exceeds nz_E$(component[ii]) = $(nz_Ez).")); end
                        elseif (~use_2D_TM && ~(ndims(data) == 2 || ndims(data) == 4)) || (use_2D_TM && ~(ndims(data) == 2 || ndims(data) == 3))
                            throw(ArgumentError("C[$ii].data[$jj] must be a 2D or 4D numeric array for 3D systems, or a 2D or 3D numeric array for 2D systems, when C[$ii].pos[$jj] is given."))
                        end

                        if ndims(data) == 2
                            M_ii = size(data, 2) # number of inputs
                        else
                            if ~use_2D_TM
                                M_ii = size(data, 4) # number of inputs
                            else
                                M_ii = size(data, 3) # number of inputs
                            end
                        end
                        if use_iv_pairs
                            # convert to linear indices
                            n_list = repeat((n1:n2), 1, pos[5], pos[6])
                            m_list = repeat(transpose(m1:m2), pos[4], 1, pos[6])
                            l_list = repeat(reshape((l1:l2),1,1,:), pos[4], pos[5], 1)
                            #ind = LinearIndices((nx_list[ii], ny_list[ii], nz_list[ii]))[CartesianIndex.(n_list, m_list, l_list)]
                            ind = Base._sub2ind((nx_list[ii], ny_list[ii], nz_list[ii]), n_list, m_list, l_list)
                        end
                    else
                        # C_struct.ind[jj] specifies linear indices of [n, m, l]
                        ind = C_struct.ind[jj]
                        if ~(minimum(ind)>0 && maximum(ind)<= nt_list[ii] && length(ind)==length(unique(ind)))
                            throw(ArgumentError("C[$ii][$jj].ind must be a vector with no repeated elements, where every element is an integer between 1 and nx_E$(component[ii])*ny_E$(component[ii])*nz_E$(component[ii]) = $(nt_list[ii])."))
                        elseif size(data, 1) != length(ind)
                            throw(ArgumentError("size(C[$ii][$jj].data, 1) = %d does not match length(C[$ii][$jj].ind) = $(length(ind))"))
                        end
                        nxyz_data = length(ind) # number of linear indices
                        M_ii = size(data, 2) # number of inputs
                    end
                    if use_iv_pairs
                        # Build index-value pairs: (a_list, ind_list, v_list)
                        N_ii = nxyz_data*M_ii; # number of nonzero elements in the jj-th part of matrix C_ii
                        ind_temp = N .+ (1:N_ii)
                        ind_list[ind_temp] = repeat(ind, M_ii) # spatial index
                        a_list[ind_temp] = reshape(repeat(M.+(1:M_ii), nxyz_data), N_ii) # input index
                        v_list[ind_temp] = reshape(data, N_ii)
                        N = N + N_ii # number of nonzero elements in matrix C up to the ii-th part
                        M = M + M_ii # number of columns in matrix C up to the ii-th part
                    else                
                        # If the cuboid of C_struct.pos is a single vertical slice or if it spans the full height of ny, we can stack reshaped C_struct[jj].data with zeros to avoid the overhead of building C with index-value pairs. But then C must be built incrementally, which is slow when length(C_struct) is large.
                        #nxyz_before = LinearIndices((nx_list[ii], ny_list[ii], nz_list[ii]))[CartesianIndex(n1, m1, l1)] - 1
                        nxyz_before = Base._sub2ind((nx_list[ii], ny_list[ii], nz_list[ii]), n1, m1, l1) - 1
                        #nxyz_after = nt_list[ii] - LinearIndices((nx_list[ii], ny_list[ii], nz_list[ii]))[CartesianIndex(n2, m2, l2)]
                        nxyz_after = nt_list[ii] - Base._sub2ind((nx_list[ii], ny_list[ii], nz_list[ii]), n2, m2, l2)
                        #C_ii[ii] = [C_ii[ii]; [spzeros(M_ii, nxyz_before) transpose(reshape(data, nxyz_data, M_ii)) spzeros(M_ii, nxyz_after)]]
                        C_ii[ii] = permutedims([transpose(C_ii[ii]) [spzeros(nxyz_before, M_ii); reshape(data, nxyz_data, M_ii); spzeros(nxyz_after, M_ii)]],(2,1))
                    end                
                end
                if use_iv_pairs
                    C_ii[ii] = sparse(a_list, ind_list, v_list, M, nt_list[ii])
                    if opts.clear_BC 
                        data = nothing; C_struct = nothing; m_list = nothing; n_list = nothing; l_list = nothing
                        ind = nothing; ind_list = nothing; a_list = nothing; v_list = nothing; ind_temp = nothing
                        GC.gc()
                    end
                else
                    M = size(C_ii[ii],1)
                end
                if opts.clear_BC
                    data = nothing; C_struct = nothing
                    GC.gc()
                end
            end
        end      
        if ~isassigned(C_ii)
            throw(ArgumentError("C cannot be no source for all x, y, and z components."))            
        end
        for ii = 1:length(C_ii)
            if ~isassigned(C_ii, ii)  
                C_ii[ii] = spzeros(M, nt_list[ii])
            end
        end
        if ~use_2D_TM
            if size(C_ii[1],1) != size(C_ii[2],1) || size(C_ii[1],1) != size(C_ii[3],1)
                throw(ArgumentError("C cannot have different number of inputs for x, y, or z components."))                     
            end
            matrices.C = [C_ii[1] C_ii[2] C_ii[3]] # C = [Cx Cy Cz]
        else
            matrices.C = C_ii[1] # C = Cx for 2D TM field
        end
        if opts.clear_BC
            C_ii = nothing
            GC.gc()
        end
    elseif ~opts.return_field_profile && use_transpose_B
        matrices.C = "transpose(B)"
    end
    
    # Check matrix sizes
    (sz_B_1, sz_B_2) = size(matrices.B)
    # nt_list = [nt_Ex, nt_Ey, nt_Ez] for 3D
    # nt_list = nt_Ex for 2D TM
    if sz_B_1 != sum(nt_list)
        if ~use_2D_TM
            throw(ArgumentError("size(matrices.B,1) must equal nt_Ex+nt_Ey+nt_Ez for 3D systems; size(matrices.B,1) = $(sz_B_1), nt_Ex+nt_Ey+nt_Ez = $(nt_Ex+nt_Ey+nt_Ez).")) 
        else
            throw(ArgumentError("size(matrices.B,1) must equal nt_Ex for 2D systems; size(matrices.B,1) = $(sz_B_1), nt_Ex = $(nt_Ex)."))
        end
    end
    if ~opts.return_field_profile && ~use_transpose_B
        (sz_C_1, sz_C_2) = size(matrices.C)
        # nt_list = [nt_Ex, nt_Ey, nt_Ez] for 3D
        # nt_list = nt_Ex for 2D TM
        if sz_C_2 != sum(nt_list)
            if ~use_2D_TM
                throw(ArgumentError("size(matrices.C,2) must equal nt_Ex+nt_Ey+nt_Ez for 3D systems; size(matrices.C,2) = $(sz_C_2), nt_Ex+nt_Ey+nt_Ez = $(nt_Ex+nt_Ey+nt_Ez)."))
            else
                throw(ArgumentError("size(matrices.C,2) must equal nt_Ex for 2D systems; size(matrices.C,2) = $(sz_C_2), nt_Ex = $(nt_Ex)."))
            end
        end
    elseif ~isa(D, Nothing)
        sz_C_1 = sz_B_2 # to-be used before subtracting D, taking C = transpose(B)
    end
                        
    t2 = time() 
    timing_build_BC = t2 - t1
    if opts.verbal
        @printf("elapsed time: %7.3f secs\nBuilding A  ... ", timing_build_BC)
    end

    ## Part 2.2: Build matrix A by calling mesti_build_fdfd_matrix()
    if (sz_B_2 == 0 || (~opts.return_field_profile && ~use_transpose_B && sz_C_1 == 0)) && ~(isdefined(opts, :store_ordering) && opts.store_ordering)
        # No need to build A if length(S) = 0 and we don't need to keep the ordering
        matrices.A = spzeros(sum(nt_list), sum(nt_list))
        is_symmetric_A = true
    else
        # Build the finite-difference differential operator
        k0dx = (2*pi/syst.wavelength)*(syst.dx)
        if ~use_2D_TM
            if include_off_diagonal_epsilon
                (matrices.A, is_symmetric_A, xPML, yPML, zPML) = mesti_build_fdfd_matrix(syst.epsilon_xx, syst.epsilon_xy, syst.epsilon_xz, syst.epsilon_yx, syst.epsilon_yy, syst.epsilon_yz, syst.epsilon_zx, syst.epsilon_zy, syst.epsilon_zz, k0dx, xBC, yBC, zBC, xPML, yPML, zPML, use_UPML)
            else
                (matrices.A, is_symmetric_A, xPML, yPML, zPML) = mesti_build_fdfd_matrix(syst.epsilon_xx, syst.epsilon_yy, syst.epsilon_zz, k0dx, xBC, yBC, zBC, xPML, yPML, zPML, use_UPML)
            end
        else 
            # 2D TM
            (matrices.A, is_symmetric_A, yPML, zPML) = mesti_build_fdfd_matrix(syst.epsilon_xx, k0dx, yBC, zBC, yPML, zPML, use_UPML)
        end
    end
    if opts.clear_memory
        GC.gc()    
    end
    if opts.clear_syst
        syst = nothing
        GC.gc()
    end
    opts.is_symmetric_A = is_symmetric_A
                        
    t3 = time() 
    timing_build_A = t3 - t2
    if opts.verbal
        @printf("elapsed time: %7.3f secs\n", timing_build_A)
    end
                        
    ## Part 3: Compute C*inv(A)*B or inv(A)*B by calling mesti_matrix_solver!()
    # This is where most of the computation is done
    # Note that A, B, C in this workspace may be cleared after calling mesti_matrix_solver!() if opts.clear_memory = true (which is the default), so we should no longer use A, B, C below
    (S, info) = mesti_matrix_solver!(matrices, opts)
    if opts.clear_memory
        GC.gc()
    end
    
    info.timing_init = info.timing_init + timing_init
    info.timing_build = info.timing_build + timing_build_BC + timing_build_A # combine with build time for K
    if ~use_2D_TM
        if xPML[1].npixels != 0 && xPML[2].npixels != 0        
            info.xPML = xPML
        end
    end
    if yPML[1].npixels != 0 && yPML[2].npixels != 0        
        info.yPML = yPML
    end
    if zPML[1].npixels != 0 && zPML[2].npixels != 0        
        info.zPML = zPML
    end

    t1 = time()

    if ~info.opts.analysis_only
        # Include the prefactor
        S = (opts.prefactor)*S
        
        # subtract D
        if ~isa(D, Nothing)
            if opts.return_field_profile 
                throw(ArgumentError("Input argument D must be nothing when C is nothing."))
            end
            (sz_D_1, sz_D_2) = size(D)
            if sz_D_1 !=sz_C_1 
                throw(ArgumentError("size(D,1) must equal size(C,1); size(D,1) = $(sz_D_1), size(C,1) = $(sz_C_1)."))
            end
            if sz_D_2 !=sz_B_2 
                throw(ArgumentError("size(D,2) must equal size(B,2); size(D,2) = $(sz_D_2), size(B,2) = $(sz_B_2).")) 
            end
            S = S - D
        end
    end
    t2 = time()
    info.timing_solve = info.timing_solve + t2 - t1 # Add the little bit of post-processing time

    if ~opts.return_field_profile
        info.timing_total = t2 - t0
        if opts.verbal && ~called_from_mesti2s 
            @printf("          Total elapsed time: %7.3f secs\n", info.timing_total) 
        end   
    end
    
    return S, info
end

# The following are mesti functions to take different number of input arguments, but all of them will
# call the mesti main function.


# When syst and B are specified; return the components of the field profiles.
function mesti(syst::Syst, B::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2},Vector{Source_struct}})
    opts = Opts() 
    
    # Check if 2D TM fields are required
    if ndims(syst.epsilon_xx) == 2
        use_2D_TM = true
    else
        use_2D_TM = false
    end
    
    if ~use_2D_TM
        (Ex, Ey, Ez, info) = mesti(syst, B, opts)
        return Ex, Ey, Ez, info
    else
        (Ex, info) = mesti(syst, B, opts)
        return Ex, info
    end
end

# When syst, B, and opts are specified; return the components of the field profiles.
function mesti(syst::Syst, B::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2},Vector{Source_struct}}, opts::Opts)
    
    (S, info) = mesti(syst, B, nothing, nothing, opts)
    
    # Check if 2D TM fields are required
    if ndims(syst.epsilon_xx) == 2
        use_2D_TM = true
    else
        use_2D_TM = false
    end

    t1 = time()

    # Number of inputs
    M = size(S,2)
    
    if ~use_2D_TM
        # Number of grid points in x, y, and z for Ex, Ey, and Ez
        (nx_Ex, ny_Ex, nz_Ex) = size(syst.epsilon_xx)
        (nx_Ey, ny_Ey, nz_Ey) = size(syst.epsilon_yy)
        (nx_Ez, ny_Ez, nz_Ez) = size(syst.epsilon_zz)  

        # Total number of grid points for Ex, Ey, and Ez    
        nt_Ex = nx_Ex*ny_Ex*nz_Ex
        nt_Ey = nx_Ey*ny_Ey*nz_Ey
        nt_Ez = nx_Ez*ny_Ez*nz_Ez        
    else
        # Number of grid points in y and z for Ex
        (ny_Ex, nz_Ex) = size(syst.epsilon_xx)
        
        # Total number of grid points for Ex
        nt_Ex = ny_Ex*nz_Ex
    end
    
    if info.opts.exclude_PML_in_field_profiles
        # Exclude the PML pixels from the returned field profiles
        if ~use_2D_TM
            n_start = 1; n_Ex_end = nx; n_Ey_end = nx; 
            Ex = reshape(S[1:nt_Ex, :], nx_Ex, ny_Ex, nz_Ex, M)[(info.xPML[1].npixels+1):(nx_Ex-info.xPML[2].npixels),(info.yPML[1].npixels+1):(ny_Ex-info.yPML[2].npixels),(info.zPML[1].npixels+1):(nz_Ex-info.zPML[2].npixels),:]
            Ey = reshape(S[nt_Ex+1:nt_Ex+nt_Ey, :], nx_Ey, ny_Ey, nz_Ey, M)[(info.xPML[1].npixels+1):(nx_Ey-info.xPML[2].npixels),(info.yPML[1].npixels+1):(ny_Ey-info.yPML[2].npixels),(info.zPML[1].npixels+1):(nz_Ey-info.zPML[2].npixels),:]
            Ez = reshape(S[nt_Ex+nt_Ey+1:nt_Ex+nt_Ey+nt_Ez, :], nx_Ez, ny_Ez, nz_Ez, M)[(info.xPML[1].npixels+1):(nx_Ez-info.xPML[2].npixels),(info.yPML[1].npixels+1):(ny_Ez-info.yPML[2].npixels),(info.zPML[1].npixels+1):(nz_Ez-info.zPML[2].npixels),:]
        else
            Ex = reshape(S[1:nt_Ex, :], ny_Ex, nz_Ex, M)[(info.yPML[1].npixels+1):(ny_Ex-info.yPML[2].npixels),(info.zPML[1].npixels+1):(nz_Ex-info.zPML[2].npixels),:]
        end
    else
        if ~use_2D_TM
            Ex = reshape(S[1:nt_Ex, :], nx_Ex, ny_Ex, nz_Ex, M)
            Ey = reshape(S[nt_Ex+1:nt_Ex+nt_Ey, :], nx_Ey, ny_Ey, nz_Ey, M)
            Ez = reshape(S[nt_Ex+nt_Ey+1:nt_Ex+nt_Ey+nt_Ez, :], nx_Ez, ny_Ez, nz_Ez, M)
        else
            Ex = reshape(S[1:nt_Ex, :], ny_Ex, nz_Ex, M)
        end
    end
    
    t2 = time()
    
    info.timing_solve = info.timing_solve + t2 - t1 # Add the little bit of post-processing time

    info.timing_total = info.timing_total + t2 - t1 # Add the little bit of post-processing time
    if info.opts.verbal
        @printf("          Total elapsed time: %7.3f secs\n", info.timing_total) 
    end
    if ~use_2D_TM
        return Ex, Ey, Ez, info
    else
        return Ex, info
    end
end

# When syst, B, and C are specified; return the generalized scattering matrix.
function mesti(syst::Syst, B::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64, 2},Vector{Source_struct}}, C::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2}, Vector{Source_struct},String})
    return mesti(syst, B, C, nothing, nothing)
end


# When syst, B, and C, and opts are specified; return the generalized scattering matrix.
function mesti(syst::Syst, B::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64, 2},Vector{Source_struct}}, C::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2}, Vector{Source_struct},String}, opts::Opts)
    return mesti(syst, B, C, nothing, opts)
end


# When syst, B, and C, and D are specified; return the generalized scattering matrices
function mesti(syst::Syst, B::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64, 2},Vector{Source_struct}}, C::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2}, Vector{Source_struct},String}, D::Union{SparseMatrixCSC{Int64,Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64,Int64},Array{Int64,2},Array{Float64,2},Array{ComplexF64,2},Array{Int32,2},Array{Float32,2},Array{ComplexF32,2}})
    return mesti(syst, B, C, D, nothing)
end
