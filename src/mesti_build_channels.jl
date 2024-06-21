# Export composite data types
export Channels_one_sided
export Channels_two_sided

# Export a function mesti_build_channels()
export mesti_build_channels

abstract type Channels end

mutable struct Channels_two_sided <: Channels
    # A composite data type to store the items on two sides, transverse functions, and wave vectors
    # See also: mesti_build_channels
    f_x_n::Function
    f_y_n::Function
    f_z_n::Function
    df_z_n::Function  
    f_x_m::Function
    f_y_m::Function
    f_z_m::Function
    df_z_m::Function
    kxdx_all::Union{Vector{Float64},Nothing}
    kydx_all::Vector{Float64}
    low::Side
    high::Side
    Channels_two_sided()=new()
end

mutable struct Channels_one_sided <: Channels
    # A composite data type to store the items on one side, transverse functions, and wave vectors
    # See also: mesti_build_channels        
    f_x_n::Function
    f_y_n::Function
    f_z_n::Function
    df_z_n::Function      
    f_x_m::Function
    f_y_m::Function
    f_z_m::Function
    df_z_m::Function    
    kxdx_all::Union{Vector{Float64},Nothing}
    kydx_all::Vector{Float64}
    N_prop::Integer
    kzdx_all::Vector{ComplexF64}
    ind_prop::Vector{Int64}
    kxdx_prop::Union{Vector{Float64},Nothing}
    kydx_prop::Vector{Float64}
    kzdx_prop::Vector{Float64}
    sqrt_nu_prop::Vector{Float64}
    ind_prop_conj::Vector{Int64}
    Channels_one_sided()=new()
end 

"""
    MESTI_BUILD_CHANNELS Set up properties of channels in the homogeneous space. 
        MESTI_BUILD_CHANNELS(nx_Ex, nx_Ey, xBC, ny_Ex, ny_Ey, yBC, k0dx, epsilon_low, epsilon_high, 
        use_continuous_dispersion, n0, m0) returns a structure containing properties of the propagating 
        and evanescent channels in a homogeneous space in the transverse direction, boundary
        condition xBC along x and yBC along y, background relative permittivity epsilon_low, epsilon_high, 
        and dimensionless frequency k0dx = (2*pi/vacuum_wavelength)*dx. 

        MESTI_BUILD_CHANNELS(syst) does the same but with the parameters extracted
        from structure syst; see ''? mesti2s'' for the fields required for
        structure "syst".

        === Input Arguments === 
        nx_Ex (positive integer scalar; required in 3D):
            Number of grid points in the x-direction for Ex. 
        nx_Ey (positive integer scalar; required in 3D):
            Number of grid points in the x-direction for Ey. 
        xBC (string or scalar number; required in 3D):
            Boundary condition in the x-direction. 
            When xBC is a character vector, available choices are (case-insensitive): 
                "periodic" - periodic BC on both sides
                "PEC"      - PEC on both sides
                "PMC"      - PMC on both sides
                "PECPMC"   - PEC on low side and PMC on high side
                "PMCPEC"   - PMC on low side and PEC on high side
            When xBC is a scalar number, the Bloch periodic boundary condition is
            used with p(m+ny) = p(m)*exp(1i*xBC); in other words, xBC = kx_B*nx_Ex*dx =
            kx_B*nx_Ey*dx = kx_B*p where ky_B is the Bloch wave number and p = nx*dx is the
            periodicity in x. 
        ny_Ex (positive integer scalar; required in 3D):
            Number of grid points in the y-direction for Ex. 
        ny_Ey (positive integer scalar; required ):
            Number of grid points in the y-direction for Ey. 
        yBC (string or scalar number; required):
            Boundary condition in y-direction, analogous to xBC. 
        k0dx (numeric scalar, real or complex; required):
            Dimensionless frequency, k0*dx = (2*pi/vacuum_wavelength)*dx. 
        epsilon_low (numeric scalar, real or complex; required):
            Relative permittivity of the homogeneous space on the low. 
        epsilon_high (numeric scalar, real or complex, or nothing; optional):
            Relative permittivity of the homogeneous space on the high. Only the
            low side will be considered if epsilon_high is not given or is nothing. 
        use_continuous_dispersion (logical scalar; optional, defaults to false):
            Whether to use the dispersion equation of the continuous wave equation
            when building the input/output channels. Defaults to false, in which case
            the finite-difference dispersion is used. 
        n0 (real numeric scalar, optional, defaults to 0):
            Center of the 1D transverse mode profile along x-direction with periodic or Bloch periodic
            boundary condition (nx = nx_Ex = nx_Ey), f_{n,a} = exp(i*kx(a)*dx*(n-n0))/sqrt(nx_Ex), 
            where kx(a) = kx_B + a*(2*pi/nx*dx).
        m0 (real numeric scalar, optional, defaults to 0):
            Center of the 1D transverse mode profile along x-direction with periodic or Bloch periodic
            boundary condition, analogous to n0. 

        === Output Arguments === 
        channels (scalar structure):
            channels.kxdx_all (1-by-nx_Ex+delta_(xBC,"Dirichlet") real row vector):
                Dimensionless transverse wave number kx*dx for all nx channels,
                including both propagating and evanescent ones. They are real-valued
                and are ordered from small to large. 
            channels.kydx_all (1-by-ny_Ey+delta_(yBC,"Dirichlet") real row vector):
                Dimensionless transverse wave number ky*dx for all ny channels,
                including both propagating and evanescent ones. They are real-valued
                and are ordered from small to large. 
            channels.f_x_n (function):
                A 1D transverse function along x-direction for Ex 
            channels.f_y_n (function):
                A 1D transverse function along x-direction for Ey 
            channels.f_z_n (function):
                A 1D transverse function along x-direction for Ez 
            channels.df_z_n (function):
                A derivative of 1D transverse function along x-direction for Ez 
            channels.f_x_m (function):
                A 1D transverse function along y-direction for Ex 
            channels.f_y_m (function):
                A 1D transverse function along y-direction for Ey 
            channels.f_z_m (function):
                A 1D transverse function along y-direction for Ez 
            channels.df_z_m (function):
                A derivative of 1D transverse function along y-direction for Ez 
            channels.low (scalar structure):
                When epsilon_low and epsilon_high are both given (i.e., epsilon_high is given
                and is not nothing), the properties specific to the low and high sides
                are returned in channels.low and channels.high; channels.low and channels.high
                are both scalar structures, and their fields are described below.
                    When only epsilon_low = epsilon_bg is given; there will be no
                channels.low and channels.high; instead, the fields that would have been
                assigned to channels.low will be assigned to channels directly. For
                example, channels.low.N_prop will be channels.N_prop instead.
            channels.low.N_prop (integer scalar):
                Number of propagating channels. 
            channels.low.kzdx_all (1-by-nx_Ex+delta_(xBC,"Dirichlet")*ny_Ey+delta_(yBC,"Dirichlet") complex row vector):
                Dimensionless longitudinal wave number kz*dx for all channels,
                including both propagating and evanescent ones. Due to the
                discretization, kzdx is equivalent to kzdx + 2*pi. Whenever kzdx is a
                solution, -kzdx is also a solution. Here, we choose the sign of kzdx
                such that:
                When k0dx is real, we have 
                    - Propagating channels: 0 < Re(kzdx) < pi, Im(kzdx) = 0. 
                    - Evanescent channels: Im(kzdx) >= 0, mod(Re(kzdx),pi) = 0. 
                When k0dx is complex, we analytically continue the above choice onto
                the complex-frequency plane. Specifically, we pick the kzdx that is
                continuously connected to one with real k0dx through a vertical line
                in the complex-(k0dx^2) plane. 
            channels.low.ind_prop (1-by-N_prop integer row vector):
                Indices of the N_prop propagating channels among all channels.low.kzdx_all. 
            channels.low.kxdx_prop (1-by-N_prop real row vector):
                Dimensionless longitudinal wave number kx*dx for the N_prop propagating
                channels, equal to channels.kxdx_all(channels.low.ind_prop). 
            channels.low.kydx_prop (1-by-N_prop real row vector):
                Dimensionless transverse wave number ky*dx for the N_prop propagating
                channels, equal to channels.kydx_all(channels.low.ind_prop). 
            channels.low.kzdx_prop (1-by-N_prop real row vector):
                Dimensionless transverse wave number kz*dx for the N_prop propagating
                channels, equal to channels.low.kzdx_all(channels.low.ind_prop).  
            channels.low.sqrt_nu_prop (1-by-N_prop row vector):
                Square root of the normalized longitudinal group velocity of the
                propagating channels, sqrt_nu_prop = sqrt(sin(kzdx)). The longitudinal
                group velocity is v_g = (sin(kzdx)/k0dx)*(c/epsilon_low). 
            channels.low.ind_prop_conj (1-by-N_prop integer row vector; optional):
                A permutation vector that switches one propagating channel with one
                having a complex-conjugated transverse profile. It only works for 2D now.
            channels.high (scalar structure; optional):
                Structure containing properties specific to the high side,
                similar to channels.high; only provided when epsilon_high is given. 
"""
function mesti_build_channels(nx_Ex::Union{Int,Nothing}, nx_Ey::Union{Int,Nothing}, xBC::Union{String,Real,Complex,Nothing}, ny_Ex::Int, ny_Ey::Union{Int,Nothing}, yBC::Union{String,Real,Complex}, k0dx::Union{Real,Complex}, epsilon_low::Union{Real,Complex}, epsilon_high::Union{Real,Complex, Nothing}=nothing, use_continuous_dispersion::Bool=false, n0::Union{Real,Nothing}=0, m0::Union{Real}=0)
    
    # Check whehter 2D TM case or not
    if nx_Ex == nothing && nx_Ey == nothing && ny_Ex != nothing && ny_Ey == nothing
       # 2D TM case
       use_2D_TM = true
       if xBC != nothing
            @warn "Only yBC is required for 2D TM fields Ex(y,z). xBC will be ignored."
       end
       if n0 != nothing
            @warn "Only m0 is required for 2D TM fields Ex(y,z). n0 will be ignored."
       end
    else
        # 3D case
        use_2D_TM = false
    end
    
    # Check input parameters
    if nx_Ex != nothing && ~(nx_Ex > 0)
        throw(ArgumentError("Input argument nx_Ex must be a positive integer scalar if given."))
    elseif nx_Ey != nothing && ~(nx_Ey > 0) 
        throw(ArgumentError("Input argument nx_Ey must be a positive integer scalar if given."))
    elseif ~(ny_Ex > 0 ) 
        throw(ArgumentError("Input argument ny_Ex must be a positive integer scalar."))       
    elseif ny_Ey != nothing && ~(ny_Ey > 0) 
        throw(ArgumentError("Input argument ny_Ey must be a positive integer scalar if given."))
    end

    if isa(epsilon_high, Nothing)
        two_sided = false    
    else
        two_sided = true
    end

    # Convert BC to take care of lowercase or uppercase
    if ~use_2D_TM
        xBC = convert_BC(xBC, "x")
    end
    yBC = convert_BC(yBC, "y")    
    
    if ~use_2D_TM
        # Check number of grid points with the boundary conditions   
        check_BC_and_grid(xBC, nx_Ex, nx_Ey, "x")
        check_BC_and_grid(yBC, ny_Ex, ny_Ey, "y")
    end

    # Convert PEC/PMC to Dirichlet/Neumann based on component and direction
    # BC_j_k means boundary condition for Ej (j = x, y, z) component along k-direction (k = x, y)
    if ~use_2D_TM
        BC_x_x = convert_BC_to_transverse(xBC, "x", "x")
        BC_y_x = convert_BC_to_transverse(xBC, "y", "x")    
        BC_y_y = convert_BC_to_transverse(yBC, "y", "y")
        BC_z_x = convert_BC_to_transverse(xBC, "z", "x")
        BC_z_y = convert_BC_to_transverse(yBC, "z", "y")
    end
    # This line also applies to 2D TM case
    BC_x_y = convert_BC_to_transverse(yBC, "x", "y")

    # These are used only for periodic and Bloch periodic boundary conditions; otherwise they stay nothing
    kLambda_x = nothing # kx_B*Lambda_x
    kLambda_y = nothing # ky_B*Lambda_y
    ind_zero_kx = nothing
    ind_zero_ky = nothing

    # Handle periodic and Bloch periodic boundary conditions
    if ~use_2D_TM
        if isa(xBC, Number)
            kLambda_x = xBC
            xBC = "Bloch"
            # kLambda_x must be real for f_x_n(kxdx_x) and f_y_n(kxdx_y) to be unitary
            if ~isa(kLambda_x, Real)
                @warn("kx_B*Lambda_x = $(real(kLambda_x)) + 1im*$(imag(kLambda_x)) is a complex number; must be real for a complete orthonormal transverse basis.")
            end
        elseif lowercase(xBC) == "bloch"
            throw(ArgumentError("To use Bloch periodic boundary condition in mesti_build_channels(), set the input argument xBC to kx_B*Lambda_x where kx_B is the Bloch wave number and Lambda_x is the periodicity along x-direction."))
        elseif xBC == "periodic"
            kLambda_x = 0
            xBC = "Bloch"
        end
    end

    if isa(yBC, Number)
        kLambda_y = yBC
        yBC = "Bloch"
        # kLambda_y must be real for f_x_m(kydx_x) and f_y_m(kydx_y) to be unitary
        if ~isa(kLambda_y, Real)
            @warn("ky_B*Lambda_y = $(real(kLambda_y)) + 1im*$(imag(kLambda_y)) is a complex number; must be real for a complete orthonormal transverse basis.")
        end
    elseif lowercase(yBC) == "bloch"
        throw(ArgumentError("To use Bloch periodic boundary condition in mesti_build_channels(), set the input argument yBC to ky_B*Lambda_y where ky_B is the Bloch wave number and Lambda_y is the periodicity along y-direction."))
    elseif yBC == "periodic"
        kLambda_y = 0
        yBC = "Bloch"
    end

    # p = [p(1), ..., p(nx)].'
    # For periodic and Bloch periodic boundary, we order channels.kxdx_all (channels.kydx_all) such that it increases monotonically from negative to positive
    # For other boundary conditions, kx >= 0 (ky >= 0), and we order channels.kxdx_all (channels.kydx_all) such that it increases monotonically from smallest to largest
    if two_sided
        channels = Channels_two_sided()
    else
        channels = Channels_one_sided()
    end
    
    if ~use_2D_TM
        (channels.f_x_n, channels.kxdx_all) = mesti_build_transverse_function(nx_Ex, BC_x_x, n0, true)
        (channels.f_x_m, _)                 = mesti_build_transverse_function(ny_Ex, BC_x_y, n0)
        (channels.f_y_n, _)                 = mesti_build_transverse_function(nx_Ey, BC_y_x, n0)    
        (channels.f_y_m, channels.kydx_all) = mesti_build_transverse_function(ny_Ey, BC_y_y, m0, true)

        channels.f_z_n = channels.f_y_n # f_z_n and f_y_n are same transverse function  
        channels.f_z_m = channels.f_x_m # f_z_m and f_x_m are same transverse function  

        channels.df_z_n = mesti_build_transverse_function_derivative(nx_Ey, BC_z_x, n0, 1)
        channels.df_z_m = mesti_build_transverse_function_derivative(ny_Ex, BC_z_y, m0, 1)
    else
        (channels.f_x_m, channels.kydx_all) = mesti_build_transverse_function(ny_Ex, BC_x_y, m0)
        channels.kxdx_all = nothing
    end
    
    if ~use_2D_TM && xBC == "Bloch"
        if mod(nx_Ex,2) == 1
            ind_zero_kx = Int(round((nx_Ex+1)/2))
        else
            ind_zero_kx = Int(round(nx_Ex/2))
        end
    end

    if yBC == "Bloch"
        if mod(ny_Ex,2) == 1
            ind_zero_ky = Int(round((ny_Ex+1)/2))
        else
            ind_zero_ky = Int(round(ny_Ex/2))
        end
    end

    # Properties for the homogeneous space on the low (kzdx, sqrt_nu_prop, number of propagating channels, etc; depends on epsilon_low/high)
    side = mesti_setup_longitudinal(k0dx, epsilon_low, channels.kxdx_all, channels.kydx_all, kLambda_x, kLambda_y, ind_zero_kx, ind_zero_ky, use_continuous_dispersion)

    if two_sided
        channels.low = side
        # Homogeneous space on the high
        if epsilon_high == epsilon_low
            channels.high = side
        elseif ~isnan(epsilon_high)
            channels.high = mesti_setup_longitudinal(k0dx, epsilon_high, channels.kxdx_all, channels.kydx_all, kLambda_x, kLambda_y, ind_zero_kx, ind_zero_ky, use_continuous_dispersion)
        end
    else
        # Add the fields of "side" to "channels"
        for fnames in fieldnames(typeof(side))
            setfield!(channels, fnames, getfield(side, fnames))
        end
    end
    return channels
end

"""
    MESTI_BUILD_CHANNELS(syst) does the same but with the parameters extracted
    from structure syst; see '? mesti2s' for the fields required for
    structure 'syst'.
"""
function mesti_build_channels(syst::Syst)
    # Check whehter 2D TM case or not    
    if ndims(syst.epsilon_xx) == 2
        # 2D TM case
        use_2D_TM = true
        if ((isdefined(syst, :epsilon_yy) && syst.epsilon_yy != nothing) || (isdefined(syst, :epsilon_zz) && syst.epsilon_zz != nothing))
            @warn "Only syst.epsilon_xx is required for 2D TM fields Ex(y,z). Other components will be ignored."
        end
    else
        # 3D case
        use_2D_TM = false
    end
    
    if ~use_2D_TM
        (nx_Ex, ny_Ex, _) = size(syst.epsilon_xx)
        (nx_Ey, ny_Ey, _) = size(syst.epsilon_yy)
    else
        (ny_Ex, _) = size(syst.epsilon_xx)
    end

    if ~isdefined(syst, :epsilon_low); throw(ArgumentError("Input argument syst must have field \"epsilon_low\".")); end
    epsilon_low = syst.epsilon_low

    if ~isdefined(syst, :wavelength); throw(ArgumentError("Input argument syst must have field \"wavelength\".")); end
    if ~isdefined(syst, :dx); throw(ArgumentError("Input argument syst must have field \"dx\".")); end
    if ~(syst.dx > 0); throw(ArgumentError("syst.dx must be a positive scalar.")); end
    k0dx = (2*pi/syst.wavelength)*(syst.dx)

   # Check boundary condition in x    
    if ~use_2D_TM
        if isdefined(syst, :kx_B) 
            if isdefined(syst, :xBC) && lowercase(syst.xBC) != "bloch"
                throw(ArgumentError("When syst.kx_B is given, syst.xBC must be \"Bloch\" if specified."))
            end
            syst.xBC = "Bloch"
            # mesti_build_channels() uses kx_B*periodicity as the input arguments xBC for Bloch BC
            xBC = (syst.kx_B)*(nx_Ex*syst.dx) # dimensionless
        else
            # Defaults to Dirichlet boundary condition unless syst.kx_B is given
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
        # Defaults to Dirichlet boundary condition unless syst.ky_B is given
        if ~isdefined(syst, :yBC)
            throw(ArgumentError("Input argument syst must have non-empty field \"yBC\" when syst.ky_B is not given."))       
        elseif ~(lowercase(syst.yBC) in ["bloch", "periodic", "pec", "pmc", "pecpmc", "pmcpec"])
            throw(ArgumentError("syst.yBC = \"$(syst.yBC)\" is not a supported option; type ''? mesti2s'' for supported options."))
        elseif lowercase(syst.yBC) == "bloch"
            throw(ArgumentError("syst.yBC = \"Bloch\" but syst.ky_B is not given."))
        end
        yBC = syst.yBC
    end

    if ~isdefined(syst, :epsilon_high)
        epsilon_high = nothing
    else
        epsilon_high = syst.epsilon_high
    end
    use_continuous_dispersion = false
    n0 = 0
    m0 = 0

    if ~use_2D_TM
        channels = mesti_build_channels(nx_Ex, nx_Ey, xBC, ny_Ex, ny_Ey, yBC, k0dx, epsilon_low, epsilon_high, use_continuous_dispersion, n0, m0)
    else
        channels = mesti_build_channels(nothing, nothing, nothing, ny_Ex, nothing, yBC, k0dx, epsilon_low, epsilon_high, use_continuous_dispersion, nothing, m0)
    end

    return channels
end


"""
    MESTI_BUILD_CHANNELS(ny_Ex, yBC, k0dx, epsilon_low, epsilon_high, use_continuous_dispersion, m0) set up  
        properties of channels in the homogeneous space for 2D TM fields and returns a structure containing  
        properties of the propagating and evanescent channels in a homogeneous space with ny pixels in the
        transverse (y) direction, boundary condition yBC along y, background relative permittivity epsilon_low, 
        epsilon_high, and dimensionless frequency k0dx = (2*pi/vacuum_wavelength)*dx where dx is the discretization grid size.
"""
function mesti_build_channels(ny_Ex::Int, yBC::Union{String,Real,Complex}, k0dx::Union{Real,Complex}, epsilon_low::Union{Real,Complex}, epsilon_high::Union{Real,Complex, Nothing}=nothing, use_continuous_dispersion::Bool=false, m0::Union{Real}=0)
    return mesti_build_channels(nothing, nothing, nothing, ny_Ex, nothing, yBC, k0dx, epsilon_low, epsilon_high, use_continuous_dispersion, nothing, m0)
end


""" 
    CONVERT_BC_TO_TRANSVERSE is a helper function to convert PEC, PMC, and etc. to Dirichlet, Neumann, and etc. according to the component and direction
"""
function convert_BC_to_transverse(BC::Union{String,Real,Complex},component::String,direction::String)
    if isa(BC, Number)
        return BC
    elseif BC == "PEC"
        if direction == component
            return "Neumann"            
        else
            return "Dirichlet"
        end
    elseif BC == "PMC"
        if direction == component
            return "Dirichlet"            
        else            
            return "Neumann"
        end            
    elseif BC == "PECPMC"
        if direction == component
            return "NeumannDirichlet"            
        else                            
            return "DirichletNeumann"
        end            
    elseif BC == "PMCPEC"
        if direction == component
            return "DirichletNeumann"            
        else                                                                    
            return "NeumannDirichlet"
        end            
    elseif BC == "periodic"
        return "periodic"
    elseif lowercase(BC) == "bloch"
        throw(ArgumentError("To use Bloch periodic boundary condition in $(direction)-direction, set the input argument $(direction)BC to k$(direction)_B*p_$(direction) where k$(direction)_B is the Bloch wave number and p_$(direction) is the periodicity along $(direction)-direction."))        
    else
        throw(ArgumentError("Input argument $(direction)BC = \"$(BC)\" is not a supported option."))        
    end
end
