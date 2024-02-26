# Export a composite data type
export Side

mutable struct Side
    # A composite data type to store the items on a side
    # See also: mesti_setup_longitudinal
    N_prop::Integer
    kzdx_all::Vector{ComplexF64}
    ind_prop::Vector{Int64}
    kxdx_prop::Union{Vector{Float64},Nothing}
    kydx_prop::Vector{Float64}
    kzdx_prop::Vector{Float64}
    sqrt_nu_prop::Vector{Float64}
    ind_prop_conj::Vector{Int64}
    Side() = new() 
end

"""
    MESTI_SETUP_LONGITUDINAL sets up a structure "Side" for one component in the homogeneous space.
        === Input Arguments === 
        k0dx (numeric scalar, real or complex; required):
            Dimensionless frequency, k0*dx = (2*pi/vacuum_wavelength)*dx. 
        epsilon_bg (numeric scalar, real or complex; required):
            Relative permittivity of the homogeneous space.
        kxdx_all (1-by-nx_Ex+delta_(xBC,"Dirichlet") real row vector):
            Dimensionless transverse wave number kx*dx for all nx channels,
            including both propagating and evanescent ones. They are real-valued
            and are ordered from small to large. 
        kydx_all (1-by-ny_Ey+delta_(yBC,"Dirichlet") real row vector):
            Dimensionless transverse wave number ky*dx for all ny channels,
            including both propagating and evanescent ones. They are real-valued
            and are ordered from small to large. 
        kLambda_x (numeric scalar, real or complex; optional):
            kx_B*Lambda_x, where Lambda_x is the periodicity along x-direction.
            This is used only for periodic and Bloch periodic boundary.
        kLambda_y (numeric scalar, real or complex; optional):
            ky_B*Lambda_y, where Lambda_y is the periodicity along y-direction.
            This is used only for periodic and Bloch periodic boundary.
        ind_zero_kx (numeric scalar; optional):
            The transverse mode index where kxdx = kx_B*Lambda_x/nx.
            This is used only for periodic and Bloch periodic boundary.
        ind_zero_ky (numeric scalar; optional):
            The transverse mode index where kydx = ky_B*Lambda_y/ny.
            This is used only for periodic and Bloch periodic boundary.
        use_continuous_dispersion (logical scalar; optional, defaults to false):
            Whether to use the dispersion equation of the continuous wave equation
            when building the input/output channels. Defaults to false, in which case
            the finite-difference dispersion is used.
        
        === Output Arguments === 
        side (scalar structure):
            side.N_prop (integer scalar):
                Number of propagating channels. 
            side.kzdx_all (1-by-nx_Ex+delta_(xBC,"Dirichlet")*ny_Ey+delta_(yBC,"Dirichlet") complex row vector):
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
            side.ind_prop (1-by-N_prop integer row vector):
                Indices of the N_prop propagating channels among all channels.low.kzdx_all. 
            side.kxdx_prop (1-by-N_prop real row vector):
                Dimensionless longitudinal wave number kx*dx for the N_prop propagating
                channels, equal to channels.kxdx_all(channels.low.ind_prop). 
            side.kydx_prop (1-by-N_prop real row vector):
                Dimensionless transverse wave number ky*dx for the N_prop propagating
                channels, equal to channels.kydx_all(channels.low.ind_prop). 
            side.kzdx_prop (1-by-N_prop real row vector):
                Dimensionless transverse wave number kz*dx for the N_prop propagating
                channels, equal to channels.low.kzdx_all(channels.low.ind_prop).  
            side.sqrt_nu_prop (1-by-N_prop row vector):
                Square root of the normalized longitudinal group velocity of the
                propagating channels, sqrt_nu_prop = sqrt(sin(kzdx)). The longitudinal
                group velocity is v_g = (sin(kzdx)/k0dx)*(c/epsilon_low). 
            side.ind_prop_conj (1-by-N_prop integer row vector; optional):
                A permutation vector that switches one propagating channel with one
                having a complex-conjugated transverse profile. Now it only uses in 2D systems.
"""
function mesti_setup_longitudinal(k0dx::Union{Real,Complex}, epsilon_bg::Union{Real,Complex}, kxdx_all::Union{StepRangeLen{<:Real}, Vector{<:Real}, Nothing}, kydx_all::Union{StepRangeLen{<:Real}, Vector{<:Real}}, kLambda_x::Union{Real,Complex,Nothing}=nothing, kLambda_y::Union{Real,Complex,Nothing}=nothing, ind_zero_kx::Union{Int,Nothing}=nothing, ind_zero_ky::Union{Int,Nothing}=nothing, use_continuous_dispersion::Bool=false)
    
    if kxdx_all == nothing && kLambda_x == nothing && ind_zero_kx == nothing
        use_2D_TM = true 
    else
        use_2D_TM = false
    end
    
    k0dx2_epsilon = (k0dx^2)*epsilon_bg
    if ~use_2D_TM
        nx = size(kxdx_all,1)
    end
    ny = size(kydx_all,1)

    side = Side()

    if ~use_continuous_dispersion
        # use the finite-difference dispersion for homogeneous space: 
        # 2D: k0dx2_epsilon =                   4*sin^2(kydx/2) + 4*sin^2(kzdx/2)
        # 3D: k0dx2_epsilon = 4*sin^2(kxdx/2) + 4*sin^2(kydx/2) + 4*sin^2(kzdx/2)
        # sin_kzdx_over_two_sq = sin^2(kzdx/2)
        if use_2D_TM
            sin_kzdx_over_two_sq = 0.25*k0dx2_epsilon .- sin.(kydx_all/2).^2
        else
            sin_kzdx_over_two_sq = reshape(0.25*k0dx2_epsilon .- sin.(kxdx_all/2).^2 .- sin.(reshape(vcat(kydx_all),1,:)/2).^2,nx*ny)
        end
        #sin_kzdx_over_two_sq = reshape(0.25*k0dx2_epsilon .- sin.(repeat(kxdx_all,1,ny)/2).^2 .- sin.(repeat(transpose(kydx_all),nx,1)/2).^2, nx*ny)
        # Dimensionless longitudinal wave number
        # asin(sqrt(z)) has two branch points (at z=0 and z=1) and with the branch cuts going outward on the real-z axis; we will address the branch choice below
        # Note kzdx is only defined up to modulo 2*pi (ie, kxdx is equivalent to kzdx + 2*pi, kzdx - 2*pi, etc) because sin(kzdx) and exp(1i*kzdx) are both invariant under 2pi shifts.
        side.kzdx_all = 2*asin.(sqrt.(Complex.((sin_kzdx_over_two_sq))))
        # Indices of the propagating channels
        # When k0dx2_epsilon is real, these are indicies of the channels with real-valued kzdx
        # When k0dx2_epsilon is complex, these are indices of the channels we consider "propagating-like"; they have complex kzdx with 0 < real(kzdx) < pi. When k0dx2_epsilon is tuned to a real number continuously, this set continuously becomes that at real k0dx2_epsilon.
        side.ind_prop = findall(x-> (real(x) > 0 && real(x) < 1), sin_kzdx_over_two_sq)  
        # Here we address the sign choice of kzdx, namely its branch
        # When k0dx2_epsilon is real, we choose the sign of kzdx such that:
        # 1) 0 < kzdx < pi for propagating channels (where kzdx is real)
        # 2) Im(kzdx) >= 0 for evanescent channels
        # Using the correct sign is important when we build the retarded Green's function of the semi-infinite homogeneous space.
        # The default branch choice of asin(sqrt(z)) returns the correct sign for the most part, except when z > 1. We need to flip the sign of those (which can only occur if k0dx2_epsilon > 6).
        # When k0dx2_epsilon is complex-valued, it is not always possible to unambiguously choose the sign that is "physical", because kzdx will be complex-valued, and the sign we "want" for real(kzdx) and the sign we want for imag(kzdx) may be incompatible.
        # What we do with complex k0dx2_epsilon is that we choose the sign for the (complex-valued) kxdx such that when k0dx2_epsilon is tuned to a real number continuously by fixing Re(k0dx2_epsilon) and varying Im(k0dx2_epsilon), the kzdx we choose continuously becomes the "correct" one at real k0dx2_epsilon without crossing any branch cut. To do so, we rotate the two branch cuts of asin(sqrt(z)) by 90 degrees to the lower part of the complex-z plane (ie, the lower part of the complex-k0dx2_epsilon plane), and we pick the branch with the correct sign when k0dx2_epsilon is real. This is implemented by flipping the sign of kzdx for the ones that require flipping.
        # Note that we will get a discontinuity whenever k0dx2_epsilon crosses one of those vertical-pointing branch cuts. That is unavoidable.
        # The following few lines implement the "flipping".
        if ~(isa(k0dx2_epsilon, Real)) || (~use_2D_TM && isa(k0dx2_epsilon, Real) && k0dx2_epsilon > 6) || (use_2D_TM && isa(k0dx2_epsilon, Real) && k0dx2_epsilon > 4)
            # Note that when imag(sin_kzdx_over_two_sq)=0, flipping is needed for sin_kzdx_over_two_sq>1 but not needed for sin_kzdx_over_two_sq<0.
            ind_flip = findall(x->(real(x)<0 && imag(x)<0) || (real(x)>1 && imag(x)<=0), sin_kzdx_over_two_sq);  
            side.kzdx_all[ind_flip] = -side.kzdx_all[ind_flip]
        end
    else
        # use the continuous dispersion for homogeneous space:
        # 2D: k0dx2_epsilon =          kydx^2 + kzdx^2
        # 3D: k0dx2_epsilon = kxdx^2 + kydx^2 + kzdx^2
        if use_2D_TM
            kzdx2 = k0dx2_epsilon .- kydx_all.^2
        else
            kzdx2 = reshape(k0dx2_epsilon .- kxdx_all.^2 .- reshape(vcat(kydx_all),1,:).^2,nx*ny)
        end
        
        side.kzdx_all = sqrt.(Complex.(kzdx2))
        side.ind_prop = findall(x-> real(x) > 0, kzdx2)
        if ~isa(k0dx2_epsilon, Real)
            ind_flip = findall(x-> (real(x) < 0 && imag(x) < 0), kzdx2)
            side.kzdx_all[ind_flip] = -side.kzdx_all[ind_flip]
        end
    end        
    # Number of propagating channels
    side.N_prop = length(side.ind_prop)

    # Wave numbers of the propagating channels
    side.kzdx_prop = side.kzdx_all[side.ind_prop]
    if use_2D_TM
        side.kxdx_prop = nothing
        side.kydx_prop = kydx_all[side.ind_prop]
    else
        side.kxdx_prop = kxdx_all[((side.ind_prop).%nx) .+ nx*((side.ind_prop).%nx .== 0)]
        side.kydx_prop = kydx_all[Int.(ceil.((side.ind_prop)./nx))]
    end
    
    # Square root of the normalized longitudinal group velocity, sqrt(sin(kzdx)), for the propagating channels
    # nu = sin(kzdx)
    # When k0dx2_epsilon is real, sqrt_nu_prop is also real. When k0dx2_epsilon is complex, sqrt_nu_prop is also complex.
    side.sqrt_nu_prop = sqrt.(sin.(side.kzdx_prop))
    
    if ~use_2D_TM
        # TODO: think how to do and implement side.ind_prop_conj in 3D
        side.ind_prop_conj = 1:side.N_prop
    else
        # Permutation that switches one propagating channel with one having a complex-conjugated transverse profile.
        if isa(kLambda_y, Nothing)
            # For Dirichlet and Neumann boundaries, u_x_m is real, so no permutation needed
            side.ind_prop_conj = 1:side.N_prop
        elseif kLambda_y == 0
            # For periodic boundary condition, complex conjugation switches ky and -ky
            if (ind_zero_ky in side.ind_prop) || (mod(side.N_prop,2) == 0)
                # Simply flip the ordering
                side.ind_prop_conj = side.N_prop:-1:1
            else
                # The last channel has -ky equal to ky due to aliasing so should not be flipped
                side.ind_prop_conj = vcat((side.N_prop-1):-1:1,side.N_prop)
            end
            # TODO: implement side.ind_prop_conj when kLambda_y == pi
        end
    end
    
    return side
end
