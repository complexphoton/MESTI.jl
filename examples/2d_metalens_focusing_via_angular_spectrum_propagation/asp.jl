using FFTW # Perform fast Fourier transform

"""
ASP use angular spectrum propagation (ASP) to propagate a scalar field
   from f0 = f(x=0,y) to f(x,y).

   === Input Arguments ===
   f0 (numeric column vector or matrix):
       Initial field profile at x = 0 plane, as a column vector.
       If f0 has more than one column, each column is treated as a
       distinct initial field profile. size(f0,1) = ny.
   x (numeric scalar or row vector):
       Distance to propagate to. When multiple distances are of interest,
       they can be given as a row vector; the profiles f(x,y) at these
       distances will be returned as the different columns of f. When f0
       has more than one column, x must be a scalar (i.e., cannot
       propagate to multiple distances when multiple initial fields are
       given).
   kx_prop (N_prop-by-1 column vector):
       Longitudinal wave number kx. It can either (1) include all of the
       ny_tot wave numbers (both propagating and evanescent ones) in which
       case N_prop = length(kx_prop) must equal ny_tot and all of them
       will be propagated, or (2) include a subset of the wave numbers up
       to some cutoff in ky (e.g., the propagating ones) in which case
       N_prop = length(kx_prop) must be an odd number and only these
       components will be propagated.
   ny_tot (integer scalar; optional):
       Total number of grid points in y direction, including padded
       zeros. Must be no smaller than ny and N_prop. Defaults to ny.
   ny_pad_low (integer scalar; optional):
       Number of zeros to pad on the low side. Defaults to
       round((ny_tot-ny)/2).

   === Output Arguments ===
   f (numeric column vector or matrix):
       Field profile(s) f(x,y). Different columns correspond to
       different initial profiles (if f0 has more than one column) or
       different distance x (if x is a row vector).
"""

function asp(f0::Union{Matrix{Int32}, Matrix{Int64}, Matrix{Float32}, Matrix{Float64}, Matrix{ComplexF32}, Matrix{ComplexF64}, Vector{Int32}, Vector{Int64}, Vector{Float32}, Vector{Float64}, Vector{ComplexF32},  Vector{ComplexF64}}, x::Union{Int64, Float64, Matrix{Int64}, Matrix{Float64}}, kx_prop::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}, ny_tot::Union{Int64, Nothing}, ny_pad_low::Union{Int64, Nothing})
    
    # Check if f0 has more than one column and enforce x to be a scalar
    if size(f0, 2) > 1 && !(isa(x, Int64) || isa(x, Float64))
        error("x must be a scalar when f0 has more than one column.")
    end

    # Ensure x is a row vector when it's a matrix
    if isa(x, Matrix) && size(x, 1) != 1
        error("x must be a row vector.")
    end
    
    ny = size(f0,1)

    N_prop = length(kx_prop)

    # no zero-padding unless ny_tot is given
    if isa(ny_tot, Nothing)
        ny_tot = ny
    elseif ny_tot < ny
        throw(ArgumentError("ny_tot, when given, must be no smaller than size(f0,1) = $(ny)."))
    elseif ny_tot < N_prop
        throw(ArgumentError("ny_tot, when given, must be no smaller than length(kx_prop) = $(N_prop)."))
    end

    # pad zeros symmetrically by default
    if isa(ny_pad_low, Nothing)
        ny_pad_low = Int(round((ny_tot-ny)/2))
    elseif ny_pad_low + ny > ny_tot
        throw(ArgumentError("ny_pad_low + ny must be no greater than size(f0,1) = $(ny)."))
    end

    # Fourier transform f(x=0,y) to f(x=0,ky), as in Eq. (S43) of the APF paper.
    # To get a finer spacing in ky, zeros are padded below and above f0.
    f0_fft = exp.((-2im*pi*ny_pad_low/ny_tot).*collect(0:(ny_tot-1))).*fft([f0; zeros(ny_tot-ny,size(f0,2))],(1,))

    # Remove the evanescent components of f(x=0,ky), propagate it to f(x,ky),
    # and ifft back to f(x,y), as in Eqs. (S41-S42) of the APF paper.
    if N_prop == ny_tot
        f = ifft(exp.(1im.*kx_prop.*x).*f0_fft, (1,))
    else
        if mod(N_prop,2) != 1
            throw(ArgumentError("length(kx_prop) = $(N_prop) must be an odd number when it is not ny_tot."))
        end
        a_max = Int(round((N_prop-1)/2))
        ind_prop = vcat(collect(1:(a_max+1)), collect((ny_tot-a_max+1):ny_tot))
        f_fft_prop = exp.(1im.*kx_prop.*x).*f0_fft[ind_prop,:]
        f = exp.((-2im*pi*a_max/ny_tot).*collect(0:(ny_tot-1))).*ifft([circshift(f_fft_prop, a_max); zeros(ny_tot-N_prop,size(f_fft_prop,2))], (1,))
    end
    
    return f

end

# The following are asp functions to take different number of input arguments, but all of them will
# call the asp main function.

function asp(f0::Union{Matrix{Int32}, Matrix{Int64}, Matrix{Float32}, Matrix{Float64}, Matrix{ComplexF32}, Matrix{ComplexF64}, Vector{Int32}, Vector{Int64}, Vector{Float32}, Vector{Float64}, Vector{ComplexF32},  Vector{ComplexF64}}, x::Union{Int64, Float64, Matrix{Int64}, Matrix{Float64}}, kx_prop::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})
    return asp(f0, x, kx_prop, nothing, nothing)
end

function asp(f0::Union{Matrix{Int32}, Matrix{Int64}, Matrix{Float32}, Matrix{Float64}, Matrix{ComplexF32}, Matrix{ComplexF64}, Vector{Int32}, Vector{Int64}, Vector{Float32}, Vector{Float64}, Vector{ComplexF32},  Vector{ComplexF64}}, x::Union{Int64, Float64, Matrix{Int64}, Matrix{Float64}}, kx_prop::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}, ny_tot::Union{Int64, Nothing})
    return asp(f0, x, kx_prop, ny_tot, nothing)
end
