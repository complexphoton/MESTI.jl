# Metalens Focusing via Angular Spectrum Propagation

This example considers a metalens with a smaller input aperture than the output aperture, enabling diffraction-limited focusing across a 60-deg field of view (FOV). For detailed design information,
see the [referenced paper](https://doi.org/10.1364/OPTICA.514907).
We show how to use mesti() to compute the field response of a metalens for various incident angles, and use the angular spectrum propagation (ASP) to obtain the field profile before/after the metalens in free space.

```julia
# Include the function to perform the angular spectrum propagation (ASP)
include("asp.jl")

# Whether to compute the field inside a metalens (defaults to false)
# For focal plane field calculations, propagating the transmitted field immediately after the metalens is sufficient
# and more cost-effective than computing the entire field inside the metalens.
# Enable this option if the internal metalens field is required.
calculate_field_inside_metalens = false

# Load necessary packages
using MESTI, Interpolations, FFTW, Printf, LinearAlgebra, MAT, Plots, LaTeXStrings
```

## System parameters

```julia
n_air       = 1.0            # Refractive index of air on the right
n_sub       = 1.0            # Refractive index of substrate on the left
n_struct    = 2.0            # Refractive index of silicon nitride
wavelength  = 1              # Vacuum wavelength
dx          = wavelength/40  # Discretization grid size [wavelength]
FOV         = 60             # FOV of the metalens in the air [deg]

h            = 5*wavelength           # Thickness of the metalens [wavelength]
D_out        = 50*wavelength          # Output diameter of the metalens [wavelength]
D_in         = 25*wavelength          # Input diameter of the metalens [wavelength]
NA           = 0.9                    # Numerical aperture NA
focal_length = D_out/2/tan(asin(NA))  # Focal length [wavelength]

# Parameters for the input and output
W_out = D_out + 4*wavelength  # Width where we sample the transmitted field [wavelength]
                              # (a larger W_out ensures all transmitted light is captured)

# W_out > D_out, so we will pad extra pixels
ny_R_extra_half = Int(round((W_out-D_out)/dx/2))

# Number of pixels in the y direction for the source (on the left) and the projection (on the right)
ny = Int(ceil(D_out/dx))
ny_L = Int(ceil(D_in/dx))
ny_R = ny + 2*ny_R_extra_half
# Number of pixels of the metasurface in the z direction
nz = Int(ceil(h/dx))

nPML = 20      # Number of pixels of PMLs
# Add PMLs on the z direction, add one pixel of free space for source and projection
nz_extra_left = 1 + nPML
nz_extra_right = nz_extra_left
# Add PMLs on the y direction
ny_extra_low = ny_R_extra_half + nPML
ny_extra_high = ny_extra_low
# Number of pixels for the whole system with PMLs
ny_tot = ny + ny_extra_low + ny_extra_high
nz_tot = nz + nz_extra_left + nz_extra_right

k0dx = 2*pi/wavelength*dx   # Dimensionless frequency k0*dx
epsilon_L = n_sub^2         # Relative permittivity on the left
epsilon_R = n_air^2         # Relative permittivity on the right and on the top & bottom

# Obtain properties of propagating channels on the input side of metalens.
# We will use them to build the truncated plane-wave sources
BC = "periodic"  # Periodic boundary condition means the propagating channels are plane waves
use_continuous_dispersion = true  # Use continuous dispersion relation for (kx,ky), i.e. kx^2 + ky^2 = k0^2 in free space 
channels_L = mesti_build_channels(ny, BC, k0dx, epsilon_L, nothing, use_continuous_dispersion)

# List of incident angles of interest in air [degree]
theta_in_list = collect(-FOV/2:1:FOV/2)

# For the incident plane waves, we only take channels within the desired FOV:
# |ky| <= n_air*k0*sin(FOV/2)
kydx_FOV = k0dx*sind.(theta_in_list)   # kydx within the FOV at 1-deg increments
kxdx_FOV = sqrt.((k0dx)^2 .- (kydx_FOV).^2)      # kxdx within the FOV
N_L = length(kydx_FOV)   # Number of incident angles of interest

# Build the propagating channels on the left of metalens
B_basis = channels_L.f_x_m(channels_L.kydx_prop)       

# Multiply the flux-normalized prefactor sqrt(nu)
sqrt_nu_L_basis = reshape(channels_L.sqrt_nu_prop,1,:)  # row vector

# Build the truncated incident plane waves at various angles
ind_source_out = findall(x-> (abs(x) > D_in/2), (collect(0.5:1:ny) .- ny/2)*dx)  # y coordinates outside of the input aperture
B_trunc = channels_L.f_x_m(kydx_FOV)*sqrt(ny/ny_L)
B_trunc[ind_source_out,:] .= 0      # Block light outside of the input aperture

# We project each truncated input plane wave, with |y|<=D_in/2, onto propagating channels at the metalens input side,
# and use their combinations as sources, such that the inputs satisfy the free-space wave equation.
B_L = zeros(ComplexF64,ny,N_L)
for ii = 1:N_L
	for jj = 1:channels_L.N_prop
		Psi_in = sum(B_trunc[:,ii].*conj(B_basis[:,jj]),dims=1)
		B_L[:,ii] .= B_L[:,ii] .+ Psi_in*sqrt_nu_L_basis[jj]*exp((-1im*1/2)*channels_L.kzdx_prop[jj]).*B_basis[:,jj]
	end
end

# Load the permittivity profile of the wide-field-of-view metalens
file = matopen("permittivity_of_metalens.mat")
epsilon_metalens = read(file, "epsilon_metalens")
close(file)
epsilon_syst = ones(ny_tot,nz_tot)
epsilon_syst[ny_extra_low.+(1:ny),nz_extra_left.+(1:nz)] = epsilon_metalens # Add pixels for PMLs

# Specify the system for mesti()
syst = Syst()
syst.epsilon_xx = epsilon_syst
syst.wavelength = wavelength
syst.dx = dx
pml = PML(nPML)
pml.direction = "all"
syst.PML = [pml]  # PML on all four sides

# In mesti(), B_struct.pos = [m1, n1, h, w] specifies the position of a
# block source, where (m1, n1) is the index of the smaller-(y,z) corner,
# and (h, w) is the height and width of the block.
# Here, we put line sources (w=1) on the metalens's left surface (n1=n_L) and centered around it.
n_L = nz_extra_left                    # z pixel immediately before the metalens (where the source is placed)
m1_L = ny_extra_low + 1                # first y pixel of the metalens
B_struct = Source_struct()
B_struct.pos = [[m1_L, n_L, ny, 1]]
B_struct.data = [B_L]

# Use MUMPS with single-precision arithmetic to improve the computational efficiency
opts = Opts()
opts.prefactor = -2im
opts.solver = "MUMPS"
opts.use_single_precision_MUMPS = true
```

## Angular Spectrum Propagation (ASP) parameters

```julia
# System width used for ASP to remove periodic wrapping artifact.
W_ASP_min = 2*D_out  # Minimal ASP window [wavelength]

# dx = wavelength/40 is not necessary. Down-sample to a coarser resolution for ASP.
dy_ASP = 3*dx    # ASP grid size [wavelength]

# If dy_ASP/dx is not a integer, we round it to a integer.
if round(dy_ASP/dx) != dy_ASP/dx
    @warn "dy_ASP/dx = $(dy_ASP/dx) should be a positive integer to ensure down-sampling is possible; rounding it to $(maximum([1, Int(round(dy_ASP/dx))]))."
    dy_ASP = maximum([1, round(dy_ASP/dx)])*dx
end
# If ny_R is an even number and dy_ASP/dx is also even number, 
# we reduce dy_ASP/dx to dy_ASP/dx-1, making it odd number.
if mod(ny_R, 2) == 0 && mod(dy_ASP/dx, 2) == 0
    @warn "When ny_R is an even number, dy_ASP/dx should be an odd number so down-sampling can be symmetric; reducing it to $(Int((dy_ASP/dx)-1))."
    dy_ASP = ((dy_ASP/dx)-1)*dx
end

# fft is more efficient when the length is a power of 2, so we make ny_ASP a power of 2.
ny_ASP = nextpow(2,Int(round(W_ASP_min/dy_ASP))) # Number of pixels for ASP
W_ASP = ny_ASP*dy_ASP                            # Actual ASP window [micron]

# y index of the points we down-sample for ASP (if any).
ind_ASP = Int.(collect(round.(1:(dy_ASP/dx):ny_R)))

# Make the sampling points symmetric around the middle.
if ind_ASP[end] != ny_R
    # Make sure that ny_R - ind_ASP[end] is an even number.
    if mod(ny_R - ind_ASP[end], 2) != 0
        ind_ASP = ind_ASP[1:(end-1)]
    end
    ind_ASP = Int.(ind_ASP .+ (ny_R .- ind_ASP[end])./2)
end

ny_ASP_pad = ny_ASP - length(ind_ASP)         # Total number of zeros to pad
ny_ASP_pad_low = Int(round(ny_ASP_pad/2))     # Number of zeros to pad on the low side
ny_ASP_pad_high = ny_ASP_pad - ny_ASP_pad_low # Number of zeros to pad on the high side

# y coordinates of the ASP sampling points, including the padded zeros 
y_ASP = (collect(0.5:ny_ASP) .- 0.5*(ny_ASP + ny_ASP_pad_low - ny_ASP_pad_high))*dy_ASP

ny_ASP_half = Int(round(ny_ASP/2))   # recall that ny_ASP is an even number

# List of (kx,ky) in ASP
# Note that after fft, ky_ASP will be (2*pi/W_ASP)*(0:(ny_ASP-1)), where the
# latter half needs to be unwrapped to negative ky.
ky_ASP = (2*pi/W_ASP).*[collect(0:(ny_ASP_half-1)); collect(-ny_ASP_half:-1)]  # [1/wavelength]
kx_ASP = sqrt.(Complex.((n_air*2*pi/wavelength)^2 .- ky_ASP.^2))  # [1/wavelength]

# Propagating components in ASP
ky_ASP_prop = ky_ASP[findall(x-> (abs(x) < (n_air*2*pi/wavelength)), ky_ASP)] # must be a column vector per asp() syntax
kx_ASP_prop = sqrt.((n_air*2*pi/wavelength)^2 .- (ky_ASP_prop).^2)
        
# Focal spot position for each incident angle
focal_spot_list = focal_length.*tand.(theta_in_list)
interp_linear = linear_interpolation(y_ASP,collect(1:ny_ASP))
ind_plot = Int.(round.(interp_linear(focal_spot_list)))
```

## Field propagation to the focal plane in free space using ASP

```julia
if !calculate_field_inside_metalens  
    # Extract the field right after the metalens using matrix C_R
    # avoiding calculating and storing the entire field inside the metalens
    n_R = n_L + nz + 1                     # z pixel immediately after the metalens
    m1_R = nPML + 1                        # first y pixel of the output projection window
            
    C_R = zeros(ny_R,length(ind_ASP))
    C_R[CartesianIndex.(ind_ASP, 1:length(ind_ASP))] .= 1
    
    C_struct = Source_struct()
    C_struct.pos = [[m1_R, n_R, ny_R, 1]]
    C_struct.data = [C_R]
            
    # Use mesti() to compute the field immediately after the metalens.
    field_right_after_metalens, _ = mesti(syst, [B_struct], [C_struct], opts)

    # Select 5 incident angles from the FOV for demonstration
    ind_incident_angle = indexin([0; 7; 15; 22; 30], Int.(round.(theta_in_list)))
    num_angle_plot = length(ind_incident_angle)
    x_plot = focal_length 
    field_at_focal_plane = zeros(ComplexF64, ny_ASP, num_angle_plot)
    for ii = 1:num_angle_plot
        # Pick the field profile right after the metalens
        Ex0_ASP = field_right_after_metalens[:,ind_incident_angle[ii]] 

        # Propagate Ex0_ASP in the air to the right using ASP
        # Here we only consider progagating waves as evanescent waves decay over long distances.
        field_at_focal_plane[:,ii] = asp(Ex0_ASP, x_plot, kx_ASP_prop, ny_ASP)
    end

    # Plot the intensity profile at the focal plane for various incident angles
    p = plot(layout = (1, 5), legend = false, size = (900, 300),
        titlefontsize = 10, xtickfontsize = 10, ytickfontsize = 10, framestyle=:box,
        left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

    colors = palette(:darkrainbow, num_angle_plot)

    for ii = 1:num_angle_plot
        plot!(p[ii], y_ASP, abs.(field_at_focal_plane[:, ii]).^2, lw=1.5,
        xlim = (round(y_ASP[ind_plot[ind_incident_angle[ii]]];digits=1)-3, round(y_ASP[ind_plot[ind_incident_angle[ii]]];digits=1)+3), 
        ylim = (0, 0.4),
        xticks = [round(y_ASP[ind_plot[ind_incident_angle[ii]]];digits=1)-3, round(y_ASP[ind_plot[ind_incident_angle[ii]]];digits=1), round(y_ASP[ind_plot[ind_incident_angle[ii]]];digits=1)+3],
        yticks = 0:0.2:0.4,
        title = L"$\theta_\mathrm{in} =$" * string(Int(round(theta_in_list[ind_incident_angle[ii]]))) * L"$^{\circ}$",
        titlefontsize = 10,
        color = colors[ii],
        xlabel = (ii == 3 ? "Lateral position y/λ" : ""),
        ylabel = (ii == 1 ? L"$|E_x(y,z=f)|^2$" : ""))
    end

    plot!(p)
    display(p)
```

![intensity_at_focal_plane](https://github.com/complexphoton/MESTI.jl/assets/109620064/12406723-e979-4584-bd34-416d912327b0)

## Calculate the field on the left side, inside and after the metalens

```julia
else
    # Use mesti() to compute the field inside the metalens.
    field_inside_metalens, _ = mesti(syst, [B_struct], opts)
                    
    # Field propagation on the right side of the metalens using ASP                 
    x_plot_right = reshape(collect(dx:dx:round(24/dx)*dx), 1, :)   # Keep it a 1-by-N array for asp() syntax
    nx_plot = length(x_plot_right)
    field_after_metalens = zeros(ComplexF64, ny_ASP, nx_plot, N_L)
    for ii = 1:N_L
        # Pick and down-sample the field profile right after the metalens
        Ex0_ASP = field_inside_metalens[nPML .+ ind_ASP,end-nPML,ii] 

        # Propagate Ex0_ASP to the right in the air using ASP.
        # Here we consider both propagating and evanescent waves to accurately capture the near field near the metalens.
        field_after_metalens[:,:,ii] = asp(Ex0_ASP, x_plot_right, kx_ASP, ny_ASP)
    end
                        
    # Use mesti() to compute the free-space field generated by source B_L
    syst.epsilon_xx = ones(ny_tot, nz_extra_left+nz_extra_right)
    field_incident, _ = mesti(syst, [B_struct], opts)
                        
    # Incident field propagation on the left side of the metalens using ASP                 
    x_plot_left = reshape(collect(-round(7/dx)*dx:dx:-dx), 1, :)   # Keep it a 1-by-N array for asp() syntax
    nx_plot = length(x_plot_left)
    field_before_metalens_incident = zeros(ComplexF64, ny_ASP, nx_plot, N_L)
    for ii = 1:N_L
        # Pick and down-sample the incident field
        Ex0_ASP = field_incident[nPML .+ ind_ASP,nPML+1,ii] 

        # Propagate Ex0_ASP to the left in the air using ASP.
        # Here we consider propagating channels as the incident waves consist of only propagating components.
        field_before_metalens_incident[:,:,ii] = asp(Ex0_ASP, x_plot_left, kx_ASP_prop, ny_ASP)
    end
                            
    # Refelcted field propagation on the left side of the metalens using ASP                 
    x_plot_left = reshape(collect(round(7/dx)*dx:-dx:dx), 1, :)   # Keep it a 1-by-N array for asp() syntax
    nx_plot = length(x_plot_left)
    field_before_metalens_reflect = zeros(ComplexF64, ny_ASP, nx_plot, N_L)
    for ii = 1:N_L
        # Pick and down-sample the reflected field
        Ex0_ASP = field_inside_metalens[nPML .+ ind_ASP,nPML+1,ii] - field_incident[nPML .+ ind_ASP,nPML+1,ii] 

        # Propagate Ex0_ASP to the left in the air using ASP.
        # Here we consider both propagating and evanescent waves to accurately capture the near field near the metalens.
        field_before_metalens_reflect[:,:,ii] = asp(Ex0_ASP, x_plot_left, kx_ASP, ny_ASP)
    end
    
    # Adjust the scale of colormap                            
    hot_cm = palette(:hot, 256)
    scale_factor = 0.8
    hot_scaled = [RGB(c.r ^ scale_factor, c.g ^ scale_factor, c.b ^ scale_factor) for c in hot_cm]
                                
    coolwarm_cm = palette(:coolwarm, 256)
    scale_factor = 0.8
    coolwarm_scaled = [RGB(c.r ^ scale_factor, c.g ^ scale_factor, c.b ^ scale_factor) for c in coolwarm_cm]

    # Apply a phase factor to smooth phase transitions between frames without affecting focusing quality.
    _, ind_focal_plane = findmin(abs.(x_plot_right .- focal_length))
    phase_shift = zeros(1,N_L)
    for ii = 1:N_L
        phase_shift[ii] = -angle(field_after_metalens[ind_plot[ii],ind_focal_plane[2],ii])
    end
                                    
    # Animation of focusing within the FOV
    anim = Animation()  # Initialize the animation

    ind_y_plot = findall(x-> (abs(x) < W_out/2), y_ASP)
    for ii in 1:N_L
        p = plot(layout = (2,1), size = (400, 600) , framestyle=:box,
                left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

        # Intensity profile
        heatmap!(p[1], y_ASP[ind_y_plot], -round(7/dx)*dx-dx/2-h:dx:round(24/dx)*dx+dx/2, abs.(transpose(hcat(field_before_metalens_incident[ind_y_plot,:,ii] .+ field_before_metalens_reflect[ind_y_plot,:,ii], field_inside_metalens[nPML .+ ind_ASP,nPML+1:end-nPML,ii], field_after_metalens[ind_y_plot,:,ii]))).^2, 
            xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24),
            xticks = -20:10:20, yticks = [-10;0:12:24], xtickfontsize = 10, ytickfontsize = 10,
            xlabel = "Lateral position y/λ", ylabel = "Axial position z/λ",
            c=hot_scaled, clims = (0,0.2), colorbar_ticks=nothing, aspect_ratio=:equal, dpi=150)
        annotate!(p[1], 0, 30, text(L"$|E_x(y,z)|^2$", 11, :center, :top))
        annotate!(p[1], -24, 23, text(L"$\theta_\mathrm{in} =$" * string(Int(round(theta_in_list[ii]))) * L"$^{\circ}$",
            11, :left, :top, :white))
        plot!(p[1], [-D_in/2, D_in/2], -h*[1,1], linestyle = :dash, lw = 1, color = :white, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
        plot!(p[1], [-D_out/2, D_out/2], 0*[1,1], linestyle = :dash, lw = 1, color = :white, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
        plot!(p[1], [-D_in/2, -D_out/2], [-h,0], linestyle = :dash, lw = 1, color = :white, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
        plot!(p[1], [D_in/2, D_out/2], [-h,0], linestyle = :dash, lw = 1, color = :white, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))

        # Phase profile
        heatmap!(p[2], y_ASP[ind_y_plot], -round(7/dx)*dx-dx/2-h:dx:round(24/dx)*dx+dx/2, real.(transpose(hcat(field_before_metalens_incident[ind_y_plot,:,ii] .+ field_before_metalens_reflect[ind_y_plot,:,ii], field_inside_metalens[nPML .+ ind_ASP,nPML+1:end-nPML,ii], field_after_metalens[ind_y_plot,:,ii])*exp(1im*phase_shift[ii]))),
            xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24),
            xticks = -20:10:20, yticks = [-10;0:12:24], xtickfontsize = 10, ytickfontsize = 10,
            xlabel = "Lateral position y/λ", ylabel = "Axial position z/λ",
            c=coolwarm_scaled, clims = (-0.16,0.16), colorbar_ticks=nothing, aspect_ratio=:equal, dpi=150)
        annotate!(p[2], 0, 29, text(L"$\mathrm{Re}[E_x(y,z)]$", 11, :center, :top))
        annotate!(p[2], -24, 23, text(L"$\theta_\mathrm{in} =$" * string(Int(round(theta_in_list[ii]))) * L"$^{\circ}$",
            11, :left, :top))
        plot!(p[2], [-D_in/2, D_in/2], -h*[1,1], linestyle = :dash, lw = 1, color = :black, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
        plot!(p[2], [-D_out/2, D_out/2], 0*[1,1], linestyle = :dash, lw = 1, color = :black, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
        plot!(p[2], [-D_in/2, -D_out/2], [-h,0], linestyle = :dash, lw = 1, color = :black, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
        plot!(p[2], [D_in/2, D_out/2], [-h,0], linestyle = :dash, lw = 1, color = :black, legend = false, xlim = (-W_out/2,W_out/2), ylim = (-(7+h),24))
                                                        
        frame(anim, p)
    end

    gif(anim, "focusing_animation.gif", fps = 12)
end
```

![focusing_animation](https://github.com/complexphoton/MESTI.jl/assets/109620064/cfe707e5-445b-40b8-ab9a-6275041b1a7a)
