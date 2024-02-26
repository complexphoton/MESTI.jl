## Open Channel Through Disorder
# In this example, we show how to use mesti2s() to compute the transmission
# matrix of a strongly scattering disordered medium, analyze the transmission
# matrix to determine an incident wavefront that can penetrate the disorder with
# almost 100% transmission (called an "open channel"), and then use mesti2s() 
# again to compute the field profile of the open channel while comparing to that
# of a typical plane-wave input.

# Call necessary packages
using MESTI, GeometryPrimitives, Arpack, Printf

# Include the function to build epsilon_xx for the disordered
include("build_epsilon_disorder.jl")

## System parameters
# dimensions of the system, in units of the wavelength lambda_0
dx      = 1/15  # discretization grid size
W       = 360   # width of the scattering region
L       = 90    # thickness of the scattering region
L_tot   = 150   # full length of the system for plotting
r_min   = 0.2   # minimal radius of the cylindrical scatterers
r_max   = 0.4   # maximal radius of the cylindrical scatterers
min_sep = 0.05  # minimal separation between cylinders
number_density = 1.3 # number density, in units of 1/lambda_0^2
rng_seed = 0    # random number generator seed

# relative permittivity, unitless
epsilon_scat = 1.2^2  # cylindrical scatterers
epsilon_bg   = 1.0^2  # background in the scattering region
epsilon_low  = 1.0^2  # frees space on the low side
epsilon_high = 1.0^2  # frees space on the high side

yBC = "periodic" # boundary condition in y

# generate a random collection of non-overlapping cylinders
build_TM = true
(epsilon_xx, y0_list, z0_list, r0_list, y_Ex, z_Ex) = 
build_epsilon_disorder(W, L, r_min, r_max, min_sep,  
                       number_density, rng_seed, dx,
                       epsilon_scat, epsilon_bg, build_TM)

## Compute the transmission matrix
syst = Syst()
syst.epsilon_xx = epsilon_xx
syst.epsilon_low = epsilon_low
syst.epsilon_high = epsilon_high
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC

# specify the input and output
input = channel_type()
output = channel_type()
# input from the low side
input.side = "low"
# output to the high side
output.side = "high"

# put PML along z-direction
pml = mesti_optimal_pml_params(syst.wavelength/syst.dx)
pml_npixels = 15
pml.npixels = pml_npixels
syst.zPML = [pml]

# transmission matrix: input from the low side, output to the high side
t, channels, _ = mesti2s(syst, input, output)

## Compare an open channel and a plane-wave input
# The most-open channels is the singular vector of the transmission matrix with 
# the largest singular value.
(_, sigma_max, v_open), _, _, _, _ = svds(t, nsv=1)

N_prop_low = channels.low.N_prop # number of propagating channels on the low side
ind_normal = Int(round((N_prop_low+1)/2)) # index of the normal-incident plane-wave

# compare the transmission
T_avg = sum(abs.(t).^2)/N_prop_low # average over all channels
T_PW  = sum(abs.(t[:,ind_normal]).^2) # normal-incident plane-wave
T_open = sigma_max[1].^2 # open channel

println(" T_avg   = ", @sprintf("%.2f", T_avg), "\n T_PW    = ", @sprintf("%.2f", T_PW), "\n T_open  = ", @sprintf("%.2f", T_open))

# specify two input incident wavefronts:
# (1) normal-incident plane-wave
# (2) open channel
v_low = zeros(ComplexF64, N_prop_low, 2)
v_low[ind_normal, 1] = 1
v_low[:, 2] = v_open
input = wavefront()
input.v_low = v_low

# we will also get the field profile in the free spaces on the two sides, for
# plotting purpose.
opts = Opts()
nz_low  = round(Int,(L_tot-L)/2/dx)
nz_high = nz_low
opts.nz_low = nz_low
opts.nz_high = nz_high

# compute the field-profiles through mesti()
Ex, _, _ = mesti2s(syst, input, opts)

## Define the matrix B by users and compare the field results
# do the same field-profile computations through defining matrix B by users and using mesti()
syst = Syst()
ny_Ex = size(epsilon_xx,1) # total number of pixel along y for Ex grids

# in mesti(), syst.epsilon_low and syst.epsilon_high are not needed
# in mesti(), we provide the whole epsilon_xx including the scattering region, source/detection region, and PML region
syst.epsilon_xx = cat(epsilon_low*ones(ny_Ex,pml_npixels+1), epsilon_xx, epsilon_high*ones(ny_Ex,pml_npixels+1), dims=2)
# in previous mesti2s() calculation, 
# syst.epsilon_xx = epsilon_xx
# syst.epsilon_low = epsilon_low   
# syst.epsilon_high = epsilon_high
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC 

# put PML along z-direction
pml = mesti_optimal_pml_params(syst.wavelength/syst.dx)
pml.npixels = pml_npixels
pml.direction = "z" # put
syst.PML = [pml]
# in previous mesti2s() calculation, 
# syst.zPML = [pml] (and do not need to specify pml.direction = "z")

Bx = Source_struct()
u_prop_low_Ex = channels.u_x_m(channels.low.kydx_prop) # the transverse profiles, u_Ex, for the propagating ones

# If it is a single propagating channel, a line source of -2i*sqrt(nu)*u_Ex(m) at l=0 will generate an z-flux-normalized incident field of exp(i*kzdx*|l|)*u_Ex(l)/sqrt(nu), where nu = sin(kzdx)
# here, we want input wavefront fields on the low side, which is a superposition of propagating channels
# these input wavefront fields on the low side are specified by v_low, and we take superpositions of the propagating channels using the v coefficients, with the sqrt(nu)*exp(-i*kzdx*0.5) prefactors included.
# we will multiply the -2i prefactor at the end.
# the source B_Ex would generate the input wavefront fields
B_Ex_low = u_prop_low_Ex*(channels.low.sqrt_nu_prop.*exp.((-1im*0.5)*channels.low.kzdx_prop).*v_low) # 0.5 pixel backpropagation indicates that the source is half a pixel away from z = 0
Bx.data = [B_Ex_low] # the value of the source

# the position of the source specify by a rectangle.
# [m1, l1, w, h] specifies the location and the size of the
# rectangle. Here, (m1, l1) is the index of the (y,z) coordinate of
# the smaller-y, smaller-z corner of the rectangle, at the location
# of Ex(m1, l1); (w, h) is the width and height of the rectangle
# line source is put on the lower source region, where is one pixel outside PML region.
Bx.pos = [[1,pml_npixels+1,ny_Ex,1]] 

opts = Opts()
# the -2i prefactor would be multiplied by.
opts.prefactor = -2im

Ex_prime, _ = mesti(syst, [Bx], opts)

# Excluding the extra padding and PML region, compare these field profiles
println("Maximum absolute value of field difference between constructing the source matrix B through mesti2s() and constructing by users = ", 
         maximum(abs.(Ex[:,nz_low+1:end-nz_high-1,:] - Ex_prime[:,pml_npixels+1+1:end-pml_npixels-1-1,:])))

## Animate the field profiles
using Plots
# normalize the field amplitude with respect to the plane-wave-input profile
Ex = Ex/maximum(abs.(Ex[:,:,1]))

nframes_per_period = 20

# extend the x coordinate to include free spaces on the two sides
z_Ex = vcat(z_Ex[1] .- (opts.nz_low:-1:1)*dx, z_Ex, z_Ex[end] .+ (1:opts.nz_high)*dx)

# animate the field profile with plane-wave input
anim_pw = @animate for ii ∈ 0:(nframes_per_period-1)
    plt1 = (heatmap(z_Ex,collect(y_Ex),real.(Ex[:,:,1]*exp(-1im*2*π*ii/nframes_per_period)),
            xlabel = "z", ylabel = "y", c = :balance, clims=(-1, 1), aspect_ratio=:equal, dpi = 450,
            xlimits=(-25,115), ylimits=(0,360)))
    scatter!(plt1, z0_list, y0_list,markersize=r0_list*0.88, alpha=0.1, color=:black, legend=false, dpi = 450)
end
gif(anim_pw, "disorder_PW_input.gif", fps = 10)

# animate the field profile of the open channel
anim_open_ch = @animate for ii ∈ 0:(nframes_per_period-1)
    plt2 = (heatmap(z_Ex,collect(y_Ex),real.(Ex[:,:,2]*exp(-1im*2*π*ii/nframes_per_period)),
            xlabel = "z", ylabel = "y", c = :balance, clims=(-1, 1), aspect_ratio=:equal, dpi = 450,
            xlimits=(-25,115), ylimits=(0,360)))
    scatter!(plt2, z0_list, y0_list,markersize=r0_list*0.88, alpha=0.3, color=:black, legend=false, dpi = 450)
end
gif(anim_open_ch, "disorder_open_channel.gif", fps = 10)
