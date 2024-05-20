## Open Channel Through Disorder in 3D System
# In this example, we show how to use mesti2s() to compute the transmission
# matrix of a 3D scattering disordered medium, analyze the transmission
# matrix to determine an incident wavefront that can penetrate the disorder with
# almost 100% transmission (called an "open channel") in 3D system.

# Note this transmission matrix calculation in 3D system may take between one to two hours in single-core
# utilizting mutlithreading is highly recommended 

# call necessary packages
using MESTI, Arpack, Printf

# include the function to build permittivity tensor for the disordered 
# and plot the transmission eigenvalue distribution
include("build_epsilon_disorder_3d.jl")
include("plot_and_compare_distribution.jl")

## System parameters
# dimensions of the system, in units of µm
lambda_0   = 0.532       # free space wavelength (µm) 
dx         = lambda_0/10 # discretization grid size (µm)
W_x        = 4           # width of the scattering region along x-direction (µm)
W_y        = 4           # width of the scattering region along y-direction (µm)
L          = 2           # thickness of the scattering region (µm)
r_min      = 0.1         # minimal radius of the cylindrical scatterers (µm)
r_max      = 0.1         # maximal radius of the cylindrical scatterers (µm)
min_sep    = 0.05        # minimal separation between cylinders (µm)
rng_seed   = 0           # random number generator seed
number_density = 28      # number density, in units of 1/lambda_0^2

# relative permittivity, unitless
epsilon_scat = 2.54^2  # cylindrical scatterers
epsilon_bg   = 1.00^2  # background in the scattering region
epsilon_low  = 1.00^2  # frees space on the low side
epsilon_high = 1.00^2  # frees space on the high side

# edge of the scattering region (where we put scatterers)
x1 = 0  
x2 = ceil(Int, W_x/dx) * dx
y1 = 0
y2 = ceil(Int, W_y/dx) * dx
z1 = r_min
z2 = (ceil(Int, L/dx)+1) * dx - r_min

# generate a random collection of non-overlapping cylinders
(epsilon_xx, epsilon_xy, epsilon_xz, 
 epsilon_yx, epsilon_yy, epsilon_yz, 
 epsilon_zx, epsilon_zy, epsilon_zz) = build_epsilon_disorder_3d(W_x, W_y, L, r_min, 
                                                                r_max, min_sep, number_density, 
                                                                rng_seed, dx, epsilon_scat, epsilon_bg, 
                                                                x1, x2, y1, y2, z1, z2)

## Compute the transmission matrix
syst= Syst();
syst.epsilon_low = epsilon_low
syst.epsilon_high = epsilon_high
syst.length_unit = "µm"
syst.wavelength = lambda_0
syst.dx = dx
syst.xBC = "periodic"
syst.yBC = "periodic"
syst.epsilon_xx = epsilon_xx
syst.epsilon_xy = epsilon_xy
syst.epsilon_xz = epsilon_xz
syst.epsilon_yx = epsilon_yx
syst.epsilon_yy = epsilon_yy
syst.epsilon_yz = epsilon_yz
syst.epsilon_zx = epsilon_zx
syst.epsilon_zy = epsilon_zy
syst.epsilon_zz = epsilon_zz

# specify inputs and output
input = channel_type()
output = channel_type()
# input from the low side with both s-polarization and p-polarization
input.side = "low"
input.polarization = "both"
# output to the high side with both s-polarization and p-polarization
output.side = "high"
output.polarization = "both"

# put PML along z-direction
pml = mesti_optimal_pml_params(syst.wavelength/syst.dx)
pml_npixels = 20
pml.npixels = pml_npixels
syst.zPML = [pml]

opts = Opts()
# clear variables to reduce peak memory usage
opts.clear_memory = true
opts.clear_BC = true

# note this transmission matrix calculation may take between one to two hours in single-core
# utilizting mutlithreading is highly recommended 
t, channels, _ = mesti2s(syst, input, output, opts)

# the most-open channels is the singular vector of the transmission matrix with 
# the largest singular value.
(_, sigma_max, v_open), _, _, _, _ = svds(t, nsv=1)

N_prop_low = channels.low.N_prop # number of propagating channels on the low side
ind_normal = Int(round((N_prop_low+1)/2)) # index of the normal-incident plane-wave

# compare the transmission
T_avg = sum(abs.(t).^2)/(2*N_prop_low) # average over all channels
T_PW  = sum(abs.(t[:,ind_normal]).^2) # normal-incident plane-wave
T_open = sigma_max[1].^2 # open channel

println(" T_avg   = ", @sprintf("%.4f", T_avg), "\n T_PW    = ", @sprintf("%.4f", T_PW), "\n T_open  = ", @sprintf("%.4f", T_open))


## Plot the the transmission eigenvalue distribution
# plot the transmission eigenvalue distribution from transmission matrix t
# and compare it with the analytic distribution (bimodal DMPK distribution)
plot_and_compare_distribution(t)