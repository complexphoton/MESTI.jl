# Open Channel Through Disorder in 3D System

In this example, we show how to use mesti2s() to compute the transmission matrix of a 3D scattering disordered medium, analyze the transmission matrix to determine an incident wavefront that can penetrate the disorder with almost 100% transmission (called an "open channel") in 3D system.

```julia
# load necessary packages
using MESTI, LinearAlgebra, Statistics, Printf

# include the function to build permittivity tensor for the disordered 
# and plot the transmission eigenvalue distribution
include("build_epsilon_disorder_3d.jl")
include("plot_and_compare_distribution.jl");
```

# System parameters

```julia
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
                                                                rng_seed, dx, epsilon_scat,                                                                           epsilon_bg, 
                                                                x1, x2, y1, y2, z1, z2);
```

# Compute the transmission matrix 

```julia
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
pml_npixels = 20
pml = PML(pml_npixels)
syst.zPML = [pml]

opts = Opts()
# clear variables to reduce peak memory usage
opts.clear_memory = true
opts.clear_BC = true

# note this transmission matrix calculation may take between one to two hours in single-core,
# but a few minutes in mutlithreading calculation.
# utilizting mutlithreading is highly recommended 
t, channels, _ = mesti2s(syst, input, output, opts)

(_, sigma_max, v_open), _, _, _, _ = svds(t, nsv=1);
```
```text:Output
===System size===
nx_Ex = 76, ny_Ex = 76; nz_Ex = 38 => 80
nx_Ey = 76, ny_Ey = 76; nz_Ey = 38 => 80
nx_Ez = 76, ny_Ez = 76; nz_Ez = 39 => 81
[N_prop_low, N_prop_high] = [185, 185] per polarization
xBC = periodic; yBC = periodic; zBC = [PML, PML]
Building B,C... elapsed time:   7.597 secs
            ... elapsed time:   2.131 secs
Building A  ... elapsed time:   8.962 secs
< Method: APF using MUMPS in single precision with METIS ordering >
Building K  ... elapsed time:   0.368 secs
Analyzing   ... elapsed time:  14.966 secs
Factorizing ... elapsed time: 421.012 secs
          Total elapsed time: 463.109 secs
```

# Compare an open channel and a plane-wave input

```julia
# perform the singular value decomposition (SVD) on the transmission matrix
_, sigma, v = svd(t)

# transmission eigenvalue (eigenvalue of t^(dag)*t) is the square of the signular value of t 
tau = sigma.^2

# The most-open channels is the singular vector of the transmission matrix with 
# the largest singular value.
v_open = v[:, 1]
sigma_open = sigma[1]

# The most-closed channels is the singular vector of the transmission matrix with 
# the smallest singular value.
v_closed = v[:, end]
sigma_closed = sigma[end]

N_prop_low_per_pol = channels.low.N_prop # number of propagating channels per polarization on the low side
ind_normal = Int(round((N_prop_low_per_pol+1)/2)) # index of the normal-incident plane-wave

# compare the transmission
T_avg = sum(abs.(t).^2)/(2*N_prop_low_per_pol) # average over all channels
T_PW  = sum(abs.(t[:,ind_normal]).^2) # normal-incident plane-wave
T_closed = sigma_closed.^2 # closed channel
T_open = sigma_open.^2 # open channel

println(" T_avg   = ", @sprintf("%.4f", T_avg), "\n T_PW    = ", @sprintf("%.4f", T_PW), "\n T_closed= ", @sprintf("%.4f", T_closed), "\n T_open  = ", @sprintf("%.4f", T_open))
```
```text:Output
 T_avg   = 0.1449
 T_PW    = 0.1501
 T_closed= 0.0000
 T_open  = 0.9997
```

# Plot the the transmission eigenvalue distribution

```julia
# transmission eigenvalue (eigenvalue of t^(dag)*t) is the square of the signular value of t 
tau = sigma.^2

# plot the transmission eigenvalue distribution from transmission matrix t
# and compare it with the analytic distribution (bimodal DMPK distribution)
using Plots
plot_and_compare_distribution(tau);
```
<img src="https://github.com/complexphoton/MESTI.jl/assets/44913081/f9426518-ec70-4410-81f6-a08440f02739" width="800" height="600">

# Compute field profiles for a closed channel, an open channel and a plane-wave input

```julia
# specify three input incident wavefronts:
# (1) closed channel
# (2) open channel
# (3) normal-incident plane-wave

# specify inputs and output
input = wavefront()
output = nothing

input.v_low_s = zeros(ComplexF64, N_prop_low_per_pol, 2)
input.v_low_p = zeros(ComplexF64, N_prop_low_per_pol, 3)
# wavefront for closed channel
input.v_low_s[:,1] = v_closed[1:N_prop_low_per_pol,1]
input.v_low_p[:,1] = v_closed[N_prop_low_per_pol+1:N_prop_low_per_pol*2,1]
# wavefront for open channel
input.v_low_s[:,2] = v_open[1:N_prop_low_per_pol,1]
input.v_low_p[:,2] = v_open[N_prop_low_per_pol+1:N_prop_low_per_pol*2,1]
# wavefront for normal-incident plane-wave
input.v_low_p[Int(ceil(N_prop_low_per_pol/2)),3] = 1

opts = Opts()
# clear variables to reduce peak memory usage
opts.clear_memory = true
opts.clear_BC = true
opts.use_METIS = true
# use the disk to store the LU factors to save memory usage in core
# It will require approximately 43 GiB of disk storage to store the LU factors.
opts.write_LU_factor_to_disk = true 
ENV["MUMPS_OOC_TMPDIR"] = "." # write the out-of-core (LU factors) files in current folder 

# note this field-profile calculation may take between one to two hours in single-core,
# but a few minutes in mutlithreading calculation.
# utilizting mutlithreading is highly recommended 
(Ex, Ey, Ez, channels, info)= mesti2s(syst, input, output, opts);
# Ex, Ey, and Ez are 4D array. 
# For example, Ex(:,:,:,i) is the field profile Ex given the i-th input source profile.
```

```text:Output
===System size===
nx_Ex = 76, ny_Ex = 76; nz_Ex = 38 => 80
nx_Ey = 76, ny_Ey = 76; nz_Ey = 38 => 80
nx_Ez = 76, ny_Ez = 76; nz_Ez = 39 => 81
[N_prop_low, N_prop_high] = [185, 185] per polarization
xBC = periodic; yBC = periodic; zBC = [PML, PML]
Building B,C... elapsed time:   1.203 secs
            ... elapsed time:   0.618 secs
Building A  ... elapsed time:   2.757 secs
< Method: factorize_and_solve using MUMPS in single precision with METIS ordering >
Analyzing   ... elapsed time:  14.699 secs
Factorizing ... elapsed time: 488.263 secs
Solving     ... elapsed time: 246.227 secs
          Total elapsed time: 755.135 secs
```

# Plot the field profiles

```julia
# resulting field for Closed/Open from s-polarization input wavefront combined with p-polarization input wavefront
Ex_closed = Ex[:,:,:,1]+Ex[:,:,:,3]
Ex_open = Ex[:,:,:,2]+Ex[:,:,:,4]
# resulting field for normal-incident plane-wave 
Ex_pw = Ex[:,:,:,5]

# normalize the field amplitude with respect to the open-channel profile
Ex_open = Ex_open ./maximum(abs.(Ex_open[:,:,:]))
Ex_closed = Ex_closed ./maximum(abs.(Ex_open[:,:,:]))
Ex_pw = Ex_pw ./maximum(abs.(Ex_open[:,:,:]))

# compare and plot the field profiles (real(Ex)) in the plane x = 2 µm
heatmap1 = heatmap(real.(Ex_closed[38,:,:]), c = :balance, clims=(-1, 1), 
        legend = :none, aspect_ratio=:equal, dpi = 450, 
        ticks = false, framestyle = :none, title = "Closed channel")
heatmap2 = heatmap(real.(Ex_open[38,:,:]), c = :balance, clims=(-0.3, 0.3), 
        legend = :none, aspect_ratio=:equal, dpi = 450, 
        ticks = false, framestyle = :none, title = "Open channel")        
heatmap3 = heatmap(real.(Ex_pw[38,:,:]), c = :balance, clims=(-1, 1), 
        legend = :none, aspect_ratio=:equal, dpi = 450, 
        ticks = false, framestyle = :none, title = "Plane wave")
plot(heatmap1, heatmap2, heatmap3, layout = (1, 3))
```

<img src="https://github.com/complexphoton/MESTI.jl/assets/44913081/7b26770a-0846-4b7e-afc7-17de8665b197" width="900" height="600">

