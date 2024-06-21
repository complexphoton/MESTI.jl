# Focusing Inside Disorder With Phase Conjugation

In this example, we show how to use mesti() to project field generated from a point source inside disorder onto propagating channels through APF method, do phase conjugation to determine an incident wavefront that can focus inside the disorder, and then use mesti2s() again to compute the field profile to show its focus.

```julia
# Call necessary packages
using MESTI, GeometryPrimitives, LinearAlgebra, Statistics, printf

# Include the function to build epsilon_xx for the disordered
include("build_epsilon_disorder.jl")
```

# System parameters

```julia
# dimensions of the system, in units of the wavelength lambda_0
dx      = 1/15  # discretization grid size
W       = 360   # width of the scattering region
L       = 90    # thickness of the scattering region
L_tot   = 270   # full length of the system for plotting
r_min   = 0.2   # minimal radius of the cylindrical scatterers
r_max   = 0.4   # maximal radius of the cylindrical scatterers
min_sep = 0.05  # minimal separation between cylinders
number_density = 1.3  # number density, in units of 1/lambda_0^2
rng_seed = 0   # random number generator seed

# relative permittivity, unitless
epsilon_scat = 1.2^2  # cylindrical scatterers
epsilon_bg   = 1.0^2  # background in the scattering region
epsilon_low  = 1.0^2  # free space on the low side
epsilon_high = 1.0^2  # free space on the high side

yBC = "periodic" # boundary condition in y

# generate a random collection of non-overlapping cylinders
build_TM = true
no_scatterered_center = true
(epsilon, y0_list, z0_list, r0_list, y_Ex, z_Ex) =
build_epsilon_disorder(W, L, r_min, r_max, min_sep,  
                       number_density, rng_seed, dx,
                       epsilon_scat, epsilon_bg, build_TM; 
                       no_scatterer_center = true)
```

# Projecting field generated from a point source onto propagating channels through APF

```julia
syst = Syst()
pml_npixels = 20
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC
# specify the permittivity profile of the simulation domain including the low side, scattering region and the high side.
syst.epsilon_xx = cat(epsilon_low*ones(size(epsilon,1),pml_npixels+1), epsilon, epsilon_high*ones(size(epsilon,1),pml_npixels+1), dims=2)

# specify the input (point source in the middle of the disordered)
m0_focus = Int((W/dx)/2)
l0_focus = Int((L/dx)/2)
Bx = Source_struct()
Bx.pos = [[m0_focus,l0_focus+pml_npixels+1,1,1]]
Bx.data = [ones(1,1)]

# put PML along z-direction
pml = PML(pml_npixels)
pml.direction = "z"
syst.PML = [pml]

# build channels for the low side (air)
channels_low = mesti_build_channels(Int(W/dx), yBC, 2*pi*syst.wavelength*dx, epsilon_low)
N_prop_low = channels_low.N_prop # number of propagating channels on the low side 

# build projection profiles C on the low side
Cx = Source_struct()
ny_Ex = Int(W/dx)
Cx.pos = [[1,pml_npixels+1,ny_Ex,1]]
# sqrt(nu)*conj(f_x(m)) serves as a projection profile for a propagating channel, where nu = sin(kzdx)
# here we project the result onto all propagating channels as a output basis 
C_low = (conj(channels_low.f_x_m(channels_low.kydx_prop))).*reshape(channels_low.sqrt_nu_prop,1,:).*reshape(exp.((-1im*0.5)*channels_low.kzdx_prop),1,:) # 0.5 pixel backpropagation indicates that the projection(detect) region is half a pixel away from z = 0
Cx.data = [C_low]

# Calculate output channel amplitude coefficients (w) for the propagating channels from the point source through APF
w, _ = mesti(syst, [Bx], [Cx])

# We can also compute the field profile Ex_field and project the field on the detect region onto the projection profiles. 
# It would generate the same output channel amplitude coefficients
# Ex_field, _ = mesti(syst, [Bx])
# C = transpose(C_low)
# w = C_low*Ex_field[:,pml_npixels+1]
```
```text:Output
===System size===
ny_Ex = 5400; nz_Ex = 1391 for Ex(y,z)
UPML on -z +z sides; ; yBC = periodic; zBC = PEC
Building B,C... elapsed time:   1.664 secs
Building A  ... elapsed time:   8.106 secs
< Method: APF using MUMPS in single precision with AMD ordering >
Building K  ... elapsed time:   2.822 secs
false
Analyzing   ... elapsed time:   5.076 secs
Factorizing ... elapsed time:  76.478 secs
          Total elapsed time:  99.805 secs
```

# Compute the regular focusing and phase-conjugated focusing profiles

```julia
# Specify the system for mesti2s() and mesti_build_channels()
syst = Syst()
syst.epsilon_xx = epsilon
syst.length_unit  = "lambda_0"
syst.wavelength = 1
syst.dx = dx
syst.yBC = yBC
syst.epsilon_low = epsilon_low
syst.epsilon_high = epsilon_high
syst.zPML = [pml]

# equivalent average epsilon for this disordered system
epsilon_ave = mean(epsilon)

# build channels for the equivalent average epsilon
channels_ave_epsilon = mesti_build_channels(Int(W/dx), yBC, 2*pi*syst.wavelength*dx, epsilon_ave)
N_prop_ave_epsilon = channels_ave_epsilon.N_prop # number of propagating channels on the equivalent average epsilon 

# regular focus wavefront
wf_reg_focus = exp.(-1im*channels_ave_epsilon.kydx_prop*(m0_focus)) .* exp.(-1im*channels_ave_epsilon.kzdx_prop*(l0_focus-0.5)) # 0.5 pixel indicates that the source is half a pixel away from z = 0

# specify two input incident wavefronts:
# (1) regular focusing wavefront
# (2) phase-conjugated focusing wavefront
input = wavefront()
input.v_low = zeros(ComplexF64, N_prop_low, 2)
input.v_low[:, 1] = wf_reg_focus[Int((N_prop_ave_epsilon-N_prop_low)/2+1):Int(end-(N_prop_ave_epsilon-N_prop_low)/2)]/norm(wf_reg_focus[Int((N_prop_ave_epsilon-N_prop_low)/2+1):Int(end-(N_prop_ave_epsilon-N_prop_low)/2)])

# for the phased conjugated input: 
# conj(coefficient*u) = conj(coefficient)*conj(u) = conj(coefficient)*perm(u(ky)) = perm(conj(coefficient))*u(ky)
# perm() means permute a vector that switches one propagating channel with one
# having a complex-conjugated transverse profile
# for the periodic boundary, this flips the sign of ky. 
input.v_low[:, 2] = conj(w)[channels_low.ind_prop_conj]/norm(w)

# we will also get the field profile in the free spaces on the two sides, for
# plotting purpose.
opts = Opts()
opts.nz_low = round((L_tot-L)/2/dx)
opts.nz_high = opts.nz_low

# for field-profile computations
Ex, _, _ = mesti2s(syst, input, opts)
```
```text:Output
===System size===
ny_Ex = 5400; nz_Ex = 1349 => 1391 for Ex(y,z)
[N_prop_low, N_prop_high] = [725, 725] per polarization
yBC = periodic; zBC = [PML, PML]
Building B,C... elapsed time:   0.801 secs
            ... elapsed time:   0.807 secs
Building A  ... elapsed time:   5.320 secs
< Method: factorize_and_solve using MUMPS in single precision with AMD ordering >
Analyzing   ... elapsed time:   3.084 secs
Factorizing ... elapsed time:  75.164 secs
Solving     ... elapsed time:   7.529 secs
            ... elapsed time:  25.039 secs
          Total elapsed time: 122.167 secs
```

# Animate the field profiles and compare the intensity profiles

```julia
using Plots
# normalize the field amplitude with respect to the phase-conjugated-input profile
Ex = Ex/maximum(abs.(Ex[:,:,2]))

nframes_per_period = 10

# extend the x coordinate to include free spaces on the two sides
z_Ex = vcat(z_Ex[1] .- (opts.nz_low:-1:1)*dx, z_Ex, z_Ex[end] .+ (1:opts.nz_high)*dx)

# animate the field profile with the regular focusing input
anim_regular_focusing = @animate for ii ∈ 0:(nframes_per_period-1)
    plt1 = (heatmap(z_Ex, collect(y_Ex), real.(Ex[:,:,1]*exp(-1im*2*π*ii/nframes_per_period)),
             c = :balance, clims=(-1, 1), legend = :none,
             aspect_ratio=:equal, dpi = 450, 
             ticks = false, framestyle = :none,
             xlimits=(-90,180), ylimits=(0,360)))
    scatter!(plt1, z0_list, y0_list, markersize=r0_list, alpha=0.3, 
             color=:black, legend=false, dpi = 450)
end
gif(anim_regular_focusing, "regular_focusing.gif", fps = 5)
```

![regular_focusing.gif](https://github.com/complexphoton/MESTI.jl/assets/44913081/f2078aef-ff6d-42d3-a056-86dc72aa753b)

```julia
# animate the field profile of the phase-conjugated focusing input
anim_phase_congjuation = @animate for ii ∈ 0:(nframes_per_period-1)
    plt2 = (heatmap(z_Ex, collect(y_Ex), real.(Ex[:,:,2]*exp(-1im*2*π*ii/nframes_per_period)),
             c = :balance, clims=(-1, 1), legend = :none,
             aspect_ratio=:equal, dpi = 450, 
             ticks = false, framestyle = :none,
             xlimits=(-90,180), ylimits=(0,360)))
    scatter!(plt2, z0_list, y0_list,markersize=r0_list, alpha=0.3, 
             color=:black, legend=false, dpi = 450)   
end
gif(anim_phase_congjuation_focusing, "phase_conjugated_focusing.gif", fps = 5)
```

![phase_conjugated_focusing.gif](https://github.com/complexphoton/MESTI.jl/assets/44913081/a1ddcb6a-f43a-44c1-a0b0-6a611e57bd98)

```julia
# plot the intensity profiles and compare them
# limit the plotting to the small region around the center of the focusing region between y ∈ [175, 185] and z ∈ [40, 50]
y_Ex_ind_focusing_range = searchsortedfirst(y_Ex, 175)-1:searchsortedfirst(y_Ex, 185)
z_Ex_ind_focusing_range  = searchsortedfirst(z_Ex, 40)-1:searchsortedfirst(z_Ex, 50)

# find the index for the scatters within the plotting range
scatterer_ind_focusing_range = (y0_list .+ r0_list) .>= 175 .&& (y0_list .- r0_list) .<= 185 .&& (z0_list .+ r0_list) .>= 40 .&& (z0_list .- r0_list) .<= 50
y0_list_focusing_range = y0_list[scatterer_ind_focusing_range]
z0_list_focusing_range = z0_list[scatterer_ind_focusing_range]
r0_list_focusing_range = r0_list[scatterer_ind_focusing_range]

theta = range(0, stop=2π, length=100)

plt3 = heatmap(z_Ex[z_Ex_ind_focusing_range], y_Ex[y_Ex_ind_focusing_range], 
               abs.(Ex[y_Ex_ind_focusing_range, z_Ex_ind_focusing_range, 1]).^2,
               xlabel="z", ylabel="y", title="Regular focusing", 
               c=cgrad(:copper, rev=false), clims=(0, 1), aspect_ratio=:equal, dpi=600)

               for i in 1:length(r0_list_focusing_range)
    y0 = y0_list_focusing_range[i]
    z0 = z0_list_focusing_range[i]
    r0 = r0_list_focusing_range[i]

    z0_circle = z0 .+ r0 * cos.(theta)
    y0_circle = y0 .+ r0 * sin.(theta)

    plot!(plt3, z0_circle, y0_circle, lw=0.5, color=:white, legend=false, 
          xlims=(40, 50), ylims=(175, 185))
end
            
plt4 = heatmap(z_Ex[z_Ex_ind_focusing_range], y_Ex[y_Ex_ind_focusing_range], 
                 abs.(Ex[y_Ex_ind_focusing_range, z_Ex_ind_focusing_range, 2]).^2,
                 xlabel="z", ylabel="y", title="Phase conjugated focusing", 
                 c=cgrad(:copper, rev=false), clims=(0, 1), aspect_ratio=:equal, dpi=600)

for i in 1:length(r0_list_focusing_range)
    y0 = y0_list_focusing_range[i]
    z0 = z0_list_focusing_range[i]
    r0 = r0_list_focusing_range[i]

    z0_circle = z0 .+ r0 * cos.(theta)
    y0_circle = y0 .+ r0 * sin.(theta)

    plot!(plt4, z0_circle, y0_circle, lw=0.5, color=:white, legend=false, 
          xlims=(40, 50), ylims=(175, 185))
end

intensity_plot = plot(plt3, plt4, layout = @layout([a b]), size=(800, 400))
display(intensity_plot)
```

![intensity_comparison](https://github.com/complexphoton/MESTI.jl/assets/44913081/e98d33c5-bd5a-43a6-b146-60fdc7ee95c0)

```julia
# compare the ratio of intensity on the focus point
println("I_phase_congugation(y_0,z_0)/I_reg(y_0,z_0) = ", @sprintf("%d", round(abs.(Ex[m0_focus,opts.nz_low+l0_focus,2]).^2/abs.(Ex[m0_focus,opts.nz_low+l0_focus,1]).^2, digits=-2)))
```
```text:Output
I_phase_congugation(y_0,z_0)/I_reg(y_0,z_0) = 1700
```
