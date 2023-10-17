# This code test the unitarity of scattering matrices S upon random permittivity profiles.
# We check S'*S ≈ I, an identity matrix.

# Specify parameters of the system
syst = Syst()
syst.xBC = "periodic"
syst.yBC = "periodic"
syst.dx = 1
syst.wavelength = 5
syst.epsilon_low = 1
syst.epsilon_high = 1
epsilon_max = 4
epsilon_min = 1

# Define the size of the scattering region (1 wavelength by 1 wavelength by 1 wavelength)
nx, ny, nz = 2, 5, 2
nx_Ex = nx; ny_Ex = ny; nz_Ex = nz -1
nx_Ey = nx; ny_Ey = ny; nz_Ey = nz -1
nx_Ez = nx; ny_Ez = ny; nz_Ez = nz

# Use optimized PML parameters for this resolution to reduce error
zpml = get_optimal_PML(syst.wavelength/syst.dx)
zpml.npixels = 25
syst.zPML = [zpml]

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from both sides with both s-polarization and p-polarization
input.side = "both"
input.polarization = "both"
# Output to both sides with both s-polarization and p-polarization
output.side = "both"
output.polarization = "both"

# Test the functionality in a test set
@testset "unitarity(diagonal ε):  " begin
    for i ∈ 1:3
	# Random permittivity profiles, whose value is between 1 and 4
        # The condition for the scattering region to be lossless is the permittivity tensor ε to be hermitian: adjoint(ε) = ε
        syst.epsilon_xx = rand(nx_Ex,ny_Ex,nz_Ex)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_yy = rand(nx_Ey,ny_Ey,nz_Ey)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_zz = rand(nx_Ez,ny_Ez,nz_Ez)* (epsilon_max-epsilon_min) .+ epsilon_min

	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.
        (S, _, _)= mesti2s(syst, input, output, opts)
        @test maximum(abs.((S'*S) - I(size(S, 1)))) ≤ 1e-6 # Check the unitarity of the scattering matrix
    end
end
                    
@testset "unitarity(general ε):   " begin
    for i ∈ 1:3
	# Random permittivity profiles, whose value is between 1 and 4
        # The condition for the scattering region to be lossless is the permittivity tensor ε to be hermitian: adjoint(ε) = ε
        syst.epsilon_xx = rand(nx_Ex,ny_Ex,nz_Ex)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_yy = rand(nx_Ey,ny_Ey,nz_Ey)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_zz = rand(nx_Ez,ny_Ez,nz_Ez)* (epsilon_max-epsilon_min) .+ epsilon_min
        syst.epsilon_xy = rand(nx_Ez,ny_Ez,nz_Ex) * (epsilon_max-epsilon_min) .+ epsilon_min + 1im*(rand(nx_Ez,ny_Ez,nz_Ex) * (epsilon_max-epsilon_min) .+ epsilon_min)
        syst.epsilon_yx = conj(syst.epsilon_xy)
        syst.epsilon_yz = rand(nx_Ey,ny_Ex,nz_Ex) * (epsilon_max-epsilon_min) .+ epsilon_min + 1im*(rand(nx_Ez,ny_Ez,nz_Ex) * (epsilon_max-epsilon_min) .+ epsilon_min)
        syst.epsilon_zy = conj(syst.epsilon_yz)
        syst.epsilon_xz = rand(nx_Ey,ny_Ex,nz_Ey) * (epsilon_max-epsilon_min) .+ epsilon_min + 1im*(rand(nx_Ez,ny_Ez,nz_Ex) * (epsilon_max-epsilon_min) .+ epsilon_min)
        syst.epsilon_zx = conj(syst.epsilon_xz)

	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.

        (S, _, _)= mesti2s(syst, input, output, opts)
	@test maximum(abs.((S'*S) - I(size(S, 1)))) ≤ 1e-6 # Check the unitarity of the scattering matrix
    end
end
