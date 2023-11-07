# This code test the reflection coefficient r and transmission coefficient t of a 1D single interface system,
# where material 1 is on the low side and material 2 is on the high side.
# Their refractive index is n1 and n2, respectively.
# We compare the numerical result and the analytic restuls.

# Specify system
syst = Syst()
syst.yBC = "periodic" # 1D system at normal incidence has periodic boundary in y
resolution = 10
syst.wavelength = 1 # vacuum wavelength  
syst.dx = 1/resolution # grid size

# Use optimized PML parameters for this resolution to reduce error
zpml = get_optimal_PML(syst.wavelength/syst.dx)
zpml.npixels = 30
syst.zPML = [zpml] 
k0dx = 2*pi/syst.wavelength*syst.dx # Dimensionless frequency k0*dx

# Specify inputs and output
input = channel_type()
output = channel_type()
# Input from the low side
input.side = "low"
# Output to the high side
output.side = "both"

# Test the functionality in a test set
@testset "1D r and t coefficients:" begin
    for i ∈ 1:6
	# Make an interface for two materials whose refractive indices are n1 and n2
        n1 = 1+1.5*rand() # n1 is a random number between 1 and 2.5
        n2 = 1+1.5*rand() # n2 is a random number between 1 and 2.5
	
        syst.epsilon_low = n1^2 # relative permittivity on the low side
	syst.epsilon_high = n2^2 # relative permittivity on the high side

        # For 1D system, we only 1 pixel in transverse direction (y) with a periodic boundary, 
        # so the relative permittivity profile is translationally invariant in y
        # In longitudnal direction (z), 0 pixel means no scattering region.
        # So the system is a single interface system with two homogenous spaces whose refrative indices are n1 and n2
        syst.epsilon_xx = ones(1,0)  

	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.
        
	(S, channels, _)= mesti2s(syst, input, output, opts) # Calculate a scattering matrix, S = [r;t]

        kzdx_1 = channels.low.kzdx_prop[1] # Dimensionless longitudinal wave number  for material 1
	kzdx_2 = channels.high.kzdx_prop[1] # Dimensionless longitudinal wave number for material 2
	# Analytic expression for reflection coefficient r of a 1D single interface system
        r = (exp(-1im*kzdx_1/2)*exp(1im*kzdx_2)-exp(1im*kzdx_1/2))/(exp(-1im*kzdx_1/2)-exp(1im*kzdx_1/2)*exp(1im*kzdx_2))
	# Analytic expression for transmission coefficient t of a 1D single interface system
	t = (sqrt(sin(kzdx_2))/sqrt(sin(kzdx_1)))*(exp(1im*kzdx_1*3/2)-exp(-1im*kzdx_1/2))/(exp(1im*kzdx_2/2)*exp(1im*kzdx_1)-exp(-1im*kzdx_2/2))
	@test abs(S[1,1]-r) ≤ 1e-4 # Test absolute error of r 
	@test abs(S[2,1]-t) ≤ 1e-4 # Test absolute error of t
    end
end
