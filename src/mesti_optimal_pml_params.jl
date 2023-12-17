# Export a function mesti_optimal_pml_params()
export mesti_optimal_pml_params

function mesti_optimal_pml_params(wavelength_over_dx)
    # Return an optimal PML parameters from double-log FOM given the mesh resolution wavelength_over_dx
    # wavelength_over_dx should be in the range of [5, 500].

    if wavelength_over_dx < 5 || wavelength_over_dx > 500
	@warn "The resolution $(wavelength_over_dx), is beyond the range of 5 to 500. Since we fit our optimal PML parameters within this specified range, the obtained PML parameters may not be optimal."
    end

    # Initialize mutable struct PML  
    pml = PML(10)

    # Get sigma_max_over_omega
    a = -26.33; b = 26.65; c = 2.781e-2
    pml.sigma_max_over_omega = exp(a + b * wavelength_over_dx^c)

    # Get power_sigma
    a = 1.9; b = 200; c = 45; d = 1.4
    pml.power_sigma = a + b / (c + wavelength_over_dx^d)

    # Get kappa_max
    a = -6.667; b = 7.285; c = 8.308e-2
    pml.kappa_max = exp(a + b * wavelength_over_dx^c)

    # Get power_kappa
    a = 2.0; b = 2.4e4; c = 2.3e3; d = 2.9
    pml.power_kappa = a + b / (c + wavelength_over_dx^d)

    # Set alpha to zero
    pml.alpha_max_over_omega = 0
    pml.power_alpha = 1

    return pml
end
