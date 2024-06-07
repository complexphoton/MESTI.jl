"""
    MESTI_SET_PML_PARAMS sets up perfectly matched layer (PML) parameters.
        If users already specifies (partial) PML parameters, we check and keep them.
        If some PML parameters are not defined, these parameters are set to 
        optimized values based on resolution and the background refractive index.
        More details are included in the upcoming paper Z. Wang et al (in preparation).
        === Input Arguments ===
        pml (a vector of PML structure; required):
            Parameters for perfectly matched layer (PML).
            A vector contains two PML structure and each PML structure represents 
            PML parameters on minus and plus background side in given direction. 
            PML is used to simulate an open boundary, which attenuates outgoing 
            waves with minimal reflection. 
            Each PML scalar structure has the following fields:    
                npixels (positive integer scalar; required): Number of PML pixels.
                    This number of pixels is added in addition to the
                    scattering region.
                power_sigma (non-negative scalar; optional): 
                    Power of the polynomial grading for the conductivity sigma.
                sigma_max_over_omega (non-negative scalar; optional):
                    Conductivity at the end of the PML. This is used to attenuate propagating waves.
                power_kappa (non-negative scalar; optional): 
                    Power of the polynomial grading for the real-coordinate-stretching 
                    factor kappa.
                kappa_max (real scalar no smaller than 1; optional):
                    Real-coordinate-stretching factor at the end of the PML. This is used to 
                    accelerate the attenuation of evanescent waves. kappa_max = 1 means no 
                    real-coordinate stretching.
                power_alpha (non-negative scalar; optional): 
                    Power of the polynomial grading for the CFS alpha factor.
                alpha_max_over_omega (non-negative scalar; optional): 
                    Complex-frequency-shifting (CFS) factor at the beginning 
                    of the PML. This is typically used in time-domain simulations 
                    to suppress late-time (low-frequency) reflections. 
                    We don't use it by default (alpha_max_over_omega = 0) 
                    since we are in frequency domain.
            We use the following PML coordinate-stretching factor:
                s(p) = kappa(p) + sigma(p)./(alpha(p) - i*omega)
            with
                sigma(p)/omega = sigma_max_over_omega*(p.^power_sigma),
                kappa(p) = 1 + (kappa_max-1)*(p.^power_kappa),
                alpha(p)/omega = alpha_max_over_omega*((1-p).^power_alpha),
            where omega is frequency, and p goes linearly from 0 at the beginning
            of the PML to 1 at the end of the PML. 
                With real-coordinate stretching, PML can attenuate evanescent waves
            more efficiently than free space, so there is no need to place free
            space in front of PML.
                The PML thickness should be chosen based on the acceptable level of
            reflectivity given the discretization resolution and the range of wave
            numbers (i.e., angles) involved; more PML pixels gives lower
            reflectivity. Typically 10-40 pixels are sufficient.
        k0dx (numeric scalar, real or complex; required):
            Normalized frequency k0*dx = (2*pi/vacuum_wavelength)*dx.
        epsilon_bg (numeric scalar, real or complex; required):
            Relative permittivity of the background space.
        direction (string; required):
            The specific direction where the PMLs are padded. 
        === Output Arguments ===
        pml (a vector of PML structure; required)
            Parameters for perfectly matched layer (PML).
            A vector contains two PML structure and each PML structure represents 
            PML parameters on minus and plus background side in given direction.
"""

function mesti_set_PML_params(pml::Vector{PML}, k0dx::Union{Real,Complex}, epsilon_bg::Union{Vector{Int64},Vector{Float64}}, direction::String)
    
    if ~(length(pml) == 2)
        throw(ArgumentError("$(direction)PML must be a vector of PML containing two elements."))
    end

    # No PML on the both sides
    if (pml[1].npixels == 0 && pml[2].npixels == 0)
        return pml
    end

    # loop over PML parameters on two sides
    for ii = 1:2
        wavelength_over_dx = ((2*pi)/k0dx)/sqrt(epsilon_bg[ii])
	# Number of PML pixels
        if ~isdefined(pml[ii], :npixels)
            throw(ArgumentError("$(direction)PML[$(ii)] must contain field npixels."))
        else
            temp = pml[ii].npixels
            if temp < 0
                throw(ArgumentError("$(direction)PML[$(ii)].npixels must be a non-negative scalar."))             
            end
            if temp == 0
                continue
            end
        end

        # Power of polynomial grading for the conductivity sigma
        if ~isdefined(pml[ii], :power_sigma)
            pml[ii].power_sigma = 3
        else
            temp = pml[ii].power_sigma
            if temp < 0
                throw(ArgumentError("$(direction)pml[$(ii)].power_sigma must be a non-negative scalar."))
            end
        end

        # Conductivity at the end of the PML
        if ~isdefined(pml[ii], :sigma_max_over_omega)
            a = -3.0138; b = 0.9303; c = 1.0128
            pml[ii].sigma_max_over_omega = a + b * wavelength_over_dx^c
        else
            temp = pml[ii].sigma_max_over_omega            
            if temp < 0
                throw(ArgumentError("$(direction)pml[$(ii)].sigma_max_over_omega must be a non-negative scalar."))
            end
        end

        # Maximal coordinate stretching factor
        if ~isdefined(pml[ii], :kappa_max)
            a = -2.0944; b = 0.6617; c = 1.0467
            pml[ii].kappa_max = a + b * wavelength_over_dx^c
        else
            temp = pml[ii].kappa_max
            if temp < 1
                throw(ArgumentError("$(direction)pml[$ii].kappa_max must be a real scalar that equals or is larger than 1."))
            end
        end

        # Power of polynomial grading for the coordinate stretching factor kappa
        if ~isdefined(pml[ii], :power_kappa)
            pml[ii].power_kappa = 3
        else
            temp = pml[ii].power_kappa
            if temp < 0
                throw(ArgumentError("$(direction)pml[$ii].power_kappa must be a non-negative scalar."))    
            end
        end

        # Maximal alpha factor for complex frequency shifting (CFS)
        if ~isdefined(pml[ii], :alpha_max_over_omega)
            # CFS is meant to suppress reflection of low-frequency components for time-domain simulations; it is not necessary for frequency-domain simulations
            pml[ii].alpha_max_over_omega = 0
        else
            temp = pml[ii].alpha_max_over_omega
            if temp < 0
                throw(ArgumentError("$(direction)pml[$ii].alpha_max_over_omega must be a non-negative scalar."))
            end
        end

        # Power of polynomial grading for the alpha factor
        if ~isdefined(pml[ii], :power_alpha)
            # not relevant when alpha_max_over_omega = 0
            pml[ii].power_alpha = 1
        else
            temp = pml[ii].power_alpha
            if temp < 0
                throw(ArgumentError("$(direction)pml[$ii].power_alpha must be a non-negative scalar."))   
            end
        end         

    end
    
    return pml       
end
