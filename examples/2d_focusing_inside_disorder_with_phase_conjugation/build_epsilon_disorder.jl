using Random

function build_epsilon_disorder(W, L, r_min, r_max, min_sep, number_density, rng_seed, dx, epsilon_scat, epsilon_bg, build_TM, build_TE = false, yBC = "periodic", y1 = 0, y2 = W, z1 = 0, z2 = L; no_scatterer_center = false)
    # Generate a random collection of cylinders and build the relative permittivity profile.
    # The full system spans y in [0, W], z in [0, L].
    # The scattering region spans [z1, z2], [y1, y2].
    # Note that subpixel smoothing is not used.
    #   
    #   === Input Arguments ===
    #   W       = width of the full system
    #   L       = thickness of the full system
    #   r_min   = minimal radius of the cylinders
    #   r_max   = maximal radius of the cylinders
    #   min_sep = minimal separation between cylinders
    #   number_density = number density of the cylinders
    #   rng_seed = random number generator seed
    #   dx       = discretization grid size
    #   epsilon_scat = relative permittivity of the cylinders
    #   epsilon_bg   = relative permittivity of the background
    #   build_TM = whether to build epsilon for TM polarization
    #   build_TE = whether to build 1/epsilon for TE polarization
    #   yBC = boundary condition in y (required only when build_TE = true)
    #   y1 (optional) = bottom edge of the scattering region (defaults to 0)
    #   y2 (optional) = top edge of the scattering region (defaults to W)
    #   z1 (optional) = left edge of the scattering region (defaults to 0)
    #   z2 (optional) = right edge of the scattering region (defaults to L)
    #   no_scatterer_center (optional) = whether to avoid put scatterer on the center of the system
    #   === Output Arguments ===
    #   epsilon_xx (ny_Ex-by-nz_Ex matrix): discretized epsilon_zz
    #   inv_epsilon_yy (ny_Hx-by-nz_Ex matrix): discretized (1/epsilon)_yy
    #   inv_epsilon_zz (ny_Ex-by-nz_Hx matrix): discretized (1/epsilon)_zz
    #   inv_epsilon_yz (ny_Ex-by-nz_Ex matrix): discretized (1/epsilon)_yz
    #   y0_list (column vector): y coordinates of the cylinder centers
    #   z0_list (column vector): z coordinates of the cylinder centers
    #   r0_list (column vector): radii of the cylinders
    #   y_Ex (row vector): y coordinates of the centers of the Ex pizels
    #   z_Ex (row vector): z coordinates of the centers of the Ex pixels
    #   y_Hx (row vector): y coordinates of the centers of the Hx pixels
    #   z_Hx (row vector): z coordinates of the centers of the Hx pixels

    # round the system size
    W = round(Int, W/dx) * dx
    L = round(Int, L/dx) * dx

    if y1 < 0 || y2 > W || z1 < 0 || z2 > L
        error("The ranges [y1, y2] and [z1, z2] are invalid.")
    end

    N_scatterers = round(Int, number_density * (y2 - y1) * (z2 - z1))
    if no_scatterer_center
        N_scatterers = N_scatterers + 1
    end

    if N_scatterers > 3000
        # store the indices of the existing cylinders by their approximate locations,
        # to facilitate overlap checking when N_scatterers is large
        use_cell = true
        cell_size = 2 * r_max + min_sep
        nz_cell = ceil(Int, (z2 - z1) / cell_size)
        ny_cell = ceil(Int, (y2 - y1) / cell_size)
        ind_scatterers = Array{Any}(nothing, ny_cell, nz_cell)
    else
        use_cell = false
    end

    # pick the radii and locations of N_scatterers non-overlapping cylinders,
    # and then generate the discretized permittivity profile
    Random.seed!(rng_seed)  # set random number generator seed
    y0_list, z0_list, r0_list = zeros(N_scatterers), zeros(N_scatterers), zeros(N_scatterers)

    for ii in 1:N_scatterers
        # randomly pick the radius of the new cylinder
        r0 = r_min + rand() * (r_max - r_min)

        if use_cell
            # distance for overlap checking
            dist_check = r0 + min_sep + r_max
        end

        # min/max values for the coordinates of the cylinder center
        y_min = y1 + r0
        y_max = y2 - r0
        z_min = z1 + r0
        z_max = z2 - r0

        # keep randomly picking the cylinder location until we get one that doesn't overlap with the existing cylinders
        found = false

        while !found
            global y0
            global z0
            if no_scatterer_center && ii == 1
                y0 = (y1+y2)/2
                z0 = (z1+z2)/2
                found = true
            else
                # center of the new cylinder
                y0 = y_min + rand() * (y_max - y_min)
                z0 = z_min + rand() * (z_max - z_min)

                if use_cell
                    # indices of cells we need to check for overlap
                    l1 = max(1,       round(Int, 0.5 + (z0 - z1 - dist_check) / cell_size))
                    l2 = min(nz_cell, round(Int, 0.5 + (z0 - z1 + dist_check) / cell_size))
                    m1 = max(1,       round(Int, 0.5 + (y0 - y1 - dist_check) / cell_size))
                    m2 = min(ny_cell, round(Int, 0.5 + (y0 - y1 + dist_check) / cell_size))

                    # indices of cylinders we need to check for overlap
                    ind_check = ind_scatterers[m1:m2, l1:l2][ind_scatterers[m1:m2, l1:l2] .!= nothing]
                    if ~isempty(ind_check)
		    	ind_check = reduce(vcat, ind_check)
	            end

                    # check those cylinders for overlap
                    if isempty(ind_check) || minimum((z0_list[ind_check] .- z0).^2 .+ (y0_list[ind_check] .- y0).^2 .- (r0_list[ind_check] .+ (r0 + min_sep)).^2) > 0
                        found = true
                    end
                else
                    # check all existing cylinders for overlap
                    if ii == 1 || minimum((y0_list[1:(ii-1)] .- y0).^2 .+ (z0_list[1:(ii-1)] .- z0).^2 .- (r0_list[1:(ii-1)] .+ (r0 + min_sep)).^2) > 0
                        found = true
                    end
                end
            end
        end
        
        y0_list[ii] = y0
        z0_list[ii] = z0
        r0_list[ii] = r0

        if use_cell
            # store the index of this cylinder by its approximate location
            m0 = min(ny_cell, max(1, round(Int, 0.5 + (y0 - y1) / cell_size)))
            l0 = min(nz_cell, max(1, round(Int, 0.5 + (z0 - z1) / cell_size)))
            if isa(ind_scatterers[m0, l0], Nothing)
                ind_scatterers[m0, l0] = [ii]
            else
                push!(ind_scatterers[m0, l0], ii)
            end
        end
    end

    if no_scatterer_center
        y0_list = y0_list[2:end]
        z0_list = z0_list[2:end]
        r0_list = r0_list[2:end]    
    end
    
    # Prepare inputs for mesti_subpixel_smoothing()
    domain = Cuboid([W/2,L/2], [W,L]) # domain for subpixel smoothing cetering at (W/2, L/2) and with witdth W and thickness L
    domain_epsilon = epsilon_bg # epsilon of the domain for subpixel smoothing
    object_list::Vector{Shape}= [] # object list
    object_epsilon_list::Vector{Float64}=[] # object epsilon list

    # Put the scatters in the object_list and object_epsilon_list 
    for jj = 1:length(r0_list)
        push!(object_list, Ball([y0_list[jj], z0_list[jj]], r0_list[jj])) 
        push!(object_epsilon_list, epsilon_scat)
    end

    # Call mesti_subpixel_smoothing() and return the results
    if build_TM && build_TE
        epsilon_xx, inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz = mesti_subpixel_smoothing(dx, domain, domain_epsilon, object_list, object_epsilon_list, yBC, "PEC", build_TM, build_TE)
        ny_Ex, nz_Ex = size(epsilon_xx); ny_Hx = size(inv_epsilon_yy, 1); nz_Hx = size(inv_epsilon_yy, 2)
        y_Ex = dx*(0:(ny_Ex-1)); z_Ex = dx*(1:nz_Ex); y_Hx = dx*(0.5:ny_Hx); z_Hx = dx*(0.5:nz_Hx)
        return epsilon_xx, inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz, y0_list, z0_list, r0_list, y_Ex, z_Ex, y_Hx, z_Hx
    end

    if build_TM && ~build_TE
        epsilon_xx = mesti_subpixel_smoothing(dx, domain, domain_epsilon, object_list, object_epsilon_list, yBC, "PEC", build_TM, build_TE)
        ny_Ex, nz_Ex = size(epsilon_xx)
        y_Ex = dx*(0:(ny_Ex-1)); z_Ex = dx*(1:nz_Ex)        
        return epsilon_xx, y0_list, z0_list, r0_list, y_Ex, z_Ex
    end

    if ~build_TM && build_TE
        inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz  = mesti_subpixel_smoothing(dx, domain, domain_epsilon, object_list, object_epsilon_list, yBC, "PEC", build_TM, build_TE)
        ny_Hx = size(inv_epsilon_yy, 1); nz_Hx = size(inv_epsilon_yy, 2)
        return inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz, y0_list, z0_list, r0_list, y_Hx, z_Hx
        y_Hx = dx*(0.5:ny_Hx); z_Hx = dx*(0.5:nz_Hx)
    end

end
