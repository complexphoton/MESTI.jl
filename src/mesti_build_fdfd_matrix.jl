# Export a composite data type PML
export PML

# Export a function mesti_build_fdfd_matrix()
export mesti_build_fdfd_matrix

mutable struct PML
    npixels::Integer
    sigma_max_over_omega::Real
    power_sigma::Real       
    alpha_max_over_omega::Real 
    power_alpha::Real
    kappa_max::Real
    power_kappa::Real 
    
    # Below are only used in mesti() and mesti2s() 
    direction::String
    side::String
    
    # Below is only used in mesti2s() 
    npixels_spacer::Union{Integer,Nothing}
    
    # Construct default parameters from Table 7.1 of Taflove & Hagness's 2005 FDTD book
    # Without specifying pml.npixels
    PML() = (pml = new(); pml.power_sigma=3.0; pml.alpha_max_over_omega=0.0; 
              pml.power_alpha=1.0; pml.kappa_max=15.0; pml.power_kappa=3.0; return pml)
    # With specifying pml.npixels   
    PML(n) = (pml = new(); pml.npixels=n; pml.power_sigma=3.0; pml.alpha_max_over_omega=0.0; 
              pml.power_alpha=1.0; pml.kappa_max=15.0; pml.power_kappa=3.0; return pml) 
end

"""
    MESTI_BUILD_FDFD_MATRIX The finite-difference frequency-domain operator in 3D.
       A = mesti_build_fdfd_matrix(epsilon_xx, epsilon_xy, epsilon_xz, epsilon_yx, epsilon_yy, epsilon_yz, epsilon_zx, epsilon_zy, epsilon_zz, k0dx, xBC, yBC, zBC, xPML, yPML, zPML, use_UPML) returns A as a sparse matrix representing 
       wave operator [curl(curl) - k0^2*epsilon_r(x,y,z)]*(dx^2)
       discretized on a square grid with grid size dx through center difference.
       for the electric field: (Ex, Ey, Ez).
    
       === Input Arguments ===
       epsilon_xx (array; required):
          epsilon_xx is a nx_Ex-by-ny_Ex-by-nz_Ex numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_xx on Ex-sites. 
          The notation ni_Ej is the total number of grids for Ej in i-direction.
       epsilon_xy (array; optional):
          epsilon_xy is a nx_Ez-by-ny_Ez-by-nz_Ex numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_xy on lower corner of the Yee lattice.
       epsilon_xz (array; optional):
          epsilon_xz is a nx_Ey-by-ny_Ex-by-nz_Ey numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_xz on lower corner of the Yee lattice. 
       epsilon_yx (array; optional):
          epsilon_yx is a nx_Ez-by-ny_Ez-by-nz_Ey numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_yx on lower corner of the Yee lattice.
       epsilon_yy (array; required):
          epsilon_yy is a nx_Ey-by-ny_Ey-by-nz_Ey numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_yy on Ey-sites. 
       epsilon_yz (array; optional):
          epsilon_yz is a nx_Ey-by-ny_Ex-by-nz_Ex numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_yz on lower corner of the Yee lattice. 
       epsilon_zx (array; optional):
          epsilon_zx is a nx_Ey-by-ny_Ez-by-nz_Ey numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_zx on lower corner of the Yee lattice.
       epsilon_zy (array; optional):
          epsilon_zy is a nx_Ez-by-ny_Ex-by-nz_Ex numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_zy on lower corner of the Yee lattice.
       epsilon_zz (array; required):
          epsilon_zz is a nx_Ez-by-ny_Ez-by-nz_Ez numeric array (real or complex)
          discretizing the relative permittivity profile epsilon_zz on Ez-sites.
       k0dx (numeric scalar, real or complex; required):
          Normalized frequency k0*dx = (2*pi/vacuum_wavelength)*dx.
       xBC (character vector or numeric scalar; required):
             Boundary condition (BC) at the two ends in x direction, effectively
             specifying Ex(n,m,l) at n=0 and n=nx_Ex+1,
                        Ey(n,m,l) at n=0 and n=nx_Ey+1,
                        Ez(n,m,l) at n=0 and n=nx_Ez+1.    
             one pixel beyond the computation domain. Available choices are:
               "periodic" - Ex(n+nx_Ex,m,l) = Ex(n,m,l),
                            Ey(n+nx_Ey,m,l) = Ey(n,m,l),
                            Ez(n+nx_Ez,m,l) = Ez(n,m,l).   
               "PEC"      - Ex(0,m,l) = Ex(1,m,l); Ex(nx_Ex+1,m,l) = Ex(nx_Ex,m,l),
                            Ey(0,m,l) = Ey(nx_Ey+1,m,l) = 0,
                            Ez(0,m,l) = Ez(nx_Ez+1,m,l) = 0.   
               "PMC"      - Ex(0,m,l) = Ex(nx_Ex+1,m,l) = 0,
                            Ey(0,m,l) = Ey(1,m,l); Ey(nx_Ey+1,m,l) = Ey(nx_Ey,m,l),
                            Ez(0,m,l) = Ez(1,m,l); Ez(nx_Ez+1,m,l) = Ez(nx_Ez,m,l).    
               "PECPMC"   - Ex(0,m,l) = Ex(1,m,l); Ex(nx_Ex+1,m,l) = 0,
                            Ey(0,m,l) = 0; Ey(nx_Ey+1,m,l) = Ey(nx_Ey,m,l),
                            Ez(0,m,l) = 0; Ez(nx_Ez+1,m,l) = Ez(nx_Ez,m,l),    
               "PMCPEC"   - Ex(0,m,l) = 0; Ex(nx_Ex+1,m,l) = Ex(nx_Ex,m,l),
                            Ey(0,m,l) = Ey(1,m,l); Ey(nx_Ey+1,m,l) = 0,
                            Ez(0,m,l) = Ez(1,m,l); Ez(nx_Ez+1,m,l) = 0.
             where PEC stands for perfect electric conductor and PMC stands for perfect
             magnetic conductor.
             When xBC is a numeric scalar, the Bloch periodic boundary condition is
          used with f(n+nx,m,l) = f(n,m,l)*exp(1i*xBC), f=Ex,Ey,Ez; in other words, 
          xBC = kx_B*nx*dx = kx_B*Lambda where kx_B is the Bloch wave number and 
          Lambda = nx*dx is the periodicity in x. In this case, nx = nx_Ex = nx_Ey = nx_Ez.
       yBC (character vector or numeric scalar; required):
          Boundary condition in y direction, analogous to xBC.
       zBC (character vector or numeric scalar; required):
          Boundary condition in z direction, analogous to xBC.    
       xPML (a two-element PML vector; optional):
          Parameters for perfectly matched layer (PML) in x direction.
          xPML = [PML_left, PML_right] 
          If users do not want to put PML on either side, just set the field "npixels" = 0.
          For example, PML_left.npixels = 0: no PML on the left.
          In each case, PML is a PML structure with the following fields:
             npixels (non-negative integer scalar; required): Number of PML pixels.
                Note this is within syst.epsilon or syst.inv_epsilon.
             power_sigma (non-negative scalar; optional): Power of the polynomial 
                grading for the conductivity sigma; defaults to 3.
             sigma_max_over_omega (non-negative scalar; optional):
                Conductivity at the end of the PML; defaults to
                   0.8*(power_sigma+1)/(k0dx*sqrt(epsilon_bg)).
                where epsilon_bg is the average relative permittivity along the
                last slice of the PML. This is used to attenuate propagating waves.
             power_kappa (non-negative scalar; optional): Power of the polynomial
                grading for the real-coordinate-stretching factor kappa; defaults
                to 3.
             kappa_max (real scalar no smaller than 1; optional):
                Real-coordinate-stretching factor at the end of the PML; defaults
                to 15. This is used to accelerate the attenuation of evanescent 
                waves. kappa_max = 1 means no real-coordinate stretching.
             power_alpha (non-negative scalar; optional): Power of the polynomial
                grading for the CFS alpha factor; defaults to 1.
             alpha_max_over_omega (non-negative scalar; optional): Complex-
                frequency-shifting (CFS) factor at the beginning of the PML. This
                is typically used in time-domain simulations to suppress late-time
                (low-frequency) reflections. We do not use it by default 
                (alpha_max_over_omega = 0) since we are in frequency domain.
          We use the following PML coordinate-stretching factor:
             s(p) = kappa(p) + sigma(p)./(alpha(p) - i*omega)
          with
             sigma(p)/omega = sigma_max_over_omega*(p.^power_sigma),
             kappa(p) = 1 + (kappa_max-1)*(p.^power_kappa),
             alpha(p)/omega = alpha_max_over_omega*((1-p).^power_alpha),
          where omega is frequency, and p goes linearly from 0 at the beginning of
          the PML to 1 at the end of the PML. 
       yPML (a two-element vector; optional):
          Parameters for PML in y direction, analogous to xPML.
       zPML (a two-element vector; optional):
          Parameters for PML in z direction, analogous to xPML.
       use_UPML (logical scalar; optional, defaults to true):
          Whether to use uniaxial PML (UPML) or not. If not, stretched-coordinate
          PML (SC-PML) will be used.
    
       === Output Arguments ===
       A (sparse matrix):
          (nt_Ex+nt_Ey+nt_Ez)-by-(nt_Ex+nt_Ey+nt_Ez) sparse matrix representing 
          the 3D FDFD operator. nt_Ex/nt_Ey/nt_Ez is the total number of grids 
          for Ex/Ey/Ez.
       is_symmetric_A (logical scalar):
          Whether matrix A is symmetric or not.
       xPML (two-element vector):
          PML parameters used on the low and high sides of x direction, if any.
       yPML (two-element vector):
          PML parameters used on the low and high sides of y direction, if any.
       zPML (two-element vector):
          PML parameters used on the low and high sides of z direction, if any.
"""
function mesti_build_fdfd_matrix(epsilon_xx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Matrix{Int64},Matrix{Float64},Matrix{ComplexF64}}, epsilon_xy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_xz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_yx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_yy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_yz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_zx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_zy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing},     epsilon_zz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, k0dx::Union{Real,Complex}, xBC::Union{String,Real,Complex,Nothing}, yBC::Union{String,Real,Complex}, zBC::Union{String,Real,Complex}, xPML::Union{Vector{PML},Nothing} = [PML(0), PML(0)], yPML::Vector{PML} = [PML(0), PML(0)], zPML::Vector{PML} = [PML(0), PML(0)], use_UPML::Bool=true)    
    # Make deepcopy of them to avoid mutating input argument 
    xPML = deepcopy(xPML); yPML = deepcopy(yPML); zPML = deepcopy(zPML)
    
    # Take care of the 2D TM case
    if ndims(epsilon_xx) == 2
        use_2D_TM = true
        if ~isa(epsilon_yy, Nothing) || ~isa(epsilon_zz, Nothing) ||  ~isa(epsilon_xy, Nothing) || ~isa(epsilon_xz, Nothing) || ~isa(epsilon_yx, Nothing) || ~isa(epsilon_yz, Nothing) || ~isa(epsilon_zx, Nothing) || ~isa(epsilon_zy, Nothing)
            throw(ArgumentError("Only epsilon_xx is required for 2D TM fields Ex(y,z), but other components should not be given or should be nothing"))
        end
        if xBC != nothing
            @warn "Only yBC and zBC are required for 2D TM fields Ex(y,z). xBC will be ignored."
        end
        if xPML != nothing
            @warn "Only yPML and zPML are required for 2D TM fields Ex(y,z). xPML will be ignored."
        end
    else
        use_2D_TM = false
    end
 
    # Check the presence off-diagonal part in epsilon
    if ~isa(epsilon_xy, Nothing) || ~isa(epsilon_xz, Nothing) || ~isa(epsilon_yx, Nothing) || ~isa(epsilon_yz, Nothing) || ~isa(epsilon_zx, Nothing) || ~isa(epsilon_zy, Nothing)
        include_off_diagonal_epsilon = true
    else
        include_off_diagonal_epsilon = false        
    end
    
    if use_2D_TM # 2D TM field
        # Number of sites in y and z for Ex
        (ny_Ex, nz_Ex) = size(epsilon_xx)
        
        # Convert BC to take care of lowercase or uppercase
        yBC = convert_BC(yBC, "y")
        zBC = convert_BC(zBC, "z")
        
        # Total number of grid points for Ex
        nt_Ex = ny_Ex*nz_Ex

        # Checking the whether the permittivity is lossless or not 
        if (maximum(imag(epsilon_xx)) >= sqrt(eps()))
            @warn("syst.epsilon_xx contains imaginary part; the permittivity profile is not lossless.")
        end
        
    else # 3D case
        # Number of grid points in x, y, and z for Ex, Ey, and Ez
        (nx_Ex, ny_Ex, nz_Ex) = size(epsilon_xx)
        (nx_Ey, ny_Ey, nz_Ey) = size(epsilon_yy)
        (nx_Ez, ny_Ez, nz_Ez) = size(epsilon_zz)

        # Convert BC to take care of lowercase or uppercase
        xBC = convert_BC(xBC, "x")
        yBC = convert_BC(yBC, "y")
        zBC = convert_BC(zBC, "z")

        # Check number of grid points with the boundary conditions
        if nx_Ey != nx_Ez; throw(ArgumentError("Number of grids along x provided by epsilon_yy and epsilon_zz should be same.")); end
        if ny_Ex != ny_Ez; throw(ArgumentError("Number of grids along y provided by epsilon_xx and epsilon_zz should be same.")); end
        if nz_Ex != nz_Ey; throw(ArgumentError("Number of grids along z provided by epsilon_xx and epsilon_yy should be same.")); end
        check_BC_and_grid(xBC, nx_Ex, nx_Ey, nx_Ez, "x")
        check_BC_and_grid(yBC, ny_Ex, ny_Ey, ny_Ez, "y")
        check_BC_and_grid(zBC, nz_Ex, nz_Ey, nz_Ez, "z")
        if (~isa(epsilon_xy, Nothing) && ~(size(epsilon_xy) == (nx_Ez, ny_Ez, nz_Ex)))
            throw(ArgumentError("The size of epsilon_xy should be should be (size(epsilon_zz, 1), size(epsilon_zz, 2), size(epsilon_xx, 3)) = ($(size(epsilon_zz, 1)), $(size(epsilon_zz, 2)), $(size(epsilon_xx, 3)))."))
        end
        if (~isa(epsilon_xz, Nothing) && ~(size(epsilon_xz) == (nx_Ey, ny_Ex, nz_Ey)))
            throw(ArgumentError("The size of epsilon_xz should be should be (size(epsilon_yy, 1), size(epsilon_xx, 2), size(epsilon_yy, 3)) = ($(size(epsilon_yy, 1)), $(size(epsilon_xx, 2)), $(size(epsilon_yy, 3)))."))
        end
        if (~isa(epsilon_yx, Nothing) && ~(size(epsilon_yx) == (nx_Ez, ny_Ez, nz_Ey)))
            throw(ArgumentError("The size of epsilon_yx should be should be (size(epsilon_zz, 1), size(epsilon_zz, 2), size(epsilon_yy, 3)) = ($(size(epsilon_zz, 1)), $(size(epsilon_zz, 2)), $(size(epsilon_yy, 3)))."))
        end
        if (~isa(epsilon_yz, Nothing) && ~(size(epsilon_yz) == (nx_Ey, ny_Ex, nz_Ex)))
            throw(ArgumentError("The size of epsilon_yz should be should be (size(epsilon_yy, 1), size(epsilon_xx, 2), size(epsilon_xx, 3)) = ($(size(epsilon_yy, 1)), $(size(epsilon_xx, 2)), $(size(epsilon_xx, 3)))."))
        end
        if (~isa(epsilon_zx, Nothing) && ~(size(epsilon_zx) == (nx_Ey, ny_Ez, nz_Ey)))
            throw(ArgumentError("The size of epsilon_zx should be should be (size(epsilon_yy, 1), size(epsilon_zz, 2), size(epsilon_yy, 3)) = ($(size(epsilon_yy, 1)), $(size(epsilon_zz, 2)), $(size(epsilon_yy, 3)))."))
        end
        if (~isa(epsilon_zy, Nothing) && ~(size(epsilon_zy) == (nx_Ez, ny_Ex, nz_Ex)))
            throw(ArgumentError("The size of epsilon_zy should be should be (size(epsilon_zz, 1), size(epsilon_xx, 2), size(epsilon_xx, 3)) = ($(size(epsilon_zz, 1)), $(size(epsilon_xx, 2)), $(size(epsilon_xx, 3)))."))
        end
        
        # Checking the whether the permittivity is lossless or not 
        if (maximum(imag(epsilon_xx)) >= sqrt(eps()))
            @warn("syst.epsilon_xx contains imaginary part; the permittivity profile is not lossless.")
        end
        if (maximum(imag(epsilon_yy)) >= sqrt(eps()))
            @warn("syst.epsilon_yy contains imaginary part; the permittivity profile is not lossless.")
        end
        if (maximum(imag(epsilon_zz)) >= sqrt(eps()))
            @warn("syst.epsilon_zz contains imaginary part; the permittivity profile is not lossless.")
        end        
        if (~isa(epsilon_xy, Nothing) && ~isa(epsilon_yx, Nothing) && maximum(abs.(epsilon_xy - conj(epsilon_yx))) >= sqrt(eps()))
            @warn("syst.epsilon_xy != conj(syst.epsilon_yx); the permittivity profile is not lossless.")
        end
        if (~isa(epsilon_xz, Nothing) && ~isa(epsilon_zx, Nothing) && maximum(abs.(epsilon_xz - conj(epsilon_zx))) >= sqrt(eps()))
            @warn("syst.epsilon_xz != conj(syst.epsilon_zx); the permittivity profile is not lossless.")
        end
        if (~isa(epsilon_yz, Nothing) && ~isa(epsilon_zy, Nothing) && maximum(abs.(epsilon_yz - conj(epsilon_zy))) >= sqrt(eps()))
            @warn("syst.epsilon_yz != conj(syst.epsilon_zy); the permittivity profile is not lossless.")
        end
        
        # Total number of grid points for Ex, Ey, and Ez    
        nt_Ex = nx_Ex*ny_Ex*nz_Ex
        nt_Ey = nx_Ey*ny_Ey*nz_Ey
        nt_Ez = nx_Ez*ny_Ez*nz_Ez
    end

    # Estimate background permittivity, used to assign default sigma_max_over_omega for PML
    if use_2D_TM # 2D TM field
        epsilon_bg_y = [1, 1]
        epsilon_bg_z = [1, 1]
        if yPML[1].npixels != 0 || yPML[2].npixels != 0
            epsilon_bg_y = real(vcat(mean(epsilon_xx[1,:]), mean(epsilon_xx[end,:])))
            # Make sure that the two sides are the same when the system is periodic; this ensures that the s-factor is continuous across the periodic boundary.
            if isa(yBC, Number) || yBC == "periodic"
                epsilon_bg_y = mean(epsilon_bg_y)*[1, 1]
            end
        end
        if zPML[1].npixels != 0 || zPML[2].npixels != 0
            epsilon_bg_z = real(vcat(mean(epsilon_xx[:,1]), mean(epsilon_xx[:,end])))
            if isa(zBC, Number) || zBC == "periodic"
                epsilon_bg_z = mean(epsilon_bg_z)*[1, 1]
            end
        end
        
        # Set default values for PML parameters
        yPML = set_PML_params(yPML, k0dx, epsilon_bg_y, "y")  
        zPML = set_PML_params(zPML, k0dx, epsilon_bg_z, "z")
    else # 3D case
        epsilon_bg_x_Ex = [1, 1]
        epsilon_bg_y_Ey = [1, 1]
        epsilon_bg_z_Ez = [1, 1]
        if xPML[1].npixels != 0 || xPML[2].npixels != 0
            epsilon_bg_x_Ex = real(vcat(mean(epsilon_xx[1,:,:]), mean(epsilon_xx[end,:,:])))
            # Make sure that the two sides are the same when the system is periodic; this ensures that the s-factor is continuous across the periodic boundary.
            if isa(xBC, Number) || xBC == "periodic"
                epsilon_bg_x_Ex = mean(epsilon_bg_x_Ex)*[1, 1]
            end
        end
        if yPML[1].npixels != 0 || yPML[2].npixels != 0
            epsilon_bg_y_Ey = real(vcat(mean(epsilon_yy[:,1,:]), mean(epsilon_yy[:,end,:])))    
            if isa(yBC, Number) || yBC == "periodic"
                epsilon_bg_y_Ey = mean(epsilon_bg_y_Ey)*[1, 1]
            end    
        end
        if zPML[1].npixels != 0 || zPML[2].npixels != 0
            epsilon_bg_z_Ez = real(vcat(mean(epsilon_zz[:,:,1]), mean(epsilon_zz[:,:,end])))
            if isa(zBC, Number) || zBC == "periodic"
                epsilon_bg_z_Ez = mean(epsilon_bg_z_Ez)*[1, 1]
            end
        end    
    
        # Set default values for PML parameters
        xPML = set_PML_params(xPML, k0dx, epsilon_bg_x_Ex, "x")  
        yPML = set_PML_params(yPML, k0dx, epsilon_bg_y_Ey, "y")
        zPML = set_PML_params(zPML, k0dx, epsilon_bg_z_Ez, "z")  
    end
    
    # Build the first derivative and use Kronecker outer product to go from 1D to 2D or from 1D to 3D
    if use_2D_TM
        # Build the first derivative and average matrices on E, such that
        #   ddx_E*f = df/dx
        #   avg_x_E*f = average of f among two neighboring pixels
        # where f = Ex is a 1d vector.
        # Note that sx_E and sx_H are column vectors.
        # Later we will use Kronecker outer product to go from 1D to 2D.
        (ddy_E, avg_y_E, sy_E, sy_H, ind_yPML_E) = build_ddx_E(ny_Ex, yBC, yPML, "y") # ddy_E operates on Ex or Ez
        (ddz_E, avg_z_E, sz_E, sz_H, ind_zPML_E) = build_ddx_E(nz_Ex, zBC, zPML, "z") # ddz_E operates on Ex or Ey
        
        # The derivative matrices on H
        # We need sparse() to force datatype conversion.
        ddy_H = sparse(-ddy_E') # ddy_H operates on Hx or Hz
        ddz_H = sparse(-ddz_E') # ddz_H operates on Hx or Hy
        
        # Number of grid points in y and z for Hx
        ny_Hx = size(ddy_E, 1) # ny_Hx = ny_Hz = ny_Ey
        nz_Hx = size(ddz_E, 1) # nz_Hx = nz_Hy = nz_Ez
        
        # Include coordinate streching into the first derivatives
        # So far we only implement for 2D TM fields Ex(y,z); will add 2D TE fields later
        ddy_E = spdiagm(ny_Hx, ny_Hx, 1 ./sy_H)*ddy_E
        ddz_E = spdiagm(nz_Hx, nz_Hx, 1 ./sz_H)*ddz_E
        
        if ~use_UPML
            ddy_H = spdiagm(ny_Ex, ny_Ex, 1 ./sy_E)*ddy_H
            ddz_H = spdiagm(nz_Ex, nz_Ex, 1 ./sz_E)*ddz_H
        end
        
        # Build the operators; use Kronecker outer product to go from 1D to 2D.
        # TM: A = [- (d/dy)^2 - (d/dz)^2 - k0^2*epsilon(y,z)]*(dx^2)
        if ~use_UPML
            # Stretched-coordinate PML (SC-PML); here, both 1/s_E and 1/s_H have already been multiplied onto the first derivatives.
            A = - kron(ddz_H*ddz_E, sparse(I,ny_Ex,ny_Ex)) - kron(sparse(I,nz_Ex,nz_Ex), ddy_H*ddy_E) - spdiagm(nt_Ex, nt_Ex, (k0dx^2)*epsilon_xx[:])
        else
            # Uniaxial PML (UPML); here, the first 1/s has already been multiplied onto the first derivatives.
            # A_UPML = (sy_E*sz_E)*A_SCPML for TM fields
            syz_E = sy_E.*reshape(sz_E, 1, :)
            A = - kron(ddz_H*ddz_E, spdiagm(ny_Ex, ny_Ex, sy_E)) - kron(spdiagm(nz_Ex, nz_Ex, sz_E), ddy_H*ddy_E) - spdiagm(nt_Ex, nt_Ex, (k0dx^2)*(syz_E[:].*epsilon_xx[:]))
        end
    else               
        # Build the first derivatives in 1D, such that ddx*f = df/dx, ddy*f = df/dy, ddz*f = df/dz
        # Later we will use Kronecker outer product to go from 1D to 3D. 
        # The operator notations follow
        #   ddx_HzEy: x-derivative operates on Ey producing Hz
        #   avg_x_Ey: average of Ey among two neighboring pixels along x-direction
        #   sx_Ey: coordinate-stretching factor for Ey along x-direction
        #   sx_Hz: coordinate-stretching factor for Hz along x-direction
        # ddx_HyEz = ddx_HzEy, sx_Ez = sx_Ey, sx_Hy = sx_Hz, because Ez(Hy) and Ey(Hz)
        # share the same x-coordinate inside same Yee-cell. The same reason applies to others. 
        # We defer three of the derivatives to later
        (ddx_HzEy, avg_x_Ey, sx_Ey, sx_Hz, ind_xPML_E) = build_ddx_E(nx_Ey, xBC, xPML, "x") 
        #ddx_HyEz = ddx_HzEy
        (ddy_HxEz, avg_y_Ez, sy_Ez, sy_Hx, ind_yPML_E) = build_ddx_E(ny_Ez, yBC, yPML, "y")
        #ddy_HzEx = ddy_HxEz
        (ddz_HyEx, avg_z_Ex, sz_Ex, sz_Hy, ind_zPML_E) = build_ddx_E(nz_Ex, zBC, zPML, "z")
        #ddz_HxEy = ddz_HyEx

        avg_x_Ex = build_ave_x_Ex(nx_Ex, xBC, "x") 
        avg_y_Ey = build_ave_x_Ex(ny_Ey, yBC, "y")
        avg_z_Ez = build_ave_x_Ex(nz_Ez, zBC, "z")
        
        avg_x_Ez = avg_x_Ey
        avg_y_Ex = avg_y_Ez
        avg_z_Ey = avg_z_Ex
        
        # The derivative matrices on H
        # We need sparse() to force datatype conversion.
        ddx_EyHz = sparse(-ddx_HzEy'); #ddx_EzHy = sparse(-ddx_HyEz')
        ddy_EzHx = sparse(-ddy_HxEz'); #ddy_ExHz = sparse(-ddy_HzEx')
        ddz_ExHy = sparse(-ddz_HyEx'); #ddz_EyHx = sparse(-ddz_HxEy')

        # Number of grid points in x, y, and z for Hx, Hy, and Hz
        nx_Hz = size(ddx_HzEy, 1); nx_Hy = nx_Hz; nx_Hx = nx_Ey
        ny_Hx = size(ddy_HxEz, 1); ny_Hz = ny_Hx; ny_Hy = ny_Ez
        nz_Hy = size(ddz_HyEx, 1); nz_Hx = nz_Hy; nz_Hz = nz_Ex

        # Total number of grid points for Hx, Hy, and Hz    
        nt_Hx = nx_Hx*ny_Hx*nz_Hx
        nt_Hy = nx_Hy*ny_Hy*nz_Hy
        nt_Hz = nx_Hz*ny_Hz*nz_Hz    

        # Start from stretched-coordinate PML (SC-PML); here, 1/s and 1/s_d are both multiplied onto the gradient. 
        # Convert coordinate-stretching factor vector to diagonal matrices    
        inv_sx_Hz = spdiagm(nx_Hz, nx_Hz, 1 ./sx_Hz); #inv_sx_Hy = inv_sx_Hz
        inv_sy_Hx = spdiagm(ny_Hx, ny_Hx, 1 ./sy_Hx); #inv_sy_Hz = inv_sy_Hx
        inv_sz_Hy = spdiagm(nz_Hy, nz_Hy, 1 ./sz_Hy); #inv_sz_Hx = inv_sz_Hy
        inv_sx_Ey = spdiagm(nx_Ey, nx_Ey, 1 ./sx_Ey); #inv_sx_Ez = inv_sx_Ey
        inv_sy_Ez = spdiagm(ny_Ez, ny_Ez, 1 ./sy_Ez); #inv_sy_Ex = inv_sy_Ez
        inv_sz_Ex = spdiagm(nz_Ex, nz_Ex, 1 ./sz_Ex); #inv_sz_Ey = inv_sz_Ex

        # Include coordinate stretching into the first derivatives    
        ddx_HzEy = inv_sx_Hz*ddx_HzEy; ddx_HyEz = ddx_HzEy #ddx_HyEz = inv_sx_Hy*ddx_HyEz
        ddy_HxEz = inv_sy_Hx*ddy_HxEz; ddy_HzEx = ddy_HxEz #ddy_HzEx = inv_sy_Hz*ddy_HzEx
        ddz_HyEx = inv_sz_Hy*ddz_HyEx; ddz_HxEy = ddz_HyEx #ddz_HxEy = inv_sz_Hx*ddz_HxEy
        ddx_EyHz = inv_sx_Ey*ddx_EyHz; ddx_EzHy = ddx_EyHz #ddx_EzHy = inv_sx_Ez*ddx_EzHy
        ddy_EzHx = inv_sy_Ez*ddy_EzHx; ddy_ExHz = ddy_EzHx #ddy_ExHz = inv_sy_Ex*ddy_ExHz
        ddz_ExHy = inv_sz_Ex*ddz_ExHy; ddz_EyHx = ddz_ExHy #ddz_EyHx = inv_sz_Ey*ddz_EyHx

        # Expand derivative matrices to 3D by Kronecker outer product
        # Note that even though ddx_HzEy = ddx_HyEz, Dx_HzEy and Dx_HyEz are not the same 
        # because Ez and Ey can have different number of points along y and z
        Dx_HzEy = kron(sparse(I,nz_Ey,nz_Ey), kron(sparse(I,ny_Ey,ny_Ey), ddx_HzEy))
        Dx_HyEz = kron(sparse(I,nz_Ez,nz_Ez), kron(sparse(I,ny_Ez,ny_Ez), ddx_HyEz))
        Dy_HxEz = kron(sparse(I,nz_Ez,nz_Ez), kron(ddy_HxEz, sparse(I,nx_Ez,nx_Ez))) 
        Dy_HzEx = kron(sparse(I,nz_Ex,nz_Ex), kron(ddy_HzEx, sparse(I,nx_Ex,nx_Ex)))
        Dz_HyEx = kron(ddz_HyEx, kron(sparse(I,ny_Ex,ny_Ex), sparse(I,nx_Ex,nx_Ex)))
        Dz_HxEy = kron(ddz_HxEy, kron(sparse(I,ny_Ey,ny_Ey), sparse(I,nx_Ey,nx_Ey)))

        Dx_EyHz = kron(sparse(I,nz_Ey,nz_Ey), kron(sparse(I,ny_Ey,ny_Ey), ddx_EyHz))
        Dx_EzHy = kron(sparse(I,nz_Ez,nz_Ez), kron(sparse(I,ny_Ez,ny_Ez), ddx_EzHy))
        Dy_EzHx = kron(sparse(I,nz_Ez,nz_Ez), kron(ddy_EzHx, sparse(I,nx_Ez,nx_Ez))) 
        Dy_ExHz = kron(sparse(I,nz_Ex,nz_Ex), kron(ddy_ExHz, sparse(I,nx_Ex,nx_Ex)))
        Dz_ExHy = kron(ddz_ExHy, kron(sparse(I,ny_Ex,ny_Ex), sparse(I,nx_Ex,nx_Ex)))
        Dz_EyHx = kron(ddz_EyHx, kron(sparse(I,ny_Ey,ny_Ey), sparse(I,nx_Ey,nx_Ey)))

        # Construct curl_E and curl_H   
        # curl_E: curl operator matrix operates on E-field producing H-field.  
        # curl_H: curl operator matrix operates on H-field producing E-field.      
        curl_E = vcat(hcat(spzeros(nt_Hx, nt_Ex),     -Dz_HxEy,               Dy_HxEz),
                      hcat(     Dz_HyEx,          spzeros(nt_Hy, nt_Ey),     -Dx_HyEz),
                      hcat(    -Dy_HzEx,               Dx_HzEy,          spzeros(nt_Hz, nt_Ez)))   
        curl_H = vcat(hcat(spzeros(nt_Ex, nt_Hx),     -Dz_ExHy,               Dy_ExHz),
                      hcat(     Dz_EyHx,          spzeros(nt_Ey, nt_Hy),     -Dx_EyHz),
                      hcat(    -Dy_EzHx,               Dx_EzHy,          spzeros(nt_Ez, nt_Hz)))    

        # Construct the vectorial Maxwell matrix A
        epsilon_diagonal = spdiagm(nt_Ex+nt_Ey+nt_Ez, nt_Ex+nt_Ey+nt_Ez, vcat(epsilon_xx[:], epsilon_yy[:], epsilon_zz[:]))
        A = curl_H*curl_E-(k0dx)^2*epsilon_diagonal
        
        # Construct the off-diagonal part from the the relative permittivity tensor epsilon_ij, when i does not equal j
        if include_off_diagonal_epsilon
            A_off_diagonal_epsilon = spzeros(ComplexF64, nt_Ex+nt_Ey+nt_Ez, nt_Ex+nt_Ey+nt_Ez)

            # Following Oskooi et al, Optics Letters 34, 2778 (2009), we average two points of Ey along y, multiply by epsilon_xy to get Ex, and then average two such points along direction x. To summarize: avg_x*epsilon_xy*avg_y*Ey. Similarly for the other terms.
            if ~isa(epsilon_xy, Nothing)
                matrix_epsilon_xy = spdiagm(nx_Ez*ny_Ez*nz_Ex, nx_Ez*ny_Ez*nz_Ex, epsilon_xy[:])
                A_off_diagonal_epsilon[1:nt_Ex, (nt_Ex+1):(nt_Ex+nt_Ey)] = A_off_diagonal_epsilon[1:nt_Ex, (nt_Ex+1):(nt_Ex+nt_Ey)] - (k0dx)^2 * kron(sparse(I,nz_Ex,nz_Ex), kron(sparse(I,ny_Ex,ny_Ex), avg_x_Ex'))*matrix_epsilon_xy*kron(sparse(I,nz_Ey,nz_Ey), kron(avg_y_Ey, sparse(I,nx_Ey,nx_Ey))) # avg_x*epsilon_xy*avg_y*Ey
            end
            if ~isa(epsilon_xz, Nothing)
                matrix_epsilon_xz = spdiagm(nx_Ey*ny_Ex*nz_Ey, nx_Ey*ny_Ex*nz_Ey, epsilon_xz[:])
                A_off_diagonal_epsilon[1:nt_Ex, (nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez)] = A_off_diagonal_epsilon[1:nt_Ex, (nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez)] - (k0dx)^2 * kron(sparse(I,nz_Ex,nz_Ex), kron(sparse(I,ny_Ex,ny_Ex), avg_x_Ex'))*matrix_epsilon_xz*kron(avg_z_Ez, kron(sparse(I,ny_Ez,ny_Ez), sparse(I,nx_Ez,nx_Ez))) # avg_x*epsilon_xz*avg_z*Ez
            end
            if ~isa(epsilon_yx, Nothing)
                matrix_epsilon_yx = spdiagm(nx_Ez*ny_Ez*nz_Ey, nx_Ez*ny_Ez*nz_Ey, epsilon_yx[:])
                A_off_diagonal_epsilon[(nt_Ex+1):(nt_Ex+nt_Ey), 1:nt_Ex] = A_off_diagonal_epsilon[(nt_Ex+1):(nt_Ex+nt_Ey), 1:nt_Ex] - (k0dx)^2 * kron(sparse(I,nz_Ey,nz_Ey), kron(avg_y_Ey', sparse(I,nx_Ey,nx_Ey)))*matrix_epsilon_yx*kron(sparse(I,nz_Ex,nz_Ex), kron(sparse(I,ny_Ex,ny_Ex), avg_x_Ex)) # avg_y*epsilon_yx*avg_x*Ex
            end
            if ~isa(epsilon_yz, Nothing)
                matrix_epsilon_yz = spdiagm(nx_Ey*ny_Ex*nz_Ex, nx_Ey*ny_Ex*nz_Ex, epsilon_yz[:])
                A_off_diagonal_epsilon[(nt_Ex+1):(nt_Ex+nt_Ey), (nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez)] = A_off_diagonal_epsilon[(nt_Ex+1):(nt_Ex+nt_Ey), (nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez)] - (k0dx)^2 * kron(sparse(I,nz_Ey,nz_Ey), kron(avg_y_Ey', sparse(I,nx_Ey,nx_Ey)))*matrix_epsilon_yz*kron(avg_z_Ez, kron(sparse(I,ny_Ez,ny_Ez), sparse(I,nx_Ez,nx_Ez))) # avg_y*epsilon_yz*avg_z*Ez
            end
            if ~isa(epsilon_zx, Nothing)
                matrix_epsilon_zx = spdiagm(nx_Ey*ny_Ez*nz_Ey, nx_Ey*ny_Ez*nz_Ey, epsilon_zx[:])
                A_off_diagonal_epsilon[(nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez), 1:nt_Ex] = A_off_diagonal_epsilon[(nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez), 1:nt_Ex] - (k0dx)^2 * kron(avg_z_Ez', kron(sparse(I,ny_Ez,ny_Ez), sparse(I,nx_Ez,nx_Ez)))*matrix_epsilon_zx*kron(sparse(I,nz_Ex,nz_Ex), kron(sparse(I,ny_Ex,ny_Ex), avg_x_Ex)) # avg_z*epsilon_zx*avg_x*Ex
            end
            if ~isa(epsilon_zy, Nothing)
                matrix_epsilon_zy = spdiagm(nx_Ez*ny_Ex*nz_Ex, nx_Ez*ny_Ex*nz_Ex, epsilon_zy[:])
                A_off_diagonal_epsilon[(nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez), (nt_Ex+1):(nt_Ex+nt_Ey)] = A_off_diagonal_epsilon[(nt_Ex+nt_Ey+1):(nt_Ex+nt_Ey+nt_Ez), (nt_Ex+1):(nt_Ex+nt_Ey)] - (k0dx)^2 * kron(avg_z_Ez', kron(sparse(I,ny_Ez,ny_Ez), sparse(I,nx_Ez,nx_Ez)))*matrix_epsilon_zy*kron(sparse(I,nz_Ey,nz_Ey), kron(avg_y_Ey, sparse(I,nx_Ey,nx_Ey))) # avg_z*epsilon_zy*avg_y*Ey
            end
            A = A + A_off_diagonal_epsilon
        end
        
        if use_UPML
            # sx_Ex = sx_Hz, because Ex and Hz share the same x-coordinate (and sx) inside same Yee-cell.
            # The same reason applies to others.        
            sx_Ex = sx_Hz
            sy_Ey = sy_Hx
            sz_Ez = sz_Hy

            sy_Ex = sy_Ez
            sz_Ey = sz_Ex
            sx_Ez = sx_Ey

            # Construct a 3D S_E to transform between A_SCPML to A_UPML
            (sx_Ex_3d,sy_Ex_3d,sz_Ex_3d) = ndgrid(sx_Ex,sy_Ex,sz_Ex)
            (sx_Ey_3d,sy_Ey_3d,sz_Ey_3d) = ndgrid(sx_Ey,sy_Ey,sz_Ey)
            (sx_Ez_3d,sy_Ez_3d,sz_Ez_3d) = ndgrid(sx_Ez,sy_Ez,sz_Ez)
            S_E = spdiagm(nt_Ex+nt_Ey+nt_Ez, nt_Ex+nt_Ey+nt_Ez, vcat((sy_Ex_3d.*sz_Ex_3d./sx_Ex_3d)[:], (sx_Ey_3d.*sz_Ey_3d./sy_Ey_3d)[:], (sx_Ez_3d.*sy_Ez_3d./sz_Ez_3d)[:]))
            A = S_E*A            
        end
    end
    
    # determine the symmetry of matrix A, assuming no spatial symmetry in epsilon or inv_epsilon
    # 2D TE fields will be added later
    is_symmetric_A = true
    
    if (~use_2D_TM && isa(xBC, Number) && xBC != 0 && xBC != pi && (nx_Ex > 1 || ny_Ex > 1 || nz_Ex > 1)) || 
       (isa(yBC, Number) && yBC != 0 && yBC != pi && ((use_2D_TM && ny_Ex > 1) || (~use_2D_TM && (nx_Ey > 1 || ny_Ey > 1 || nz_Ey > 1)))) ||
       (isa(zBC, Number) && zBC != 0 && zBC != pi && ((use_2D_TM && nz_Ex > 1) || (~use_2D_TM && (nx_Ez > 1 || ny_Ez > 1 || nz_Ez > 1))))   
       # Bloch periodic boundary condition with k_B*periodicity != 0 or pi breaks the symmetry of A because its ddx is complex-valued, except when there is only one pixel.
        is_symmetric_A = false
    elseif (~use_2D_TM && (xPML[1].npixels != 0 || xPML[2].npixels != 0 || yPML[1].npixels != 0 || yPML[2].npixels != 0 || zPML[1].npixels != 0 || zPML[2].npixels != 0)) || (use_2D_TM && ~use_UPML && (yPML[1].npixels != 0 || yPML[2].npixels != 0 || zPML[1].npixels != 0 || zPML[2].npixels != 0))
        is_symmetric_A = false
    end
    
    if ~use_2D_TM
        return (A, is_symmetric_A, xPML, yPML, zPML)
    else
        return (A, is_symmetric_A, yPML, zPML)
    end
end

"""
    MESTI_BUILD_FDFD_MATRIX The finite-difference frequency-domain operator for 2D TM waves.
       A = mesti_build_fdfd_matrix(epsilon_xx, k0dx, yBC, zBC, yPML, zPML, use_UPML)
           returns A as a sparse matrix representing wave operator [- (d/dy)^2 - (d/dz)^2 - k0^2*epsilon(y,z)]*(dx^2)
           for the Ex(y,z) component of transverse-magnetic (TM) fields (Ex, Hy, Hz),
           discretized on a square grid with grid size dx through center difference.
           Matrix A has size [ny_Ex*nz_Ex, ny_Ex*nz_Ex].   
        Here are the differences between 3D and 2D TM versions
           epsilon_xx (matrix; required):
              epsilon_xx is a ny_Ex-by-nz_Ex matrix (real or complex)
              discretizing the relative permittivity profile on Ex-sites. ny_Ex/nz_Ex 
              is the total number of grids for Ex in y/z-direction.
           A (sparse matrix):
              nt_Ex-by-nt_Ex sparse matrix representing the 2D FDFD operator. 
              nt_Ex = ny_Ex*nz_Ex is the total number of grids for Ex.

       === Input Arguments ===
       epsilon_xx (matrix; required):
          epsilon_xx is a ny_Ex-by-nz_Ex matrix (real or complex)
          discretizing the relative permittivity profile epsilon_xx on Ex-sites.
          The notation ni_Ej is the total number of grids for Ej in i-direction.
       k0dx (numeric scalar, real or complex; required):
          Normalized frequency k0*dx = (2*pi/vacuum_wavelength)*dx.
       yBC (character vector or numeric scalar; required):
             Boundary condition (BC) at the two ends in y direction, effectively
                        Ex(n,m,l) at m=0 and m=ny_Ex+1.    
             one pixel beyond the computation domain. Available choices are:
               "periodic" - Ex(m+ny_Ex,l) = Ex(m,l),
               "PEC"      - Ex(0,l) = Ex(ny_Ex,l) = 0,
               "PMC"      - Ex(0,l) = Ex(1,l); Ex(ny_Ex+1,l) = Ex(ny_Ex,l),
               "PECPMC"   - Ex(0,l) = 0; Ex(ny_Ex+1,m,l) = Ex(ny_Ex,m,l),
               "PMCPEC"   - Ex(0,l) = Ex(1,l); Ex(ny_Ex+1,m,l) = 0,
             where PEC stands for perfect electric conductor and PMC stands for perfect
             magnetic conductor.
             When yBC is a numeric scalar, the Bloch periodic boundary condition is
          used with Ex(m+ny_Ex,l) = Ex(m,l)*exp(1i*yBC); in other words, 
          yBC = ky_B*ny_Ex*dx = ky_B*Lambda where ky_B is the Bloch wave number and 
          Lambda = ny_Ex*dx is the periodicity in y.
       zBC (character vector or numeric scalar; required):
          Boundary condition in z direction, analogous to yBC.    
       yPML (a two-element PML vector; optional):
          Parameters for perfectly matched layer (PML) in y direction.
          xPML = [PML_left, PML_right] 
          If users do not want to put PML on either side, just set the field "npixels" = 0.
          For example, PML_left.npixels = 0: no PML on the left.
          In each case, PML is a PML structure with the following fields:
             npixels (non-negative integer scalar; required): Number of PML pixels.
                Note this is within syst.epsilon or syst.inv_epsilon.
             power_sigma (non-negative scalar; optional): Power of the polynomial 
                grading for the conductivity sigma; defaults to 3.
             sigma_max_over_omega (non-negative scalar; optional):
                Conductivity at the end of the PML; defaults to
                   0.8*(power_sigma+1)/(k0dx*sqrt(epsilon_bg)).
                where epsilon_bg is the average relative permittivity along the
                last slice of the PML. This is used to attenuate propagating waves.
             power_kappa (non-negative scalar; optional): Power of the polynomial
                grading for the real-coordinate-stretching factor kappa; defaults
                to 3.
             kappa_max (real scalar no smaller than 1; optional):
                Real-coordinate-stretching factor at the end of the PML; defaults
                to 15. This is used to accelerate the attenuation of evanescent 
                waves. kappa_max = 1 means no real-coordinate stretching.
             power_alpha (non-negative scalar; optional): Power of the polynomial
                grading for the CFS alpha factor; defaults to 1.
             alpha_max_over_omega (non-negative scalar; optional): Complex-
                frequency-shifting (CFS) factor at the beginning of the PML. This
                is typically used in time-domain simulations to suppress late-time
                (low-frequency) reflections. We do not use it by default 
                (alpha_max_over_omega = 0) since we are in frequency domain.
          We use the following PML coordinate-stretching factor:
             s(p) = kappa(p) + sigma(p)./(alpha(p) - i*omega)
          with
             sigma(p)/omega = sigma_max_over_omega*(p.^power_sigma),
             kappa(p) = 1 + (kappa_max-1)*(p.^power_kappa),
             alpha(p)/omega = alpha_max_over_omega*((1-p).^power_alpha),
          where omega is frequency, and p goes linearly from 0 at the beginning of
          the PML to 1 at the end of the PML. 
       zPML (a two-element vector; optional):
          Parameters for PML in z direction, analogous to yPML.
       use_UPML (logical scalar; optional, defaults to true):
          Whether to use uniaxial PML (UPML) or not. If not, stretched-coordinate
          PML (SC-PML) will be used.

       === Output Arguments ===
       A (sparse matrix):
          nt_Ex-by-nt_Ex sparse matrix representing the 2D FDFD operator. 
          nt_Ex = ny_Ex*nz_Ex is the total number of grids for Ex.
       is_symmetric_A (logical scalar):
          Whether matrix A is symmetric or not.
       yPML (two-element vector):
          PML parameters used on the low and high sides of y direction, if any.
       zPML (two-element vector):
          PML parameters used on the low and high sides of z direction, if any.
"""
function mesti_build_fdfd_matrix(epsilon_xx::Union{Matrix{Int64},Matrix{Float64},Matrix{ComplexF64}}, k0dx::Union{Real,Complex}, yBC::Union{String,Real,Complex}, zBC::Union{String,Real,Complex}, yPML::Vector{PML} = [PML(0), PML(0)], zPML::Vector{PML} = [PML(0), PML(0)], use_UPML::Bool=true)
    return mesti_build_fdfd_matrix(epsilon_xx, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, k0dx, nothing, yBC, zBC, nothing, yPML, zPML, use_UPML)
end


# When only the diagonal terms of epsilon (i.e. epsilon_xx, epsilon_yy, epsilon_zz) are specified in the 3D case.
function mesti_build_fdfd_matrix(epsilon_xx::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Matrix{Int64},Matrix{Float64},Matrix{ComplexF64}}, epsilon_yy::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, epsilon_zz::Union{Array{Int64,3},Array{Float64,3},Array{ComplexF64,3},Nothing}, k0dx::Union{Real,Complex}, xBC::Union{String,Real,Complex,Nothing}, yBC::Union{String,Real,Complex}, zBC::Union{String,Real,Complex}, xPML::Union{Vector{PML},Nothing} = [PML(0), PML(0)], yPML::Vector{PML} = [PML(0), PML(0)], zPML::Vector{PML} = [PML(0), PML(0)], use_UPML::Bool=true)
    return mesti_build_fdfd_matrix(epsilon_xx, nothing, nothing, nothing, epsilon_yy, nothing, nothing, nothing, epsilon_zz, k0dx, xBC, yBC, zBC, xPML, yPML, zPML, use_UPML)
end

"""
    CHECK_BC_AND_GRID(BC, n_Ex, n_Ey, n_Ez, direction) is a helper function and checks the consistency of boundary condition and number of grids in n_Ex, n_Ey, and n_Ez. 
"""
function check_BC_and_grid(BC::Union{String,Real,Complex}, n_Ex::Int, n_Ey::Int, n_Ez::Int, direction::String)
    if isa(BC, Number)
        if (n_Ex != n_Ez || n_Ex != n_Ey)
            throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx, epsilon_yy, and epsilon_zz should be the same for Bloch periodic boundary condition in $(direction)-direction."))
        end      
    elseif (n_Ex != n_Ez || n_Ex != n_Ey) && (BC == "periodic" ||  BC == "PECPMC" || BC == "PMCPEC")
        throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx, epsilon_yy, and epsilon_zz should be the same for $(BC) boundary condition in $(direction)-direction."))     
    else
        if direction == "x"
            if (n_Ex != n_Ez+1 || n_Ex != n_Ey+1) && (BC == "PEC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx should be 1 grid more than epsilon_yy (epsilon_zz) for $(BC) boundary condition in $(direction)-direction.")); end
            if (n_Ex != n_Ez-1 || n_Ex != n_Ey-1) && (BC == "PMC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx should be 1 grid less than epsilon_yy (epsilon_zz) for $(BC) boundary condition in $(direction)-direction.")); end
        elseif direction == "y"
            if (n_Ey != n_Ex+1 || n_Ey != n_Ez+1) && (BC == "PEC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_yy should be 1 grid more than epsilon_xx (epsilon_zz) for $(BC) boundary condition in $(direction)-direction.")); end
            if (n_Ey != n_Ex-1 || n_Ey != n_Ez-1) && (BC == "PMC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_yy should be 1 grid less than epsilon_xx (epsilon_zz) for $(BC) boundary condition in $(direction)-direction.")); end            
        elseif direction == "z"
            if (n_Ez != n_Ey+1 || n_Ez != n_Ex+1) && (BC == "PEC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_zz should be 1 grid more than epsilon_xx (epsilon_yy) for $(BC) boundary condition in $(direction)-direction.")); end
            if (n_Ez != n_Ey-1 || n_Ez != n_Ex-1) && (BC == "PMC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_zz should be 1 grid less than epsilon_xx (epsilon_yy) for $(BC) boundary condition in $(direction)-direction.")); end              
        end
    end
end

"""
    CHECK_BC_AND_GRID(BC, n_Ex, n_Ey, direction) is a helper function and checks the consistency of boundary condition and number of grids in n_Ex, and n_Ey. 
"""
function check_BC_and_grid(BC::Union{String,Real,Complex}, n_Ex::Int, n_Ey::Int, direction::String)
    if isa(BC, Number)
        if (n_Ex != n_Ey)
            throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx, epsilon_yy, should be the same for Bloch periodic boundary condition in $(direction)-direction."))
        end      
    elseif (n_Ex != n_Ey) && (BC == "periodic" ||  BC == "PECPMC" || BC == "PMCPEC")
        throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx, epsilon_yy, should be the same for $(BC) boundary condition in $(direction)-direction."))     
    else
        if direction == "x"
            if (n_Ex != n_Ey+1) && (BC == "PEC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx should be 1 grid more than epsilon_yy for $(BC) boundary condition in $(direction)-direction.")); end
            if (n_Ex != n_Ey-1) && (BC == "PMC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_xx should be 1 grid less than epsilon_yy for $(BC) boundary condition in $(direction)-direction.")); end
        elseif direction == "y"
            if (n_Ey != n_Ex+1) && (BC == "PEC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_yy should be 1 grid more than epsilon_xx for $(BC) boundary condition in $(direction)-direction.")); end
            if (n_Ey != n_Ex-1) && (BC == "PMC"); throw(ArgumentError("Number of grids along $(direction) provided by epsilon_yy should be 1 grid less than epsilon_xx for $(BC) boundary condition in $(direction)-direction.")); end            
        end
    end
end

"""
    CONVERT_BC is a helper function to handles the case of string in BC
"""
function convert_BC(BC::Union{String,Real,Complex},direction::String)
    if isa(BC, Number)
        return BC
    elseif lowercase(BC) == "pec"
        return "PEC"
    elseif lowercase(BC) == "pmc"
        return "PMC"
    elseif lowercase(BC) == "pecpmc"
        return "PECPMC"
    elseif lowercase(BC) == "pmcpec"
        return "PMCPEC"
    elseif lowercase(BC) == "periodic"
        return "periodic"
    elseif lowercase(BC) == "bloch"
        throw(ArgumentError("To use Bloch periodic boundary condition in $(direction)-direction, set the input argument $(direction)BC to k$(direction)_B*p_$(direction) where k$(direction)_B is the Bloch wave number and Lambda_$(direction) is the periodicity along $(direction)-direction."))        
    else
        ethrow(ArgumentError("Input argument $(direction)BC = \"$(BC)\" is not a supported option."))         
    end
end

"""
    BUIlD_DDX_E builds the first-derivative matrix
        The following documentation all take direction = x for example and can be analogous to other direction.
        ddx: x-direction first-derivative matrix (acting on Ey or Ez)
        avg: average matrix (acting on Ey, or Ez)
        s_E: x-coordinate-stretching factor for Ey, or Ez (on integer sites)
        s_H: x-coordinate-stretching factor for Hy, or Hz (on half-integer sites)
        n_E: number of sites in x for Ey or Ez
        ind_PML_E: a vector containing indices of the PML pixels for Ey or Ez (on integer sites)
"""
function build_ddx_E(n_E::Int, BC::Union{String,Real,Complex}, pml::Vector{PML}, direction::String)
    
    # Handle periodic and Bloch periodic boundary (adding the feature later) conditions
    if isa(BC, Number)
        kLambda = BC
        BC = "Bloch"
        if ~isa(kLambda, Real)
            @warn "k$(direction)_B*periodicity = $(real(kLambda)) + 1im*$(imag(kLambda)) is a complex number."
        end
    elseif BC == "periodic"
        kLambda = 0
        BC = "Bloch"
    end

    # Build the first-derivative matrix and the average matrix on E
    # f = [f(1), ..., f(n_E)].' with f = Ey or Ez, on integer sites
    # df = (df/dx)*dx, proportional to Hy or Hz, on half-integer sites
    # avg_f = average of f between two neighboring sites, on half-integer sites
    if BC == "Bloch"
        # f(n_E+1) = f(1)*exp(1i*kLambda); f(0) = f(n_E)*exp(-1i*kLambda)
        # ddx*f = [df(1.5), ..., df(n_E+0.5)].'
        ddx = spdiagm(n_E, n_E, +1 => ones(n_E-1), 0 => -ones(n_E), 1-n_E => exp(1im*kLambda)*ones(1))
        # avg*f = [avg_f(1.5), ..., avg_f(n_E+0.5)].'
        avg = spdiagm(n_E, n_E, +1 => ones(n_E-1)/2, 0 => ones(n_E)/2, 1-n_E => exp(1im*kLambda)*ones(1)/2)
    elseif BC == "PEC" # PEC on both sides
        # f(0) = f(n_E+1) = 0
        # ddx*f = df = [df(0.5), ..., df(n+0.5)].'
        ddx = spdiagm(n_E+1, n_E, 0 => ones(n_E), -1 => -ones(n_E))
        # avg*f = [avg_f(0.5), ..., avg_f(n_E+0.5)].'
        avg = spdiagm(n_E+1, n_E, 0 => ones(n_E)/2, -1 => ones(n_E)/2)        
    elseif BC == "PMC" # PMC on both sides
        # f(0) = f(1); f(n_E+1) = f(n_E)
        # ddx*f = [df(1.5), ..., df(n_E-0.5)].'; exclude df(0.5) and df(n_E+0.5) because they are zero
        ddx = spdiagm(n_E-1, n_E, +1 => ones(n_E-1), 0 => -ones(n_E-1))
        # avg*f = [avg_f(1.5), ..., avg_f(n_E-0.5)].'; exclude avg_f(0.5) and avg_f(n_E+0.5)
        avg = spdiagm(n_E-1, n_E, +1 => ones(n_E-1)/2, 0 => ones(n_E-1)/2)
    elseif BC == "PECPMC" # PEC on the low side, PMC on the high side
        # f(0) = 0; f(n_E+1) = f(n_E)
        # ddx*f = [df(0.5), ..., df(n_E-0.5)].'; exclude df(n_E+0.5) because it is zero
        ddx = spdiagm(n_E, n_E, 0 => ones(n_E), -1 => -ones(n_E-1))
        # avg*f = [avg_f(0.5), ..., avg_f(n_E-0.5)].'; we exclude avg_f(n_E+0.5)
        avg = spdiagm(n_E, n_E, 0 => ones(n_E)/2, -1 => ones(n_E-1)/2)
    elseif BC == "PMCPEC" # PMC on the low side, PEC on the high side
        # f(0) = f(1); f(n_E+1) = 0
        # ddx*f = [df(1.5), ..., df(n_E+0.5)].'; exclude df(0.5) because it is zero
        ddx = spdiagm(n_E, n_E, +1 => ones(n_E-1), 0 => -ones(n_E))
        # avg*f = [avg_f(1.5), ..., avg_f(n_E+0.5)].'; exclude avg_f(0.5)
        avg = spdiagm(n_E, n_E, +1 => ones(n_E-1)/2, 0 => ones(n_E)/2)
    else
        throw(ArgumentError("Input argument $(direction)BC = \"$(BC)\" is not a supported option."))              
    end

    n_H = size(ddx,1) # number of sites for df/dx (ie, Hy or Hz)

    # Coordinate-stretching factor: s(x) = kappa(x) + sigma(x)/(alpha(u) - i*omega)
    s_E = ones(ComplexF64, n_E) # s-factor for Ey or Ez (on integer sites)
    s_H = ones(ComplexF64, n_H) # s-factor for Hy or Hz (on half-integer sites)

    # no coordinate stretching if no PML is specified   
    if (length(pml) == 2 && pml[1].npixels == 0 && pml[2].npixels == 0)
        ind_PML_E = [nothing, nothing]
        return ddx, avg, s_E, s_H, ind_PML_E
    end

    # Number of PML pixels on the low and high sides
    npixels = [pml[1].npixels, pml[2].npixels]

    n = n_E  # Also works for 2D TM fields
             # Will add 2D TE case later
    # Cannot have more PML pixels than the number of pixels
    # In our implementation, the conductivity goes to zero at one pixel before PML. So there needs to be at least one more pixel in addition to PML.
    if sum(npixels) >= n
        throw(ArgumentError("Total number of pixels = $(n) in $(direction) direction must be greater than the number of PML pixels = $(npixels[1]) + $(npixels[2]) = $(sum(npixels)) but is not."))
    end

    # Below, p(x) is a function that goes linearly from p(x)=0 at one site before PML to p(x)=1 at the "end of PML".
    # Note that where the "end of PML" is depends on the boundary condition.
    # Let index i=0 be one site before PML, so p(i=0)=0.
    # Then, i=1 is the first site of PML and i=npixels is the last site of PML we explicitly simulate. But we do not set p(i=npixels)=1.
    # For Dirichlet BC, the BC is such that f=0 at i=npixels+1, so we let the end of PML be i=npixels+1, with an effective PML thickness of npixels+1 pixels.
    # For Neumann BC, the BC is such that df=0 at i=npixels+0.5, so we let the end of PML be i=npixels+0.5, with an effective PML thickness of npixels+0.5 pixels.
    # For periodic and Bloch periodic BC, we let p(x) be symmetric on the two sides of f with p(0.5)=p(n+0.5)=1; note that f(0.5) and f(n+0.5) are the same site (with a possible Bloch phase difference). If the PML parameters on the two sides are the same, the s-factor will be continuous across the periodic boundary.

    # Construct p(x) and their corresponding indices for df/dx
    # Also works for 2D TM fields; will add 2D TE case later
    if BC == "Bloch"
        if npixels[1]*npixels[2] == 0 && (npixels[1] != 0 || npixels[2] != 0)
            @warn "Bloch periodic boundary condition is used with a single-sided PML in $(direction) direction; transmission through PML will only undergo single-pass attenuation."
        end
        npixels_effective = [npixels[1]+0.5, npixels[2]+0.5]        
        # f = Ez; ddx*f = [df(1.5), ..., df(n_E+0.5)].'
        # no s-factor for df(0.5) since we only consider df(n_E+0.5)
        p_PML_2 = [((1:npixels[1]).-0.5)/npixels_effective[1], ((1:(npixels[2]+1)).-0.5)/npixels_effective[2]]
        ind_PML_2 = [reverse(1:npixels[1]), (n+1).-reverse(1:(npixels[2]+1))]                
    elseif BC == "PEC" # PEC on both sides    
        npixels_effective = [npixels[1]+1, npixels[2]+1]        
        p_PML_2 = [((1:(npixels[1]+1)).-0.5)/npixels_effective[1], ((1:(npixels[2]+1)).-0.5)/npixels_effective[2]]
        ind_PML_2 = [reverse(1:(npixels[1]+1)), (n+2).-reverse(1:(npixels[2]+1))]                
    elseif BC == "PMC" # PMC on both sides
        npixels_effective = [npixels[1]+0.5, npixels[2]+0.5]
        p_PML_2 = [((1:npixels[1]).-0.5)/npixels_effective[1], ((1:(npixels[2])).-0.5)/npixels_effective[2]]
        ind_PML_2 = [reverse(1:npixels[1]), n.-reverse(1:(npixels[2]))]        
    elseif BC == "PECPMC" # PEC on the low side, PMC on the high side
        npixels_effective = [npixels[1]+1, npixels[2]+0.5]
        p_PML_2 = [((1:npixels[1]+1).-0.5)/npixels_effective[1], ((1:(npixels[2])).-0.5)/npixels_effective[2]]
        ind_PML_2 = [reverse(1:npixels[1]+1), (n+1).-reverse(1:(npixels[2]))]        
    elseif BC == "PMCPEC" # PMC on the low side, PEC on the high side
        npixels_effective = [npixels[1]+0.5, npixels[2]+1]
        p_PML_2 = [((1:npixels[1]).-0.5)/npixels_effective[1], ((1:(npixels[2]+1)).-0.5)/npixels_effective[2]]
        ind_PML_2 = [reverse(1:npixels[1]), (n+1).-reverse(1:(npixels[2]+1))] 
    else
        throw(ArgumentError("Input argument $(direction)BC = \"$(BC)\" is not a supported option."))                        
    end        

    # Construct p(x) and their corresponding indices for f and d^2f/dx^2    
    p_PML_1 = [(1:npixels[1])/npixels_effective[1], (1:npixels[2])/npixels_effective[2]]
    ind_PML_1 = [reverse(1:npixels[1]), (n+1).-reverse(1:npixels[2])]   
    
    # Also works for 2D TM modes; will add 2D TE case later
    p_PML_E = p_PML_1
    ind_PML_E = ind_PML_1
    p_PML_H = p_PML_2
    ind_PML_H = ind_PML_2    
    
    # Loop over PML on the two sides
    for ii = 1:2
        if npixels[ii] > 0
            # coordinate-stretching factor s(p) from Eq 7.73 of Taflove & Hagness's 2005 FDTD book
            # sigma is the conductivity, equivalent to imag-coordinate stretching, used to attenuate propagating waves.
            # kappa is real-coordinate stretching, used to accelerate the attenuation of evanescent waves.
            # alpha is used for complex frequency shifting (CFS) to suppress reflection of low-frequency components for time-domain simulations.
            # In general, kappa, sigma, alpha can all be arbitrary functions of position.
            # To minimize discretization-induced reflection, kappa should start from 1, and sigma should start from 0.
            # We use a polynomial grading for all of them: Eqs 7.60 and 7.79 of Taflove & Hagness's 2005 FDTD book.
            kappa(p) = 1 .+ (pml[ii].kappa_max-1)*(p.^(pml[ii].power_kappa))
            sigma_over_omega(p) = pml[ii].sigma_max_over_omega*(p.^(pml[ii].power_sigma))
            alpha_over_omega(p) = pml[ii].alpha_max_over_omega*((1 .- p).^(pml[ii].power_alpha))
            func_s(p) = kappa(p) .+ sigma_over_omega(p)./(alpha_over_omega(p) .- 1im)

            # Evaluate s(p) on integer and half-integer sites
            # Note s_E and s_H should be sampled directly from the polynomial functions; if we use the polynomial for one and linear interpolation for the other, the performance will be much worse.
            s_E[Int.(ind_PML_E[ii])] = func_s(p_PML_E[ii]) # column vector
            s_H[Int.(ind_PML_H[ii])] = func_s(p_PML_H[ii]) # column vector   
        end
    end        

    return ddx, avg, s_E, s_H, ind_PML_E
end


"""
    BUIlD_AVE_X_EX builds the average matrix that acts on Ex along the x-direction and can be analogous to other direction.
        avg: average matrix (acting on Ex along x-direction)
"""
function build_ave_x_Ex(n_E::Int, BC::Union{String,Real,Complex}, direction::String)
    
    # Handle periodic and Bloch periodic boundary (adding the feature later) conditions
    if isa(BC, Number)
        kLambda = BC
        BC = "Bloch"
        if ~isa(kLambda, Real)
            @warn "k$(direction)_B*periodicity = $(real(kLambda)) + 1im*$(imag(kLambda)) is a complex number."
        end
    elseif BC == "periodic"
        kLambda = 0
        BC = "Bloch"
    end

    # Build the first-derivative matrix and the average matrix on E
    # f = [f(1), ..., f(n_E)].' with f = Ey or Ez, on integer sites
    # df = (df/dx)*dx, proportional to Hy or Hz, on half-integer sites
    # avg_f = average of f between two neighboring sites, on half-integer sites
    if BC == "Bloch"
        # f(n_E+1) = f(1)*exp(1i*kLambda); f(0) = f(n_E)*exp(-1i*kLambda)
        # avg*f = [avg_f(1.5), ..., avg_f(n_E+0.5)].'
        avg = spdiagm(n_E, n_E, +1 => ones(n_E-1)/2, 0 => ones(n_E)/2, 1-n_E => exp(1im*kLambda)*ones(1)/2)
    elseif BC == "PEC" # PEC on both sides
        # f(0) = f(1); f(n_E+1) = f(n_E)
        # avg*f = [avg_f(1.5), ..., avg_f(n_E-0.5)].'; exclude avg_f(0.5) and avg_f(n_E+0.5)
        avg = spdiagm(n_E-1, n_E, +1 => ones(n_E-1)/2, 0 => ones(n_E-1)/2)
    elseif BC == "PMC" # PMC on both sides
        # f(0) = f(n_E+1) = 0
        # avg*f = [avg_f(0.5), ..., avg_f(n_E+0.5)].'
        avg = spdiagm(n_E+1, n_E, 0 => ones(n_E)/2, -1 => ones(n_E)/2)        
    elseif BC == "PECPMC" # PEC on the low side, PMC on the high side
        # f(0) = f(1); f(n_E+1) = 0
        # avg*f = [avg_f(1.5), ..., avg_f(n_E+0.5)].'; exclude avg_f(0.5)
        avg = spdiagm(n_E, n_E, +1 => ones(n_E-1)/2, 0 => ones(n_E)/2)
    elseif BC == "PMCPEC" # PMC on the low side, PEC on the high side
        # f(0) = 0; f(n_E+1) = f(n_E)
        # avg*f = [avg_f(0.5), ..., avg_f(n_E-0.5)].'; we exclude avg_f(n_E+0.5)
        avg = spdiagm(n_E, n_E, 0 => ones(n_E)/2, -1 => ones(n_E-1)/2)
    else
        throw(ArgumentError("Input argument $(direction)BC = \"$(BC)\" is not a supported option."))              
    end

    return avg
end



"""
    SET_PML_PARAMS sets default values for PML parameters
"""
function set_PML_params(pml::Vector{PML}, k0dx::Union{Real,Complex}, epsilon_bg::Union{Vector{Int64},Vector{Float64}}, direction::String)
    
    if ~(length(pml) == 2)
        throw(ArgumentError("$(direction)PML must be a vector of PML containing two elements."))
    end

    # No PML on the both sides
    if (pml[1].npixels == 0 && pml[2].npixels == 0)
        return pml
    end

    # loop over PML parameters on two sides
    for ii = 1:2
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
            # From Table 7.1 of Taflove & Hagness's 2005 FDTD book
            pml[ii].power_sigma = 3
        else
            temp = pml[ii].power_sigma
            if temp < 0
                throw(ArgumentError("$(direction)pml[$(ii)].power_sigma must be a non-negative scalar."))
            end
        end

        # Conductivity at the end of the PML
        if ~isdefined(pml[ii], :sigma_max_over_omega)
            # Eq 7.67 of Taflove & Hagness's 2005 FDTD book
            pml[ii].sigma_max_over_omega = 0.8*(pml[ii].power_sigma+1)/(k0dx*sqrt(epsilon_bg[ii])) 
        else
            temp = pml[ii].sigma_max_over_omega            
            if temp < 0
                throw(ArgumentError("$(direction)pml[$(ii)].sigma_max_over_omega must be a non-negative scalar."))
            end
        end

        # Maximal coordinate stretching factor
        if ~isdefined(pml[ii], :kappa_max)
            # From Table 7.1 of Taflove & Hagness's 2005 FDTD book
            pml[ii].kappa_max = 15
        else
            temp = pml[ii].kappa_max
            if temp < 1
                throw(ArgumentError("$(direction)pml[$ii].kappa_max must be a real scalar that equals or is larger than 1."))
            end
        end

        # Power of polynomial grading for the coordinate stretching factor kappa
        if ~isdefined(pml[ii], :power_kappa)
            # From Table 7.1 of Taflove & Hagness's 2005 FDTD book
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
