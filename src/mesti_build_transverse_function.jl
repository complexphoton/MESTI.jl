"""
    MESTI_BUILD_TRANSVERSE_FUNCTION sets up 1D transverse function and wave number.

        === Input Arguments ===
        nx (positive integer scalar; required):
            Number of grid points in the transverse direction. 
        xBC (string or scalar number; required):
            Boundary condition in the transverse direction. 
            When xBC is a character vector, available choices are (case-insensitive): 
                "periodic"              - p(n+nx) = p(n) 
                "Dirichlet"             - p(0) = p(nx+1,n) = 0 
                "Neumann"               - p(0) = p(1); p(nx+1) = p(nx) 
                "DirichletNeumann"      - p(0) = 0; p(nx+1) = p(nx) 
                "NeumannDirichlet"      - p(0) = p(1); p(nx+1) = 0 
            When xBC is a scalar number, the Bloch periodic boundary condition is
            used with p(n+nx) = p(n)*exp(1i*xBC); in other words, xBC = kx_B*nx*dx =
            kx_B*p where kx_B is the Bloch wave number and p = nx*dx is the
            periodicity in the transverse direction. 
        n0 (real numeric scalar, optional, defaults to 0):
            Center of the transverse mode along x with periodic or Bloch periodic
            boundary condition, u_{n,a} = exp(i*kx[a]*dx*(n-n0))/sqrt(nx), where
            kx(a) = kx_B + a*(2*pi/nx*dx).
        offset (logical scalar; optional, defaults to false):
            Whether to use the offset transverse function. Due to the staggered properties 
            of the Yee grid, it should be true for the x-direction for transverse function Ex.
            It only applies to xBC is a scalar number or xBC == "periodic", since we have 
            already included the offset in other BCs.

        === Output Arguments ===
        fun_f_1d (function):
            A function that, given one element of kxdx_all as the input, 
            returns its normalized transverse field profile as an nx vector;
            when the input is a vector, it returns a matrix where each column
            is the respective transverse profile. The transverse modes form a
            complete and orthonormal set, so the nx-by-nx matrix
            fun_f_1d(kxdx_all) is unitary for periodic (Bloch), Neumann, 
            DirichletNeumann, and NeumannDirichlet, but it is not the case for Dirichlet.
            For Dirichlet, we intend to maintain the trivial solution in this transverse function,
            since it would give us non-trivial solution in another component when we build up 
            2D transverse function.
        kxdx_all (1-by-nx_Ex+delta_(xBC,"Dirichlet") real row vector):
            Dimensionless transverse wave number kx*dx for all ny channels,
            including both propagating and evanescent ones. They are real-valued
            and are ordered from small to large. 

"""
function mesti_build_transverse_function(nx::Int, xBC::Union{String,Real}, n0::Real=0, offset::Bool=false)    
    # Check input parameters        
    if ~(nx>=0)
        throw(ArgumentError("Input argument nx must be a natural number."))
    end

    # Convert BC to take care of lowercase or uppercase
    xBC = convert_BC_1d(xBC, "x")
        
    # Handle periodic and Bloch periodic boundary conditions
    if isa(xBC, Number)
        ka_x = xBC
        xBC = "Bloch"
        # ka_x must be real for u_x(kxdx_all) to be unitary
        if ~isa(ka_x, Real)
            @warn("kx_B*a = $(real(ka_x)) + 1im*$(imag(ka_x)) is a complex number; must be real for a complete orthonormal transverse basis.")
        end
    elseif lowercase(xBC) == "periodic"
        ka_x = 0
        xBC = "Bloch"
    end

    # f = [p(1), ..., p(nx)].'
    # For periodic and Bloch periodic boundary, we order kxdx_all such that it increases monotonically from negative to positive
    # For other boundary conditions, kx >= 0, and we order kxdx_all such that it increases monotonically from smallest to largest
    # Transverse modes in x (form a complete basis in x)
    if xBC == "Bloch"
        # p(nx+1) = p(1)*exp(1i*ka_x); p(0) = p(nx)*exp(-1i*ka_x)
        # The transverse mode index where kxdx = ka_x/nx
        if mod(nx,2) == 1
            ind_zero_kx = round((nx+1)/2)
        else
            ind_zero_kx = round(nx/2)
        end
        kxdx_all = (ka_x/nx) .+ ((1:nx).-ind_zero_kx)*(2*pi/nx) 
    elseif xBC == "Dirichlet"
        # p(0) = p(nx+1) = 0
        # Note that here kxdx_all starts from 0
        # Even it gives us trivial solution in this transverse function,
        # but this kxdx = 0 would give us non-trivial solution in another component. 
        kxdx_all = (0:nx)*(pi/(nx+1))        
    elseif xBC == "Neumann"
        # p(0) = p(1), p(nx+1) = p(nx)
        kxdx_all = ((1:nx).-1)*(pi/(nx))        
    elseif xBC == "DirichletNeumann"
        # p(0) = 0, p(nx+1) = p(nx)
        kxdx_all = ((0.5:nx))*(pi/(nx+0.5))                
    elseif xBC == "NeumannDirichlet"
        # p(0) = p(1), p(ny+1) = 0
        kxdx_all = ((0.5:nx))*(pi/(nx+0.5))                        
    else
        error("Input argument xBC = $(xBC) is not a supported option.")        
    end

    # Build up the tranverse function with n = 1...nx   
    fun_f_1d(kxdx) =  if xBC == "Bloch"
                            # Dimensionless transverse mode profile: u_{n,a} = exp(1im*(n-n0)*kxdx(a))/sqrt(nx)
                            exp.(((1:nx).+0.5*(offset).-n0)*(1im*reshape(vcat(kxdx),1,:)))/sqrt(nx)
                      elseif xBC == "Dirichlet"
                            # Dimensionless transverse mode profile: u_{n,a} = i*sin(n*kxdx(a))*sqrt(2/(nx+1))
                            1im*sin.(((1:nx))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+1))
                      elseif xBC == "Neumann"
                            # Dimensionless transverse mode profile:
                            # When kxdx == 0: u_{n,a} = sqrt(1/nx)
                            # When kxdx != 0: u_{n,a} = cos((n-0.5)*kxdx(a))*sqrt(2/nx)
                            # So we subtract reshape(vcat(kxdx.== 0),1,:)*(1-sqrt(1/2)) from the cos() which is nonzero only when kxdx=0
                            (cos.(((0.5:nx))*reshape(vcat(kxdx),1,:)).-reshape(vcat(kxdx.== 0),1,:)*(1-sqrt(1/2)))*sqrt(2/(nx))
                      elseif xBC == "DirichletNeumann"
                            # Dimensionless transverse mode profile: u_{n,a} = i*sin(n*kxdx(a))*sqrt(2/(nx+0.5))
                            1im*sin.(((1:nx))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+0.5))
                      elseif xBC == "NeumannDirichlet"
                            # Dimensionless transverse mode profile: u_{n,a} = cos((n-0.5)*kxdx(a))*sqrt(2/(nx+0.5))        
                            cos.(((0.5:nx))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+0.5)) 
                      end 
    return fun_f_1d, kxdx_all
end
"""
    MESTI_BUILD_TRANSVERSE_FUNCTION_DERIVATIVE sets up derivative of a transverse function.

        === Input Arguments ===
        nx (positive integer scalar; required):
            Number of grid points in the transverse direction. 
        xBC (string or scalar number; required):
            Boundary condition in the transverse direction. 
            When xBC is a character vector, available choices are (case-insensitive): 
                "periodic"              - p(n+nx) = p(n) 
                "Dirichlet"             - p(0) = p(nx+1,n) = 0 
                "Neumann"               - p(0) = p(1); p(nx+1) = p(nx) 
                "DirichletNeumann"      - p(0) = 0; p(nx+1) = p(nx) 
                "NeumannDirichlet"      - p(0) = p(1); p(nx+1) = 0 
            When xBC is a scalar number, the Bloch periodic boundary condition is
            used with p(n+nx) = p(n)*exp(1i*xBC); in other words, xBC = kx_B*nx*dx =
            kx_B*p where kx_B is the Bloch wave number and p = nx*dx is the
            periodicity in x. 
        n0 (real numeric scalar, optional, defaults to 0):
            Center of the transverse mode along x with periodic or Bloch periodic
            boundary condition, u_{n,a} = exp(i*kx(a)*dx*(n-n0))/sqrt(nx), where
            kx(a) = kx_B + a*(2*pi/nx*dx).
        changegrid (real numeric scalar, optional, defaults to 0):
            When the derivative tranverse function want to add or subtract,
            another tranverse function, it may have dimension mismatch.
            This input primarily handle this scenario. 
        
        === Output Arguments ===
            fun_df_1d (function_handle):
                A function that, given one element of (kxdx_all) as the input, 
                returns its normalized derivative transverse field profile as an nx vector;
                when the input is a vector, it returns a matrix where each column
                is the respective derivative transverse profile. 
"""
function mesti_build_transverse_function_derivative(nx::Int, xBC::Union{String,Real}, n0::Real=0, changegrid::Real=0)    
    # Check input parameters        
    if ~(nx>=0)
        error("Input argument nx must be a natural number.")
    end

    # Convert BC to take care of lowercase or uppercase
    xBC = convert_BC_1d(xBC, "x")
        
    # Handle periodic and Bloch periodic boundary conditions
    if isa(xBC, Number)
        ka_x = xBC
        xBC = "Bloch"
        # ka_x must be real for u_x(kxdx_all) to be unitary
        if ~isa(ka_x, Real)
            @warn("kx_B*a = $(real(ka_x)) + 1im*$(imag(ka_x)) is a complex number; must be real for a complete orthonormal transverse basis.")
        end
    elseif lowercase(xBC) == "periodic"
        ka_x = 0
        xBC = "Bloch"
    end
    
    # Forward difference of dimensionless transverse mode profile
    fun_df_1d(kxdx) = if xBC == "Bloch"
                            exp.(((2:(nx+1)).-n0)*(1im*reshape(vcat(kxdx),1,:)))/sqrt(nx)-
                            exp.(((1:nx).-n0)*(1im*reshape(vcat(kxdx),1,:)))/sqrt(nx)
                      elseif xBC == "Dirichlet"
                            1im*sin.((((2-changegrid):(nx+1)))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+1)) -        
                            1im*sin.((((1-changegrid):nx))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+1))
                      elseif xBC == "Neumann"
                            (cos.(((1.5:(nx+1-changegrid)))*reshape(vcat(kxdx),1,:)).-reshape(vcat(kxdx.== 0),1,:)*(1-sqrt(1/2)))*sqrt(2/(nx))-
                            (cos.(((0.5:(nx-changegrid)))*reshape(vcat(kxdx),1,:)).-reshape(vcat(kxdx.== 0),1,:)*(1-sqrt(1/2)))*sqrt(2/(nx))
                      elseif xBC == "DirichletNeumann"
                            1im*sin.((((2-changegrid):(nx+1-changegrid)))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+0.5))-        
                            1im*sin.((((1-changegrid):(nx-changegrid)))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+0.5))
                      elseif xBC == "NeumannDirichlet"
                            cos.(((1.5:(nx+1)))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+0.5))-        
                            cos.(((0.5:nx))*reshape(vcat(kxdx),1,:))*sqrt(2/(nx+0.5)) 
                      end 
    return fun_df_1d
end

"""
    CONVERT_BC_1D is a helper function to handle the case of string in BC
"""
function convert_BC_1d(BC::Union{String,Real,Complex},direction::String)
    if isa(BC, Number)
        return BC
    elseif lowercase(BC) == "dirichlet"
        return "Dirichlet"
    elseif lowercase(BC) == "neumann"
        return "Neumann"
    elseif lowercase(BC) == "dirichletneumann"
        return "DirichletNeumann"
    elseif lowercase(BC) == "neumanndirichlet"
        return "NeumannDirichlet"
    elseif lowercase(BC) == "periodic"
        return "periodic"
    elseif lowercase(BC) == "bloch"
        error("To use Bloch periodic boundary condition in $(direction)-direction, set the input argument $(direction)BC to k$(direction)_B*p_$(direction) where k$(direction)_B is the Bloch wave number and p_$(direction) is the periodicity along $(direction)-direction.")        
    else
        error("Input argument $(direction)BC = \"$(BC)\" is not a supported option.")         
    end
end
