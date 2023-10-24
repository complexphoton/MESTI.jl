# Export composite data types
export Matrices
export Opts
export Info

# Export a function mesti_matrix_solver!()
export mesti_matrix_solver!

mutable struct Matrices
    # A composite data type to store the matrices, A, B, and C
    # After construct matrix K, we will set them nothing to reduce memory usage
    # See also: mesti_matrix_solver!
    A::Union{SparseMatrixCSC{Int64, Int64},SparseMatrixCSC{Complex{Int64}, Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64, Int64},Nothing}
    B::Union{SparseMatrixCSC{Int64, Int64},SparseMatrixCSC{Complex{Int64}, Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64, Int64},Array{Int64, 2},Array{Float64, 2},Array{ComplexF64, 2},Nothing}
    C::Union{SparseMatrixCSC{Int64, Int64},SparseMatrixCSC{Complex{Int64}, Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64, Int64},Array{Int64, 2},Array{Float64, 2},Array{ComplexF64, 2},Nothing,String}
    Matrices() = new()
end

mutable struct Opts
    # A composite data type to store options
    # See also: mesti_matrix_solver!, mesti2s, and mesti 
    is_symmetric_A::Union{Integer, Nothing}
    verbal::Integer    
    prefactor::Union{Number, Nothing}
    solver::Union{String, Nothing}     
    method::Union{String, Nothing}
    clear_BC::Integer
    clear_syst::Integer    
    clear_memory::Integer     
    verbal_solver::Integer 
    use_single_precision_MUMPS::Integer
    use_METIS::Integer    
    nrhs::Integer 
    store_ordering::Integer
    ordering    
    analysis_only::Integer    
    nthreads_OMP::Integer    
    iterative_refinement::Integer

    exclude_PML_in_field_profiles::Integer
    return_field_profile::Integer
    use_given_ordering::Integer

    n0::Union{Real,Nothing}
    m0::Union{Real,Nothing}
    use_continuous_dispersion::Union{Integer,Nothing}    
    symmetrize_K::Union{Integer, Nothing}
    nz_low::Union{Integer, Nothing}
    nz_high::Union{Integer, Nothing}

    # the following four are for block low-rank
    use_BLR::Integer
    threshold_BLR::Real
    icntl_36::Integer
    icntl_38::Integer
    
    # this is only for MUMPS solver
    parallel_dependency_graph::Integer

    Opts() = new()
end

mutable struct Info
    # A composite data type to store information 
    # for users to look up    
    opts::Opts
    timing_total::Real
    timing_init::Real
    timing_build::Real
    timing_analyze::Real
    timing_factorize::Real
    timing_solve::Real
    ordering_method # This datatype changes during the code
    ordering # TODO:: Test and check this data type
    itr_ref_nsteps::Integer
    itr_ref_omega_1::Array{Real}
    itr_ref_omega_2::Array{Real}
    xPML::Vector{PML}
    yPML::Vector{PML}
    zPML::Vector{PML}

    # Below are used in mesti2s() only    
    ind_in_trivial_ch::Vector{Int64}
    ind_out_trivial_ch::Vector{Int64}
    ind_in_nontrivial_ch::Vector{Int64}
    ind_out_nontrivial_ch::Vector{Int64}
    
    Info() = new()
end

"""
    MESTI_MATRIX_SOLVER! Computes matrices.C*inv(matrices.A)*matrices.B or inv(matrices.A)*matrices.B.
        (X, info) = MESTI_MATRIX_SOLVER!(matrices), when matrices.A != nothing, matrices.B != nothing, 
        and matrices.C = nothing, returns X = inv(matrices.A)*matrices.B for sparse matrix matrices.A 
        and (sparse or dense) matrix matrices.B, with the information of the computation 
        returned in structure "info".
        
        (S, info) = MESTI_MATRIX_SOLVER!(matrices), when matrices.A != nothing, matrices.B != nothing, 
        and matrices.C != nothing, returns S = matrices.C*inv(matrices.A)*matrices.B 
        where matrix matrices.C is either sparse or dense. When the MUMPS3 is available, this is done by 
        computing the Schur complement of an augmented matrix K = [matrices.A,matrices.B;matrices.C,0] 
        through a partial factorization.
        
        (X, info) = MESTI_MATRIX_SOLVER!(matrices, opts), when matrices.A != nothing, matrices.B != nothing, 
        and matrices.C = nothing, and
        (S, info) = MESTI_MATRIX_SOLVER!(matrices, opts), when matrices.A != nothing, matrices.B != nothing, 
        and matrices.C != nothing, 
        allow detailed options to be specified with structure "opts" of the input arguments.
        
        === Input Arguments ===
        matrices (scalar structure; required):
            A structure specifies matrices. It contain the following fields:
            A (sparse matrix; required):
                Matrix A in the C*inv(A)*B or inv(A)*B returned.
            B (numeric matrix; required):
                Matrix B in the C*inv(A)*B or inv(A)*B returned.
            C (numeric matrix or "transpose(B)" or nothing; optional):
                Matrix C in the C*inv(A)*B returned.
                    If C = transpose(B), the user can set C = "transpose(B)" as a
                character vector, which will be replaced by transpose(B) in the code. If
                matrix A is symmetric, C = "transpose(B)", and opts.method = "APF", the
                matrix K = [A,B;C,0] will be treated as symmetric when computing its
                Schur complement to lower computing time and memory usage.
                    To compute X = inv(A)*B, the user can simply omit C from the input
                argument if there is no need to change the default opts. If opts is
                needed, the user can set C = nothing here.
        opts (scalar structure; optional):
            A structure that specifies the options of computation; defaults to an
            undefined Opts() structure. It can contain the following fields (all optional):
            opts.verbal (logical scalar; optional, defaults to true):
                Whether to print system information and timing to the standard output.
            opts.is_symmetric_A (logical scalar; optional):
                Whether matrix A is symmetric or not. This is only used when
                opts.solver = "MUMPS", in which case opts.is_symmetric_A will be
                determined by the issymmetric(A) command if not specified by users.
            opts.solver (character vector; optional):
                The solver used for sparse matrix factorization. Available choices are
                (case-insensitive):
                    "MUMPS"  - (default when MUMPS is available) Use MUMPS. Its JULIA 
                            interface MUMPS3.jl must be installed.
                    "JULIA" -  (default when MUMPS is not available) Uses the built-in 
                            lu() function in JULIA, which uses UMFPACK. 
                MUMPS is faster and uses less memory than lu(), and is required for
                the APF method.
            opts.method (character vector; optional):
                The solution method. Available choices are (case-insensitive):
                    "APF" - Augmented partial factorization. C*inv(A)*B is obtained
                            through the Schur complement of an augmented matrix
                            K = [A,B;C,0] using a partial factorization. Must have
                            opts.solver = "MUMPS". This is the most efficient method,
                            but it cannot be used for computing X=inv(A)*B or with
                            iterative refinement.
                    "FG"  - Factorize and group. Factorize A=L*U, and obtain C*inv(A)*B
                            through C*inv(U)*inv(L)*B with optimized grouping. Must
                            have opts.solver = "JULIA". This is slightly better than
                            "FS" when MUMPS is not available, but it cannot be used for
                            computing X=inv(A)*B.
                    "FS"  - Factorize and solve. Factorize A=L*U, solve for X=inv(A)*B
                            with forward and backward substitutions, and project with
                            C as C*inv(A)*B = C*X. Here, opts.solver can be either
                            "MUMPS" or "JULIA", and it can be used for computing
                            X=inv(A)*B or with iterative refinement.
                    "C*inv(U)*inv(L)*B"   - Same as "FG".    
                    "factorize_and_solve" - Same as "FS".
                By default, if C is given and opts.iterative_refinement = false, then
                "APF" is used when opts.solver = "MUMPS", and "C*inv(U)*inv(L)*B" is
                used when opts.solver = "JULIA". Otherwise, "factorize_and_solve" is
                used.
            opts.verbal_solver (logical scalar; optional, defaults to false):
                Whether to have the solver print detailed information to the standard
                output. Note the behavior of output from MUMPS depends on compiler.
            opts.clear_memory (logical scalar; optional, defaults to false):
                Whether or not to clear variables to reduce peak memory usage. When
                opts.clear_memory = true, the following variables may be cleared in
                the caller's workspace if they exist: A, B, C. Some other variables
                inside mesti_matrix_solver() will be cleared too. However, currently,
                we are trying to figure out how to do it in JULIA language.
            opts.use_single_precision_MUMPS (boolean scalar; optional, defaults to true):
                Whether to use single precision version of MUMPS; used only when 
                opts.solver = "MUMPS". Using single precision version of MUMPS can 
                reduce memory usage and computing time.
            opts.use_METIS (logical scalar; optional, defaults to false):
                Whether to use METIS (instead of the default AMD) to compute the
                ordering in MUMPS. Using METIS can sometimes reduce memory usage
                and/or factorization and solve time, but it typically takes longer at
                the analysis (i.e., ordering) stage stage in 2D. In 3D METIS is general 
                better than AMD.
            opts.nrhs (positive integer scalar; optional):
                The number of right-hand sides (number of columns of matrix B) to
                consider simultaneously, used only when opts.method =
                "factorize_and_solve" and C is given. Defaults to 1 if
                opts.iterative_refinement = true, 10 if opts.solver = "MUMPS" with
                opts.iterative_refinement = false, 4 otherwise.
            opts.store_ordering (logical scalar; optional, defaults to false):
                Whether to store the ordering sequence (permutation) for matrix A or
                matrix K; only possible when opts.solver = "MUMPS". If
                opts.store_ordering = true, the ordering will be returned in
                info.ordering.
            opts.ordering (positive integer vector; optional):
                A user-specified ordering sequence for matrix A or matrix K, used only
                when opts.solver = "MUMPS". Using the ordering from a previous
                computation can speed up the analysis stage, but the matrix size must
                be the same.
            opts.analysis_only (logical scalar; optional, defaults to false):
                When opts.analysis_only = true, the factorization and solution steps
                will be skipped, and S = nothing will be returned. The user can use
                opts.analysis_only = true with opts.store_ordering = true to return
                the ordering for A or K; only possible when opts.solver = 'MUMPS'.
            opts.nthreads_OMP (positive integer scalar; optional):
                Number of OpenMP threads used in MUMPS; overwrites the OMP_NUM_THREADS
                environment variable.
            opts.parallel_dependency_graph (logical scalar; optional):
                If MUMPS is multithread, whether to use parallel dependency graph in MUMPS.
                This typically improve the time performance, but marginally increase 
                the memory usage.
            opts.iterative_refinement (logical scalar; optional, defaults to false):
                Whether to use iterative refinement in MUMPS to lower round-off
                errors. Iterative refinement can only be used when opts.solver =
                "MUMPS" and opts.method = "factorize_and_solve" and C is given, in
                case opts.nrhs must equal 1. When iterative refinement is used, the
                relevant information will be returned in info.itr_ref_nsteps,
                info.itr_ref_omega_1, and info.itr_ref_omega_2.
            use_BLR (logical scalar; optional, defaults to false):
                Whether to use block low-rank approximation in MUMPS to possibly lower computational
                cost (but in most case it does not). It can only be used when opts.solver = "MUMPS".
            threshold_BLR (positive real scalar; optional):
                The dropping parameter controls the accuracy of the block low-rank approximations. 
                It can only be used when opts.solver = "MUMPS" and opts.use_BLR = true.
                Please refer to the section of BLR API in MUMPS userguide.
            icntl_36 (positive integer scalar; optional):
                It controls the choice of the BLR factorization variant. 
                It can only be used when opts.solver = "MUMPS" and opts.use_BLR = true.
                Please refer to the section of BLR API in MUMPS userguide.
            icntl_38 (positive integer scalar; optional):
                It estimated compression rate of LU factors.
                It can only be used when opts.solver = "MUMPS" and opts.use_BLR = true.
                Please refer to the section of BLR API in MUMPS userguide.

        === Output Arguments ===
        S (full numeric matrix):
            C*inv(A)*B or inv(A)*B.
        info (scalar structure):
            A structure that contains the following fields:
            info.opts (scalar structure):
                The final "opts" used, excluding the user-specified matrix ordering.
            info.timing_init (numeric scalar):
                Timing of the initial stages, in seconds
            info.timing_build (numeric scalar):
                Timing of the building stages, in seconds
            info.timing_analyze (numeric scalar):
                Timing of the analyzing stages, in seconds
            info.timing_factorize (numeric scalar):
                Timing of the factorizing stages, in seconds
            info.timing_solve (numeric scalar):
                Timing of the solving stages, in seconds    
            info.timing_total (numeric scalar):
                Timing of the total, in seconds    
            info.ordering_method (character vector; optional):
                Ordering method used in MUMPS.
            info.ordering (positive integer vector; optional):
                Ordering sequence returned by MUMPS when opts.store_ordering = true.
            info.itr_ref_nsteps (integer vector; optional):
                Number of steps of iterative refinement for each input, if
                opts.iterative_refinement = true; 0 means no iterative refinement.
            info.itr_ref_omega_1 (real vector; optional):
                Scaled residual omega_1 at the end of iterative refinement for each
                input; see MUMPS user guide section 3.3.2 for definition.
            info.itr_ref_omega_2 (real vector; optional):
                Scaled residual omega_2 at the end of iterative refinement for each
                input; see MUMPS user guide section 3.3.2 for definition.
"""
function mesti_matrix_solver!(matrices::Matrices, opts::Union{Opts,Nothing}=nothing)
    opts = deepcopy(opts)
    
    ## Part 1: Initialization
    ## Check validity & consistency of input arguments and assign default values
    t0 = time() 
        
    if ~isdefined(matrices, :C)
         matrices.C = nothing
    end
    
    if isa(matrices.C, Nothing)  
        # return_X = isa(matrices.C, Nothing), but is easier to understand the meaning when we use return_X
        return_X = true
    elseif isa(matrices.C, AbstractArray{<:Number})
        return_X = false
        use_transpose_B = false
    elseif matrices.C == "transpose(B)"
        return_X = false
        use_transpose_B = true
    else
        throw(ArgumentError("Input argument matrices.C must be a numeric matrix or \"transpose(B)\" or nothing, if given."))
    end

    if ~(isa(opts, Nothing) || isa(opts, Opts))
        throw(ArgumentError("Input argument opts must be a Opts structure or nothing, if given."))
    end

    if isa(opts, Nothing)
        opts = Opts()
    end

    # Check that the user did not accidentally use options only in mesti2s()
    if isdefined(opts, :symmetrize_K) && !isa(opts.symmetrize_K, Nothing)
        throw(ArgumentError("opts.symmetrize_K is not used in mesti_matrix_solver(); to symmetrize matrix K = [A,B;C,0], set matrices.C = \"transpose(B)\", make sure matrix A is symmetric, set opts.solver = \"MUMPS\", and set opts.method = \"APF\"."))
    end

    # Turn on verbal output by default
    if ~isdefined(opts, :verbal)
        opts.verbal = true
    elseif ~isa(opts.verbal, Bool)
        throw(ArgumentError("opts.verbal must be a boolean, if given."))
    end
    
    # Use MUMPS for opts.solver when it is available
    MUMPS_available = @isdefined(Mumps)
    if ~isdefined(opts, :solver) || isa(opts.solver, Nothing)
        if MUMPS_available
            opts.solver = "MUMPS"
        else
            opts.solver = "JULIA"
        end
    else
        opts.solver = uppercase(opts.solver)
        if ~(opts.solver == "MUMPS") && ~(opts.solver == "JULIA") 
            throw(ArgumentError("opts.solver = \"$(opts.solver)\" is not a supported option; use \"MUMPS\" or \"JULIA\"."))
        elseif opts.solver == "MUMPS" && ~MUMPS_available
            throw(ArgumentError("opts.solver = \"$(opts.solver)\" but package MUMPS3.jl is not found."))
        end
    end
    
    # No iterative refinement by default; only used in factorize_and_solve when computing S=C*inv(A)*B with MUMPS
    str_itr_ref = nothing
    if ~isdefined(opts, :iterative_refinement)
        opts.iterative_refinement = false
    elseif ~isa(opts.iterative_refinement, Bool)
        throw(ArgumentError("opts.iterative_refinement must be a boolean, if given."))
    elseif opts.iterative_refinement
        throw(ArgumentError("opts.iterative_refinement = true is not supported for current version, but it will be supported for the future version."))
        str_itr_ref = "with iterative refinement"
    end

    # By default, if C is given and opts.iterative_refinement = false, then "APF" is used when opts.solver = "MUMPS", and "C*inv(U)*inv(L)*B" is used when opts.solver = "JULIA". Otherwise, "factorize_and_solve" is used.
    if ~isdefined(opts, :method) || isa(opts.method, Nothing)
        if return_X || opts.iterative_refinement
            opts.method = "factorize_and_solve"
        else
            if uppercase(opts.solver) == "MUMPS"
                opts.method = "APF"
            else
                opts.method = "C*inv(U)*inv(L)*B"
            end                        
        end
    elseif ~((lowercase(opts.method) in ["apf", "fs", "factorize_and_solve", "fg", "c*inv(u)*inv(l)*b"]))
        throw(ArgumentError("opts.method = \"$(opts.method)\" is not a supported option; use \"APF\", \"factorize_and_solve\", or \"C*inv(U)*inv(L)*B\"."))
    elseif return_X && ~((lowercase(opts.method) in ["fs", "factorize_and_solve"]))
        throw(ArgumentError("opts.method = \"$(opts.method)\" cannot be used when matrices.C = nothing; use opts.method = \"factorize_and_solve\" instead."))
    elseif lowercase(opts.method) == "fs" || lowercase(opts.method) == "factorize_and_solve"
        opts.method = "factorize_and_solve" # opts.method = "FS" is short for opts.method = "factorize_and_solve"
    elseif lowercase(opts.method) == "apf"
        opts.method = "APF"
    elseif lowercase(opts.method) == "fg" || lowercase(opts.method) == "c*inv(u)*inv(l)*b"
        opts.method = "C*inv(U)*inv(L)*B"        
    end

    # When opts.method = "APF" and opts.solver = "JULIA", the solution method is not actually APF, so we throw an error to the users
    str_method = opts.method
    if opts.method == "APF" && opts.solver == "JULIA"
        throw(ArgumentError("opts.method = \"APF\" requires opts.solver = \"MUMPS\"."))
    end
    if opts.method == "C*inv(U)*inv(L)*B" && opts.solver == "MUMPS"
        throw(ArgumentError("opts.method = \"C*inv(U)*inv(L)*B\" requires opts.solver = \"JULIA\"."))
    end
    if opts.iterative_refinement && ~(~return_X && opts.method == "factorize_and_solve" && opts.solver == "MUMPS")
        throw(ArgumentError("To use opts.iterative_refinement = true, input argument matrices.C must not be nothing, opts.method must be \"factorize_and_solve\", and opts.solver must be \"MUMPS\".\nHere isa(matrices.C, Nothing) = $(isa(matrices.C, Nothing)), opts.method = \"$(opts.method)\", opts.solver = \"$(opts.solver)\"."))
    end

    # Turn off solver's verbal output by default
    if ~isdefined(opts, :verbal_solver)
        opts.verbal_solver = false
    elseif ~isa(opts.verbal_solver, Bool)
        throw(ArgumentError("opts.verbal_solver must be a boolean, if given."))
    end

    # Determine whether matrix matrices.C will be used
    if return_X
        use_C = false
    elseif use_transpose_B
        if opts.method == "APF"
            # In this case, we keep matrices.C = "transpose(B)" here and use transpose(B) later so the memory of transpose(B) can be automatically cleared after use
            use_C = false
        else
            # In other cases, we may as well allocate the memory for matrices.C now
            matrices.C = permutedims(matrices.B, (2,1))
            use_C = true
        end
    else
        use_C = true
    end

    # At this point, there are two possibilities for which use_C = false:
    # (1) matrices.C = nothing, return_X = true, opts.method = "factorize_and_solve"
    # (2) matrices.C = "transpose(B)", return_X = false, use_transpose_B = true, opts.method = "APF", opts.solver = "MUMPS"

    # Check matrix sizes
    (sz_A_1, sz_A_2) = size(matrices.A)
    (sz_B_1, sz_B_2) = size(matrices.B)
    if isa(matrices.C, Nothing) || matrices.C == "transpose(B)"
        sz_C_1=0; sz_C_2=0
    else
        (sz_C_1, sz_C_2) = size(matrices.C)
    end
    if sz_A_1 != sz_A_2; throw(ArgumentError("Input argument A must be a square matrix; size(A) = [$(sz_A_1), $(sz_A_2)].")); end
    if sz_A_2 != sz_B_1; throw(ArgumentError("size(matrices.A,2) must equal size(matrices.B,1); size(matrices.A,2) = $(sz_A_2), size(matrices.B,1) = $(sz_B_1).")); end
    if (sz_C_2 != sz_A_1 && use_C); throw(ArgumentError("size(matrices.C,2) must equal size(matrices.A,1); size(matrices.C,2) = $(sz_C_2), size(matrices.A,1) = $(sz_A_1).")); end

    # By default, we don't clear variables unless specified by users
    if ~isdefined(opts, :clear_memory)
        opts.clear_memory = false
    elseif ~isa(opts.clear_memory, Bool)
        throw(ArgumentError("opts.clear_memory must be a logical scalar, if given."))
    end
    
    # Go over options only relevant when MUMPS is used
    opts.use_given_ordering = false
    str_ordering = nothing
    str_sym_K = nothing
    str_MUMPS_precision = nothing
    if opts.solver == "MUMPS"
        # Determine the symmetry of matrix matrices.A if not specified
        # To skip this step (which can be slow), the user should specify opts.is_symmetric_A
        if ~isdefined(opts, :is_symmetric_A)
            opts.is_symmetric_A = issymmetric(matrices.A)
        elseif ~isa(opts.is_symmetric_A, Bool)
            throw(ArgumentError("opts.is_symmetric_A must be a boolean, if given."))
        end

        # Whether matrix K = [A,B;C,0] will be treated as symmetric
        if opts.method == "APF" && use_transpose_B && opts.is_symmetric_A
            str_sym_K = " (symmetric K)"
        end

        # Use AMD by default because it does not require additional efforts to compile METIS
	if ~isdefined(opts, :use_METIS)
	    opts.use_METIS = false
	elseif ~isa(opts.use_METIS, Bool)
	    throw(ArgumentError("opts.use_METIS must be a boolean, if given."))   
	end                 

        if opts.use_METIS
            str_ordering = " with METIS ordering"
        else
            str_ordering = " with AMD ordering"
        end
        
        # Use single-precision MUMPS by default
        if ~isdefined(opts, :use_single_precision_MUMPS)
            opts.use_single_precision_MUMPS = true
        elseif ~isa(opts.use_single_precision_MUMPS, Bool)
            throw(ArgumentError("opts.use_single_precision_MUMPS must be a boolean, if given."))   
        end                            
        if opts.use_single_precision_MUMPS
            str_MUMPS_precision = " in single precision"
        else
            str_MUMPS_precision = " in double precision"
        end        
        
        # We don't use KEEP(401) = 1 by default        
        if ~isdefined(opts, :parallel_dependency_graph)
            opts.parallel_dependency_graph = false
        elseif ~isa(opts.parallel_dependency_graph, Bool)
            throw(ArgumentError("opts.parallel_dependency_graph must be a boolean, if given."))    
        end

        # We don't use BLR by default
        if ~isdefined(opts, :use_BLR)
            opts.use_BLR = false
        elseif ~isa(opts.use_BLR, Bool)
            throw(ArgumentError("opts.use_BLR must be a boolean, if given."))
        end
        
        # If BLR is activated, we check or pick the threshold of BLR.
        if ~opts.use_BLR && isdefined(opts, :threshold_BLR)
            throw(ArgumentError("opts.threshold_BLR should not be given when opts.use_BLR = false."))    
        elseif opts.use_BLR && isdefined(opts, :threshold_BLR) && opts.threshold_BLR < 0
            throw(ArgumentError("opts.threshold_BLR should be positive, if given"))   
        elseif opts.use_BLR && ~isdefined(opts, :threshold_BLR)
            opts.threshold_BLR = 1e-7
        end
        
        # We don't store matrix ordering by default
        if ~isdefined(opts, :store_ordering)
            opts.store_ordering = false
        elseif ~isa(opts.store_ordering, Bool)
            throw(ArgumentError("opts.store_ordering must be a boolean, if given."))
        elseif opts.store_ordering          
            throw(ArgumentError("opts.store_ordering = true is not supported for current version, but it will be supported for the future version."))
        end

        # Use the user-specified ordering, if given
        if isdefined(opts, :ordering)
            opts.use_given_ordering = true
            if opts.use_given_ordering
                throw(ArgumentError("opts.ordering is not supported for current version, but it will be supported for the future version."))        
            end
            str_ordering = " with user-specified ordering"
            if opts.use_METIS
                throw(ArgumentError("opts.use_METIS cannote be true when opts.ordering is given."))
            end
        end

        # Whether to skip the factorization and solution steps
        if ~isdefined(opts, :analysis_only)
            opts.analysis_only = false
        elseif ~isa(opts.analysis_only, Bool)
            throw(ArgumentError("opts.analysis_only must be a logical scalar, if given."))
        elseif opts.analysis_only && ~opts.store_ordering
            throw(ArgumentError("When opts.analysis_only = true, opts.store_ordering must also be true."))
        end
        
        # Number of openMP threads in MUMPS; leave empty if not specified
        if isdefined(opts, :nthreads_OMP)
            if ~(opts.nthreads_OMP>0)
                throw(ArgumentError("opts.nthreads_OMP must be a positive integer, if given."))
            end
        end
    else
        if isdefined(opts, :is_symmetric_A) && ~isa(opts.is_symmetric_A, Nothing)
            opts.is_symmetric_A = nothing	    
        end

        #if isdefined(opts, :use_METIS)
        #    @warn("opts.use_METIS is only used when opts.solver = \"MUMPS\"; will be ignored.")
        #end

        if isdefined(opts, :use_single_precision_MUMPS)
            @warn("opts.use_single_precision_MUMPS is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end
        
        if isdefined(opts, :store_ordering) && opts.store_ordering == true
            throw(ArgumentError("opts.store_ordering = true can only be used when opts.solver = \"MUMPS\"."))
        end

        if isdefined(opts, :ordering)
            @warn("opts.ordering is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end

        if isdefined(opts, :analysis_only) && opts.analysis_only == true
            throw(ArgumentError("opts.analysis_only = true can only be used when opts.solver = \"MUMPS\"."))
        end
        opts.analysis_only = false
                
        if isdefined(opts, :nthreads_OMP)
            @warn("opts.nthreads_OMP is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end
        
        if isdefined(opts, :use_BLR)
            @warn("opts.use_BLR is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end

        if isdefined(opts, :threshold_BLR)
            @warn("opts.threshold_BLR is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end
        
        if isdefined(opts, :icntl_36)
            @warn("opts.icntl_36 is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end

        if isdefined(opts, :icntl_38)
            @warn("opts.icntl_38 is only used when opts.solver = \"MUMPS\"; will be ignored.")
        end        
    end

    # Number of columns to solve for simultaneously; only used in factorize_and_solve when computing S=C*inv(A)*B
    str_nrhs = nothing
    if opts.method == "factorize_and_solve" && ~return_X
        if ~isdefined(opts, :nrhs)
            if opts.solver == "MUMPS"
                if opts.iterative_refinement
                    opts.nrhs = 1 # iterative refinement requires nrhs = 1
                else
                    opts.nrhs = min(10, sz_B_2)
                end
            else
                opts.nrhs = min(4, sz_B_2)
            end
        elseif ~(opts.nrhs>0)
                throw(ArgumentError("opts.nrhs must be a positive integer, if given."))
        elseif opts.iterative_refinement && opts.nrhs != 1
            throw(ArgumentError("When opts.iterative_refinement = true, opts.nrhs must be 1, if given."))
        end
        str_nrhs = @sprintf(" with nrhs = %d", opts.nrhs)
    elseif isdefined(opts, :nrhs)
        @warn("opts.nrhs is not used when opts.method = \"$(opts.method)\" and isempty(C) = $(isa(C, Nothing)); will be ignored.")
    end

    t2 = time() 
    timing_init = t2-t0

    ## Computation Part
    # No need to compute if length(S) = 0 and we don't need to keep the ordering
    if (sz_B_2 == 0 || (sz_C_1 == 0 && use_C)) && ~opts.store_ordering
        opts.method = "None"
        opts.solver = "None"
        if opts.verbal; @printf("No computation needed\n"); end
    elseif opts.verbal
        #@printf("< Method: %s using %s%s%s%s%s >\n", opts.method, opts.solver, str_nrhs, str_ordering, str_itr_ref, str_sym_K)
        @printf("< Method: %s using %s", opts.method, opts.solver)
        if ~isa(str_MUMPS_precision, Nothing); @printf("%s", str_MUMPS_precision); end        
        if ~isa(str_nrhs, Nothing); @printf("%s", str_nrhs); end
        if ~isa(str_ordering, Nothing); @printf("%s", str_ordering); end
        if ~isa(str_itr_ref, Nothing); @printf("%s", str_itr_ref); end
        if ~isa(str_sym_K, Nothing); @printf("%s", str_sym_K); end
        @printf(" >\n")
    end

    if opts.method == "None"
        S = zeros(sz_C_1, sz_B_2)
        info.timing_build = 0
        info.timing_analyze = 0
        info.timing_factorize = 0
        info.timing_solve = 0
    elseif opts.method == "APF"
    ## Compute S=C*inv(A)*B with APF (augmented partial factorization)
        # Build matrix K=[A,B;C,0] and use MUMPS to compute its Schur complement -C*inv(A)*B with the LU factors discarded.
        t1 = time()
        if opts.verbal; @printf("Building K  ... "); end
        N = size(matrices.A,1)
        is_symmetric_K = opts.is_symmetric_A
        if use_transpose_B
            M_tot = size(matrices.B, 2)
            D = spzeros(M_tot, M_tot) # zero matrix
            K = hcat(vcat(matrices.A,transpose(matrices.B)), vcat(matrices.B,D)) # matrices.C = transpose(B)
        else
            # pad zeros so that size(B,2) = size(C,1)
            M_in = size(matrices.B,2)
            M_out = size(matrices.C,1)
            M_tot = max(M_in, M_out)
            if M_tot > M_in
                # pad M_tot-M_in columns of zeros to B
                matrices.B = hcat(matrices.B, spzeros(N, M_tot-M_in))
            elseif M_tot > M_out
                # pad M_tot-M_out rows of zeros to C
                matrices.C = vcat(matrices.C, spzeros(M_tot-M_out, N))
            end
            D = spzeros(M_tot, M_tot) # zero matrix
            K = hcat(vcat(matrices.A,matrices.C), vcat(matrices.B,D))
            is_symmetric_K = false # even if A is symmetric, generally C won't equal transpose(B); we will not check whether C equals B.' or not; the user should set C = "transpose(B)" if C=B.'
        end        
        if opts.clear_memory
            matrices.A = nothing; matrices.B = nothing; matrices.C = nothing; D = nothing;
            GC.gc()
        end
        ind_schur = N .+ (1:M_tot) # indices for the Schur variables; must be a row vector
        
        t2 = time(); timing_build = t2-t1
        if opts.verbal; @printf("elapsed time: %7.3f secs\n", timing_build); end
        
        # Call MUMPS to analyze and compute the Schur complement (using a partial factorization)
        # This is typically the most memory-consuming part of the whole simulation
        (id, info) = MUMPS_analyze_and_factorize(K, opts, is_symmetric_K, ind_schur)
        
        info.timing_build = timing_build  # the build time for A, B, C will be added in addition to this
        t1 = time()

        if opts.analysis_only
            S = nothing
        else
            # Retrieve C*inv(A)*B = -H = -K/A, stored as a dense matrix
            S = -get_schur_complement(id)

            # Remove the padded zeros
            if ~use_transpose_B
                if M_tot > M_in
                    S = S[:, 1:M_in]
                elseif M_tot > M_out
                    S = S[1:M_out, :]
                end
            end
        end

        # Destroy the MUMPS instance and deallocate memory
        finalize!(id)

        #MPI.Finalize() # We cannot finalize the MPI, since we may do another computation

        if opts.clear_memory
            GC.gc()
        end        

        t2 = time()

        # When opts.analysis_only = true, the factorization time could be
        # nonzero due to the instance termination time.            
        info.timing_factorize = info.timing_factorize + t2-t1
        info.timing_solve = 0
    elseif opts.method == "C*inv(U)*inv(L)*B"               
    ## Compute C*inv(U)*inv(L)*B where A=LU, with the order of multiplication based on matrix nnz
         # Factorize as P*inv(R)*A*Q = L*U where R is diagonal, L and U are lower and upper triangular, and P and Q are permutation matrices
         # For simplicity, we refer to this as A = L*U below
        (L, U, P, Q, R, info) = JULIA_factorize(matrices.A, opts)
        info.timing_build = 0  # the build time for A, B, C will be added in addition to this
        if opts.clear_memory
            matrices.A = nothing
            GC.gc()
        end

        # Here, we evaluate C*inv(A)*B, not necessarily as C*[inv(A)*B], but more generally as C*inv(U)*inv(L)*B.
        # The full expression is C*inv(A)*B = C*Q*inv(U)*inv(L)*P*inv(R)*B.
        # There are a few ways to group the mldivide or mrdivide operations. Like matrix multiplications, it is generally faster and more memory efficient to group such that we operate onto the side with fewer elements first.
        if opts.verbal; @printf("Solving     ... "); end
        t1 = time()
        if ~isa(matrices.B, SparseMatrixCSC); matrices.B = sparse(matrices.B); end
        if ~isa(matrices.C, SparseMatrixCSC); matrices.C = sparse(matrices.C); end
        nnz_B = nnz(matrices.B)
        nnz_C = nnz(matrices.C)
        if nnz_B <= nnz_C
            # Operate onto B first
            inv_L_B = sparse(L\(P*(R\matrices.B))) # inv(L)*B
            if opts.clear_memory
                L = nothing; P = nothing; R = nothing; B = nothing
                GC.gc()
            end
            if nnz_C < nnz(inv_L_B)
                S = Matrix(((matrices.C*Q)/U)*inv_L_B) # [C*inv(U)]*[inv(L)*B]
            else
                # This version essentially is the same as factorize_and_solve except that here we project with C after the whole X is computed
                S = (matrices.C*Q)*Matrix(U\inv_L_B)   # C*[inv(U)*[inv(L)*B]]
                if issparse(S); S = Matrix(S); end
            end
        else
            # Operate onto C first
            C_inv_U = sparse((matrices.C*Q)/U) # C*inv(U)
            if opts.clear_memory
                U = nothing; Q = nothing; C = nothing
                GC.gc()
            end
            if nnz_B < nnz(C_inv_U)
                S = Matrix(C_inv_U*(L\(P*(R\matrices.B)))) # [C*inv(U)]*[inv(L)*B]
            else
                S = Matrix(C_inv_U/L)*(P*(R\matrices.B))   # [[C*inv(U)]*inv(L)]*B
            end
        end
        t2 = time(); info.timing_solve = t2-t1;
        if opts.verbal; @printf("elapsed time: %7.3f secs\n", info.timing_solve); end
    elseif opts.method == "factorize_and_solve"
    ## Compute S=C*inv(A)*B or X=inv(A)*B by factorizing A and solving for X column by column
        # Factorize A = L*U where L and U are upper and lower triangular, with permutations
        if opts.solver == "MUMPS"
            (id, info) = MUMPS_analyze_and_factorize(matrices.A, opts, opts.is_symmetric_A)
        else
            (L, U, P, Q, R, info) = JULIA_factorize(matrices.A, opts)
            if opts.clear_memory
                matrices.A = nothing
                GC.gc()
            end
        end
        info.timing_build = 0  # the build time for A, B, C will be added in addition to this

        if opts.analysis_only
            # Skip the solve stage
            S = nothing
            info.timing.solve = 0
            # Destroy the MUMPS instance and deallocate memory
            set_job!(id,-2) # what to do: terminate the instance
            finalize!(id)
            #MPI.Finalize() # We cannot finalize the MPI, since we may do another computation
        else       
            # Solve stage (forward and backward substitutions)
            if opts.verbal; @printf("Solving     ... "); end
            t1 = time()
            if return_X # Compute X=inv(A)*B; we call X as S here since S is what mesti_matrix_solver() returns
                if opts.solver == "MUMPS"
                    set_job!(id,3) # what to do: solve
                    # Note that we need to specify that the RHS is sparse first, and then provide RHS
                    set_icntl!(id,20,1;displaylevel=0) # tell MUMPS that the RHS is sparse
                    if opts.use_single_precision_MUMPS
                        # Convert the double-precision to single-precision
                        if eltype(matrices.B) == ComplexF64
                            provide_rhs!(id,convert(SparseMatrixCSC{ComplexF32, Int64}, matrices.B))
                        elseif eltype(matrices.B) == Float64
                            provide_rhs!(id,convert(SparseMatrixCSC{Float32, Int64}, matrices.B))
                        elseif eltype(matrices.B) == Complex{Int64}
                            provide_rhs!(id,convert(SparseMatrixCSC{Complex{Int32}, Int64}, matrices.B))
                        else eltype(matrices.B) == Int64
                            provide_rhs!(id,convert(SparseMatrixCSC{Int32, Int64}, matrices.B))
                        end
                    else
                        provide_rhs!(id,matrices.B) # no need to loop since we keep everything     
                    end
                    invoke_mumps!(id)  # perform the solve
                    if id.infog[1] < 0; error("$(MUMPS_error_message(id.infog))"); end # check for errors    
                    S = get_sol(id) # X = id.SOL
                    # Destroy the MUMPS instance and deallocate memory
                    finalize!(id)
                    #MPI.Finalize() # We cannot finalize the MPI, since we may do another computation
                else
                    # Forward and backward substitutions + undo scaling and ordering
                    # X = Q*U\(L\(P*(R\B)))
                    # Do it in two steps so we can clear L to reduce peak memory usage
                    inv_L_B = L\(P*(R\matrices.B)) # inv(L)*B
                    if opts.clear_memory
                        L = nothing; P = nothing; R = nothing; matrices.B = nothing
                        GC.gc()
                    end
                    S = Q*Matrix(U\inv_L_B)
                end
            else # Compute S=C*inv(A)*B
                M_in = size(matrices.B, 2)
                M_out = size(matrices.C, 1)
                S = zeros(ComplexF64, M_out, M_in)
                # Storing the whole X=inv(A)*B wastes memory, so we solve for opts.nrhs columns of X each time and only keep its projection onto C.
                if opts.solver == "MUMPS"
                    if opts.iterative_refinement
                        info.itr_ref_nsteps = zeros(M_in,1)
                        info.itr_ref_omega_1 = zeros(M_in,1)
                        info.itr_ref_omega_2 = zeros(M_in,1)
                    end
                    for k = 1:opts.nrhs:M_in
                        in_list = k:min(k+opts.nrhs-1, M_in)
                        set_job!(id,3)  # what to do: solve
                        # Note that we need to specify that the RHS is sparse first, and then provide RHS
                        set_icntl!(id,20,1;displaylevel=0) # tell MUMPS that the RHS is sparse
                        if opts.use_single_precision_MUMPS
                            # Convert the double-precision to single-precision
                            if eltype(matrices.B) == ComplexF64
                                provide_rhs!(id,convert(SparseMatrixCSC{ComplexF32, Int64}, matrices.B[:,in_list]))
                            elseif eltype(matrices.B) == Float64
                                provide_rhs!(id,convert(SparseMatrixCSC{Float32, Int64}, matrices.B[:,in_list]))
                            elseif eltype(matrices.B) == Complex{Int64}
                                provide_rhs!(id,convert(SparseMatrixCSC{Complex{Int32}, Int64}, matrices.B[:,in_list]))
                            else eltype(matrices.B) == Int64
                                provide_rhs!(id,convert(SparseMatrixCSC{Int32, Int64}, matrices.B[:,in_list]))
                            end
                        else
                            provide_rhs!(id,matrices.B[:,in_list])
                        end
                        invoke_mumps!(id)  # perform the solve
                        if id.infog[1] < 0; error(MUMPS_error_message(id.infog)); end # check for errors
                        S[:,in_list] = matrices.C*get_sol(id) # X = id.SOL
                        if opts.iterative_refinement  # we must have opts.nrhs = 1 in this case
                            info.itr_ref_nsteps[k] = id.infog[15] # number of steps of iterative refinement
                            info.itr_ref_omega_1[k] = id.rinfog[7] # scaled residual 1; see MUMPS user guide section 3.3.2
                            info.itr_ref_omega_2[k] = id.rinfog[8] # scaled residual 2; see MUMPS user guide section 3.3.2
                        end
                    end
                    finalize!(id)
                    #MPI.Finalize() # We cannot finalize the MPI, since we may do another computation
                else
                    CQ = matrices.C*Q
                    if opts.clear_memory
                        matrices.C = nothing
                        GC.gc()
                    end
                    for k = 1:opts.nrhs:M_in
                        in_list = k:min(k+opts.nrhs-1, M_in)
                        # Forward and backward substitutions + undo scaling and ordering
                        S[:,in_list] = CQ*Matrix(U\(L\(P*(R\matrices.B[:,in_list]))))
                    end
                end
            end
            if issparse(S); S = Matrix(S); end
            t2 = time(); info.timing_solve = t2-t1
            if opts.verbal; @printf("elapsed time: %7.3f secs\n", info.timing_solve); end
        end
    else
        throw(ArgumentError("opts.method = \"$(opts.method)\" is not a supported option."))
    end
    
    # Convert info.ordering_method from integer to character array, per MUMPS definition of icntl(7)
    if opts.solver == "MUMPS"
        if info.ordering_method == 0
            str_ordering = "AMD"
        elseif info.ordering_method == 1
            str_ordering = "user specified"    
        elseif info.ordering_method == 2
            str_ordering = "AMF"             
        elseif info.ordering_method == 3
            str_ordering = "SCOTCH"
        elseif info.ordering_method == 4
            str_ordering = "PORD"
        elseif info.ordering_method == 5
            str_ordering = "METIS"                
        elseif info.ordering_method == 6
            str_ordering = "QAMD"                
        else                
            str_ordering = "N/A"
        end
        # Check if METIS ordering was actually used in MUMPS
        if opts.use_METIS && info.ordering_method != 5
            @warn("opts.use_METIS = true, but $(str_ordering) ordering was actually used in MUMPS.")
        end
        info.ordering_method = str_ordering        
    end

    ## Need to thinik more how to clear memory in Julia                
    #if opts.use_given_ordering; opts = rmfield(opts, 'ordering'); end # We don't return the user-specified ordering again since it can be large
    info.opts = opts # Return the parameters used for user's reference
    info.timing_init = timing_init # Initialization time
    t2 = time(); info.timing_total = t2-t0 # Total computing time
    
    return (S, info)
end

"""
    JULIA_FACTORIZE calls JULIA's lu() to factorize matrix A
"""
function JULIA_factorize(A::Union{SparseMatrixCSC{Int64, Int64},SparseMatrixCSC{Complex{Int64}, Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64, Int64}}, opts::Opts)

    if opts.verbal; @printf("Factorizing ... "); end
    t1 = time()

    # P*inv(R)*A*Q = L*U where R is diagonal, L and U are lower and upper triangular, and P and Q are permutation matrices    
    F = lu(A)
    L = F.L; U=F.U; Q = F.q; P = F.p; R = F.Rs
    if opts.clear_memory; F = nothing; end    
    Q = hcat(Q,(1:length(Q)))
    P = hcat(P,(1:length(P)))
    Q = sparse(Q[:,1],Q[:,2],1)
    P = sparse(P[:,2],P[:,1],1)
    R = spdiagm(size(R,1),size(R,1), 1 ./R)
    
    info = Info()
    t2 = time()
    info.timing_factorize = t2-t1
    info.timing_analyze = 0 # the analysis time is already counted in the factorization time 
    if opts.verbal; @printf("elapsed time: %7.3f secs\n", info.timing_factorize); end
    return (ComplexF64.(L), ComplexF64.(U), ComplexF64.(P), ComplexF64.(Q), ComplexF64.(R), info)
end

"""
    MUMPS_ANALYZE_AND_FACTORIZE calls MUMPS to analyze and factorize matrix A (if ind_schur is not given) or to compute its Schur complement (if ind_schur is given)
"""
function MUMPS_analyze_and_factorize(A::Union{SparseMatrixCSC{Int64, Int64},SparseMatrixCSC{Complex{Int64}, Int64},SparseMatrixCSC{Float64, Int64},SparseMatrixCSC{ComplexF64, Int64}}, opts::Opts, is_symmetric::Bool, ind_schur::Union{UnitRange{Int64},Nothing} = nothing, par = 1::Int64)

    ## Initialize MUMPS
    N = size(A,1)
    if is_symmetric
        sym = 2  # specify that matrix A is symmetric but not positive definite
    else
        sym = 0  # specify that matrix A is not symmetric
    end
    
    if opts.use_single_precision_MUMPS
        # Convert the element data type of A to the single-precision     
        if eltype(A) == ComplexF64
            A = convert(SparseMatrixCSC{ComplexF32, Int64}, A)
        elseif eltype(A) == Float64
            A = convert(SparseMatrixCSC{Float32, Int64}, A)
        elseif eltype(A) == Complex{Int64}
            A = convert(SparseMatrixCSC{Complex{Int32}, Int64}, A)
        else eltype(A) == Int64
            A = convert(SparseMatrixCSC{Int32, Int64}, A)
        end
    end
        
    MPI.Initialized() ? nothing : MPI.Init()
    id = Mumps(A, sym=sym, par=par) # get the default parameters    

    set_icntl!(id,4,0;displaylevel=0); # Turn off diagnostic messages 
    if opts.verbal_solver
        # Output to standard output stream, which is labeled by 6 in fortran
        # Note that the output behavior depends on the compiler used to compile MUMPS:
        # - ifortran will let output go to the standard output
        # - other compilers like gfortran will sometimes write output to a file fort.6, sometimes give nothing
        set_icntl!(id,3,6;displaylevel=0)
    else
        set_icntl!(id,3,0;displaylevel=0)  # turn off output
    end

    # Set the number of OpenMP threads, if given (only avaiable after MUMPS 5.2.0)
    if isdefined(opts, :nthreads_OMP)
        set_icntl!(id,16,opts.nthreads_OMP;displaylevel=0)
    end 
    
    # Activate the BLR, if given
    if opts.use_BLR
        set_icntl!(id,35,2;displaylevel=0)
        set_cntl!(id,7,opts.threshold_BLR;displaylevel=0)
        if isdefined(opts, :icntl_36)
            set_icntl!(id,36,opts.icntl_36;displaylevel=0)
        end
        if isdefined(opts, :icntl_38)
            set_icntl!(id,38,opts.icntl_38;displaylevel=0)            
        end        
    end
    
    # Specify where the Schur block is
    # We should allow ind_schur to be an empty vector (for which the Schur complement is an empty matrix).
    if ~isa(ind_schur, Nothing)
        if ~(all(x-> (isa(x, Integer) &&  x > 0 && x <= size(A,1)), ind_schur) && size(ind_schur,2) == 1)
            throw(ArgumentError("ind_schur must be a row vector of positive integers not exceeding size(A,1) = $(N)."))
        end
        # This line should be changed in the future to use the distributed by columns.
        set_schur_centralized_by_column!(id, ind_schur)

        # Discard factors, since all we need is the Schur complement
        set_icntl!(id,31,1;displaylevel=0)
    elseif opts.iterative_refinement
        set_cntl!(id,10,1000;displaylevel=0)# Enable iterative refinement and set maximum number of iterations; usually 1-2 iterations is sufficient
        # Lower the stopping criterion (omega_1 + omega_2 < id.cntl[2]) for iterative refinement to machine precision.
        # Note that iterative refinement will also stop when omega_1 + omega_2 does not decrease by at least a factor of 5.
        #@printf("Changing stopping criterion of iterative refinement from id.cntl[2] = %g to %g\n", id.cntl[2], 1e-16);        
        set_cntl!(id,2,1e-16;displaylevel=0)      
    end

    ## Analysis stage
    if opts.verbal; @printf("Analyzing   ... "); end
    t1 = time()
    if opts.parallel_dependency_graph
        # Split dependency graph and processed independently by OpenMP threads.
        # This typically improve the time performance, but marginally increase the memory usage in full multithread.
        set_keep!(id,401,1)
    end
    set_job!(id,1) # what to do: analysis  
    if opts.use_given_ordering
        # Use a user-specified ordering, if given
        # We don't check the validity of opts.ordering other than its size because it can be slow to check a long vector
        if length(opts.ordering) != N
            throw(ArgumentError("length(opts.ordering) = $(length(opts.ordering)) does not equal matrix size = $(N)."))
        end
        set_icntl!(id,7,1) # use the ordering in id.PERM_IN
        id.perm_in = opts.ordering # specify the ordering (a permutation vector)
    else
        if opts.use_METIS
            # Use the METIS package to compute ordering.
            # This typically gives a more sparse factorization, so the memory usage is lower and the factorization and solve stages will be faster. But the analyze stage will be much slower.
            set_icntl!(id,7,5)
        else
            # use AMD to computer ordering by default
            set_icntl!(id,7,0)
        end
    end
    invoke_mumps!(id) # run the analysis
    if id.infog[1] < 0; error("$(MUMPS_error_message(id.infog))"); end # check for errors

    info = Info()
    # Store information on the ordering used
    info.ordering_method = id.infog[7] # the ordering methosd that was actually used
    if opts.store_ordering
        # Store the ordering (a permutation vector) that was computed
        # The following line should be checked or we should define what is the ordering in the interface.
        info.ordering = zeros(1,N)
        info.ordering[id.SYM_PERM] = 1:N
    end

    t2 = time(); info.timing_analyze = t2-t1
    if opts.verbal; @printf("elapsed time: %7.3f secs\nFactorizing ... ", info.timing_analyze); end   

    ## Factorize stage
    t1 = time()
    set_job!(id,2) # what to do: factorize
    invoke_mumps!(id) # run the factorization
    
    if id.infog[1] < 0 error("$(MUMPS_error_message(id.infog))") end # check for errors

    t2 = time(); info.timing_factorize = t2-t1
    if opts.verbal; @printf("elapsed time: %7.3f secs\n", info.timing_factorize); end
    
    return id, info
end


"""
    MUMPS_ERROR_MESSAGE interprets some of the error messages from MUMPS
"""
function MUMPS_error_message(infog)

    @printf("\n")
    for nn = 1:length(infog)
        @printf("infog[%d] = %d\n", nn, infog[nn])
    end

    # Interpret some of the error values; look at MUMPS user guide for complete listing
    if infog[1] == -1
        err_msg = @sprintf("An error occurred on processor %d", infog[2])        
    elseif infog[1] == -2
        err_msg = @sprintf("NNZ (or NZ) = %d is out of range", infog[2])
    elseif infog[1] == -3
        err_msg = "MUMPS was called with an invalid value for JOB"        
    elseif infog[1] == -4
        err_msg = @sprintf("Error in user-provided permutation array PERM_IN at position %d", infog[2])        
    elseif infog[1] == -5
        err_msg = @sprintf("Not enough memory for real workspace allocation during analysis; infog(2) = %d", infog[2])
    elseif infog[1] == -6
        err_msg = @sprintf("Matrix is singular in structure; structural rank %d", infog[2])        
    elseif infog[1] == -7
        err_msg = @sprintf("Not enough memory for integer workspace allocation during analysis; infog[2] = %d", infog[2])        
    elseif infog[1] == -8
        err_msg = "Integer workarray IS too small for factorization; should increase ICNTL(14) and call again"
    elseif infog[1] == -9
        err_msg = "Real/complex workarray S too small for factorization; should increase ICNTL(14) and call again"
    elseif infog[1] == -10
        err_msg = "Matrix is numerically singular"       
    elseif infog[1] == -13
        err_msg = @sprintf("Not enough memory for real/complex workspace allocation during factorization; infog[2] = %d; estimated memory usage = %d MB; actual memory allocation = %d MB", infog[2], infog[17], infog[19])        
    elseif infog[1] == -14
        err_msg = "Integer workarray IS too small for solution; should increase icntl[14] and call again"
    elseif infog[1] == -15
        err_msg = "Integer workarray IS too small for iterative refinement and/or error analysis; should increase icntl(14) and call again"       
    elseif infog[1] == -16
         err_msg = @sprintf("N = %d is out of range", infog[2])        
    elseif infog[1] == -17
        err_msg = "Internal send buffer too small; should increase icntl[14] and call again"
    elseif infog[1] == -20
        err_msg = "Internal reception buffer too small; should increase ICNTL[14] and call again"      
    elseif infog[1] == -44
        err_msg = @sprintf("The solve stage cannot be performed because the factors or part of the factors are not available; icntl[31] = %d", infog[2])        
    else        
        msg = @sprintf("MUMPS failed with infog[1]=%d, infog[2]=%d: %s.", infog[1], infog[2], err_msg)
    end

    return msg
end
