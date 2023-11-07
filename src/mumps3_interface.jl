# this file provides the low-level interface with the MUMPS 5.1.2 library
# by controlling access to the pointers contained in the Mumps object.
# Many functions are unsafe and are marked as such.
# None of these functions change the JOB parameter.

export invoke_mumps!,
set_icntl!, set_job!,
provide_matrix!,
provide_rhs!,
get_rhs!, get_rhs,
get_sol!, get_sol,
set_schur_centralized_by_column!,
get_schur_complement!, get_schur_complement


"""
    invoke_mumps_unsafe!(mumps)

Call the appropriate mumps C-library, passing to it the Mumps object `mumps`

This is a low-level function, meaning that you have complete control over what
operations are done, based on the MUMPS manual.

Be warned, a direct call can crash Julia if `mumps` is not appropriately
initialized.

See also: [`invoke_mumps!`](@ref)
"""
@inline function invoke_mumps_unsafe!(mumps::Mumps{TC,TR}) where {TC,TR}
    global MUMPS_LIB_S = joinpath(ENV["MUMPS_PREFIX"],"libsmumps")
    global MUMPS_LIB_D = joinpath(ENV["MUMPS_PREFIX"],"libdmumps")
    global MUMPS_LIB_C = joinpath(ENV["MUMPS_PREFIX"],"libcmumps")
    global MUMPS_LIB_Z = joinpath(ENV["MUMPS_PREFIX"],"libzmumps")

    global LIB_S = dlopen(MUMPS_LIB_S)
    global LIB_D = dlopen(MUMPS_LIB_D)
    global LIB_C = dlopen(MUMPS_LIB_C)
    global LIB_Z = dlopen(MUMPS_LIB_Z)

    @assert MPI.Initialized() "must call MPI.Init() exactly once before calling mumps"
    if TC==Float32
        cfun = :smumps_c
        LIB = LIB_S
    elseif TC==Float64
        cfun = :dmumps_c
        LIB = LIB_D
    elseif TC==ComplexF32
        cfun = :cmumps_c
        LIB = LIB_C
    elseif TC==ComplexF64
        cfun = :zmumps_c
        LIB = LIB_Z
    end
    sym = dlsym(LIB,cfun)
    ccall(sym,Cvoid,(Ref{Mumps{TC,TR}},), mumps)
    return nothing
end
"""
    invoke_mumps!(mumps)

Call the appropriate mumps C-library, passing to it the Mumps object `mumps`,
but checking to make sure `mumps` has been initialized first, so that it's safe.

This is a low-level function, meaning that you have complete control over what
operations are done, based on the MUMPS manual.

Be warned, a direct call can crash Julia if `mumps` is not appropriately
initialized.

See also: [`invoke_mumps_unsafe!`](@ref)
"""
function invoke_mumps!(mumps::Mumps)
    check_finalized(mumps)
    invoke_mumps_unsafe!(mumps)
end
function check_finalized(mumps::Mumps)
    if mumps._finalized
        throw(ErrorException("Mumps object already finalized"))
    end
end


"""
    set_icntl!(mumps,i,val; [displaylevel=1])

Set the integer control parameters according to ICNTL[i]=val

See also: [`display_icntl`](@ref)
"""
function set_icntl!(mumps::Mumps,i::Int,val::Int; displaylevel=mumps.icntl[4]-1)
    icntl = mumps.icntl
    mumps.icntl = (icntl[1:i-1]...,convert(MUMPS_INT,val),icntl[i+1:end]...)
    displaylevel>0 ? display_icntl(stdout,mumps.icntl,i,val) : nothing
    return nothing
end


"""
    set_cntl!(mumps,i,val; [displaylevel=1])

Set the real/complex control parameters according to CNTL[i]=val

See also: [`display_cntl`](@ref)
"""
function set_cntl!(mumps::Mumps{TC,TR},i::Int,val::Float64; displaylevel=mumps.icntl[4]-1) where {TC,TR}
    cntl = mumps.cntl
    mumps.cntl = (cntl[1:i-1]...,convert(TR,val),cntl[i+1:end]...)
    displaylevel>0 ? display_cntl(stdout,mumps.cntl,i,val) : nothing
    return nothing
end


"""
    set_job!(mumps,job)

Set the phase to `job`. See MUMPS manual for options.
"""
function set_job!(mumps::Mumps,i)
    mumps.job=i
    return nothing
end


"""
    set_save_dir!(mumps,dir)

set name of directory in which to store out-of-core files.
"""
function set_save_dir!(mumps,dir::String)
    @assert length(dir≤255) "directory name has $(length(dir)) characters, must be ≤255"
    i = length(dir+1)
    save_dir = mumps.save_dir
    mumps.save_dir = (dir...,'\0',save_dir[i+2:end]...)
    return nothing
end


"""
    set_save_prefix!(mumps,dir)

prefix for out-of-core files.
"""
function set_save_prefix!(mumps,prefix::String)
    @assert length(prefix) "prefix name has $(length(prefix)) characters, must be ≤255"
    i = length(prefix+1)
    save_prefix = mumps.save_prefix
    mumps.save_prefix = (prefix...,'\0',save_prefix[i+2:end]...)
    return nothing
end


"""
    provide_matrix!(mumps,A)

Provide a square matrix `A` to a `mumps` object. It internally converts
`A` to be consistent with the ICNTL[5] setting.

If needed, it tries to convert element type of `A` to be consistent with
type of `mumps`, throwing a warning in this case.

See also: [`provide_rhs!`](@ref)
"""
function provide_matrix!(mumps::Mumps{T},A::AbstractArray{TA}) where {T,TA}
    @assert size(A,1)==size(A,2) "matrix must be square, but it is $(size(A,1))×$(size(A,2))"
    if !(T==TA)
        @warn "matrix with element type $TA will attempt be converted to Mumps type $T"
    end
    if is_matrix_assembled(mumps)
        typeof(A)<:SparseMatrixCSC ? nothing : (@warn "matrix is dense, but ICNTL[5]=$(mumps.icntl[5]) indicates assembled. attempting to convert matrix to sparse")
        _provide_matrix_assembled!(mumps,convert(SparseMatrixCSC,A))
    else
        !(typeof(A)<:SparseMatrixCSC) ? nothing : (@warn "matrix is sparse, but ICNTL[5]=$(mumps.icntl[5]) indicates elemental. attempting to convert matrix to dense")
        _provide_matrix_elemental!(mumps,convert(Matrix,A))
    end
end
function _provide_matrix_assembled!(mumps::Mumps,A::SparseMatrixCSC{T}) where T
    if is_matrix_distributed(mumps)
        _provide_matrix_assembled_distributed!(mumps,A)
    else
        _provide_matrix_assembled_centralized!(mumps,A)
    end
end
function _provide_matrix_assembled_centralized!(mumps::Mumps{T},A::SparseMatrixCSC) where T
    if is_symmetric(mumps)
        I,J,V = findnz(triu(A))
    else
        I,J,V = findnz(A)
    end
    irn, jcn, a = convert.((Array{MUMPS_INT},Array{MUMPS_INT},Array{T}),(I,J,V))
    mumps.irn, mumps.jcn, mumps.a = pointer.((irn,jcn,a))
    mumps.n = A.n
    mumps.nnz = length(V)
    append!(mumps._gc_haven,[Ref(irn),Ref(jcn),Ref(a)])
    return nothing
end
function _provide_matrix_assembled_distributed!(mumps::Mumps{T},A::SparseMatrixCSC) where T
    throw(ErrorException("not written yet."))
    return nothing
end
function _provide_matrix_elemental!(mumps::Mumps{T},A::Array) where T
    mumps.n = size(A,1)
    mumps.nelt = 1
    eltptr = convert.(MUMPS_INT,[1,mumps.n+1])
    eltvar = convert.(MUMPS_INT,collect(1:mumps.n))
    mumps.eltptr, mumps.eltvar = pointer.(eltptr, eltvar)
    if is_symmetric(mumps)
        a_elt = convert.(T,[A[i,j] for i ∈ 1:mumps.n for j ∈ 1:i])
    else
        a_elt = convert.(T,A[:])
    end
    mumps.a_elt = pointer(a_elt)
    append!(mumps._gc_haven,[Ref(eltptr),Ref(eltvar),Ref(a_elt)])
    return nothing
end


"""
    provide_rhs!(mumps,y)

Provide a RHS matrix or vector `rhs` to a `mumps` object. It internally converts
`rhs` to be consistent with the ICNTL[20] setting, and additionally allocates
`mumps.rhs` according to the ICNTL[21] setting.

If needed, it tries to convert element type of `rhs` to be consistent with
type of `mumps`, throwing a warning in this case.

See also: [`provide_matrix!`](@ref)
"""
function provide_rhs!(mumps::Mumps,rhs::AbstractMatrix)
    if is_rhs_dense(mumps)
        provide_rhs_dense!(mumps,rhs)
    else
        provide_rhs_sparse!(mumps,rhs)
    end
end
function provide_rhs_sparse!(mumps::Mumps{T},rhs::AbstractMatrix) where T
    rhs = convert(SparseMatrixCSC,rhs)

    mumps.nz_rhs = length(rhs.nzval)
    mumps.nrhs = size(rhs,2)

    rhs_sparse = convert.(T,rhs.nzval)
    irhs_sparse = convert.(MUMPS_INT,rhs.rowval)
    irhs_ptr = convert.(MUMPS_INT,rhs.colptr)

    mumps.rhs_sparse = pointer(rhs_sparse)
    mumps.irhs_sparse = pointer(irhs_sparse)
    mumps.irhs_ptr = pointer(irhs_ptr)
    if is_sol_central(mumps)
        y = fill(convert(T,NaN),prod(size(rhs)))
        mumps.rhs = pointer(y)
        mumps.lrhs = size(rhs,1)
        push!(mumps._gc_haven,Ref(y))
    end
    append!(mumps._gc_haven,[Ref(rhs_sparse),Ref(irhs_sparse),Ref(irhs_ptr)])
    return nothing
end
function provide_rhs_dense!(mumps::Mumps{T},rhs::AbstractMatrix) where T
    y = convert(Matrix{T},rhs)[:]
    mumps.rhs = pointer(y)
    mumps.lrhs = size(rhs,1)
    mumps.nrhs = size(rhs,2)
    push!(mumps._gc_haven,Ref(y))
    return nothing
end
provide_rhs!(mumps::Mumps,rhs::AbstractVector) = provide_rhs!(mumps,repeat(rhs,1,1))


"""
    get_rhs!(y,mumps)

Retrieve right hand side from `mumps`, storing it in pre-allocated `x`

See also: [`get_rhs`](@ref), [`get_sol!`](@ref), [`get_sol`](@ref)
"""
function get_rhs!(x,mumps::Mumps)
    check_finalized(mumps)
    if has_rhs(mumps)
        get_rhs_unsafe!(x,mumps)
    else
        @warn "mumps has no rhs"
    end
    return nothing
end
function get_rhs_unsafe!(x::SparseMatrixCSC,mumps::Mumps)
    @assert !is_rhs_dense(mumps) "rhs is dense, target is sparse. try with dense target"
    for i ∈ LinearIndices(x.colptr)
        x.colptr[i] = unsafe_load(mumps.irhs_ptr,i)
    end
    for i ∈ LinearIndices(x.rowval)
        x.rowval[i] = unsafe_load(mumps.irhs_sparse,i)
        x.nzval[i] = unsafe_load(mumps.rhs_sparse,i)
    end
    return nothing
end
function get_rhs_unsafe!(x::Union{SubArray,Array},mumps::Mumps)
    @assert is_rhs_dense(mumps) "rhs is sparse, target is dense. try with sparse target"
    for i ∈ LinearIndices(x)
        x[i] = unsafe_load(mumps.s,i)
    end
    return nothing
end
"""
    get_rhs(mumps) -> y

Retrieve right hand side from `mumps`

See also: [`get_rhs!`](@ref), [`get_sol!`](@ref), [`get_sol`](@ref)
"""
function get_rhs(mumps::Mumps{T}) where T
    n = mumps.nrhs
    if !is_rhs_dense(mumps)
        m = mumps.n
        colptr = ones(MUMPS_INT,mumps.nrhs+1)
        rowval = ones(MUMPS_INT,mumps.nz_rhs)
        nzval = Array{T}(undef,mumps.nz_rhs)
        x = SparseMatrixCSC(m,n,colptr,rowval,nzval)
    else
        m = mumps.lrhs
        x = Array{T}(undef,m,n)
    end
    get_rhs!(x,mumps)
    return x
end


"""
    get_sol!(x,mumps)

Retrieve solution `x` from `mumps` into pre-allocated array.

See also: [`get_rhs!`](@ref), [`get_rhs`](@ref), [`get_sol`](@ref)
"""
function get_sol!(x::Union{SubArray,Array},mumps::Mumps)
    check_finalized(mumps)
    if mumps.job ∉ [3,5,6]
        @warn "mumps has not passed through a solution phase"
    end
    if has_rhs(mumps)
        get_sol_unsafe!(x,mumps)
    else
        @warn "mumps has no rhs"
    end
    return nothing
end
function get_sol_unsafe!(x::Union{SubArray,Array},mumps::Mumps)
    for i ∈ LinearIndices(x)
        x[i] = unsafe_load(mumps.rhs,i)
    end
    return nothing
end
"""
    get_sol(mumps) -> x

Retrieve solution from `mumps`

See also: [`get_rhs!`](@ref), [`get_rhs`](@ref), [`get_sol!`](@ref)
"""
function get_sol(mumps::Mumps{T}) where T
    if mumps.job ∉ [3,5,6]
        @warn "mumps has not passed through a solution phase"
    end
    x = Array{T}(undef,mumps.lrhs,mumps.nrhs)
    get_sol!(x,mumps)
    return x
end


"""
    get_schur_complement!(S,mumps)

Retrieve Schur complement matrix from `mumps` into pre-allocated `S`

See also: [`get_schur_complement`](@ref), [`mumps_schur!`](@ref), [`mumps_schur`](@ref)
"""
function get_schur_complement!(S,mumps::Mumps)
    @assert has_schur(mumps) "schur complement not yet allocated."
    get_schur_complement_unsafe!(S,mumps)
end
function get_schur_complement_unsafe!(S,mumps::Mumps)
    for i ∈ LinearIndices(S)
        S[i] = unsafe_load(mumps.schur,i)
    end
    return nothing
end
"""
    get_schur_complement(mumps) -> S

Retrieve Schur complement matrix `S` from `mumps`

See also: [`get_schur_complement!`](@ref), [`mumps_schur!`](@ref), [`mumps_schur`](@ref)
"""
function get_schur_complement(mumps::Mumps{T}) where T
    S = Array{T}(undef,mumps.size_schur,mumps.size_schur)
    get_schur_complement!(S,mumps)
    return S
end


"""
    set_schur_centralized_by_column!(mumps,schur_inds)

Set up Schur complement matrix calculation for the "centralized by column"
method suggested in the MUMPS manual

See also: [`mumps_schur!`](@ref), [`mumps_schur`](@ref)
"""
function set_schur_centralized_by_column!(mumps::Mumps{T},schur_inds::AbstractArray{Int}) where T
    mumps.size_schur = length(schur_inds)
    listvar_schur = convert.(MUMPS_INT,schur_inds)
    mumps.listvar_schur = pointer(listvar_schur)
    mumps.nprow = 1
    mumps.npcol = 1
    mumps.mblock = 100
    mumps.nblock = 100
    mumps.schur_lld = mumps.size_schur
    schur = Array{T}(undef,mumps.size_schur^2)
    mumps.schur = pointer(schur)
    set_icntl!(mumps,19,3)
    append!(mumps._gc_haven,[Ref(listvar_schur),Ref(schur)])
end


function LinearAlgebra.det(mumps::Mumps{T}) where T
    if !has_det(mumps)
        @warn "ICNTL[33]=0, determinant not computed"
        d = T<:Complex ? complex(NaN,NaN) : NaN
    else
        if T<:Complex
            d = complex(mumps.rinfog[12],mumps.rinfog[13])*2^mumps.infog[34]
        else
            d = mumps.rinfog[12]*2^mumps.infog[34]
        end
    end
    return convert(T,d)
end


####################################################################
### auxilliary functions
####################################################################
is_matrix_assembled(mumps::Mumps) = !(mumps.icntl[5] ∈ [1])

is_matrix_distributed(mumps::Mumps) = mumps.icntl[18] ∈ [1,2,3]

is_rhs_dense(mumps::Mumps) = mumps.icntl[20] ∉ [1,2,3]

is_sol_central(mumps::Mumps) = mumps.icntl[21] ∉ [1]

has_det(mumps::Mumps) = mumps.icntl[33] ∉ [0]

is_symmetric(mumps::Mumps) = mumps.sym ∈ [1,2]

is_posdef(mumps::Mumps) = mumps.sym ∈ [1]

has_matrix(mumps::Mumps) = mumps.n > 0

has_rhs(mumps::Mumps) = mumps.nrhs*mumps.lrhs>0 || mumps.nz_rhs>0

has_schur(mumps::Mumps) = mumps.size_schur > 0

LinearAlgebra.issymmetric(mumps::Mumps) = is_symmetric(mumps)
LinearAlgebra.isposdef(mumps::Mumps) = is_posdef(mumps)
