# this file contains the main high-level functionality of MUMPS3:
# functions to solve A*x=y, LU factorization of A, compute the determinant
# of A, get a Schur complement matrix associated with A, compute the inverse
# of said Schur complemet matrix, compute select elements of the
# inverse of A
#
# Any function which alters the JOB state lives here.
const DEFAULT_FORTRAN_COMMUNICATOR = -987654

export mumps_solve!, mumps_solve,
mumps_factorize!, mumps_factorize,
mumps_det!, mumps_det,
mumps_schur_complement!, mumps_schur_complement,
mumps_select_inv!, mumps_select_inv,
initialize!, finalize!


"""
    Mumps{T}(; sym=0, par=1) -> mumps
    Mumps(A; sym=0, par=1) -> mumps
    Mumps(A, rhs; sym=0, par=1) -> mumps

Create an instance of a `Mumps` object with the same type as matrix `A`
and simultaneously provide `A` to it.

Optional arguments `sym` and `par` correspond to the MUMPS 5.1.2 paremeters.
They are:
    `sym` = 0 for unsymmetric, 1 for symm pos def, 2 for general symm
    `par` = 0 for host participating in parallel calculation, 1 otherwise

If not arguments are passed, create an initialized but empty instance of `Mumps`
"""
Mumps{T}(;sym=0,par=1) where T = Mumps{T}(sym,par,DEFAULT_FORTRAN_COMMUNICATOR)
function Mumps(A::AbstractArray{T}; kwargs...) where T
    if !haskey(kwargs,:sym)
        if issymmetric(A)
            kwargs = (kwargs...,:sym => 2)
        else
            kwargs = (kwargs...,:sym => 0)
        end
    end
    mumps = Mumps{T}(;kwargs...)
    typeof(A)<:Array ? set_icntl!(mumps,5,1; displaylevel=0) : nothing
    provide_matrix!(mumps,A)
    return mumps
end
function Mumps(A::AbstractArray{T}, rhs::AbstractArray{TA}; kwargs...) where {T,TA}
    if !haskey(kwargs,:sym)
        if issymmetric(A)
            kwargs = (kwargs...,:sym => 2)
        else
            kwargs = (kwargs...,:sym => 0)
        end
    end
    mumps = Mumps{promote_type(T,TA)}(;kwargs...)
    typeof(A)<:Array ? set_icntl!(mumps,5,1) : nothing
    provide_matrix!(mumps,A)
    if typeof(rhs)<:AbstractSparseArray
        set_icntl!(mumps,20,1; displaylevel=0)
    elseif typeof(rhs)<:Array
        set_icntl!(mumps,20,0; displaylevel=0)
    else
        throw(ArgumentError("unrecognized array type for rhs"))
    end
    provide_rhs!(mumps,rhs)
    return mumps
end


"""
    mumps_solve!(x,A,y; kwargs...)
    mumps_solve!(x,mumps)
    mumps_solve!(x,mumps,y)

Solve `A*x=y`, saving result in pre-allocated x.
If `mumps` is given, must have previously been provided a matrix `A`.
If `y` is not given, `mumps` must have previously been provided `y`

See also: [`mumps_solve`](@ref), [`get_sol!`](@ref), [`get_sol`](@ref)
"""
function mumps_solve!(mumps::Mumps)
    @assert has_matrix(mumps) "matrix not yet provided to mumps object"
    @assert has_rhs(mumps) "rhs not yet provided to mumps object"
    if mumps.job ∈ [2,4] # if already factored, just solve
        mumps.job = 3
    elseif mumps.job ∈ [1] # if analyzed only, factorize and solve
        mumps.job=5
    elseif mumps.job ∈ [3,5,6] # is solved already, retrieve solution
        return nothing
    else # else analyze, factor, solve
        mumps.job=6
    end
    invoke_mumps!(mumps)
    return nothing
end
function mumps_solve!(x::AbstractArray,mumps::Mumps)
    mumps_solve!(mumps)
    get_sol!(x,mumps)
end
function mumps_solve!(x::AbstractArray,A::AbstractArray,rhs::AbstractArray; kwargs...)
    mumps = Mumps(A,rhs; kwargs...)
    suppress_display!(mumps)
    set_icntl!(mumps,24,1) # null pivot -- check whether this is necessary/a good idea
    mumps_solve!(x,mumps)
    finalize!(mumps)
end
function mumps_solve!(x::AbstractArray,mumps::Mumps,rhs::AbstractArray)
    provide_rhs!(mumps,rhs)
    if mumps.job ∈ [2,4] # if already factored, just solve
    elseif mumps.job ∈ [1] # if analyzed only, factorize and solve
        mumps.job=5
    elseif mumps.job ∈ [3,5,6] # is solved already, reset to solve only
        mumps.job=2
    else # else analyze, factor, solve
        mumps.job=6
    end
    mumps_solve!(x,mumps)
end
"""
    mumps_solve(A,y) -> x
    mumps_solve(mumps,y) -> x
    mumps_solve(mumps) -> x

Solve `A*x=y`
If `mumps` is given, must have previously been provided a matrix `A`.
If only input is `mumps` must also have been provided `y`.

See also: [`mumps_solve!`](@ref)
"""
function mumps_solve(mumps::Mumps)
    if mumps.job ∈ [3,5,6]
        x = get_sol(mumps)
    else
        x = get_rhs(mumps)
    end
    x = convert(Matrix,x)
    mumps_solve!(x,mumps)
    return x
end
function mumps_solve(A::AbstractArray,rhs::AbstractArray{T,N}; kwargs...) where {T,N}
    x = fill(convert(promote_type(eltype(A),eltype(rhs)),NaN),size(rhs)...)
    mumps_solve!(x,A,rhs; kwargs...)
    return x
end
function mumps_solve(mumps::Mumps,rhs::AbstractArray)
    x = copy(convert(Matrix,rhs))
    mumps_solve!(x,mumps,rhs)
    return x
end


"""
    mumps_factorize!(mumps)

LU factorize `A` previously provided to `mump`.
LU stored in `mumps`, but not in a particularly accessible way.
Useful for doing repeated solves downstream.

See also: [`mumps_factorize`](@ref)
"""
function mumps_factorize!(mumps::Mumps)
    @assert has_matrix(mumps) "matrix not yet provided to mumps object"
    if mumps.job ∈ [2,3,4,5,6] # already factored
        @warn "already factored"
        return nothing
    elseif mumps.job ∈ [1] # if analyzed only, factorize
        mumps.job=2
    else # else analyze, factor
        mumps.job=4
    end
    invoke_mumps!(mumps)
end
"""
    mumps_factorize(A) -> mumps

LU factorize `A`.
LU stored in `mumps`, but not in a particularly accessible way.
Useful for doing repeated solves downstream.

See also: [`mumps_factorize!`](@ref)
"""
function mumps_factorize(A::AbstractArray)
    mumps = Mumps(A)
    suppress_display!(mumps)
    mumps_factorize!(mumps)
    return mumps
end


"""
    mumps_det!(mumps; discard=true)

Compute determinant of `A`, which has been previously provided to
`mumps`.

Determinant can be computed from mutated `mumps` by just `det(mumps)`
[must have loaded LinearAlgebra].

Optional keyward `discard` controls whether LU factors are discarded via
ICNTL[31]. This is useful if you only care about the determinant and don't
want to do any further computation with `mumps`. Use `discard=2` to throw
away only L. See manual for [MUMPS 5.1.2](http://mumps.enseeiht.fr/doc/userguide_5.2.0.pdf#page=42).

See also: [`mumps_det`](@ref)
"""
function mumps_det!(mumps::Mumps; discard=true)
    @assert has_matrix(mumps) "matrix not yet provided to mumps object"
    if has_det(mumps) && mumps.job>1
        return nothing
    end
    set_icntl!(mumps,31,Int(discard))
    set_icntl!(mumps,33,1)
    mumps.job>0 ? mumps.job=1 : nothing
    mumps_factorize!(mumps)
end
"""
    mumps_det(A) -> det

Compute determinant of `A`.

See also: [`mumps_det!`](@ref)
"""
function mumps_det(A)
    mumps = Mumps(A)
    suppress_display!(mumps)
    mumps_det!(mumps)
    D = det(mumps)
    finalize!(mumps)
    return D
end


"""
    mumps_schur_complement!(mumps, schur_inds)
    mumps_schur_complement_complement!(mumps, x)

`schur_inds` is integer array of Schur indices.
If `x` is sparse, Schur indices determined from populated rows of `x`

See also: [`mumps_schur_complement`](@ref), [`get_schur_complement!`](@ref), [`get_schur_complement`](@ref)
"""
function mumps_schur_complement!(mumps::Mumps, schur_inds::AbstractArray{Int,1})
    @assert has_matrix(mumps) "matrix not yet provided to mumps object"
    set_schur_centralized_by_column!(mumps, schur_inds)
    if mumps.job ∈ [1] # if analyzed only, factorize
        mumps.job=2
    else # else analyze, factor
        mumps.job=4
    end
    invoke_mumps!(mumps)
end
mumps_schur_complement!(mumps::Mumps, x::SparseMatrixCSC) = mumps_schur_complement!(mumps,unique!(sort(x.rowval)))
mumps_schur_complement!(mumps::Mumps, x::SparseVector) = mumps_schur_complement!(mumps,x.nzind)
"""
    mumps_schur_complement(A,schur_inds) -> S
    mumps_schur_complement(A,x) -> S

`schur_inds` is integer array
`x` is sparse, populated rows are Schur indices
`S` is Schur complement matrix.

See also: [`mumps_schur_complement!`](@ref)
"""
function mumps_schur_complement(A::AbstractArray,x)
    mumps = Mumps(A)
    suppress_display!(mumps)
    set_icntl!(mumps,8,0) # turn scaling off, not used with schur anyway (suppresses warning message with schur)
    mumps_schur_complement!(mumps,x)
    S = get_schur_complement(mumps)
    finalize!(mumps)
    return S
end


"""
    mumps_select_inv!(x,mumps)
    mumps_select_inv!(x,A)

Compute selected elements of A⁻¹ with same sparsity pattern as `x`,
stored in `x`. If passed `mumps`, must have previously been provided
with matrix `A`.

See also: [`mumps_select_inv`](@ref), [`get_rhs!`](@ref), [`get_rhs`](@ref)
"""
function mumps_select_inv!(x::AbstractSparseArray,mumps::Mumps)
    @assert has_matrix(mumps) "matrix not yet provided to mumps object"
    set_icntl!(mumps,30,1)
    set_icntl!(mumps,20,3)
    provide_rhs!(mumps,x)
    if mumps.job ∈ [2,4] # if already factored, just solve
        mumps.job = 3
    elseif mumps.job ∈ [1] # if analyzed only, factorize and solve
        mumps.job=5
    else # else analyze, factor, solve
        mumps.job=6
    end
    invoke_mumps!(mumps)
    get_rhs!(x,mumps)
end
function mumps_select_inv!(x,A)
    @assert size(x)==size(A) "target and matrix must have same size"
    mumps = Mumps(A)
    suppress_display!(mumps)
    mumps_select_inv!(x,mumps)
    finalize!(mumps)
end
"""
    mumps_select_inv(A,x) -> A⁻¹
    mumps_select_inv(A,I,J) -> A⁻¹

Compute selected elements of A⁻¹.
If two arguments are passed, the second must be sparse, and its sparsity pattern
determines the entries of the inverse.
If three arguments are passed, the integer arrays `I` and `J` specify which
entries via `i(k),j(k) = I[k],J[k]`.

See also: [`mumps_select_inv!`](@ref)
"""
function mumps_select_inv(A,x::AbstractSparseArray)
    y = copy(x)
    mumps_select_inv!(y,A)
    return y
end
function mumps_select_inv(A,I::AbstractArray{Int},J::AbstractArray{Int})
    x = sparse(I,J,fill(one(eltype(A)),length(I)),size(A)...)
    mumps_select_inv!(x,A)
    return x
end


"""
    initialize!(mumps)

Reinitialize `mumps`, regardless of its current state

See also: [`finalize`](@ref)
"""
function initialize!(mumps::Mumps)
    suppress_display!(mumps)
    mumps._finalized=false
    mumps.job=-1
    invoke_mumps!(mumps)
end


"""
    finalize!(mumps)

Release the pointers contained in `mumps`

See also: [`initialize`](@ref)
"""
function finalize!(mumps::Mumps)
    suppress_display!(mumps)
    check_finalized(mumps)
    finalize_unsafe!(mumps)
end
function finalize_unsafe!(mumps::Mumps)
    mumps._finalized=true
    mumps.job=-2
    invoke_mumps_unsafe!(mumps)
end


function Base.:\(mumps::Mumps,y)
    suppress_display!(mumps)
    if mumps.job < 2
        mumps_factorize!(mumps)
    end
    x = mumps_solve(mumps,y)
    return x
end
function LinearAlgebra.ldiv!(mumps::Mumps,y)
    suppress_display!(mumps)
    if mumps.job < 2
        mumps_factorize!(mumps)
    else
        mumps.job=2
    end
    provide_rhs!(mumps,y)
    mumps_solve!(y,mumps)
end
function LinearAlgebra.ldiv!(x,mumps::Mumps,y)
    suppress_display!(mumps)
    if mumps.job < 2
        mumps_factorize!(mumps)
    end
    mumps_solve!(x,mumps,y)
end
function LinearAlgebra.inv(mumps::Mumps)
    suppress_display!(mumps)
    y = sparse(1:mumps.n,1:mumps.n,1:mumps.n,mumps.n,mumps.n)
    mumps_select_inv!(y,mumps)
    finalize!(mumps)
    return y
end
