export display_icntl, display_cntl, display_keep

function Base.show(io::IO,mumps::Mumps{TC,TR}) where {TC,TR}
    print(io,"Mumps{$TC,$TR}: ")
    if TC<:Float32
        println(io,"single precision real")
        lib = "smumps"
    elseif TC<:Float64
        println(io,"double precision real")
        lib = "dmumps"
    elseif TC<:ComplexF32
        println(io,"single precision complex")
        lib = "cmumps"
    elseif TC<:ComplexF64
        println(io,"double precision complex")
        lib = "zmumps"
    end
    println(io,"lib: ", lib)
    print(io,"job: ", mumps.job, " ")
    if mumps.job==-3
        println(io,"save/restore")
    elseif mumps.job==-2
        println(io,"terminate")
    elseif mumps.job==-1
        println(io,"initialize")
    elseif mumps.job==1
        println(io,"analyze")
    elseif mumps.job==2
        println(io,"factorize")
    elseif mumps.job==3
        println(io,"solve")
    elseif mumps.job==4
        println(io,"analyze + factorize")
    elseif mumps.job==5
        println(io,"factorize + solve")
    elseif mumps.job==6
        println(io,"analyze + factorize + solve")
    elseif mumps.job==7
        println(io,"save")
    elseif mumps.job==8
        println(io,"restore")
    else
        println(io,"unrecognized")
    end
    print(io,"sym: ", mumps.sym)
    if mumps.sym==1
        println(io," symmetric pos def")
    elseif mumps.sym==2
        println(io," symmetric")
    else
        println(io," unsymmetric")
    end
    print(io,"par: ", mumps.par)
    mumps.par==0 ? println(io," host not among workers") : println(io," host among workers ")
    print(io,"matrix A: ")
    if has_matrix(mumps)
        print(io,"$(mumps.n)×$(mumps.n) ")
        if is_matrix_assembled(mumps)
            println(io,"sparse matrix, with $(mumps.nnz) nonzero elements")
        else
            print(io,"elemental matrix with $(mumps.nelt) element")
            mumps.nelt>1 ? println(io,"s") : println()
        end
    else
        println(io,"uninitialized")
    end

    print(io,"rhs B:")
    rhs_type = is_rhs_dense(mumps) ? "dense" : "sparse"
    nz_rhs = is_rhs_dense(mumps) ? "" : string(",with ",mumps.nz_rhs," nonzero elements")
    if has_rhs(mumps)
        lrhs = is_rhs_dense(mumps) ?  mumps.lrhs : mumps.n
        nrhs = mumps.nrhs
        println(io," $lrhs×$nrhs ",rhs_type," matrix", nz_rhs)
    else
        println(io," uninitialized")
    end

    println(io,"ICNTL settings summary: ")
    icntl_inds = [4,9,13,19,30,33]
    for i ∈ eachindex(icntl_inds)
        print(io,"\t")
        display_icntl(io,mumps.icntl,icntl_inds[i],mumps.icntl[icntl_inds[i]])
    end
end


"""
    display_icntl(mumps)

Show the complete ICNTL integer array of `mumps`, with descriptions

See also: [`set_icntl!`](@ref)
"""
display_icntl(io::IO,mumps::Mumps) = display_icntl(io,mumps.mumps.icntl)
function display_icntl(io::IO,icntl)
    for i ∈ eachindex(icntl)
        display_icntl(io,icntl,i,icntl[i])
    end
end
function display_icntl(io::IO,icntl,i,val)
    automatic = "decided by software"
    print(io,"$i,\t$val\t")
    if i==1
        print(io,"output stream for error messages: ")
        if val≤0
            print(io,"suppressed")
        else
            print(io,"$val")
        end
    elseif i==2
        print(io,"output stream for diagnostics, statistics, warnings: ")
        if val≤0
            print(io,"suppressed")
        else
            print(io,"$val")
        end
    elseif i==3
        print(io,"output stream for global info: ")
        if val≤0
            print(io,"suppressed")
        else
            print(io,"$val")
        end
    elseif i==4
        print(io,"level of printing: ")
        if val≤0
            print(io,"error, warnings, diagnostics suppressed")
        elseif val==1
            print(io,"only error messages")
        elseif val==2
            print(io,"errors, warnings, main statistics")
        elseif val==3
            print(io,"errors, warnings, terse diagnostics")
        else
            print(io,"errors, warnings, info on input, output parameters")
        end
    elseif i==5
        print(io,"matrix input format: ")
        if val==1
            print(io,"elemental")
        else
            print(io,"assembled")
        end
    elseif i==6
        print(io,"permutation and/or scaling: ")
        if val==1
            print(io,"number of diagonal nonzeros is maximized")
        elseif val ∈ [2,3]
            print(io,"smallest diagonal value is maximized")
        elseif val==4
            print(io,"sum of diagonal values is maximized")
        elseif val ∈ [5,6]
            print(io,"product of diagonal values is maximized")
        elseif val==7
            print(automatic)
        else
            print(io,"none")
        end
    elseif i==7
        print(io,"reordering for analysis: ")
        if val==0
            print(io,"Approximate Minimum Degree (AMD)")
        elseif val==1
            print(io,"given by user via PERM_IN (see provide_perm_in)")
        elseif val==2
            print(io,"Approximate Minimum Fill (AMF)")
        elseif val==3
            print(io,"SCOTCH, if installed, else ",automatic)
        elseif val==4
            print(io,"PORD, if installed, else ",automatic)
        elseif val==5
            print(io,"METIS, if installed, else ",automatic)
        elseif val==6
            print(io,"Approximate Minimum Degree with quasi-dense row detection (QAMD)")
        else
            print(automatic)
        end
    elseif i==8
        print(io,"scaling strategy: ")
        if val==-2
            print(io,"computed during analysis")
        elseif val==-1
            print(io,"provided by user in COLSCA and ROWSCA (see provide_)")
        elseif val==1
            print(io,"diagonal, computed during factorization")
        elseif val==3
            print(io,"comlumn, computed during factorization")
        elseif val==4
            print(io,"row and column based on inf-norms, computed during factorization")
        elseif val==7
            print(io,"row and column iterative, computed during factorization")
        elseif val==8
            print(io,"row and column iterative, computed during factorization")
        elseif val==77
            print(automatic," during analysis")
        else
            print(io,"none")
        end
    elseif i==9
        print(io,"transposed: ")
        if val==1
            print(io,"false")
        else
            print(io,"true")
        end
    elseif i==10
        print(io,"iterative refinement: ")
        if val<0
            print(io,"fixed number of iterations")
        elseif val>0
            print(io,"until convergence with max number of iterations")
        else
            print(io,"none")
        end
    elseif i==11
        print(io,"statistics of error analysis: ")
        if val==1
            print(io,"all (including expensive ones)")
        elseif val==2
            print(io,"main (avoid expensive ones)")
        else
            print(io,"none")
        end
    elseif i==12
        print(io,"ordering for symmetric matrices: ")
        if val==0
            print(io,automatic)
        elseif val==2
            print(io,"on compressed graph")
        elseif val==3
            print(io,"constrained ordering, only used with AMF (see icntl 7)")
        else
            print(io,"usual ordering, nothing done")
        end
    elseif i==13
        print(io,"parallelism of root node: ")
        if val==-1
            print(io,"force splitting")
        elseif val>0
            print(io,"sequential factorization (ScaLAPACK not used) unless num workers > $val")
        else
            print(io,"parallel factorization")
        end
    elseif i==14
        print(io,"percentage increase in estimated working space: $val%")
    elseif i==18
        print(io,"distributed input matrix: ")
        if val==1
            print(io,"structure provided centralized, mumps returns mapping, user provides entries to mapping")
        elseif val==2
            print(io,"structure provided centralized at analysis, entries provided to all workers at factorization")
        elseif val==3
            print(io,"distributed matrix, pattern, entries provided")
        else
            print(io,"centralized")
        end
    elseif i==19
        print(io,"Schur complement: ")
        if val==1
            print(io,"true, Schur complement returned centralized by rows")
        elseif val ∈ [2,3]
            print(io,"true, Schur complement returned distributed by columns")
        else
            print(io,"false, complete factorization")
        end
    elseif i==20
        print(io,"rhs: ")
        0<val<4 ? print(io,"sparse, ") : nothing
        if val==1
            print(io,"sparsity-exploting acceleration of solution ", automatic)
        elseif val==2
            print(io,"sparsity not exploited in solution")
        elseif val==3
            print(io,"sparsity exploited to accelerate solution")
        else
            print(io,"dense")
        end
    elseif i==21
        print(io,"distribution of solution vectors: ")
        if val==1
            print(io,"distributed")
        else
            print(io,"assembled and stored in centralized RHS")
        end
    elseif i==22
        print(io,"out-of-core (OOC) factorization and solve: ")
        if val==0
            print(io,"false")
        elseif val==1
            print(io,"true")
        else
            @warn "not sure this is a valid setting"
        end
    elseif i==23
        print(io,"max size (in MB) of working memory per worker: ")
        if val>0
            print(io,"$val MB")
        else
            print(automatic)
        end
    elseif i==24
        print(io,"null pivot row detection: ")
        if val==1
            print(io,"true")
        else
            print(io,"false. if null pivot present, will result in error INFO(1)=-10")
        end
    elseif i==25
        print(io,"defecient matrix and null space basis: ")
        if val==-1
            print(io,"complete null space basis computed")
        elseif val==0
            print(io,"normal solution phase. if matrix found singular, one possible solution returned")
        else
            print(io,"$val-th vector of null space basis computed")
        end
    elseif i==26
        print(io,"if Schur, solution phase: ")
        if val==1
            print(io,"condense/reduce rhs on Schur")
        elseif val==2
            print(io,"expand Schur solution on complete solution variables")
        else
            print(io,"standard solution")
        end
    elseif i==27
        print(io,"blocking size for multiple rhs: ")
        if val<0
            print(automatic)
        elseif val==0
            print(io,"no blocking, same as 1")
        else
            print(io,"blocksize=min(NRHS,$val)")
        end
    elseif i==28
        print(io,"ordering computation: ")
        if val==1
            print(io,"sequential")
        elseif val==2
            print(io,"parallel")
        else
            print(automatic)
        end
    elseif i==29
        print(io,"parallel ordering tool: ")
        if val==1
            print(io,"PT-SCOTCH, if available")
        elseif val==2
            print(io,"PARMETIS, if available")
        else
            print(automatic)
        end
    elseif i==30
        print(io,"compute entries of A⁻¹: ")
        if val==1
            print(io,"true")
        else
            print(io,"false")
        end
    elseif i==31
        print(io,"discarded factors: ")
        if val==1
            print(io,"all")
        elseif val==2
            print(io,"U, for unsymmetric")
        else
            print(io,"none, except for ooc factorization of unsymmetric")
        end
    elseif i==32
        print(io,"forward elimination of rhs: ")
        if val==1
            print(io,"performed during factorization")
        else
            print(io,"not performed during factorization (standard)")
        end
    elseif i==33
        print(io,"compute determinant: ")
        if val==0
            print(io,"false")
        else
            print(io,"true")
        end
    elseif i==34
        print(io,"OOC file conservation: ")
        if val==1
            print(io,"not marked for deletion")
        else
            print(io,"marked for deletion")
        end
    elseif i==35
        print(io,"Block Low-Rank: ")
        if val==1
            print(io,"activated, options ", automatic)
        elseif val==2
            print(io,"activated during factorization and solution phases")
        elseif val==3
            print(io,"activated during factorzation only")
        else
            print(io,"not activated")
        end
    elseif i==36
        print(io,"BLR variant: ")
        if val==1
            print(io, "UCFS with low-rank updates accumulation; compression is performed earlier")
        else
            print(io, "Standard UFSC")
        end
    elseif i==38
        print(io,"Estimated compression rate of LU factors in ppt: $val")
    else
        print(io,"not used")
    end
    print(io,"\n")
    return nothing
end


"""
    display_cntl(mumps)

Show the complete CNTL real array of `mumps`, with descriptions

See also: [`set_cntl!`](@ref)
"""
display_cntl(io::IO,mumps::Mumps) = display_cntl(io,mumps.mumps.cntl)
function display_cntl(io::IO,cntl)
    for i ∈ eachindex(cntl)
        display_icntl(io,cntl,i,cntl[i])
    end
end
function display_cntl(io::IO,cntl,i,val)
    print(io,"$i,\t$val\t")
    if i==1
        print(io,"relative threshold for numerical pivoting")
    elseif i==2
        print(io,"stopping criterion for iterative refinement")
    elseif i==3
        print(io,"null pivot?")
    elseif i==4
        print(io,"threshold for state pivoting")
    elseif i==5
        print(io,"fixation for null pivots")
    elseif i==7
        print(io,"precision of dropping parameter in BLR compression")
    else
        print(io,"not used")
    end
    print(io,"\n")
    return nothing
end


"""
    display_keep(mumps)

Show the KEEP array of `mumps`

See also: [`set_keep!`](@ref)
"""
display_keep(io::IO,mumps::Mumps) = display_keep(io,mumps.mumps.keep)
function display_keep(io::IO,keep)
    for i ∈ eachindex(keep)
        display_keep(io,keep,i,icntl[i])
    end
end
function display_keep(io::IO,keep,i,val)
    if i == 401
        print(io,"keep(401): ")
        print(io,"$val")
    end
end