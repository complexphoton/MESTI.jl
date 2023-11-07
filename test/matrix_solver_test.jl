# This code test the function mesti_matrix_solver!() to solves X in the linear systems AX = B.

# Define the size of matrix and the number of RHS
N, M = 2000, 20

# Test the functionality in a test set
@testset "matrix_solver:          " begin
    for i ∈ 1:100
	matrices = Matrices()
	matrices.A = sparse(I,N,N) + sprand(N,N,1/N) # Generate a sparse matrix A
        matrices.B = sprand(N,M,1/sqrt(N*M)) # Generate a sparse random matrix for a matrix B
	
	opts = Opts()
        opts.verbal = false # Not print system information and timing to the standard output.
	
        (X, _) = mesti_matrix_solver!(matrices, opts)  # Solve the linear system using matrix solver
        @test norm(matrices.A*X-matrices.B) ≤ 1e-5 # Test the correctness of the solution
    end
end

