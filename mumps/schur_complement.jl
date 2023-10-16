# Julia Code for Testing MUMPS Solver with Sparse Matrices and Schur Complement
# This Julia script demonstrates the usage of the MUMPS solver for solving sparse linear systems,
# specifically focusing on computing the Schur complement, which MESTI would utilize in APF method.

# Set the environment variable "MUMPS_PREFIX" to the path of MUMPS library on your machine.
ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.2_par_metis/lib"

# Import necessary packages
using MUMPS3, MPI, LinearAlgebra, SparseArrays, Statistics, Test

# Check if MPI is initialized, and initialize if not
MPI.Initialized() ? nothing : MPI.Init()

# Test set for Schur complement calculation
@testset "Schur Complement (double precision): " begin
    for i ∈ 1:100
        # Define matrix parameters
        m = 50; n = 5; p = .2; T = ComplexF64

        # Generate matrices A, B, C, and D
	A = I + sprand(T,m,m,p)
        B = sprand(T,m,n,p)
        C = sprand(T,n,m,p)
        D = sprand(T,n,n,p)

        # Construct an augmented matrix K
        K = [A B; C D]

        # Initialize a Mumps object with the same type as M
        id = Mumps(K, sym=0, par=1)

        # Set Mumps options
	set_icntl!(id, 4, 0; displaylevel=0); # Turn off diagnostic messages 
	set_icntl!(id, 3, 0); # Turn off  output stream messages from MUMPS
	set_schur_centralized_by_column!(id, 51:55);  # Specify q in the matrix K
	set_job!(id, 1);  # Specify that Mumps will perform an analysis
	set_icntl!(id, 7, 5);  # Specify the Metis ordering

	# Perform the analysis
	invoke_mumps!(id)

        # Specify that Mumps will perform a factorization.
        set_job!(id, 2)

	# Perform the factorization
        invoke_mumps!(id)

	# Take the Schur Complement
	H = get_schur_complement(id)

	# Check the error of the Schur Complement.
	@test mean(abs.(((D - C*inv(Matrix(A))*B) - H)./(D - C*inv(Matrix(A))*B))) ≤ sqrt(eps(Float64))
    end
end

@testset "Schur Complement (single precision): " begin
    for i ∈ 1:100
        # Define matrix parameters
        m = 50; n = 5; p = .2; T = ComplexF32

        # Generate matrices A, B, C, and D
        A = I + sprand(T,m,m,p)
        B = sprand(T,m,n,p)
        C = sprand(T,n,m,p)
        D = sprand(T,n,n,p)

        # Construct an augmented matrix K
        K = [A B; C D]

        # Initialize a Mumps object with the same type as M
        id = Mumps(K, sym=0, par=1)

        # Set Mumps options
        set_icntl!(id, 4, 0; displaylevel=0); # Turn off diagnostic messages
        set_icntl!(id, 3, 0); # Turn off  output stream messages from MUMPS
        set_schur_centralized_by_column!(id, 51:55);  # Specify q in the matrix K
        set_job!(id, 1);  # Specify that Mumps will perform an analysis
        set_icntl!(id, 7, 5);  # Specify the Metis ordering

        # Perform the analysis
        invoke_mumps!(id)

        # Specify that Mumps will perform a factorization.
        set_job!(id, 2)

        # Perform the factorization
        invoke_mumps!(id)

        # Take the Schur Complement
        H = get_schur_complement(id)
        
	# Check the error of the Schur Complement.
	@test mean(abs.(((D - C*inv(Matrix(A))*B) - H)./(D - C*inv(Matrix(A))*B))) ≤ sqrt(eps(Float32))
    end
end
