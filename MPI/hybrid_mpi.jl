# With hybrid mpi parallalization, this code test the unitarity of scattering matrices S upon random permittivity profiles.
# We check S'*S ≈ I, an identity matrix.

# Call necessary packages
using MESTI, MPI

MPI.Initialized() ? nothing : MPI.Init()
root = 0;
comm = MPI.COMM_WORLD;

if MPI.Comm_rank(comm) != root
  # Call the worker processors when we need them (during analyzing and factorization)
  T = ComplexF32
  id = Mumps{T}(; sym=0, par=1)
  set_job!(id,1);
  invoke_mumps!(id);
  set_job!(id,2);
  invoke_mumps!(id);
  finalize!(id)
else
  # Specify parameters of the system in the main processor
  syst = Syst()
  syst.xBC = "periodic"
  syst.yBC = "periodic"
  syst.dx = 1
  syst.wavelength = 5
  syst.epsilon_low = 1
  syst.epsilon_high = 1
  epsilon_max = 4
  epsilon_min = 1

  # Define the size of the scattering region (1 wavelength by 1 wavelength by 1 wavelength)
  nx, ny, nz = 2, 5, 2
  nx_Ex = nx; ny_Ex = ny; nz_Ex = nz -1
  nx_Ey = nx; ny_Ey = ny; nz_Ey = nz -1
  nx_Ez = nx; ny_Ez = ny; nz_Ez = nz

  # Use optimized PML parameters for this resolution to reduce error
  zpml = mesti_optimal_pml_params(syst.wavelength/syst.dx)
  zpml.npixels = 25
  syst.zPML = [zpml]

  # Specify inputs and output
  input = channel_type()
  output = channel_type()
  # Input from both sides with both s-polarization and p-polarization
  input.side = "both"
  input.polarization = "both"
  # Output to both sides with both s-polarization and p-polarization
  output.side = "both"
  output.polarization = "both"

  syst.epsilon_xx = rand(nx_Ex,ny_Ex,nz_Ex)* (epsilon_max-epsilon_min) .+ epsilon_min
  syst.epsilon_yy = rand(nx_Ey,ny_Ey,nz_Ey)* (epsilon_max-epsilon_min) .+ epsilon_min
  syst.epsilon_zz = rand(nx_Ez,ny_Ez,nz_Ez)* (epsilon_max-epsilon_min) .+ epsilon_min

  opts = Opts()
  opts.verbal_solver = true # Print system information from the MUMPS solver.
  # opts.nthreads_OMP = 16 # Assign the number of threads used in MUMPS.
  (S, _, _)= mesti2s(syst, input, output, opts)
end
