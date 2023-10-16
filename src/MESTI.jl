module MESTI

using LinearAlgebra
using SparseArrays
using Statistics
using TensorCast
using LazyGrids
using Printf
# If users specfiy the ENV["MUMPS_PREFIX"], they want to utilize MUMPS, so using MUMPS3 and MPI
if haskey(ENV, "MUMPS_PREFIX")
    using MUMPS3
    using MPI
end

include("get_optimal_PML.jl")
include("build_transverse_function_1d.jl")
include("mesti_build_fdfd_matrix.jl")
include("mesti_matrix_solver.jl")
include("setup_longitudinal.jl")
include("mesti_main.jl")
include("mesti_build_channels.jl")
include("mesti2s.jl")

end
