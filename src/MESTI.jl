module MESTI

using GeometryPrimitives
using LinearAlgebra
using SparseArrays
using StaticArrays
using Statistics
using TensorCast
using LazyGrids
using Printf
using Libdl
using MPI

include("mumps3_types.jl")
include("mumps3_struc.jl")
include("mumps3_interface.jl")
include("mumps3_convenience_wrappers.jl")
include("mumps3_icntl_alibis.jl")
include("mumps3_printing.jl")

include("get_optimal_PML.jl")
include("build_transverse_function_1d.jl")
include("mesti_subpixel_smoothing.jl")
include("mesti_build_fdfd_matrix.jl")
include("mesti_matrix_solver.jl")
include("setup_longitudinal.jl")
include("mesti_main.jl")
include("mesti_build_channels.jl")
include("mesti2s.jl")

end
