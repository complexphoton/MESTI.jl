using Pkg

Pkg.add("Arpack")
Pkg.add(url="https://github.com/stevengj/GeometryPrimitives.jl") # installing the latest version of GeometryPrimitives, which fixed its dependency error 
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("Printf")
Pkg.add("Random")
Pkg.add("SparseArrays")
Pkg.add("Statistics")
Pkg.add("Test")
Pkg.add("MPI")
Pkg.add("FFTW")
Pkg.add("MAT")
Pkg.add("LaTeXStrings")
