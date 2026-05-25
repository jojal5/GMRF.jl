using GMRF
using LinearAlgebra, SparseArrays, Random
using Test

@testset "GMRF.jl" begin
    include("graphstructure_test.jl")
    include("gridstructure_test.jl")
    include("igmrf_test.jl")
    include("utils_test.jl")
end
