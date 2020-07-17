module GMRF

using LinearAlgebra, SparseArrays, StatsBase, Distributions

import Distributions.rand, Distributions.logpdf

include("graphstructure.jl")
include("gridstructure.jl")
include("igmrf.jl")

export iGMRF, rand, logpdf, fullconditionals, fullcondlogpdf, getconditional

end # module
