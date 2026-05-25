module GMRF

using LinearAlgebra, Random, SparseArrays, StatsBase, Distributions

import Distributions.rand, Distributions.logpdf
import Random.rand

include("graphstructure.jl")
include("gridstructure.jl")
include("igmrf.jl")
include("utils.jl")

export iGMRF, rand, logpdf

end # module
