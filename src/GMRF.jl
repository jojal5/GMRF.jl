module GMRF

using Distributions, LinearAlgebra, Random, SparseArrays

import Distributions.logpdf
import Random.rand

include("graphstructure.jl")
include("gridstructure.jl")
include("igmrf.jl")
include("utils.jl")

export iGMRF, rand, logpdf

end # module
