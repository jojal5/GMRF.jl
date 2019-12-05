module GMRF

using LinearAlgebra, SparseArrays, StatsBase, Distributions

import Distributions.rand, Distributions.logpdf

struct GridStructure
    m::Int64                      # Number fo grid cells
    gridSize::Tuple{Int64,Int64}    # Tuple containing the number of rows and the number of columns
    nbs::Array{Vector{Int64}}   # list of neighbors for each grid cell
    nnbs::Vector{Int64}      # Number of neighbors for each node
    condIndSubset::Array{Vector{Int64}} # Conditional independant subsets of grid cell
end

struct GMRF
    G::GridStructure
    κ::Float64                       # Precision of the field
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}       # Structure matrix minus the diagonal
end

struct iGMRF
    G::GridStructure
    rankDeficiency::Int64
    κ::Float64                       # Precision of the field
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}       # Structure matrix minus the diagonal
end

import Distributions.rand, Distributions.logpdf
include("functions.jl")

export structure_igmrf, rand, logpdf

end # module
