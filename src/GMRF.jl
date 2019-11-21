module GMRF

using LinearAlgebra, SparseArrays, StatsBase, Distributions

import Distributions.rand, Distributions.logpdf

struct GraphStructure
    order::Int               # Order of the field
    m₁::Int                  # Number of rows
    m₂::Int                  # Number of columns
    m::Int                   # Number of nodes
    nbs::Array{Vector{Int64}}   # list of neighbors
    nnbs::Vector{Int64}      # Number of neighbors for each node
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}       # Structure matrix minus the diagonal
    condIndSubset::Array{Vector{Int64}} # Conditional independant subsets
end

struct iGMRF
    κ::Real                       # Precision of the field
    G::GraphStructure
end

import Distributions.rand, Distributions.logpdf
include("functions.jl")

export structure_igmrf, rand, logpdf

end # module
