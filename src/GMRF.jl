module GMRF

using LinearAlgebra, SparseArrays, StatsBase, Distributions

import Distributions.rand, Distributions.logpdf

struct GridStructure
    m::Int                      # Number fo grid cells
    gridSize::Tuple{Int,Int}    # Tuple containing the number of rows and the number of columns
    nbs::Array{Vector{Int}}   # list of neighbors for each grid cell
    nnbs::Vector{Int}      # Number of neighbors for each node
    condIndSubset::Array{Vector{Int}} # Conditional independant subsets of grid cell
end


struct subGridStructure
    G::GridStructure
    V::Vector{Int}  # Grid cells of the subgrid
    condIndSubset::Array{Vector{Int}} # Conditional independant subsets of the grid cell
end


struct iGMRF
    G::GridStructure
    rankDeficiency::Int
    κ::Real                       # Precision of the field
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}       # Structure matrix minus the diagonal
end

import Distributions.rand, Distributions.logpdf
include("functions.jl")

export structure_igmrf, rand, logpdf

end # module
