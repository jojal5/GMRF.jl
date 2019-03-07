module GMRF

using LinearAlgebra, SparseArrays, StatsBase

include("functions.jl")

struct GraphStructure
    order::Int               # Order of the field
    m₁::Int                  # Number of rows
    m₂::Int                  # Number of columns
    m::Int                   # Number of nodes
    nbs::Array{Array{Int64}}   # list of neighbors
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    condIndSubset::Array{Array{Int64}} # Conditional independant subsets
end

struct iGMRF
    μ::Array{T} where T<:Real     # Vector of mean parameters
    κ::Real                       # Precision of the field
    G::GraphStructure
end

export structure_igmrf, rand, logpdf

end # module
