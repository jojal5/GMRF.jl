module GMRF

using LinearAlgebra, SparseArrays, StatsBase, Distributions

import Distributions.rand, Distributions.logpdf

struct GraphStructure
    m::Int64                      # Number fo grid cells
    nbs::Array{Vector{Int64}}   # list of neighbors for each grid cell
end

struct GridStructure
    gridSize::Tuple{Int64,Int64}    # Tuple containing the number of rows and the number of columns
    nbs::Array{Vector{Int64}}   # list of neighbors for each grid cell
    condIndSubset::Vector{Vector{Int64}} # Conditional independant subsets of grid cell
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}       # Structure matrix minus the diagonal
end

function showGraphStructure(io::IO, obj::GraphStructure; prefix::String = "")

    println(io, prefix, "GraphStructure")
    println(io, prefix, "m :\t", obj.m)
    println(io, prefix, "nbs :\t", typeof(obj.nbs), "[", length(obj.nbs), "]")

end

function Base.show(io::IO, obj::GraphStructure)

    showGraphStructure(io, obj)

end

function showGridStructure(io::IO, obj::GridStructure; prefix::String = "")

    println(io, prefix, "GridStructure")
    println(io, prefix, "gridSize :\t", obj.gridSize)
    println(io, prefix, "nbs :\t\t", typeof(obj.nbs), "[", length(obj.nbs), "]")

end

function Base.show(io::IO, obj::GridStructure)

    showGridStructure(io, obj)

end

import Distributions.rand, Distributions.logpdf
include("igmrf.jl")

export iGMRF, rand, logpdf, fullconditionals, fullcondlogpdf, getconditional

end # module
