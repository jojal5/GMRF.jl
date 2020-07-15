struct GridStructure
    gridSize::Tuple{Int64,Int64}    # Tuple containing the number of rows and the number of columns
    nbs::Array{Vector{Int64}}   # list of neighbors for each grid cell
    condIndSubset::Vector{Vector{Int64}} # Conditional independant subsets of grid cell
    W::SparseMatrixCSC{Int64,Int64}       # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}       # Structure matrix minus the diagonal
end

function showGridStructure(io::IO, obj::GridStructure; prefix::String = "")

    println(io, prefix, "GridStructure")
    println(io, prefix, "gridSize :\t", obj.gridSize)
    println(io, prefix, "nbs :\t\t", typeof(obj.nbs), "[", length(obj.nbs), "]")

end

function Base.show(io::IO, obj::GridStructure)

    showGridStructure(io, obj)

end
