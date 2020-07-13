struct GraphStructure
    m::Int64                      # Number fo grid cells
    nbs::Vector{Vector{Int64}}   # list of neighbors for each grid cell
end

function showGraphStructure(io::IO, obj::GraphStructure; prefix::String = "")

    println(io, prefix, "GraphStructure")
    println(io, prefix, "m :\t", obj.m)
    println(io, prefix, "nbs :\t", typeof(obj.nbs), "[", length(obj.nbs), "]")

end

function Base.show(io::IO, obj::GraphStructure)

    showGraphStructure(io, obj)

end
