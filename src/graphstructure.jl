struct GraphStructure
    m::Int64                      # Number fo grid cells
    neighbors::Vector{Vector{Int64}}   # list of neighbors for each grid cell
end

function showGraphStructure(io::IO, obj::GraphStructure; prefix::String = "")

    println(io, prefix, "GraphStructure")
    println(io, prefix, "m :\t", obj.m)
    println(io, prefix, "neighbors :\t", typeof(obj.neighbors), "[", length(obj.neighbors), "]")

end

function Base.show(io::IO, obj::GraphStructure)

    showGraphStructure(io, obj)

end
