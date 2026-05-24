

"""
    second_order_lattice_neighbors(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

Compute the second order neighbors of each point of a regular lattice of size `(m₁, m₂)`.

# Details

Returns a tuple the tuple `(nbs, W)` where
- nbs is the list of neighbors for each point
- W is the corresponding structure matrix for the lattice
"""
function second_order_lattice_neighbors(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}
    # Alternative by adding molecules. There should not be missing values in the grid.
    m = m₁ * m₂
    W = spzeros(Int64,m,m)
    pos = reshape(1:m,m₁,m₂)

    for i=1:m₁
        for j=1:m₂

            S = zeros(Int64,m₁,m₂)

            if (i-2>0)
               S[i-2:i,j] =  S[i-2:i,j] + [1, -2, 1]
            end

            if (i+2<=m₁)
               S[i:i+2,j] =  S[i:i+2,j] + [1, -2, 1]
            end

            if (j-2>0)
               S[i,j-2:j] =  S[i,j-2:j] + [1,-2, 1]
            end

            if (j+2<=m₂)
               S[i,j:j+2] =  S[i,j:j+2] + [1,-2, 1]
            end



            if (i-1>0) && (i+1<=m₁)
                S[i-1:i+1,j] = S[i-1:i+1,j] + [-2, 4, -2]
            end

            if (j-1>0) && (j+1<=m₂)
                S[i,j-1:j+1] = S[i,j-1:j+1] + [-2, 4, -2]
            end



            if (i-1>0) && (j+1<=m₂)
                S[i-1:i,j:j+1] = S[i-1:i,j:j+1] + [-2 2; 2 -2]
            end

            if (i+1<=m₁) && (j+1<=m₂)
                S[i:i+1,j:j+1] = S[i:i+1,j:j+1] + [2 -2; -2 2]
            end

            if (i-1>0) && (j-1>0)
                S[i-1:i,j-1:j] = S[i-1:i,j-1:j] + [2 -2; -2 2]
            end

            if (i+1<=m₁) && (j-1>0)
                S[i:i+1,j-1:j] = S[i:i+1,j-1:j] + [-2 2; 2 -2]
            end

            W[:,pos[i,j]] = S[:]

        end
    end

    # Compute the list of neighbors for each node
    nbs =  Array{Int64,1}[]
    for i = 1:m
        push!(nbs,findall(W[:,i] .< 0))
    end

    return (nbs, W)

end