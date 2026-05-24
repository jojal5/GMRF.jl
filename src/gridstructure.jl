struct GridStructure
    gridSize::Tuple{Int64,Int64}    # Tuple containing the number of rows and the number of columns
    nbs::Vector{Vector{Int64}}   # list of neighbors for each grid cell
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

"""
    first_order_lattice_neighbors(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

Compute the first-order neighbors of each node in a regular two-dimensional lattice of size `(m₁, m₂)`.

# Details

The lattice contains `m₁ * m₂` nodes, indexed column-wise. Two nodes are first-order
neighbors if they are horizontally or vertically adjacent on the lattice.

Returns a tuple `(nbs, W)`, where:

- `nbs` is the list of first-order neighbors for each node;
- `W` is the corresponding intrinsic CAR precision matrix, with `W[i, i]`
  equal to the number of neighbors of node `i` and `W[i, j] = -1` when
  nodes `i` and `j` are neighbors.
"""
function first_order_lattice_neighbors(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

    m₁ > 0 || throw(ArgumentError("m₁ must be positive."))
    m₂ > 0 || throw(ArgumentError("m₂ must be positive."))

    # 1-off diagonal elements
    v = ones(Int64, m₁)
    v[end] = 0
    V = repeat(v, outer = m₂)
    pop!(V)

    # m₁-off diagonal elements
    U = ones(Int64, m₁ * (m₂ - 1))

    # get the upper triangular part of the matrix
    m = m₁ * m₂
    D = sparse(1:(m - 1), 2:m, V, m, m) +
        sparse(1:(m - m₁), (m₁ + 1):m, U, m, m)

    # make D symmetric
    D = D + D'

    # compute the list of neighbors for each node
    nbs = Vector{Vector{Int64}}(undef, m)
    for i in 1:m
        nbs[i] = findall(!iszero, D[:, i])
    end

    # put the number of neighbors on the diagonal
    W = -D + spdiagm(0 => length.(nbs))

    return nbs, W
end

"""
    second_order_lattice_neighbors(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

Compute the second-order neighbors of each node in a regular two-dimensional lattice of size `(m₁, m₂)`.

# Details

The lattice contains `m₁ * m₂` nodes, indexed column-wise. The function returns a tuple
`(nbs, W)`, where:

- `nbs` is the list of second-order neighbors for each node;
- `W` is the corresponding second-order structure matrix.

The neighbor list is obtained from the negative off-diagonal entries of `W`.
"""
function second_order_lattice_neighbors(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

    m₁ > 0 || throw(ArgumentError("m₁ must be positive."))
    m₂ > 0 || throw(ArgumentError("m₂ must be positive."))

    m = m₁ * m₂
    W = spzeros(Int64, m, m)
    pos = reshape(1:m, m₁, m₂)

    for i in 1:m₁
        for j in 1:m₂

            S = zeros(Int64, m₁, m₂)

            if i - 2 > 0
                S[i-2:i, j] = S[i-2:i, j] + [1, -2, 1]
            end

            if i + 2 <= m₁
                S[i:i+2, j] = S[i:i+2, j] + [1, -2, 1]
            end

            if j - 2 > 0
                S[i, j-2:j] = S[i, j-2:j] + [1, -2, 1]
            end

            if j + 2 <= m₂
                S[i, j:j+2] = S[i, j:j+2] + [1, -2, 1]
            end

            if i - 1 > 0 && i + 1 <= m₁
                S[i-1:i+1, j] = S[i-1:i+1, j] + [-2, 4, -2]
            end

            if j - 1 > 0 && j + 1 <= m₂
                S[i, j-1:j+1] = S[i, j-1:j+1] + [-2, 4, -2]
            end

            if i - 1 > 0 && j + 1 <= m₂
                S[i-1:i, j:j+1] = S[i-1:i, j:j+1] + [-2 2; 2 -2]
            end

            if i + 1 <= m₁ && j + 1 <= m₂
                S[i:i+1, j:j+1] = S[i:i+1, j:j+1] + [2 -2; -2 2]
            end

            if i - 1 > 0 && j - 1 > 0
                S[i-1:i, j-1:j] = S[i-1:i, j-1:j] + [2 -2; -2 2]
            end

            if i + 1 <= m₁ && j - 1 > 0
                S[i:i+1, j-1:j] = S[i:i+1, j-1:j] + [-2 2; 2 -2]
            end

            W[:, pos[i, j]] = S[:]
        end
    end

    nbs = Vector{Vector{Int64}}(undef, m)
    for i in 1:m
        nbs[i] = findall(W[:, i] .< 0)
    end

    return nbs, W
end