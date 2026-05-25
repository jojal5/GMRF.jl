# struct GridStructure
#     grid_size::Tuple{Int64,Int64}             # Tuple containing the number of rows and the number of columns
#     neighbors::Vector{Vector{Int64}}          # list of neighbors for each grid cell
#     cond_ind_subset::Vector{Vector{Int64}}    # Conditional independant subsets of grid cell
#     W::SparseMatrixCSC{Int64,Int64}           # Structure matrix
#     W̄::SparseMatrixCSC{Int64,Int64}           # Structure matrix minus the diagonal
# end

struct GridStructure
    m₁::Int64   # Number of rows
    m₂::Int64   # Number of columns

    function GridStructure(m₁::Integer, m₂::Integer)
        m₁ > 0 || throw(ArgumentError("m₁ must be positive."))
        m₂ > 0 || throw(ArgumentError("m₂ must be positive."))

        return new(Int64(m₁), Int64(m₂))
    end
end

function showGridStructure(io::IO, obj::GridStructure)

    println(io, "GridStructure")
    println(io, " ", "m₁ = ", obj.m₁)
    println(io, " ", "m₂ = ", obj.m₂)

end

function Base.show(io::IO, obj::GridStructure)
    showGridStructure(io, obj)
end

"""
    first_order_igmrf_structure_matrix(G::GridStructure)::SparseMatrixCSC{Int64,Int64}

Compute the first-order iGMRF structure matrix on the regular lattice `G`.

# Details

The lattice contains `G.m₁ * G.m₂` nodes, indexed column-wise. Two nodes are
first-order neighbors if they are horizontally or vertically adjacent.

Returns the structure matrix `W`, where `W[i, i]` is the number of first-order
neighbors of node `i`, and `W[i, j] = -1` when nodes `i` and `j` are neighbors.
"""
function first_order_igmrf_structure_matrix(G::GridStructure)::SparseMatrixCSC{Int64,Int64}

    m₁ = G.m₁
    m₂ = G.m₂

    # 1-off diagonal elements
    v = ones(Int64, m₁)
    v[end] = 0
    V = repeat(v, outer = m₂)
    pop!(V)

    # m₁-off diagonal elements
    U = ones(Int64, m₁ * (m₂ - 1))

    # get the upper triangular part of the adjacency matrix
    m = m₁ * m₂
    D = sparse(1:(m - 1), 2:m, V, m, m) +
        sparse(1:(m - m₁), (m₁ + 1):m, U, m, m)

    # make D symmetric
    D = D + D'

    # iGMRF structure matrix
    W = spdiagm(0 => vec(sum(D, dims = 2))) - D

    return W
end























"""
    GridStructure(m₁::Integer, m₂::Integer; order::Integer)

Construct the neighborhood structure of a regular two-dimensional lattice of size
`(m₁, m₂)`.

The keyword argument `order` specifies the neighborhood structure. Use `order = 1`
for first-order neighbors and `order = 2` for second-order neighbors.
"""
function GridStructure(m₁::Integer, m₂::Integer; order::Integer)::GridStructure

    m₁ > 0 || throw(ArgumentError("m₁ must be positive."))
    m₂ > 0 || throw(ArgumentError("m₂ must be positive."))
    order in (1, 2) || throw(ArgumentError("order must be either 1 or 2."))

    if order == 1
        nbs, W = first_order_lattice_neighbors(m₁, m₂)
        condIndSubset = first_order_conditional_independent_subsets(m₁, m₂)
    else
        nbs, W = second_order_lattice_neighbors(m₁, m₂)
        condIndSubset = second_order_conditional_independent_subsets(m₁, m₂)
    end

    W̄ = W - spdiagm(0 => diag(W))

    return GridStructure((m₁, m₂), nbs, condIndSubset, W, W̄)
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
    first_order_conditional_independent_subsets(m₁::Integer, m₂::Integer)::Vector{Vector{Int64}}

Compute the first-order conditional independent subsets of a regular two-dimensional lattice of size `(m₁, m₂)`.

# Details

The lattice contains `m₁ * m₂` nodes, indexed column-wise. The function partitions
the nodes into two subsets such that no two first-order neighbors belong to the
same subset.
"""
function first_order_conditional_independent_subsets(m₁::Integer, m₂::Integer)::Vector{Vector{Int64}}

    m₁ > 0 || throw(ArgumentError("m₁ must be positive."))
    m₂ > 0 || throw(ArgumentError("m₂ must be positive."))

    cond_ind_subset_index = 2 * ones(Int64, m₁, m₂)
    cond_ind_subset_index[1:2:end, 1:2:end] .= 1
    cond_ind_subset_index[2:2:end, 2:2:end] .= 1

    return [findall(vec(cond_ind_subset_index) .== i) for i in 1:2]
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

"""
    second_order_conditional_independent_subsets(m₁::Integer, m₂::Integer)::Vector{Vector{Int64}}

Compute a partition of the nodes of a regular two-dimensional lattice of size `(m₁, m₂)`
into conditional independent subsets for a second-order neighborhood structure.

The lattice contains `m₁ * m₂` nodes, indexed column-wise. The function returns seven
subsets of node indices.
"""
function second_order_conditional_independent_subsets(m₁::Integer, m₂::Integer)::Vector{Vector{Int64}}

    m₁ > 0 || throw(ArgumentError("m₁ must be positive."))
    m₂ > 0 || throw(ArgumentError("m₂ must be positive."))

    cond_ind_subset_index = zeros(Int64, m₁, m₂)

    cond_ind_subset_index[1:3:end, 1:4:end] .= 1
    cond_ind_subset_index[2:3:end, 3:4:end] .= 1

    cond_ind_subset_index[1:3:end, 2:4:end] .= 2
    cond_ind_subset_index[2:3:end, 4:4:end] .= 2

    cond_ind_subset_index[1:3:end, 3:4:end] .= 3
    cond_ind_subset_index[2:3:end, 1:4:end] .= 3

    cond_ind_subset_index[1:3:end, 4:4:end] .= 4
    cond_ind_subset_index[2:3:end, 2:4:end] .= 4

    cond_ind_subset_index[3:3:end, 1:3:end] .= 5
    cond_ind_subset_index[3:3:end, 2:3:end] .= 6
    cond_ind_subset_index[3:3:end, 3:3:end] .= 7

    return [findall(vec(cond_ind_subset_index) .== i) for i in 1:7]
end