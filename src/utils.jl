"""
    fo_nbs(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

Compute the first order neighbors of each point of a regular lattice of size `(m₁, m₂)`.

# Details

Returns a tuple the tuple `(nbs, W)` where
- nbs is the list of neighbors for each point
- W is the corresponding structure matrix for the lattice
"""
function fo_nbs(m₁::Integer, m₂::Integer)::Tuple{Vector{Vector{Int64}}, SparseMatrixCSC{Int64,Int64}}

    # 1-off diagonal elements
    v = ones(Int64,m₁)
    v[end] = 0
    V = repeat(v,outer=m₂)
    pop!(V)

    # n-off diagonal elements
    U = ones(Int64,m₁*(m₂-1))

    # get the upper triangular part of the matrix
    m = m₁ * m₂
    D = sparse(1:(m-1), 2:m, V, m, m) + sparse(1:(m-m₁),(m₁+1):m, U, m, m)

    # make D symmetric
    D = D + D'

    # Compute the list of neighbors for each node
    nbs = fill(Int[], m)
    for i = 1:m
        nbs[i] = findall(!iszero, D[:,i])
    end

    # Put the number of neighbors on the diagonal
    W = -D + spdiagm(0 => length.(nbs))

    return (nbs, W)

end