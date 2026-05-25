"""
    log_pseudodet(W::SparseMatrixCSC, rank_deficiency::Integer)

Compute the logarithm of the pseudo-determinant of a symmetric sparse matrix `W`.

The argument `rank_deficiency` is the number of zero eigenvalues of `W`. The
pseudo-determinant is computed as the product of the nonzero eigenvalues, after
converting `W` to a dense symmetric matrix. For large matrix, consider another 
approach based on sparse decomposition.

# Detail

This quantity is important for computing the likelihood as defined at Eq. (3.13) by Rue & Held (2002).
"""
function log_pseudodet(W::SparseMatrixCSC, rank_deficiency::Integer)

    n = size(W, 1)

    size(W, 2) == n || throw(DimensionMismatch("W must be square."))
    0 <= rank_deficiency < n || throw(ArgumentError("rank_deficiency must be between 0 and size(W, 1) - 1."))

    λ = eigvals(Symmetric(Matrix(W)))
    sort!(λ)

    λ⁺ = λ[(rank_deficiency + 1):end]

    return sum(log, λ⁺)
end