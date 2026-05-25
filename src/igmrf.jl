struct iGMRF
    G::GridStructure
    order::Int64
    precision::Float64
    cond_ind_subset::Vector{Vector{Int64}}    # Conditional independent subsets of grid cells
    W::SparseMatrixCSC{Int64,Int64}           # Structure matrix
    W̄::SparseMatrixCSC{Int64,Int64}           # Structure matrix minus the diagonal
    log_pseudodet_W::Float64                  # Log pseudo-determinant of the structure matrix (useful for logpdf computing)
end

"""
    iGMRF(G::GridStructure; order::Integer = 1, precision::Real = 1.0)

Construct an intrinsic Gaussian Markov random field on the regular lattice `G`.

The keyword argument `order` specifies the neighborhood structure. Use `order = 1`
for a first-order iGMRF and `order = 2` for a second-order iGMRF. The argument
`precision` is the positive precision parameter.
"""
function iGMRF(G::GridStructure; order::Integer = 1, precision::Real = 1.0)::iGMRF

    order in (1, 2) || throw(ArgumentError("order must be either 1 or 2."))
    precision > 0 || throw(ArgumentError("precision must be positive."))

    order = Int64(order)
    precision = Float64(precision)

    if order == 1
        W = first_order_igmrf_structure_matrix(G)
        cond_ind_subset = first_order_igmrf_conditional_independent_subsets(G)
    else
        W = second_order_igmrf_structure_matrix(G)
        cond_ind_subset = second_order_igmrf_conditional_independent_subsets(G)
    end

    W̄ = W - spdiagm(0 => Vector(diag(W)))

    rank_deficiency = order == 1 ? 1 : 3

    log_pseudodet_W = log_pseudodet(W, rank_deficiency)

    return iGMRF(
        G,
        order,
        precision,
        cond_ind_subset,
        W,
        W̄,
        log_pseudodet_W,
    )
end

"""
    iGMRF(m₁::Integer, m₂::Integer; order::Integer = 1, precision::Real = 1.0)

Construct an intrinsic Gaussian Markov random field on a regular two-dimensional
lattice of size `(m₁, m₂)`.
"""
function iGMRF(
    m₁::Integer,
    m₂::Integer;
    order::Integer = 1,
    precision::Real = 1.0,
)::iGMRF

    G = GridStructure(m₁, m₂)

    return iGMRF(G; order = order, precision = precision)
end

function Base.show(io::IO, obj::iGMRF)

    println(io, "iGMRF")
    showGridStructure(io, obj.G; prefix = " ")
    println(io, " order = ", obj.order)
    println(io, " precision = ", obj.precision)
end

"""
    constraint_matrix(F::iGMRF)::Matrix{Float64}

Construct the constraint matrix associated with the intrinsic Gaussian Markov
random field `F`.

For a first-order iGMRF, the constraint matrix contains the constant vector. For
a second-order iGMRF, it contains the constant vector and the two coordinate
vectors.
"""
function constraint_matrix(F::iGMRF)::Matrix{Float64}

    rank_deficiency = F.order == 1 ? 1 : 3

    m₁, m₂ = F.G.m₁, F.G.m₂
    m = m₁ * m₂

    e₁ = ones(Float64, m)

    if rank_deficiency == 1
        return reshape(e₁, :, 1)
    end

    e₂ = Float64.(repeat(1:m₁, m₂))
    e₃ = Float64.(repeat(1:m₂, inner = m₁))

    return hcat(e₁, e₂, e₃)
end


"""
    rand(F::iGMRF)::Vector{Float64}
    rand(rng::AbstractRNG, F::iGMRF)::Vector{Float64}

Generate one realization from the intrinsic Gaussian Markov random field `F`.

The realization is sampled using the precision matrix `κW` and then projected
onto the constraint space associated with the rank deficiency of `F`.

Use `rand(rng, F)` with an explicit random number generator for reproducible
simulation.
"""
function rand(rng::AbstractRNG, F::iGMRF)::Vector{Float64}

    κ = F.precision

    W = F.W
    A = constraint_matrix(F)
    m = F.G.m₁ * F.G.m₂

    Q = κ * W + A * A'
    C = cholesky(Symmetric(Q))

    z = randn(rng, m)
    x = C.L' \ z

    V = C \ A
    M = A' * V
    c = A' * x

    return x - V * (M \ c)
end

rand(F::iGMRF)::Vector{Float64} = rand(Random.default_rng(), F)


"""
    logpdf(F::iGMRF, y::AbstractVector{<:Real})::Real

Compute the pseudo log-density of the intrinsic Gaussian Markov random field `F` at `y`.

The density is evaluated on the subspace of dimension `m - k`, where `m` is the
number of grid cells and `k` is the rank deficiency of the structure matrix. The
normalizing constant uses the log pseudo-determinant of the structure matrix `W`.

This implements Eq. (3.13) of Rue and Held (2002).
"""
function logpdf(F::GMRF.iGMRF, y::AbstractVector{<:Real})::Real

    κ = F.precision
    m = F.G.m₁ * F.G.m₂
    rank_deficiency = F.order == 1 ? 1 : 3
    W = F.W

    length(y) == m || throw(DimensionMismatch("length(y) must be equal to F.G.m₁ * F.G.m₂."))

    r = m - rank_deficiency

    v = W * y
    q = dot(y, v)

    return -0.5 * r * log(2π) +
            0.5 * r * log(κ) +
            0.5 * F.log_pseudodet_W -
            0.5 * κ * q
end


"""
    full_conditional_canonical_parameters(F::iGMRF, y::AbstractVector{<:Real})

Compute the canonical parameters of the full conditional distributions of the
intrinsic Gaussian Markov random field `F` at all grid cells, given the current
field values `y`.

Returns a tuple `(h, Q)`, where `h[i]` is the canonical parameter and `Q[i]` is
the precision of the full conditional distribution at grid cell `i`.
"""
function full_conditional_canonical_parameters(
    F::iGMRF,
    y::AbstractVector{<:Real},
)

    κ = F.precision

    W = F.W
    W̄ = F.W̄

    length(y) == size(W, 1) ||
        throw(DimensionMismatch("length(y) must be equal to the number of grid cells."))

    h = -κ .* (W̄ * y)

    Q = Vector(diag(W))
    Q .*= κ

    return h, Q
end


"""
    full_conditionals(F::iGMRF, y::AbstractVector{<:Real})::Vector{NormalCanon}

Compute the full conditional distributions of the intrinsic Gaussian Markov
random field `F` at all grid cells, given the current field values `y`.

The distributions are returned in canonical normal form.
"""
function full_conditionals(F::iGMRF, y::AbstractVector{<:Real})::Vector{NormalCanon}

    h, Q = full_conditional_canonical_parameters(F, y)

    return NormalCanon.(h, Q)
end

"""
    full_conditionals_logpdf(F::iGMRF, y::AbstractVector{<:Real})::Vector{Float64}

Compute the log-density of each grid-cell value under its full conditional
distribution.

For each grid cell `i`, this returns `log f(y[i] | y[-i])`, where the full
conditional distribution is represented in canonical normal form with canonical
parameter `h[i]` and precision `Q[i]`.
"""
function full_conditionals_logpdf(F::iGMRF, y::AbstractVector{<:Real})::Vector{Float64}

    h, Q = full_conditional_canonical_parameters(F, y)

    return @. h * y - 0.5 * Q * y^2 - 0.5 * log(2π) + 0.5 * log(Q) - 0.5 * h^2 / Q
end


"""
    conditional_distribution(F::iGMRF, B::AbstractVector{<:Integer}, x::AbstractVector{<:Real})::MvNormalCanon

Compute the conditional distribution of the grid cells outside `B`, given the
values `x` at the grid cells in `B`.

The vector `x` must have the same length and ordering as `B`. The returned
distribution is represented in canonical normal form.
"""
function conditional_distribution(
    F::iGMRF,
    B::AbstractVector{<:Integer},
    x::AbstractVector{<:Real}
)::MvNormalCanon

    W = F.W
    κ = F.precision
    m = F.G.m₁ * F.G.m₂

    all(1 .<= B) && all(B .<= m) || throw(ArgumentError("all indices in B must be between 1 and $m."))
    allunique(B) || throw(ArgumentError("indices in B must be unique."))
    length(x) == length(B) || throw(DimensionMismatch("length(x) must be equal to length(B)."))

    A = setdiff(1:m, B)

    W_AA = W[A, A]
    W_AB = W[A, B]

    h = -κ .* (W_AB * x)
    J = κ .* Matrix(W_AA)

    return MvNormalCanon(h, J)
end