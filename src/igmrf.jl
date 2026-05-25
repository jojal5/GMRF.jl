struct iGMRF
    G::GridStructure
    rank_deficiency::Int64
    κ::Float64                 # Precision of the field
    log_pseudodet_W::Float64   # Log pseudo-determinant of the structure matrix (useful for logpdf computing)
end

function Base.show(io::IO, obj::iGMRF)

    println(io, "iGMRF")
    println(io, "G :")
    showGridStructure(io, obj.G, prefix = "\t\t\t")
    println(io)
    println(io, "rank deficiency :\t", obj.rank_deficiency)
    println(io, "κ :\t\t\t", obj.κ)

end

"""
    iGMRF(m₁::Integer, m₂::Integer, order::Integer, κ::Real)

Construct an intrinsic Gaussian Markov random field on a regular two-dimensional
lattice of size `(m₁, m₂)`.

The argument `order` specifies the neighborhood structure. Use `order = 1` for a
first-order iGMRF and `order = 2` for a second-order iGMRF. The parameter `κ` is
the precision parameter.

The constructor builds the corresponding `GridStructure`, sets the rank deficiency,
and precomputes the log pseudo-determinant of the structure matrix.
"""
function iGMRF(m₁::Integer, m₂::Integer, order::Integer, κ::Real)::iGMRF

    order in (1, 2) || throw(ArgumentError("order must be either 1 or 2."))
    κ > 0 || throw(ArgumentError("κ must be positive."))

    G = GridStructure(m₁, m₂; order = order)

    rank_deficiency = order == 1 ? 1 : 3
    log_pseudodet_W = log_pseudodet(G.W, rank_deficiency)

    return iGMRF(G, rank_deficiency, κ, log_pseudodet_W)
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

    rank_deficiency = F.rank_deficiency
    rank_deficiency in (1, 3) ||
        throw(ArgumentError("rank_deficiency must be either 1 or 3."))

    m₁, m₂ = F.G.grid_size
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

    κ = F.κ
    κ > 0 || throw(ArgumentError("κ must be positive."))

    W = F.G.W
    A = constraint_matrix(F)
    m = prod(F.G.grid_size)

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

    κ = F.κ
    m = prod(F.G.grid_size)
    k = F.rank_deficiency
    W = F.G.W

    length(y) == m || throw(DimensionMismatch("length(y) must be equal to prod(F.G.grid_size)."))

    r = m - k

    v = W * y
    q = dot(y, v)

    return -0.5 * r * log(2π) +
            0.5 * r * log(κ) +
            0.5 * F.log_pseudodet_W -
            0.5 * κ * q
end




function fullconditionals(F::iGMRF, y::Vector{<:Real})::Vector{NormalCanon}

    κ = F.κ

    W̄ = F.G.W̄
    W = F.G.W

    Q = κ * Array(diag(F.G.W))
    h = -κ*(W̄*y)

    pd = NormalCanon.(h,Q)

    return pd

end

function fullcondlogpdf(F::iGMRF, y::Vector{<:Real})::Vector{<:Real}

    pd = fullconditionals(F::iGMRF,y::Vector{<:Real})

    clpdf = logpdf.(pd,y)

    return clpdf

end

function getconditional(F::GMRF.iGMRF, B::Vector{<:Integer}, x::Vector{<:Real})::MvNormalCanon

    W = F.G.W

    sort!(B)

    A = setdiff(1:(F.G.grid_size[1] * F.G.grid_size[2]), B)

    Waa = W[A,A]
    Wab = W[A,B]

    h = -Wab*x*F.κ

    J = Array(F.κ*Waa)

    pd = MvNormalCanon(h,J)

    return pd

end
