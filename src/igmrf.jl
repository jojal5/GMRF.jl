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

function iGMRF(m₁::Integer, m₂::Integer, order::Integer, κ::Real)::iGMRF
    order in (1, 2) || throw(ArgumentError("order must be either 1 or 2."))

    G = GridStructure(m₁, m₂, order=order)

    if order==1
        rankdef=1
    else
        rankdef = 3
    end 

    pdet = log_pseudodet(G.W, rankdef)

    return iGMRF(G, rankdef, κ, pdet)

end


function rand(F::iGMRF)::Vector{<:Real}

    @assert F.rank_deficiency == 1 || F.rank_deficiency == 3 "The rank deficiency should be either 1 or 3"

    κ = F.κ
    W = F.G.W
    m₁ = F.G.grid_size[1]
    m₂ = F.G.grid_size[2]
    m = m₁ * m₂

    if F.rank_deficiency == 1

        e₁ = ones(m,1)

        A = e₁

        Q = κ*W + e₁*e₁'

    else

        e₁ = ones(m)
        e₂ = repeat(1:m₁, m₂)
        e₃ = repeat(1:m₂,inner = m₁)

        A = hcat(e₁,e₂,e₃)

        Q = κ*W + e₁*e₁' + e₂*e₂' + e₃*e₃'

    end

    C = cholesky(Q)
    L = C.L

    z = randn(m)

    x = L'\z

#     V = zeros(m,size(A,2))
#     for ii=1:size(A,2)
#        V[:,ii] = C\A[:,ii]
#     end
    V = C\A
    W = A'*V
    U = W\(V')
    c = A' * x
    y = x - U' * c

    return y

end


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
