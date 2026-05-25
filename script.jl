using LinearAlgebra, SparseArrays, Test

using Pkg
pkg"activate ."

using GMRF

G = GMRF.GridStructure(20, 20, order =1)
F = GMRF.iGMRF(G, 1, 30.)

y = GMRF.rand(F)


@time GMRF.logpdf(F, y)

function logpdf2(F::GMRF.iGMRF, y::Array{<:Real})::Real

    κ = F.κ
    m = prod(F.G.grid_size)
    k = F.rankDeficiency

    W = F.G.W

    v = κ * (W*y)
    q = dot(y,v)

    lpdf =  .5*(m-k)*log(κ) - .5*q

    return lpdf

end

@time logpdf2(F, y)

function logpdf3(F::GMRF.iGMRF, y::AbstractVector{<:Real})::Real

    κ = F.κ
    m = prod(F.G.grid_size)
    k = F.rankDeficiency
    W = F.G.W

    κ > 0 || throw(ArgumentError("κ must be positive."))
    length(y) == m || throw(DimensionMismatch("length(y) must be equal to prod(F.G.grid_size)."))

    r = m - k

    Wy = W * y
    q = dot(y, Wy)

    return -0.5 * r * log(2π) +
            0.5 * r * log(κ) +
            0.5 * logdetW -
            0.5 * κ * q
end

@time logpdf3(F, y)



function log_pseudodet(W::SparseMatrixCSC{<:Real,<:Integer}; tol::Real = 1e-10)
    λ = eigvals(Symmetric(Matrix(W)))
    λ⁺ = λ[λ .> tol]
    return sum(log, λ⁺)
end

logdetW = log_pseudodet(W)

eigvals(F.G.W)

λ = eigvals(Symmetric(Matrix(F.G.W)))