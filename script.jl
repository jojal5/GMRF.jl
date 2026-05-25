using LinearAlgebra, SparseArrays, Test, Plots

using Pkg
pkg"activate ."

using GMRF

F = iGMRF(20, 20, order = 1, precision = 1.)

F = iGMRF(3, 2, order = 2, precision = 1.)
GMRF.second_order_igmrf_conditional_independent_subsets(F.G)
GMRF.log_pseudodet(F.W, 3)

F = iGMRF(3, 2, order = 1, precision = 1.)

GMRF.log_pseudodet(F.W, 1)

F.G

y = rand(F)

# x = reshape(y, 20, 20)
# heatmap(x)

@time logpdf(F, y)

@time pd = GMRF.full_conditionals(F,y)

@time h, Q = GMRF.full_conditional_canonical_parameters(F,y)

@time GMRF.fullcondlogpdf(F,y)

@time GMRF.full_conditionals_logpdf(F,y)











using LinearAlgebra, SparseArrays

m₁ = 1
    m₂ = 1
    
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

    D = D + D'

    sum(D, dims=1)

    W = -D + spdiagm(0 => vec(sum(D, dims=1)))