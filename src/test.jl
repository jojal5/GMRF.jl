
using GMRF, Distributions, LinearAlgebra, SparseArrays

m₁ = 15
m₂ = 15
m = m₁*m₂

order = 2

G = structure_igmrf(m₁,m₂,order)
μ_F = GMRF.iGMRF([0.0],300.0,G)
μ = 10 .+ GMRF.rand(μ_F)

e₁ = ones(Int64,m)
e₂ = repeat(1:m₁, m₂)
e₃ = repeat(1:m₂,inner = m₁)

all(isapprox.(G.W*e₁,0))
all(isapprox.(G.W*e₂,0))
all(isapprox.(G.W*e₃,0))

Q = Array{Float64}(G.W) + e₁*e₁' + e₂*e₂' + e₃*e₃'

C = cholesky(Q)
