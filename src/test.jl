
using GMRF, Distributions, LinearAlgebra, SparseArrays

m₁ = 15
m₂ = 15
m = m₁*m₂

order = 1

G = structure_igmrf(m₁,m₂,order)
F = GMRF.iGMRF(30.0,G)
X = GMRF.rand(F)


l = logpdf(F,X)
l = GMRF.fullcondlogpdf(F,X)

e₁ = ones(Int64,m)
e₂ = repeat(1:m₁, m₂)
e₃ = repeat(1:m₂,inner = m₁)

all(isapprox.(G.W*e₁,0))
all(isapprox.(G.W*e₂,0))
all(isapprox.(G.W*e₃,0))

Q = Array{Float64}(G.W) + e₁*e₁' + e₂*e₂' + e₃*e₃'

C = cholesky(Q)

C = cholesky(G.W)

W = G.W

W̄ = G.W̄

h = -W̄*X
J = 300*W

pd = NormalCanon.(h,300*G.nnbs)

l = logpdf.(pd,x)
