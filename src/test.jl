
using GMRF, Distributions, LinearAlgebra, SparseArrays

m₁ = 15
m₂ = 15
m = m₁*m₂

order = 1
G = GMRF.gridstructure_igmrf(m₁,m₂,order)
F = GMRF.iGMRF(G,1,10)
X = GMRF.rand(F)
l = logpdf(F,X)
l = GMRF.fullcondlogpdf(F,X)


order = 2
G = GMRF.gridstructure_igmrf(m₁,m₂,order)
F = GMRF.iGMRF(G,3,10)
X = GMRF.rand(F)
l = logpdf(F,X)
l = GMRF.fullcondlogpdf(F,X)
