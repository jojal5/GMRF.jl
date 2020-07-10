
using GMRF, Distributions, LinearAlgebra, SparseArrays

m₁ = 15
m₂ = 15
m = m₁*m₂

order = 1
F = iGMRF(m₁,m₂,order, 10)
X = rand(F)
l = logpdf(F,X)
l = fullcondlogpdf(F,X)


order = 2
F = iGMRF(m₁,m₂,order, 10)
X = rand(F)
l = logpdf(F,X)
l = fullcondlogpdf(F,X)
