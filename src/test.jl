
using GMRF, Distributions, LinearAlgebra, SparseArrays

m₁ = 15
m₂ = 15
m = m₁*m₂

V = collect(1:m)
A = sort(rand(1:m, 2))
B = setdiff(V, A)

order = 1
F = iGMRF(m₁,m₂,order, 10)
X = rand(F)
l = logpdf(F,X)
fc = fullconditionals(F,X)
l = fullcondlogpdf(F,X)
fc = getconditional(F, B, X[B])


order = 2
F = iGMRF(m₁,m₂,order, 10)
X = rand(F)
l = logpdf(F,X)
fc = fullconditionals(F,X)
l = fullcondlogpdf(F,X)
fc = getconditional(F, B, X[B])
