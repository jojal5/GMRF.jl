using LinearAlgebra, SparseArrays, Test, Plots

using Pkg
pkg"activate ."

using GMRF

F = iGMRF(20, 20, 1, 1.)

y = rand(F)

# x = reshape(y, 20, 20)
# heatmap(x)

@time logpdf(F, y)

@time pd = GMRF.full_conditionals(F,y)

@time h, Q = GMRF.full_conditional_canonical_parameters(F,y)

@time GMRF.fullcondlogpdf(F,y)

@time GMRF.full_conditionals_logpdf(F,y)

