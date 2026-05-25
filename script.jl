using LinearAlgebra, SparseArrays, Test, Plots

using Pkg
pkg"activate ."

using GMRF

F = iGMRF(20, 20, order = 1, precision = 1.)

y = rand(F)

# x = reshape(y, 20, 20)
# heatmap(x)

@time logpdf(F, y)

@time pd = GMRF.full_conditionals(F,y)

@time h, Q = GMRF.full_conditional_canonical_parameters(F,y)

@time l =  GMRF.full_conditionals_logpdf(F,y)

sum(l)

B = [59; 70; 100; 117; 206; 221; 338; 349; 373; 380]
x = y[B]

GMRF.conditional_distribution(F, B, x)





using LinearAlgebra, SparseArrays, Test, Plots

using Pkg
pkg"activate ."

using GMRF

F = iGMRF(3, 3, order = 1, precision = 1.)

y = rand(F)

l1 = logpdf(F, y)

l =  GMRF.full_conditionals_logpdf(F,y)

l2 = sum(l)