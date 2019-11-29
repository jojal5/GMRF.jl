
using GMRF, Distributions, LinearAlgebra, SparseArrays

m₁ = 15
m₂ = 15
m = m₁*m₂

G = GMRF.GridStructure(m,(m₁,m₂),[[1,2,3],[1,2,3]],[3,3,3],[[1,2,3],[1,2,3]])
SG = GMRF.subGridStructure(G, [2,5,10], [[1,2],[2,3]])



order = 1
F = structure_igmrf(m₁,m₂,order,30)
X = GMRF.rand(F)
l = logpdf(F,X)
l = GMRF.fullcondlogpdf(F,X)


order = 2
F = structure_igmrf(m₁,m₂,order,30)
X = GMRF.rand(F)
l = logpdf(F,X)
l = GMRF.fullcondlogpdf(F,X)
