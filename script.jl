using LinearAlgebra, SparseArrays, Test

using Pkg
pkg"activate ."

using GMRF


W_expected = sparse([
                4 -4 1 -4 2 0 1 0 0
                -4 9 -4 2 -6 2 0 1 0
                1 -4 4 0 2 -4 0 0 1
                -4 2 0 9 -6 1 -4 2 0
                2 -6 2 -6 16 -6 2 -6 2
                0 2 -4 1 -6 9 0 2 -4
                1 0 0 -4 2 0 4 -4 1
                0 1 0 2 -6 2 -4 9 -4
                0 0 1 0 2 -4 1 -4 4
            ])

G = GMRF.GridStructure(3, 3; order = 2)

G.W == W_expected
G.W̄