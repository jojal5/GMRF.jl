using SparseArrays, Test

using Pkg
pkg"activate ."

using GMRF
import GMRF.second_order_lattice_neighbors

nbs, W = second_order_lattice_neighbors(3, 3)

W_expected = sparse([
         4 -4  1 -4  2  0  1  0  0
        -4  9 -4  2 -6  2  0  1  0
         1 -4  4  0  2 -4  0  0  1
        -4  2  0  9 -6  1 -4  2  0
         2 -6  2 -6 16 -6  2 -6  2
         0  2 -4  1 -6  9  0  2 -4
         1  0  0 -4  2  0  4 -4  1
         0  1  0  2 -6  2 -4  9 -4
         0  0  1  0  2 -4  1 -4  4
    ])

W == W_expected