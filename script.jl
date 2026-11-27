using LinearAlgebra, SparseArrays, Test

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



function second_order_igmrf_structure_matrix_kron(G::GridStructure)::SparseMatrixCSC{Int64,Int64}
    function axis_matrices(n)
        Iₙ = spdiagm(0 => ones(Int64, n))
        D₁ = diff(Iₙ; dims=1)
        D₂ = diff(D₁; dims=1)
        return Iₙ, D₁' * D₁, D₂' * D₂
    end

    I₁, L₁, B₁ = axis_matrices(G.m₁)
    I₂, L₂, B₂ = axis_matrices(G.m₂)

    # Second differences along each axis, plus twice the mixed-difference penalty.
    return kron(I₂, B₁) + kron(B₂, I₁) + 2 * kron(L₂, L₁)
end


G = GridStructure(5,15)

@time expected = GMRF.second_order_igmrf_structure_matrix(G)
@test second_order_igmrf_structure_matrix_kron(G) == expected