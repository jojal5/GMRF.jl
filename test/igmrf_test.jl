
@testset "igmrf.jl" begin

    @testset "constructors" begin
        @testset "iGMRF of order 1" begin
            F = iGMRF(3, 2, order=1, precision=1.)

            @test F.G == GridStructure(3, 2)
            @test F.order == 1
            @test F.precision ≈ 1.0

            W_expected = sparse([
                2 -1 0 -1 0 0
                -1 3 -1 0 -1 0
                0 -1 2 0 0 -1
                -1 0 0 2 -1 0
                0 -1 0 -1 3 -1
                0 0 -1 0 -1 2
            ])
            @test F.W == W_expected

            W̄_expected = sparse([
                0 -1 0 -1 0 0
                -1 0 -1 0 -1 0
                0 -1 0 0 0 -1
                -1 0 0 0 -1 0
                0 -1 0 -1 0 -1
                0 0 -1 0 -1 0
            ])
            @test F.W̄ == W̄_expected

            @test F.cond_ind_subset == [[1, 3, 5], [2, 4, 6]]

            @test F.log_pseudodet_W ≈ 4.499809670330264

        end

        @testset "iGMRF of order 2" begin
            F = iGMRF(3, 2, order=2, precision=1.)

            @test F.G == GridStructure(3, 2)
            @test F.order == 2
            @test F.precision ≈ 1.0

            W_expected = sparse([
                3 -4 1 -2 2 0
                -4 8 -4 2 -4 2
                1 -4 3 0 2 -2
                -2 2 0 3 -4 1
                2 -4 2 -4 8 -4
                0 2 -2 1 -4 3
            ])
            @test F.W == W_expected

            W̄_expected = sparse([
                0 -4 1 -2 2 0
                -4 0 -4 2 -4 2
                1 -4 0 0 2 -2
                -2 2 0 0 -4 1
                2 -4 2 -4 0 -4
                0 2 -2 1 -4 0
            ])
            @test F.W̄ == W̄_expected

            @test F.cond_ind_subset == [[1], [4], [2], [5], [3], [6], []]

            @test F.log_pseudodet_W ≈ 6.068425588244111

        end

    end

    @testset "constraint_matrix" begin
        import GMRF.constraint_matrix

        @testset "first-order iGMRF" begin
            F = iGMRF(2, 3, order=1, precision=1.0)

            A = constraint_matrix(F)

            @test A == ones(Int64, 6, 1)
        end

        @testset "second-order iGMRF" begin
            F = iGMRF(3, 3, order=2, precision=1.0)

            A = constraint_matrix(F)

            A_expected = [
                1 1 1
                1 2 1
                1 3 1
                1 1 2
                1 2 2
                1 3 2
                1 1 3
                1 2 3
                1 3 3
            ]

            @test A == A_expected
        end
    end

    @testset "rand" begin
        @testset "first-order iGMRF" begin
            F = iGMRF(3, 3; order=1, precision=1.0)

            rng = MersenneTwister(1234)
            y = rand(rng, F)

            A = constraint_matrix(F)

            @test y isa Vector{Float64}
            @test length(y) == 9
            @test isapprox(A' * y, zeros(size(A, 2)); atol=1e-10) # verifies that the generated field satisfies the intrinsic constraints
        end

        @testset "second-order iGMRF" begin
            F = iGMRF(3, 3; order=2, precision=1.0)

            rng = MersenneTwister(1234)
            y = rand(rng, F)

            A = constraint_matrix(F)

            @test y isa Vector{Float64}
            @test length(y) == 9
            @test isapprox(A' * y, zeros(size(A, 2)); atol=1e-10) # verifies that the generated field satisfies the intrinsic constraints
        end
    end

    @testset "logpdf(::iGMRF)" begin

        @testset "first-order 2 x 2" begin
            F = iGMRF(2, 2; order=1, precision=2.0)

            y = [1.0, -1.0, 0.0, 0.0]
            # For this y, y'Wy = 6.

            # Eigenvalues of W are 0, 2, 2, 4.
            # Hence log_pseudodet(W) = log(16), rank deficiency is 1, and r = 3.
            expected =
                -0.5 * 3 * log(2π) +
                0.5 * 3 * log(2.0) +
                0.5 * log(16.0) -
                0.5 * 2.0 * 6.0

            @test logpdf(F, y) ≈ expected
        end

        @testset "second-order 3 x 3 finite value" begin
            F = iGMRF(3, 3; order=2, precision=2.)

            # Eigenvalues of W are 0, 0, 0, 2, 6, 6, 12, 12, 30
            # Hence log_pseudodet(W) = log(311040), rank deficiency is 3, and r = 6. 

            y = [1., 0., 1., 0., 1., 0., 1., 0., 1.]
            # For this y, y'Wy = 56.

            expected =
                -0.5 * 6 * log(2π) +
                0.5 * 6 * log(2.0) +
                0.5 * log(311040) -
                0.5 * 2.0 * 56.

            @test logpdf(F, y) ≈ expected
        end

        @testset "dimension mismatch" begin
            F = iGMRF(2, 2; order=1, precision=1.0)
            @test_throws DimensionMismatch logpdf(F, zeros(3))
        end
    end

    @testset "full_conditional_canonical_parameters" begin

        @testset "first-order 2 x 2 lattice" begin
            F = iGMRF(2, 2; order=1, precision=2.0)

            y = [1.0, 2.0, 3.0, 4.0]

            h, Q = GMRF.full_conditional_canonical_parameters(F, y)

            # For this y, W̄ * y = [-5, -5, -5, -5].
            # With κ = 2, h = -κ * W̄ * y = [10, 10, 10, 10].
            # Also Q = κ * diag(W) = [4, 4, 4, 4].

            @test h == [10.0, 10.0, 10.0, 10.0]
            @test Q == [4.0, 4.0, 4.0, 4.0]
        end

        @testset "dimension mismatch" begin
            F = iGMRF(2, 2; order=1, precision=1.0)
            @test_throws DimensionMismatch GMRF.full_conditional_canonical_parameters(F, zeros(3))
        end
    end

    @testset "full_conditionals" begin

        F = iGMRF(2, 2; order=1, precision=2.0)

        y = [1.0, 2.0, 3.0, 4.0]

        pd = full_conditionals(F, y)
        h, Q = full_conditional_canonical_parameters(F, y)

        @test length(pd) == 4
        @test all(pd .== NormalCanon.(h, Q))

    end

end

