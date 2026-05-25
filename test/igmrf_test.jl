
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

end












