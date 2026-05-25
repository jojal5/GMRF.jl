@testset "igmrf.jl" begin
    @testset "Base.show(io, obj)" begin
        # does not throw
        G = GMRF.GridStructure((1,2), [[2], [1]], [[1],[2]], spzeros(1, 2), spzeros(1, 2))
        igmrf = iGMRF(G, 1, 1.0, 0.)
        io = IOBuffer()
        @test_logs Base.show(io, igmrf)

    end

    @testset "iGMRF(m₁, m₂, order, κ, log_D)" begin
        # order not 1 or 2 throws
        @test_throws ArgumentError iGMRF(1, 1, 0, 0.)

        # Simple grid of order 1
        igmrf = iGMRF(2, 2, 1, 1.)

        @test igmrf.G.grid_size == (2, 2)
        @test igmrf.rank_deficiency == 1
        @test igmrf.κ ≈ 1.0
        @test igmrf.log_pseudodet_W ≈ 2.7725887222397803
        # TODO : Test that W̄ was derived from W

        # Simple grid of order 2
        igmrf = iGMRF(3, 3, 2, 1.)

        @test igmrf.G.grid_size == (3, 3)
        @test igmrf.rank_deficiency == 3
        @test igmrf.κ ≈ 1.0
        @test igmrf.log_pseudodet_W ≈ 12.647676800254212
        # TODO : Test that W̄ was derived from W

    end

    @testset "constraint_matrix" begin
        import GMRF.constraint_matrix

    @testset "first-order iGMRF" begin
        F = iGMRF(2, 3, 1, 1.0)

        A = constraint_matrix(F)

        @test size(A) == (6, 1)
        @test A == ones(6, 1)
    end

    @testset "second-order iGMRF" begin
        F = iGMRF(3, 3, 2, 1.0)

        A = constraint_matrix(F)

        A_expected = [
            1.0  1.0  1.0
            1.0  2.0  1.0
            1.0  3.0  1.0
            1.0  1.0  2.0
            1.0  2.0  2.0
            1.0  3.0  2.0
            1.0  1.0  3.0
            1.0  2.0  3.0
            1.0  3.0  3.0
        ]

        @test size(A) == (9, 3)
        @test A == A_expected
    end
end

    @testset "first_order_lattice_neighbors(m₁, m₂)" begin
        # Grid 1 x 1
        nbs, W = GMRF.first_order_lattice_neighbors(1, 1)

        @test nbs == [[]]
        # TODO : Test W

        # Grid 2 x 2
        nbs, W = GMRF.first_order_lattice_neighbors(2, 2)

        @test nbs == [[2, 3], [1, 4], [1, 4], [2, 3]]
        # TODO : Test W

    end

    @testset "second_order_lattice_neighbors(m₁, m₂)" begin
        # Grid 1 x 1
        nbs, W = GMRF.second_order_lattice_neighbors(1, 1)

        @test nbs == [[]]
        # TODO : Test W

        # Grid 2 x 2
        nbs, W = GMRF.second_order_lattice_neighbors(2, 2)

        # TODO : Test nbs
        # TODO : Test W

    end

    @testset "first_order_conditional_independent_subsets(m₁, m₂)" begin
        # Grid 1 x 1
        cond = GMRF.first_order_conditional_independent_subsets(1, 1)

        # TODO : Test cond

        # Grid 2 x 2
        cond = GMRF.first_order_conditional_independent_subsets(2, 2)

        # TODO : Test cond

    end

    @testset "second_order_conditional_independent_subsets(m₁, m₂)" begin
        # Grid 1 x 1
        cond = GMRF.second_order_conditional_independent_subsets(1, 1)

        # TODO : Test cond

        # Grid 2 x 2
        cond = GMRF.second_order_conditional_independent_subsets(2, 2)

        # TODO : Test cond

    end

    @testset "rand(F)" begin
        # returns plausible data (rank_deficiency == 1)
        # TODO: Test if all insupport ?

        # returns plausible data (rank_deficiency == 3)
        # TODO: Test if all insupport ?

    end

    @testset "logpdf(F, y)" begin
        # TODO : Test with known values

    end

    @testset "fullconditionals(F, y)" begin
        # TODO : Test with known values

    end

    @testset "fullcondlogpdf(F, y)" begin
        # TODO : Test with known values

    end

    @testset "getconditional(F, B, x)" begin
        # TODO : Test with known values

    end

end
