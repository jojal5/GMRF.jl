@testset "igmrf.jl" begin
    @testset "Base.show(io, obj)" begin
        # does not throw
        G = GMRF.GridStructure((1,2), [[2], [1]], [[1],[2]], spzeros(1, 2), spzeros(1, 2))
        igmrf = iGMRF(G, 1, 1.0)
        io = IOBuffer()
        @test_logs Base.show(io, igmrf)

    end

    @testset "iGMRF(m₁, m₂, order, κ)" begin
        # order not 1 or 2 throws
        @test_throws AssertionError iGMRF(1, 1, 0, 1.0)

        # Simple grid of order 1
        igmrf = iGMRF(1, 1, 1, 1.0)

        @test igmrf.G.gridSize == (1, 1)
        @test igmrf.rankDeficiency == 1
        @test igmrf.κ ≈ 1.0
        # TODO : Test that W̄ was derived from W

        # Simple grid of order 2
        igmrf = iGMRF(1, 1, 2, 1.0)

        @test igmrf.G.gridSize == (1, 1)
        @test igmrf.rankDeficiency == 3
        @test igmrf.κ ≈ 1.0
        # TODO : Test that W̄ was derived from W

    end

    @testset "fo_nbs(m₁, m₂)" begin
        # Grid 1 x 1
        nbs, W = GMRF.fo_nbs(1, 1)

        @test nbs == [[]]
        # TODO : Test W

        # Grid 2 x 2
        nbs, W = GMRF.fo_nbs(2, 2)

        @test nbs == [[2, 3], [1, 4], [1, 4], [2, 3]]
        # TODO : Test W

    end

    @testset "so_nbs(m₁, m₂)" begin
        # Grid 1 x 1
        nbs, W = GMRF.so_nbs(1, 1)

        @test nbs == [[]]
        # TODO : Test W

        # Grid 2 x 2
        nbs, W = GMRF.so_nbs(2, 2)

        # TODO : Test nbs
        # TODO : Test W

    end

    @testset "fo_condindsubsets(m₁, m₂)" begin
        # Grid 1 x 1
        cond = GMRF.fo_condindsubsets(1, 1)

        # TODO : Test cond

        # Grid 2 x 2
        cond = GMRF.fo_condindsubsets(2, 2)

        # TODO : Test cond

    end

    @testset "so_condindsubsets(m₁, m₂)" begin
        # Grid 1 x 1
        cond = GMRF.so_condindsubsets(1, 1)

        # TODO : Test cond

        # Grid 2 x 2
        cond = GMRF.so_condindsubsets(2, 2)

        # TODO : Test cond

    end

    @testset "rand(F)" begin
        # rankdeficiency != 1 or 3 throws
        G = GMRF.GridStructure((1,2), [[2], [1]], [[1],[2]], spzeros(1, 2), spzeros(1, 2))
        igmrf = iGMRF(G, 2, 1.0)
        @test_throws AssertionError rand(igmrf)

        # returns plausible data (rankdeficiency == 1)
        # TODO: Test if all insupport ?

        # returns plausible data (rankdeficiency == 3)
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
