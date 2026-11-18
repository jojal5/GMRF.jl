@testset "log_pseudodet" begin
    import GMRF.log_pseudodet

    @testset "identity matrix" begin
        W = sparse(I, 3, 3)
        @test log_pseudodet(W, 0) ≈ 0.0
    end

    @testset "diagonal matrix" begin
        W = spdiagm(0 => [1,2,3])

        @test log_pseudodet(W, 0) ≈ log(1.0) + log(2.0) + log(3.0)
        @test log_pseudodet(W, 1) ≈ log(2.0) + log(3.0)
        @test log_pseudodet(W, 2) ≈ log(3.0)
    end

    @testset "invalid arguments" begin
        W = sparse(I, 3, 3)

        @test_throws ArgumentError log_pseudodet(W, -1)
        @test_throws ArgumentError log_pseudodet(W, 3)
        @test_throws DimensionMismatch log_pseudodet(spzeros(2, 3), 0)
    end
end