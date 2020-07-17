@testset "gridstructure.jl" begin
    gs = GMRF.GridStructure((1,2), [[2], [1]], [[1],[2]], spzeros(1, 2), spzeros(1, 2))
    io = IOBuffer()

    @testset "showGridStructure(io, obj; prefix)" begin
        # does not throw
        @test_logs GMRF.showGridStructure(io, gs, prefix = "\t")

    end

    @testset "Base.show(io, obj)" begin
        # does not throw
        @test_logs Base.show(io, gs)

    end

end
