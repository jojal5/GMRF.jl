@testset "graphstructure.jl" begin
    gs = GMRF.GraphStructure(2, [[2], [1]])
    io = IOBuffer()

    @testset "showGraphStructure(io, obj; prefix)" begin
        # does not throw
        @test_logs GMRF.showGraphStructure(io, gs, prefix = "\t")

    end

    @testset "Base.show(io, obj)" begin
        # does not throw
        @test_logs Base.show(io, gs)

    end

end
