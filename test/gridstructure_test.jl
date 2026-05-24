@testset "gridstructure.jl" begin
    gs = GMRF.GridStructure((1, 2), [[2], [1]], [[1], [2]], spzeros(1, 2), spzeros(1, 2))
    io = IOBuffer()

    @testset "showGridStructure(io, obj; prefix)" begin
        # does not throw
        @test_logs GMRF.showGridStructure(io, gs, prefix="\t")

    end

    @testset "Base.show(io, obj)" begin
        # does not throw
        @test_logs Base.show(io, gs)

    end

    @testset "first_order_lattice_neighbors" begin
        import GMRF.first_order_lattice_neighbors

        @testset "2 x 2 lattice" begin
            nbs, W = first_order_lattice_neighbors(2, 2)

            @test nbs == [[2, 3], [1, 4], [1, 4], [2, 3]]

            W_expected = sparse([
                2 -1 -1 0
                -1 2 0 -1
                -1 0 2 -1
                0 -1 -1 2
            ])

            @test W == W_expected
        end

        @testset "3 x 2 lattice" begin
            nbs, W = first_order_lattice_neighbors(3, 2)

            @test nbs == [
                [2, 4],
                [1, 3, 5],
                [2, 6],
                [1, 5],
                [2, 4, 6],
                [3, 5],
            ]

            W_expected = sparse([
                2 -1 0 -1 0 0
                -1 3 -1 0 -1 0
                0 -1 2 0 0 -1
                -1 0 0 2 -1 0
                0 -1 0 -1 3 -1
                0 0 -1 0 -1 2
            ])

            @test W == W_expected
        end

        @testset "invalid dimensions" begin
            @test_throws ArgumentError first_order_lattice_neighbors(0, 3)
            @test_throws ArgumentError first_order_lattice_neighbors(3, 0)
            @test_throws ArgumentError first_order_lattice_neighbors(-1, 3)
            @test_throws ArgumentError first_order_lattice_neighbors(3, -1)
        end
    end

end
