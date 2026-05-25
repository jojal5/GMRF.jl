@testset "gridstructure.jl" begin

    @test_set "constructor" begin
        
        @test_set "valid arguments" begin
            G = GridStructure(3, 4)
            @test G.m₁ == 3
            @test G.m₂ == 4
            @test typeof(G.m₁) == Int64
            @test typeof(G.m₂) == Int64
        end

        @testset "integer conversion" begin
            G = GridStructure(Int32(3), Int32(4))

            @test G.m₁ == 3
            @test G.m₂ == 4
            @test typeof(G.m₁) == Int64
            @test typeof(G.m₂) == Int64
        end

        @testset "invalid dimensions" begin
            @test_throws ArgumentError GridStructure(0, 4)
            @test_throws ArgumentError GridStructure(3, 0)
            @test_throws ArgumentError GridStructure(-1, 4)
            @test_throws ArgumentError GridStructure(3, -1)
        end
    end

    @test_set "show function" begin
        G = GridStructure(3, 4)
        io = IOBuffer()
        # does not throw
        @test_logs GMRF.showGridStructure(io, G)
    end

    @test_set "first_order_igmrf_structure_matrix()" begin
        import GRMF.first_order_igmrf_structure_matrix
        import GMRF.GridStructure

        @test_set "line 1 x 3" begin
            W_expected = sparse([1 -1 0; -1 2 -1; 0 -1 1])
            W = first_order_igmrf_structure_matrix(GridStructure(1, 3))
            @test W == W_expected
        end

        @testset "3 x 2 lattice" begin
            W_expected = sparse([
                2 -1 0 -1 0 0
                -1 3 -1 0 -1 0
                0 -1 2 0 0 -1
                -1 0 0 2 -1 0
                0 -1 0 -1 3 -1
                0 0 -1 0 -1 2
            ])

            W = first_order_igmrf_structure_matrix(GridStructure(3, 2))
            @test W == W_expected
        end

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

    @testset "first_order_conditional_independent_subsets" begin
        subsets = GMRF.first_order_conditional_independent_subsets(3, 3)

        @test subsets == [
            [1, 3, 5, 7, 9],
            [2, 4, 6, 8],
        ]
    end

    @testset "second_order_lattice_neighbors" begin
        import GMRF.second_order_lattice_neighbors

        @testset "2 x 2 lattice" begin
            nbs, W = second_order_lattice_neighbors(2, 2)

            @test nbs == [[2, 3], [1, 4], [1, 4], [2, 3]]

            W_expected = sparse([
                2 -2 -2 2
                -2 2 2 -2
                -2 2 2 -2
                2 -2 -2 2
            ])

            @test W == W_expected
        end

        @testset "3 x 2 lattice" begin
            nbs, W = second_order_lattice_neighbors(3, 2)

            @test nbs == [
                [2, 4],
                [1, 3, 5],
                [2, 6],
                [1, 5],
                [2, 4, 6],
                [3, 5],
            ]

            W_expected = sparse([
                3 -4 1 -2 2 0
                -4 8 -4 2 -4 2
                1 -4 3 0 2 -2
                -2 2 0 3 -4 1
                2 -4 2 -4 8 -4
                0 2 -2 1 -4 3
            ])

            @test W == W_expected
        end

        @testset "3 x 3 lattice" begin
            nbs, W = second_order_lattice_neighbors(3, 3)

            @test nbs == [
                [2, 4],
                [1, 3, 5],
                [2, 6],
                [1, 5, 7],
                [2, 4, 6, 8],
                [3, 5, 9],
                [4, 8],
                [5, 7, 9],
                [6, 8],
            ]

            W_expected = sparse([
                4 -4 1 -4 2 0 1 0 0
                -4 9 -4 2 -6 2 0 1 0
                1 -4 4 0 2 -4 0 0 1
                -4 2 0 9 -6 1 -4 2 0
                2 -6 2 -6 16 -6 2 -6 2
                0 2 -4 1 -6 9 0 2 -4
                1 0 0 -4 2 0 4 -4 1
                0 1 0 2 -6 2 -4 9 -4
                0 0 1 0 2 -4 1 -4 4
            ])

            @test W == W_expected
        end

        @testset "invalid dimensions" begin
            @test_throws ArgumentError second_order_lattice_neighbors(0, 3)
            @test_throws ArgumentError second_order_lattice_neighbors(3, 0)
            @test_throws ArgumentError second_order_lattice_neighbors(-1, 3)
            @test_throws ArgumentError second_order_lattice_neighbors(3, -1)
        end
    end

    @testset "second_order_conditional_independent_subsets" begin
        subsets = GMRF.second_order_conditional_independent_subsets(3, 3)

        @test subsets == [
            [1, 8],
            [4],
            [2, 7],
            [5],
            [3],
            [6],
            [9],
        ]
    end

end



@testset "GridStructure constructor" begin
    import GMRF.GridStructure

    @testset "first-order 3 x 3 lattice" begin
        G = GridStructure(3, 3; order=1)

        @test G.grid_size == (3, 3)

        @test G.neighbors == [
            [2, 4],
            [1, 3, 5],
            [2, 6],
            [1, 5, 7],
            [2, 4, 6, 8],
            [3, 5, 9],
            [4, 8],
            [5, 7, 9],
            [6, 8],
        ]

        @test G.cond_ind_subset == [
            [1, 3, 5, 7, 9],
            [2, 4, 6, 8],
        ]

        W_expected = sparse([
            2 -1 0 -1 0 0 0 0 0
            -1 3 -1 0 -1 0 0 0 0
            0 -1 2 0 0 -1 0 0 0
            -1 0 0 3 -1 0 -1 0 0
            0 -1 0 -1 4 -1 0 -1 0
            0 0 -1 0 -1 3 0 0 -1
            0 0 0 -1 0 0 2 -1 0
            0 0 0 0 -1 0 -1 3 -1
            0 0 0 0 0 -1 0 -1 2
        ])

        @test G.W == W_expected
        @test G.W̄ == W_expected - spdiagm(0 => diag(W_expected))
    end

    @testset "second-order 3 x 3 lattice" begin
        G = GridStructure(3, 3; order=2)

        @test G.grid_size == (3, 3)

        @test G.neighbors == [
            [2, 4],
            [1, 3, 5],
            [2, 6],
            [1, 5, 7],
            [2, 4, 6, 8],
            [3, 5, 9],
            [4, 8],
            [5, 7, 9],
            [6, 8],
        ]

        @test G.cond_ind_subset == [
            [1, 8],
            [4],
            [2, 7],
            [5],
            [3],
            [6],
            [9],
        ]

        W_expected = sparse([
            4 -4 1 -4 2 0 1 0 0
            -4 9 -4 2 -6 2 0 1 0
            1 -4 4 0 2 -4 0 0 1
            -4 2 0 9 -6 1 -4 2 0
            2 -6 2 -6 16 -6 2 -6 2
            0 2 -4 1 -6 9 0 2 -4
            1 0 0 -4 2 0 4 -4 1
            0 1 0 2 -6 2 -4 9 -4
            0 0 1 0 2 -4 1 -4 4
        ])

        @test G.W == W_expected
        @test G.W̄ == W_expected - spdiagm(0 => diag(W_expected))
    end

    @testset "invalid arguments" begin
        @test_throws ArgumentError GridStructure(0, 3; order=1)
        @test_throws ArgumentError GridStructure(3, 0; order=1)
        @test_throws ArgumentError GridStructure(3, 3; order=0)
    end
end