@testset "Pair Interactions" begin

    # Test that each exchange gets mapped to the correct backend type
    function test_type_mapping()
        cryst = Sunny.diamond_crystal()
        latsize = (4, 4, 4)
        exchange_ints = diamond_test_exchanges()

        # Note: These tests rely on diamond_test_exchanges outputting
        #  the interactions in the order
        #  [heisenberg, on_site [heisen], diagonal, general]
        correct_types = [
            Sunny.HeisenbergCPU,
            Sunny.HeisenbergCPU,
            Sunny.DiagonalCouplingCPU,
            Sunny.GeneralCouplingCPU
        ]

        for (int, type) in zip(exchange_ints, correct_types)
            sys = SpinSystem(cryst, [int], latsize)
            ℋ = sys.hamiltonian
            if type == Sunny.HeisenbergCPU
                @test length(ℋ.heisenbergs) == 1
            elseif type == Sunny.DiagonalCouplingCPU
                @test length(ℋ.diag_coups) == 1
            else
                @test length(ℋ.gen_coups) == 1
            end
        end
    end

    test_type_mapping()

end
