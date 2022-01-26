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

    function Base.isapprox(int::Sunny.QuadraticInteraction, int2::Sunny.QuadraticInteraction; kwargs...)
        exchange_eq = isapprox(int.J, int2.J; kwargs...)
        bond_eq = int.bond == int2.bond
        exchange_eq && bond_eq
    end

    # Test that the output of `show` on pair interactions actually creates the correct interaction again,
    #   up to truncation error in `show`. (Does not check label on interactions, which isn't printed in `show`)
    function test_quadratic_shown_info()
        io = IOBuffer()
        context = IOContext(io)
        exchange_ints = diamond_test_exchanges()
        for int in exchange_ints
            show(context, "text/plain", int)
            output = String(take!(io))
            reparsed_int = eval(Meta.parse(output))
            @test isapprox(int, reparsed_int; atol=1e-4)
        end
    end

    test_quadratic_shown_info()

end
