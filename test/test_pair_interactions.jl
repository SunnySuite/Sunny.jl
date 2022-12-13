@testitem "Pair Interactions" begin
    include("test_shared.jl")

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
            @test isapprox(int, reparsed_int; atol=1e-3)
        end
    end

    test_quadratic_shown_info()

    # Checks that a list of bonds contains no equivalent bonds
    function no_equivalent_bonds(bonds::Vector{Bond})
        for (i, bond) in enumerate(bonds[1:end-1])
            equiv_bond = Bond(bond.j, bond.i, -bond.n)
            rest_of_bonds = bonds[i+1:end]
            if (bond ∈ rest_of_bonds) || (equiv_bond ∈ rest_of_bonds)
                return false
            end
        end
        return true
    end

    # Test that iterating over bonds on specific sublattices both only gives
    #  bonds of that sublattice, but also that we hit all bonds if we iterate
    #  over each sublattice.
    function correct_sublat_iteration(bondtable::Sunny.BondTable)
        total_iterated = 0
        for b in 1:nbasis(bondtable)
            for (bond, _) in Sunny.sublat_bonds(bondtable, b)
                if !(bond.i == b || bond.j == b)
                    return false
                end
                total_iterated += 1
            end
        end
        return total_iterated == length(bondtable)
    end

    # Test bond table
    cryst = Sunny.diamond_crystal()
    latsize = (4, 4, 4)
    exchange_ints = diamond_test_exchanges()
    sys = SpinSystem(cryst, exchange_ints, (4,4,4))
    ℋ = sys.hamiltonian
    pair_ints = [ℋ.heisenbergs..., ℋ.diag_coups..., ℋ.gen_coups...]

    for pair_int in pair_ints
        bondtable = pair_int.bondtable
        @test no_equivalent_bonds(bondtable.culled_bonds)
        @test correct_sublat_iteration(bondtable)
    end

end
