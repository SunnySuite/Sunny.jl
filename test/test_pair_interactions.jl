#=
@testitem "Pair Interactions" begin
    include("shared.jl")

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
=#