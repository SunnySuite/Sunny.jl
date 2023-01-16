
function PairExchanges()
    return PairExchanges(
        Tuple{Bool, Bond, Float64}[],
        Tuple{Bool, Bond, Mat3}[],
        Tuple{Bool, Bond, Float64}[],
    )
end


function propagate_exchange!(hamiltonian::Interactions, cryst::Crystal, J::Mat3, J_biq::Float64, bond::Bond)
    # Verify that atom indices are in range
    (1 <= bond.i <= nbasis(cryst)) || error("Atom index $(bond.i) is out of range.")
    (1 <= bond.j <= nbasis(cryst)) || error("Atom index $(bond.j) is out of range.")

    # Verify that exchange is symmetry-consistent
    if !is_coupling_valid(cryst, bond, J)
        println("Symmetry-violating exchange: $J.")
        println("Use `print_bond(crystal, $bond)` for more information.")
        error("Interaction violates symmetry.")
    end

    # Biquadratic interactions not yet supported in SU(N) mode
    N = size(hamiltonian.anisos[1].matrep, 1)
    if !iszero(J_biq) && N > 0
        error("Biquadratic interactions not yet supported in SU(N) mode.")
    end

    # Print a warning if an interaction already exists for bond
    (; heisen, quadmat, biquad) = hamiltonian.pairexch[bond.i]
    if any(x -> x[2] == bond, vcat(heisen, quadmat, biquad))
        println("Warning: Overriding exchange for bond $bond.")
    end

    isheisen = isapprox(diagm([J[1,1],J[1,1],J[1,1]]), J; atol=1e-12)

    for (i, (; heisen, quadmat, biquad)) in enumerate(hamiltonian.pairexch)
        for (bond′, J′) in zip(all_symmetry_related_couplings_for_atom(cryst, i, bond, J)...)
            # Remove any existing interactions for bond′
            matches_bond(x) = x[2] == bond′
            filter!(!matches_bond, heisen)
            filter!(!matches_bond, quadmat)
            filter!(!matches_bond, biquad)

            # The energy or force calculation only needs to see each bond once.
            # Use a trick below to select half the bonds to cull.
            bond_delta = (bond′.j - bond′.i, bond′.n...)
            @assert bond_delta != (0, 0, 0, 0)
            isculled = bond_delta > (0, 0, 0, 0)

            if isheisen
                @assert J ≈ J′
                push!(heisen, (isculled, bond′, J′[1,1]))
            else
                push!(quadmat, (isculled, bond′, J′))
            end

            if !iszero(J_biq)
                push!(biquad, (isculled, bond′, J_biq))
            end
        end

        # Sort interactions so that non-culled bonds appear first
        sort!(heisen, by=first)
        sort!(quadmat, by=first)
        sort!(biquad, by=first)
    end
end

