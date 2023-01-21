
function PairExchanges()
    return PairExchanges(
        Tuple{Bool, Bond, Float64}[],
        Tuple{Bool, Bond, Mat3}[],
        Tuple{Bool, Bond, Float64}[],
    )
end


"""
    set_exchange_with_biquadratic!(sys::System, J, J′, bond::Bond)

Sets both quadratic and biquadratic exchange interactions along `bond`, yielding
a pairwise energy ``𝐒_i⋅J 𝐒_j + J′ (𝐒_i⋅𝐒_j)²``. These interactions will be
propagated to equivalent bonds in consistency with crystal symmetry. Any
previous exchange interactions on these bonds will be overwritten.

For systems with `mode=:projected` the biquadratic interactions will
automatically be renormalized to achieve maximum consistency with the more
variationally accurate SU(_N_) mode. This renormalization introduces a
correction to the quadratic part of the exchange, which is why the two parts
must be specified concurrently.

See also [`set_exchange!`](@ref).
"""
function set_exchange_with_biquadratic!(sys::System{N}, J, J_biq, bond::Bond) where N
    if bond.i == bond.j && iszero(bond.n)
        error("Exchange interactions must connect different sites.")
    end
    J = Mat3(J isa Number ? J*I : J)

    (; crystal) = sys
    (; pairexch) = sys.interactions

    # Verify that atom indices are in range
    (1 <= bond.i <= nbasis(crystal)) || error("Atom index $(bond.i) is out of range.")
    (1 <= bond.j <= nbasis(crystal)) || error("Atom index $(bond.j) is out of range.")

    # Verify that exchange is symmetry-consistent
    if !is_coupling_valid(crystal, bond, J)
        println("Symmetry-violating exchange: $J.")
        println("Use `print_bond(crystal, $bond)` for more information.")
        error("Interaction violates symmetry.")
    end

    # Biquadratic interactions not yet supported in SU(N) mode
    if !iszero(J_biq) && sys.mode==:SUN
        error("Biquadratic interactions not yet supported in SU(N) mode.")
    end

    # Print a warning if an interaction already exists for bond
    (; heisen, quadmat, biquad) = pairexch[bond.i]
    if any(x -> x[2] == bond, vcat(heisen, quadmat, biquad))
        println("Warning: Overriding exchange for bond $bond.")
    end

    isheisen = isapprox(diagm([J[1,1],J[1,1],J[1,1]]), J; atol=1e-12)

    for (i, (; heisen, quadmat, biquad)) in enumerate(pairexch)
        for (bond′, J′) in zip(all_symmetry_related_couplings_for_atom(crystal, i, bond, J)...)
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


"""
    set_exchange!(sys::System, J, bond::Bond)

Sets a 3×3 spin-exchange matrix `J` along `bond`, yielding a pairwise
interaction energy ``𝐒_i⋅J 𝐒_j``. This interaction will be propagated to
equivalent bonds in consistency with crystal symmetry. Any previous exchange
interactions on these bonds will be overwritten. The parameter `bond` has the
form `Bond(i, j, offset)`, where `i` and `j` are atom indices within the unit
cell, and `offset` is a displacement in unit cells.

Scalar `J` implies a pure Heisenberg exchange.

As a convenience, `dmvec(D)` can be used to construct the antisymmetric part of
the exchange, where `D` is the Dzyaloshinskii-Moriya pseudo-vector. The
resulting interaction will be ``𝐃⋅(𝐒_i×𝐒_j)``.

# Examples
```julia
using Sunny, LinearAlgebra

# An explicit exchange matrix
J1 = [2 3 0;
     -3 2 0;
      0 0 2]
set_exchange!(sys, J1, bond)

# An equivalent Heisenberg + DM exchange 
J2 = 2*I + dmvec([0,0,3])
set_exchange!(sys, J2, bond)
```

See also [`set_exchange_with_biquadratic!`](@ref), [`dmvec`](@ref).
"""
function set_exchange!(sys::System{N}, J, bond::Bond) where N
    set_exchange_with_biquadratic!(sys, J, 0.0, bond)
end


"""
    dmvec(D)

Representation of the Dzyaloshinskii-Moriya interaction pseudo-vector `D` as
an antisymmetric matrix,

```
  |  0   D₃ -D₂ |
  | -D₃  0   D₁ |
  |  D₂ -D₁  0  |
```

Useful in the context of [`set_exchange!`](@ref).
"""
function dmvec(D)
   SA[ 0  D[3] -D[2];
   -D[3]     0  D[1];
    D[2] -D[1]    0]
end
