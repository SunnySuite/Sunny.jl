# Partition every nonzero bound into one of two sets
function bond_parity(bond)
    bond_delta = (bond.j - bond.i, bond.n...)
    @assert bond_delta != (0, 0, 0, 0)
    return bond_delta > (0, 0, 0, 0)
end

"""
    set_biquadratic!(sys::System, J, bond::Bond)

Sets a scalar biquadratic interaction along `bond`, yielding a pairwise energy
``J (𝐒_i⋅𝐒_j)²``. This interaction will be propagated to equivalent bonds in
consistency with crystal symmetry. Any previous biquadratic exchange
interactions on these bonds will be overwritten.

For systems restricted to dipoles, the biquadratic interactions will
automatically be renormalized to achieve maximum consistency with the more
variationally accurate SU(_N_) mode. This renormalization introduces also a
correction to the quadratic part of the exchange.

See also [`set_exchange!`](@ref).
"""
function set_biquadratic!(sys::System{N}, J, bond::Bond) where N
    # If `sys` has been reshaped, then operate first on `sys.origin`, which
    # contains full symmetry information.
    if !isnothing(sys.origin)
        set_biquadratic!(sys.origin, J, bond)
        transfer_unit_cell!(sys, sys.origin)
        return
    end

    # Biquadratic interactions not yet supported in SU(N) mode
    if !iszero(J) && sys.mode==:SUN
        error("Biquadratic interactions not yet supported in SU(N) mode.")
    end

    if !is_homogeneous(sys)
        error("Use `set_exchange_at!` for inhomogeneous systems.")
    end

    ints = interactions(sys)

    # Verify bond indices
    if bond.i == bond.j && iszero(bond.n)
        error("Bond must connect different sites.")
    end
    (1 <= bond.i <= nbasis(sys.crystal)) || error("Atom index $(bond.i) is out of range.")
    (1 <= bond.j <= nbasis(sys.crystal)) || error("Atom index $(bond.j) is out of range.")

    # Print a warning if an interaction already exists for bond
    if any(x -> x.bond == bond, ints[bond.i].biquad)
        println("Warning: Overriding biquadratic interaction for bond $bond.")
    end

    for i in 1:nbasis(sys.crystal)
        for bond′ in all_symmetry_related_bonds_for_atom(sys.crystal, i, bond)
            # Remove any existing exchange for bond′
            matches_bond(c) = c.bond == bond′
            filter!(!matches_bond, ints[i].biquad)

            # The energy or force calculation only needs to see each bond once
            isculled = bond_parity(bond′)

            # Add to list
            coupling = Coupling(isculled, bond′, J)
            push!(ints[i].biquad, coupling)
        end

        # Sort interactions so that non-culled bonds appear first
        sort!(ints[i].biquad, by=c->c.isculled)
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

See also [`set_biquadratic!`](@ref), [`dmvec`](@ref).
"""
function set_exchange!(sys::System{N}, J, bond::Bond) where N
    # If `sys` has been shaped, then operate first on `sys.origin`, which
    # contains full symmetry information.
    if !isnothing(sys.origin)
        set_exchange!(sys.origin, J, bond)
        transfer_unit_cell!(sys, sys.origin)
        return
    end

    if !is_homogeneous(sys)
        error("Use `set_exchange_at!` for inhomogeneous systems.")
    end
    ints = interactions(sys)

    # Verify bond indices
    if bond.i == bond.j && iszero(bond.n)
        error("Bond must connect different sites.")
    end
    (1 <= bond.i <= nbasis(sys.crystal)) || error("Atom index $(bond.i) is out of range.")
    (1 <= bond.j <= nbasis(sys.crystal)) || error("Atom index $(bond.j) is out of range.")

    # Convert J to Mat3
    J = Mat3(J isa Number ? J*I : J)

    # Verify that exchange is symmetry-consistent
    if !is_coupling_valid(sys.crystal, bond, J)
        println("Symmetry-violating exchange: $J.")
        println("Use `print_bond(crystal, $bond)` for more information.")
        error("Interaction violates symmetry.")
    end

    # Print a warning if an interaction already exists for bond
    if any(x -> x.bond == bond, vcat(ints[bond.i].heisen, ints[bond.i].exchange))
        println("Warning: Overriding exchange for bond $bond.")
    end

    is_heisenberg = isapprox(diagm([J[1,1],J[1,1],J[1,1]]), J; atol=1e-12)

    for i in 1:nbasis(sys.crystal)
        for (bond′, J′) in zip(all_symmetry_related_couplings_for_atom(sys.crystal, i, bond, J)...)
            # Remove any existing exchange for bond′
            matches_bond(c) = c.bond == bond′
            filter!(!matches_bond, ints[i].heisen)
            filter!(!matches_bond, ints[i].exchange)

            # The energy or force calculation only needs to see each bond once
            isculled = bond_parity(bond′)

            # Add to list
            if is_heisenberg
                coupling = Coupling(isculled, bond′, J′[1,1])
                push!(ints[i].heisen, coupling)
            else
                coupling = Coupling(isculled, bond′, J′)
                push!(ints[i].exchange, coupling)
            end
        end

        # Sort interactions so that non-culled bonds appear first
        sort!(ints[i].heisen, by=c->c.isculled)
        sort!(ints[i].exchange, by=c->c.isculled)
    end
end


"""
    dmvec(D)

Antisymmetric matrix representation of the Dzyaloshinskii-Moriya pseudo-vector,

```
  [  0    D[3] -D[2]
   -D[3]   0    D[1]
    D[2] -D[1]   0  ]
```

Useful in the context of [`set_exchange!`](@ref).
"""
function dmvec(D)
    D = Vec3(D)
    return SA[  0.0  D[3] -D[2]
              -D[3]   0.0  D[1]
               D[2] -D[1]   0.0 ]
end
