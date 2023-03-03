function validate_bond(cryst::Crystal, bond::Bond)
    # Verify bond indices
    if bond.i == bond.j && iszero(bond.n)
        error("Bond must connect different atoms.")
    end
    (1 <= bond.i <= natoms(cryst)) || error("Atom index $(bond.i) is out of range.")
    (1 <= bond.j <= natoms(cryst)) || error("Atom index $(bond.j) is out of range.")
end

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

    sys.mode==:SUN && error("Biquadratic interactions not yet supported in SU(N) mode.")
    validate_bond(sys.crystal, bond)
    is_homogeneous(sys) || error("Cannot symmetry-propagate interactions for an inhomogeneous system.")

    ints = interactions_homog(sys)


    # Print a warning if an interaction already exists for bond
    if any(x -> x.bond == bond, ints[bond.i].biquad)
        println("Warning: Overriding biquadratic interaction for bond $bond.")
    end

    for i in 1:natoms(sys.crystal)
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

    validate_bond(sys.crystal, bond)
    is_homogeneous(sys) || error("Cannot symmetry-propagate interactions for an inhomogeneous system.")

    ints = interactions_homog(sys)

    # Convert J to Mat3
    J = Mat3(J isa Number ? J*I : J)
    is_heisenberg = isapprox(diagm([J[1,1],J[1,1],J[1,1]]), J; atol=1e-12)

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

    for i in 1:natoms(sys.crystal)
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


# Given a `bond` for `cryst`, return a corresponding new bond for the reshaped
# `new_cryst`. The new bond will begin at atom `new_i`.
function transform_bond(new_cryst::Crystal, new_i::Int, cryst::Crystal, bond::Bond)
    new_ri = new_cryst.positions[new_i]

    # Verify that new_i (indexed into new_cryst) is consistent with bond.i
    # (indexed into crystal).
    @assert bond.i == position_to_index(cryst, cryst.latvecs \ new_cryst.latvecs * new_ri)

    # Positions in new fractional coordinates
    br = BondRaw(cryst, bond)
    new_rj = new_ri + new_cryst.latvecs \ cryst.latvecs * (br.rj - br.ri)

    # Construct bond using new indexing system
    new_j, new_n = position_to_index_and_offset(new_cryst, new_rj)
    return Bond(new_i, new_j, new_n)
end

# CURRENTLY UNUSED
"""
    bonded_site(sys::System, site::Site, bond::Bond)

Given a [`Site`](@ref) that acts as the first participant in a [`Bond`](@ref),
return the second [`Site`](@ref) participating in the bond. For reshaped
systems, the indices in `bond` refer to the original crystallographic unit cell.
Useful for generating inputs to [`set_exchange_at`](@ref) and
[`set_biquadratic_at`](@ref).

# Example

```julia
# Calculate site indices from a position in fractional coordinates
site1 = position_to_site(sys, r)

# Get the other site that participates in a bond
site2 = bonded_site(sys, site1, bond)

# Use both sites to set an inhomogeneous interaction
set_exchange_at!(sys, J, site1, site2)
```
"""
function bonded_site(sys::System{N}, site, bond::Bond) where N
    site = to_cartesian(site)
    bond′ = transform_bond(sys.crystal, to_atom(site), orig_crystal(sys), bond)
    cell′ = mod1(to_cell(site) .+ bond′.n, sys.latsize)
    return CartesianIndex(cell′[1], cell′[2], cell′[3], bond′.j)
end

# Internal function only
function sites_to_bond(site1::CartesianIndex{4}, site2::CartesianIndex{4})
    n = Tuple(to_cell(site2)) .- Tuple(to_cell(site1))
    return Bond(to_atom(site1), to_atom(site2), n)
end

# Internal function only
function push_coupling!(couplings, bond, J)
    isculled = bond_parity(bond)
    filter!(c -> c.bond != bond, couplings)
    push!(couplings, Coupling(isculled, bond, J))
    sort!(couplings, by=c->c.isculled)
    return
end

"""
    set_biquadratic_at!(sys::System, J, site1::Site, site2::Site)

Sets the scalar biquadratic interaction along the single bond connecting two
[`Site`](@ref)s, ignoring crystal symmetry. The system must support
inhomogeneous interactions via [`to_inhomogeneous`](@ref).

See also [`set_biquadratic!`](@ref).
"""
function set_biquadratic_at!(sys::System{N}, J, site1, site2) where N
    site1 = to_cartesian(site1)
    site2 = to_cartesian(site2)
    bond = sites_to_bond(site1, site2)

    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")
    ints = interactions_inhomog(sys)

    push_coupling!(ints[site1].biquad, bond, J)
    push_coupling!(ints[site2].biquad, reverse(bond), J')
    return
end


"""
    set_exchange_at!(sys::System, J, site1::Site, site2::Site)

Sets the exchange interaction along the single bond connecting two
[`Site`](@ref)s, ignoring crystal symmetry. The system must support
inhomogeneous interactions via [`to_inhomogeneous`](@ref).

See also [`set_exchange!`](@ref).
"""
function set_exchange_at!(sys::System{N}, J, site1, site2) where N
    site1 = to_cartesian(site1)
    site2 = to_cartesian(site2)
    bond = sites_to_bond(site1, site2)

    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")
    ints = interactions_inhomog(sys)

    # Convert J to Mat3
    J = Mat3(J isa Number ? J*I : J)
    is_heisenberg = isapprox(diagm([J[1,1],J[1,1],J[1,1]]), J; atol=1e-12)
    
    if is_heisenberg
        push_coupling!(ints[site1].heisen, bond, J[1,1])
        push_coupling!(ints[site2].heisen, reverse(bond), J[1,1]')
    else
        push_coupling!(ints[site1].exchange, bond, J)
        push_coupling!(ints[site2].exchange, reverse(bond), J')
    end
    return
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
