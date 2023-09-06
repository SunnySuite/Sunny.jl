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

# Convert J to Union{Float64, Mat3}. If J is _exactly_ the identity matrix, we
# can compactly represent it using a single float. If J is near (but not
# exactly) the identity matrix, retain the full matrix representation. This
# could hypothetically be important to preserve symmetry breaking effects. For
# example, a user might select J=diagm([a,a,a+ϵ]) for infinitesimal ϵ to favor
# the z direction.
function to_float_or_mat3(J)
    if J isa Number || J == J[1] * I
        J = Float64(first(J))
    else
        J = Mat3(J)
    end
    return J::Union{Float64, Mat3}
end

# Internal function only
function push_coupling!(sys, couplings, bond, bilin, biquad, large_S)
    # Perform renormalization derived in https://arxiv.org/abs/2304.03874
    if sys.mode == :dipole && !large_S
        # Multiplicative average is the correct result for two sites with
        # different S (derivation omitted).
        S1 = (sys.Ns[bond.i]-1)/2
        S2 = (sys.Ns[bond.j]-1)/2
        S = sqrt(S1*S2)
        r = (1 - 1/S + 1/4S^2)
    
        # J_bq (s⋅sⱼ)^2 -> J_bq (r (sᵢ⋅sⱼ)^2 - (sᵢ⋅sⱼ)/2 + S^3 + S^2/4)
        bilin = bilin - (biquad/2) * one(bilin)
        biquad = r * biquad
        
        # Drop the constant shift, J_bq (S^3 + S^2/4).
    end

    # Remove previous coupling on this bond
    filter!(c -> c.bond != bond, couplings)

    # If the new coupling is exactly zero, return early
    iszero(bilin) && iszero(biquad) && return

    # Otherwise, add the new coupling to the list
    isculled = bond_parity(bond)
    push!(couplings, PairCoupling(isculled, bond, bilin, biquad))
    sort!(couplings, by=c->c.isculled)
    
    return
end

"""
    set_exchange!(sys::System, J, bond::Bond; biquad=0, large_S=false)

Sets a 3×3 spin-exchange matrix `J` along `bond`, yielding a pairwise
interaction energy ``𝐒_i⋅J 𝐒_j``. This interaction will be propagated to
equivalent bonds in consistency with crystal symmetry. Any previous exchange
interactions on these bonds will be overwritten. The parameter `bond` has the
form `Bond(i, j, offset)`, where `i` and `j` are atom indices within the unit
cell, and `offset` is a displacement in unit cells.

The parameter `J` may be scalar or matrix-valued. As a convenience, `dmvec(D)`
can be used to construct the antisymmetric part of the exchange, where `D` is
the Dzyaloshinskii-Moriya pseudo-vector. The resulting interaction will be
``𝐃⋅(𝐒_i×𝐒_j)``.

The optional parameter `biquad` defines the strength ``b`` for scalar
biquadratic interactions of the form ``b (𝐒_i⋅𝐒_j)²`` For systems restricted
to dipoles, ``b`` will be automatically renormalized for maximum consistency
with the more variationally accurate SU(_N_) mode. Set `large_S=true` to work in
the large-``S`` limit and disable this renormalization.

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
"""
function set_exchange!(sys::System{N}, J, bond::Bond; biquad=0, large_S=false) where N
    is_homogeneous(sys) || error("Use `set_exchange_at!` for an inhomogeneous system.")
    ints = interactions_homog(sys)

    # If `sys` has been reshaped, then operate first on `sys.origin`, which
    # contains full symmetry information.
    if !isnothing(sys.origin)
        set_exchange!(sys.origin, J, bond; biquad)
        set_interactions_from_origin!(sys)
        return
    end

    validate_bond(sys.crystal, bond)

    J = to_float_or_mat3(J)

    # Verify that exchange is symmetry-consistent
    if !is_coupling_valid(sys.crystal, bond, J)
        @error """Symmetry-violating exchange: $J.
                  Use `print_bond(crystal, $bond)` for more information."""
        error("Interaction violates symmetry.")
    end

    # Print a warning if an interaction already exists for bond
    if any(x -> x.bond == bond, ints[bond.i].pair)
        warn_coupling_override("Overriding coupling for $bond.")
    end

    for i in 1:natoms(sys.crystal)
        bonds, Js = all_symmetry_related_couplings_for_atom(sys.crystal, i, bond, J)
        for (bond′, J′) in zip(bonds, Js)
            push_coupling!(sys, ints[i].pair, bond′, J′, biquad, large_S)
        end
    end
end

# Converts two sites to a bond with indices for possibly reshaped unit cell. For
# internal use only.
function sites_to_internal_bond(sys::System{N}, site1::CartesianIndex{4}, site2::CartesianIndex{4}, n_ref) where N
    (; crystal, latsize) = sys

    n0 = Tuple(to_cell(site2)) .- Tuple(to_cell(site1))

    # Try to build a bond with the provided offset n_ref
    if !isnothing(n_ref)
        if all(iszero, mod.(n_ref .- n0, latsize))
            return Bond(to_atom(site1), to_atom(site2), n_ref)
        else
            cell1 = Tuple(to_cell(site1))
            cell2 = Tuple(to_cell(site2))
            @error """Cells $cell1 and $cell2 are not compatible with the offset
                      $n_ref for a system with lattice size $latsize."""
            error("Incompatible displacement specified")
        end
    end
    
    # Otherwise, search over all possible wrappings of the bond
    ns = view([n0 .+ latsize .* (i,j,k) for i in -1:1, j in -1:1, k in -1:1], :)
    bonds = map(ns) do n
        Bond(to_atom(site1), to_atom(site2), n)
    end
    distances = global_distance.(Ref(crystal), bonds)

    # Indices of bonds, from smallest to largest
    perm = sortperm(distances)

    # If one of the bonds is much shorter than all others by some arbitrary
    # `safety` factor, then return it
    safety = 4
    if safety * distances[perm[1]] < distances[perm[2]] - 1e-12
        return bonds[perm[1]]
    else
        n1 = bonds[perm[1]].n
        n2 = bonds[perm[2]].n
        @error """Cannot find an obvious offset vector. Possibilities include $n1 and $n2.
                  Try using a bigger system size, or pass an explicit offset vector."""
        error("Ambiguous offset between sites.")
    end
end


"""
    set_exchange_at!(sys::System, J, site1::Site, site2::Site; biquad=0, large_S=false, offset=nothing)

Sets the exchange interaction along the single bond connecting two
[`Site`](@ref)s, ignoring crystal symmetry. The system must support
inhomogeneous interactions via [`to_inhomogeneous`](@ref).

If the system is relatively small, the direction of the bond can be ambiguous
due to possible periodic wrapping. Resolve this ambiguity by passing an explicit
`offset` vector, in multiples of unit cells.

See also [`set_exchange!`](@ref).
"""
function set_exchange_at!(sys::System{N}, J, site1::Site, site2::Site; biquad=0, large_S=false, offset=nothing) where N
    site1 = to_cartesian(site1)
    site2 = to_cartesian(site2)
    bond = sites_to_internal_bond(sys, site1, site2, offset)

    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")
    ints = interactions_inhomog(sys)

    J = to_float_or_mat3(J)
    push_coupling!(sys, ints[site1].pair, bond, J, biquad, large_S)
    push_coupling!(sys, ints[site2].pair, reverse(bond), J', biquad', large_S)

    return
end


"""
    remove_periodicity!(sys::System, dims)

Remove periodic interactions along the dimensions where `dims` is `true`. The
system must support inhomogeneous interactions via [`to_inhomogeneous`](@ref).

# Example

```julia
# Remove periodic boundaries along the 1st and 3rd dimensions
remove_periodicity!(sys::System, (true, false, true))
```
"""
function remove_periodicity!(sys::System{N}, dims) where N
    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")

    for site in eachsite(sys)
        ints = interactions_inhomog(sys)[site]
        filter!(ints.pair) do (; bond)
            offset_cell = Tuple(to_cell(site)) .+ bond.n

            # keep bond if it is acceptable along every dimension (either
            # `!dims` or if each cell component is within bounds)
            all(@. !dims || 1 <= offset_cell <= sys.latsize)
        end
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
