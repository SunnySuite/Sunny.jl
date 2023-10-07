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

# Perform renormalization derived in https://arxiv.org/abs/2304.03874
function renormalize_biquad(sys, bond, bilin, biquad)
    @assert sys.mode == :dipole

    # Multiplicative average is the correct result for two sites with
    # different S (derivation omitted).
    S1 = (sys.Ns[bond.i]-1)/2
    S2 = (sys.Ns[bond.j]-1)/2
    S = sqrt(S1*S2)
    r = (1 - 1/S + 1/4S^2)

    # biquad (s⋅sⱼ)^2 -> biquad (r (sᵢ⋅sⱼ)^2 - (sᵢ⋅sⱼ)/2 + S^3 + S^2/4)
    bilin = bilin - (biquad/2) * one(bilin)
    biquad = r * biquad
    
    # Drop the constant shift, `biquad (S^3 + S^2/4)`.

    return (bilin, biquad)
end

# Internal function only
function push_coupling!(couplings, bond, bilin, biquad, tensordec)
    # Remove previous coupling on this bond
    filter!(c -> c.bond != bond, couplings)

    # If the new coupling is exactly zero, return early
    iszero(bilin) && iszero(biquad) && isempty(tensordec.data) && return

    # Otherwise, add the new coupling to the list
    isculled = bond_parity(bond)
    push!(couplings, PairCoupling(isculled, bond, bilin, biquad, tensordec))

    # Sorting after each insertion will introduce quadratic scaling in length of
    # `couplings`. In typical usage, the `couplings` list will be short.
    sort!(couplings, by=c->c.isculled)
    
    return
end


function decompose_general_coupling(op, gen1, gen2; fast)
    N1 = size(gen1[1], 1)
    N2 = size(gen2[1], 1)

    if fast
        bilin = zeros(3, 3)
        for α in 1:3, β in 1:3
            v = kron(gen1[α], gen2[β])
            J = tr(v' * op) / tr(v' * v)
            @assert imag(J) < 1e-12
            bilin[α, β] = real(J)
            op = op - v * bilin[α, β]
        end
        bilin = Mat3(bilin)

        u = sum(kron(gen1[α], gen2[α]) for α in 1:3)
        if norm(op) > 1e-12 && normalize(op) ≈ normalize(u^2 + u/2)
            @info "Detected scalar biquadratic. Not yet optimized."
        end
    else
        bilin = 0.0
    end

    return bilin, TensorDecomposition(gen1, gen2, svd_tensor_expansion(op, N1, N2))
end

function Base.zero(::Type{TensorDecomposition})
    gen = spin_matrices(; N=0)
    return TensorDecomposition(gen, gen, [])
end

function transform_coupling_by_symmetry(cryst, tensordec::TensorDecomposition, symop, parity)
    (; gen1, gen2, data) = tensordec
    isempty(data) && return tensordec

    if !parity
        data = [(B, A) for (A, B) in data]
        gen2, gen1 = (gen1, gen2)
    end
    R = cryst.latvecs * symop.R * inv(cryst.latvecs)
    Q = R * det(R)
    U1 = unitary_for_rotation(Q, gen1)
    U2 = unitary_for_rotation(Q, gen2)
    # Under the symop, coherents transform as `Z -> U Z`. Then couplings must
    # transform as `A -> U A U'` so that the expected energy on a bond `⟨A⟩⟨B⟩`
    # is invariant. By analogy, spin rotates as `S -> R S` and the 3×3 exchange
    # matrix transforms as `J -> R J Rᵀ` to preserve `Sᵀ J S`.
    data = [(Hermitian(U1*A*U1'), Hermitian(U2*B*U2')) for (A, B) in data]
    return TensorDecomposition(gen1, gen2, data)
end

function Base.isapprox(op1::TensorDecomposition, op2::TensorDecomposition; kwargs...)
    isempty(op1.data) == isempty(op2.data) && return true
    op1′ = sum(kron(A, B) for (A, B) in op1.data)
    op2′ = sum(kron(A, B) for (A, B) in op2.data)
    return isapprox(op1′, op2′; kwargs...)
end

function set_pair_coupling_aux!(sys::System, bilin::Union{Float64, Mat3}, biquad::Float64, tensordec::TensorDecomposition, bond::Bond)
    # If `sys` has been reshaped, then operate first on `sys.origin`, which
    # contains full symmetry information.
    if !isnothing(sys.origin)
        set_pair_coupling_aux!(sys.origin, bilin, biquad, tensordec, bond)
        set_interactions_from_origin!(sys)
        return
    end

    # Simple checks on bond indices
    validate_bond(sys.crystal, bond)

    # Verify that couplings are symmetry-consistent
    if !is_coupling_valid(sys.crystal, bond, bilin)
        effective = isempty(tensordec.data) ? "" : " effective"
        @error """Symmetry-violating$effective exchange $bilin.
                  Use `print_bond(crystal, $bond)` for more information."""
    end
    if !is_coupling_valid(sys.crystal, bond, tensordec)
        @error """Symmetry-violating coupling. Use `print_bond(crystal, $bond)` for more information."""
        error("Interaction violates symmetry.")
    end

    # Print a warning if an interaction already exists for bond
    ints = interactions_homog(sys)
    if any(x -> x.bond == bond, ints[bond.i].pair)
        warn_coupling_override("Overriding coupling for $bond.")
    end

    # Propagate all couplings by symmetry
    for i in 1:natoms(sys.crystal)
        for bond′ in all_symmetry_related_bonds_for_atom(sys.crystal, i, bond)
            bilin′ = transform_coupling_for_bonds(sys.crystal, bond′, bond, bilin)
            tensordec′ = transform_coupling_for_bonds(sys.crystal, bond′, bond, tensordec)
            push_coupling!(ints[i].pair, bond′, bilin′, biquad, tensordec′)
        end
    end
end


"""
    set_pair_coupling!(sys::System, coupling, bond)

Sets an arbitrary `coupling` along `bond`. This coupling will be propagated to
equivalent bonds in consistency with crystal symmetry. Any previous interactions
on these bonds will be overwritten. The parameter `bond` has the form `Bond(i,
j, offset)`, where `i` and `j` are atom indices within the unit cell, and
`offset` is a displacement in unit cells. The `coupling` is a represented as a
matrix acting in the tensor product space of the two sites, and typically
originates from [`to_product_space`](@ref).

# Examples
```julia
# Add a bilinear and biquadratic exchange
S = spin_matrices(1/2)
Si, Sj = to_product_space(S, S)
set_pair_coupling!(sys, Si'*J1*Sj + (Si'*J2*Sj)^2, bond)
```
"""
function set_pair_coupling!(sys::System{N}, tensordec::Matrix{ComplexF64}, bond; fast=true) where N
    is_homogeneous(sys) || error("Use `set_pair_coupling_at!` for an inhomogeneous system.")

    gen1 = spin_operators(sys, bond.i)
    gen2 = spin_operators(sys, bond.j)
    bilin, tensordec = decompose_general_coupling(tensordec, gen1, gen2; fast)
    biquad = 0.0

    set_pair_coupling_aux!(sys, bilin, biquad, tensordec, bond)
end

"""
    set_exchange!(sys::System, J, bond::Bond; biquad=0, large_S=false)

Sets a 3×3 spin-exchange matrix `J` along `bond`, yielding a pairwise
interaction energy ``𝐒_i⋅J 𝐒_j``. This interaction will be propagated to
equivalent bonds in consistency with crystal symmetry. Any previous interactions
on these bonds will be overwritten. The parameter `bond` has the form `Bond(i,
j, offset)`, where `i` and `j` are atom indices within the unit cell, and
`offset` is a displacement in unit cells.

The parameter `J` may be scalar or matrix-valued. As a convenience, `dmvec(D)`
can be used to construct the antisymmetric part of the exchange, where `D` is
the Dzyaloshinskii-Moriya pseudo-vector. The resulting interaction will be
``𝐃⋅(𝐒_i×𝐒_j)``.

The optional parameter `biquad` defines the strength ``b`` for scalar
biquadratic interactions of the form ``b (𝐒_i⋅𝐒_j)²``. For systems restricted
to dipoles, ``b`` will be automatically renormalized for maximum consistency
with the more variationally accurate SU(_N_) mode. Set `large_S=true` to work in
the large-``S`` limit and disable this renormalization.

# Examples
```julia
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
    bilin = to_float_or_mat3(J)
    if sys.mode == :dipole && !large_S
        bilin, biquad = renormalize_biquad(sys, bond, bilin, biquad)
    end
    set_pair_coupling_aux!(sys, bilin, Float64(biquad), zero(TensorDecomposition), bond)
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

    bilin = to_float_or_mat3(J)
    if sys.mode == :dipole && !large_S
        bilin, biquad = renormalize_biquad(sys, bond, bilin, biquad)
    end

    push_coupling!(ints[site1].pair, bond, bilin, biquad, zero(TensorDecomposition))
    push_coupling!(ints[site2].pair, reverse(bond), bilin', biquad', zero(TensorDecomposition))

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
