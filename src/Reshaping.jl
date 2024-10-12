
"""
    reshape_supercell(sys::System, shape)

Maps an existing [`System`](@ref) to a new one that has the shape and
periodicity of a requested supercell. The columns of the ``3×3`` integer matrix
`shape` represent the supercell lattice vectors measured in units of the
original crystal lattice vectors. Interactions, spins, and other settings will
be inherited from `sys`.

In the special case that `shape` is a diagonal matrix, this function coincides
with [`resize_supercell`](@ref).

See also [`repeat_periodically`](@ref).
"""
function reshape_supercell(sys::System, shape)
    is_homogeneous(sys) || error("Cannot reshape system with inhomogeneous interactions.")

    orig = orig_crystal(sys)
    check_shape_commensurate(orig, shape)
    prim_cell = @something primitive_cell(orig) Mat3(I)
    shape_in_prim = prim_cell \ shape
    @assert all_integer(shape_in_prim; orig.symprec)
    shape_in_prim = round.(Int, shape_in_prim)

    # Unit cell for new system, in units of original unit cell.
    new_dims = NTuple{3, Int}(gcd.(eachcol(shape_in_prim)))
    new_shape = Mat3(shape * diagm(collect(inv.(new_dims))))
    new_cryst = reshape_crystal(orig_crystal(sys), new_shape)

    return reshape_supercell_aux(sys, new_cryst, new_dims)
end


# Transfer interactions from `src` to reshaped `sys`.
#
# Frequently `src` will refer to `sys.origin`, associated with the original
# crystal, which will make symmetry analysis available. In this case, the
# process to set a new interaction is to first modify `sys.origin`, and then to
# call this function.
function transfer_interactions!(sys::System, src::System)
    new_ints = interactions_homog(sys)

    for new_i in 1:natoms(sys.crystal)
        # Find `src` interaction either through an atom index or a site index
        if is_homogeneous(src)
            i = map_atom_to_other_crystal(sys.crystal, new_i, src.crystal)
        else
            i = map_atom_to_other_system(sys.crystal, new_i, src)
        end
        src_int = src.interactions_union[i]

        # Copy onsite couplings
        new_ints[new_i].onsite = src_int.onsite

        # Copy pair couplings
        new_pc = PairCoupling[]
        for pc in src_int.pair
            new_bond = map_bond_to_other_crystal(src.crystal, pc.bond, sys.crystal, new_i)
            push!(new_pc, PairCoupling(new_bond, pc.scalar, pc.bilin, pc.biquad, pc.general))
        end
        new_pc = sort!(new_pc, by=c->c.isculled)
        new_ints[new_i].pair = new_pc
    end
end


function reshape_supercell_aux(sys::System{N}, new_cryst::Crystal, new_dims::NTuple{3, Int}) where N
    # Allocate data for new system, but with an empty list of interactions
    new_na               = natoms(new_cryst)
    new_Ns               = zeros(Int, new_dims..., new_na)
    new_κs               = zeros(Float64, new_dims..., new_na)
    new_gs               = zeros(Mat3, new_dims..., new_na)
    new_ints             = empty_interactions(sys.mode, new_na, N)
    new_ewald            = nothing
    new_extfield         = zeros(Vec3, new_dims..., new_na)
    new_dipoles          = zeros(Vec3, new_dims..., new_na)
    new_coherents        = zeros(CVec{N}, new_dims..., new_na)
    new_dipole_buffers   = Array{Vec3, 4}[]
    new_coherent_buffers = Array{CVec{N}, 4}[]

    # Preserve origin for repeated reshapings. In other words, ensure that
    # new_sys.origin === sys.origin
    orig_sys = @something sys.origin sys

    new_sys = System(orig_sys, sys.mode, new_cryst, new_dims, new_Ns, new_κs, new_gs, new_ints, new_ewald,
        new_extfield, new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, copy(sys.rng))

    # Transfer interactions. In the case of an inhomogeneous system, the
    # interactions in `sys` have detached from `orig`, so we use the latest
    # ones.
    transfer_interactions!(new_sys, sys)

    # Copy per-site quantities
    for new_site in eachsite(new_sys)
        site = position_to_site(sys, position(new_sys, new_site))
        new_sys.Ns[new_site]        = sys.Ns[site]
        new_sys.κs[new_site]        = sys.κs[site]
        new_sys.gs[new_site]        = sys.gs[site]
        new_sys.extfield[new_site]  = sys.extfield[site]
        new_sys.dipoles[new_site]   = sys.dipoles[site]
        new_sys.coherents[new_site] = sys.coherents[site]
    end

    # Restore dipole-dipole interactions if present. This involves pre-computing
    # an interaction matrix that depends on `new_dims`.
    if !isnothing(sys.ewald)
        enable_dipole_dipole!(new_sys, sys.ewald.μ0_μB²)
    end

    return new_sys
end


# Shape of a possibly reshaped unit cell, given in multiples of the original
# unit cell.
function cell_shape(sys)
    return orig_crystal(sys).latvecs \ sys.crystal.latvecs
end

"""
    resize_supercell(sys::System, dims::NTuple{3, Int})

Creates a [`System`](@ref) with a given number of conventional unit cells in
each lattice vector direction. Interactions, spins, and other settings will be
inherited from `sys`.

Equivalent to:

```julia
reshape_supercell(sys, [dims[1] 0 0; 0 dims[2] 0; 0 0 dims[3]])
```

See also [`reshape_supercell`](@ref) and [`repeat_periodically`](@ref).
"""
function resize_supercell(sys::System, dims::NTuple{3,Int})
    return reshape_supercell(sys, diagm(Vec3(dims)))
end

"""
    repeat_periodically(sys::System, counts::NTuple{3, Int})

Creates a [`System`](@ref) identical to `sys` but repeated a given number of
times along each system axis according to the specified `counts`. This is a
special case of the more general [`reshape_supercell`](@ref).

See also [`repeat_periodically_as_spiral`](@ref), which rotates the spins
between periodic copies.
"""
function repeat_periodically(sys::System, counts::NTuple{3,Int})
    all(>=(1), counts) || error("Require at least one count in each direction.")
    return reshape_supercell_aux(sys, sys.crystal, counts .* sys.dims)
end

"""
    repeat_periodically_as_spiral(sys::System, counts::NTuple{3, Int}; k, axis)

Repeats the magnetic cell of [`System`](@ref) a number of times along each
system axis according to the specified `counts`. Spins in each system image will
be rotated according to the propagation wavevector `k` (in RLU) and the rotation
`axis` (in global Cartesian coordinates). Coincides with
[`repeat_periodically`](@ref) in the special case of `k = [0, 0, 0]`

See also [`minimize_spiral_energy!`](@ref) to find an energy-minimizing
wavevector `k` and spin dipole configuration.

# Example
```julia
k = minimize_spiral_energy!(sys, axis; k_guess=randn(3))
repeat_periodically_as_spiral(sys, counts; k, axis)
```
"""
function repeat_periodically_as_spiral(sys::System, counts::NTuple{3,Int}; k, axis)
    new_sys = repeat_periodically(sys, counts)

    # Propagation wavevector in global coordinates
    k_global = orig_crystal(sys).recipvecs * k

    # Lattice vectors for sys, the magnetic cell 
    supervecs = sys.crystal.latvecs * diagm(Vec3(sys.dims))

    # Original positions units of supervecs (components between 0 and 1)
    rs = [supervecs \ global_position(sys, site) for site in eachsite(sys)]

    # Tighted symprec suitable for scaled positions
    symprec = sys.crystal.symprec / minimum(sys.dims)

    # Copy per-site quantities
    for new_site in eachsite(new_sys)
        # Positions of new_sys in units of supervecs
        new_r = supervecs \ global_position(new_sys, new_site)

        # Find index into original sys corresponding to a periodic copy of new_r
        site = findfirst(r -> is_periodic_copy(new_r, r; symprec), rs)

        # Offset of periodic image in global coordinates
        offset = supervecs * (new_r - rs[site])
        R = axis_angle_to_matrix(axis, k_global ⋅ offset)

        # Rotate original spin of sys to new_sys
        if sys.mode == :SUN
            U = unitary_irrep_for_rotation(R; N=sys.Ns[site])
            Z = sys.coherents[site]
            set_coherent!(new_sys, U * Z, new_site)
        else
            s = sys.dipoles[site]
            set_dipole!(new_sys, R * s, new_site)
        end
    end

    return new_sys
end
