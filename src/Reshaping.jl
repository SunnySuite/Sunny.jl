
"""
    reshape_supercell(sys::System, shape)

Maps an existing [`System`](@ref) to a new one that has the shape and
periodicity of a requested supercell. The columns of the ``3×3`` integer matrix
`shape` represent the supercell lattice vectors measured in units of the
original crystal lattice vectors.
"""
function reshape_supercell(sys::System{N}, shape) where N
    is_homogeneous(sys) || error("Cannot reshape system with inhomogeneous interactions.")

    orig = orig_crystal(sys)

    supervecs = orig.latvecs * shape
    check_latvecs_commensurate(orig, supervecs)

    # Convert `shape` to multiples of primitive lattice vectors
    if !isnothing(orig.prim_latvecs)
        prim_shape = orig.prim_latvecs \ supervecs
        prim_shape′ = round.(Int, prim_shape)
        @assert prim_shape′ ≈ prim_shape
        new_latsize = NTuple{3, Int}(gcd.(eachcol(prim_shape′)))
    else
        shape′ = round.(Int, shape)
        @assert shape′ ≈ shape
        new_latsize = NTuple{3, Int}(gcd.(eachcol(shape′)))
    end

    # Unit cell for new system, in units of original unit cell.
    new_cell_shape = Mat3(shape * diagm(collect(inv.(new_latsize))))

    return reshape_supercell_aux(sys, new_latsize, new_cell_shape)
end


# Transfer interactions from `src` to reshaped `sys`.
#
# Frequently `src` will refer to `sys.origin`, associated with the original
# crystal, which will make symmetry analysis available. In this case, the
# process to set a new interaction is to first modify `sys.origin`, and then to
# call this function.
function set_interactions_from_origin!(sys::System{N}, src::System{N}) where N
    new_ints = interactions_homog(sys)

    for new_i in 1:natoms(sys.crystal)
        # Find `src` interaction either through an atom index or a site index
        if is_homogeneous(src)
            i = map_atom_to_other_crystal(sys.crystal, new_i, src.crystal)
        else
            # Use case here is to prepare for a LSWT calculation by building an
            # equivalent system with a single unit cell
            @assert sys.latsize == (1,1,1)
            @assert sys.crystal.latvecs ≈ src.crystal.latvecs * diagm(Vec3(src.latsize))
            i = map_atom_to_other_system(sys.crystal, new_i, src)
        end
        src_int = src.interactions_union[i]

        # Copy onsite couplings
        new_ints[new_i].onsite = src_int.onsite

        # Copy pair couplings
        new_pc = PairCoupling[]
        for pc in src_int.pair
            new_bond = transform_bond(sys.crystal, new_i, orig_crystal(src), pc.bond)
            isculled = bond_parity(new_bond)
            push!(new_pc, PairCoupling(isculled, new_bond, pc.scalar, pc.bilin, pc.biquad, pc.general))
        end
        new_pc = sort!(new_pc, by=c->c.isculled)
        new_ints[new_i].pair = new_pc
    end
end


function reshape_supercell_aux(sys::System{N}, new_latsize::NTuple{3, Int}, new_cell_shape::Mat3) where N

    # `origin` refers to a system with the original unit cell. For sequential
    # reshapings, `sys.origin` keeps its original meaning. Make a deep copy so
    # that the new system fully owns `origin`, and mutable updates to the
    # previous system won't affect this one.
    origin = clone_system(isnothing(sys.origin) ? sys : sys.origin)

    # If `new_cell_shape == I`, we can reuse unit cell of `origin`. Still need
    # to reallocate various fields because `new_latsize` may have changed.
    if new_cell_shape ≈ I && is_homogeneous(sys)
        new_cryst = origin.crystal

        new_na               = natoms(new_cryst)
        new_Ns               = zeros(Int, new_latsize..., new_na)
        new_κs               = zeros(Float64, new_latsize..., new_na)
        new_gs               = zeros(Mat3, new_latsize..., new_na)
        new_extfield         = zeros(Vec3, new_latsize..., new_na)
        new_dipoles          = zeros(Vec3, new_latsize..., new_na)
        new_coherents        = zeros(CVec{N}, new_latsize..., new_na)
        new_dipole_buffers   = Array{Vec3, 4}[]
        new_coherent_buffers = Array{CVec{N}, 4}[]

        # Can reuse existing interactions
        new_ints = interactions_homog(origin)
        new_ewald = nothing

        new_sys = System(nothing, origin.mode, new_cryst, new_latsize, new_Ns, new_κs, new_gs, new_ints, new_ewald,
                         new_extfield, new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))

    # Otherwise, we will need to reshape the crystal, map the interactions, and
    # keep a reference to the original system.
    else
        new_cryst = reshape_crystal(origin.crystal, new_cell_shape)

        new_na               = natoms(new_cryst)
        new_Ns               = zeros(Int, new_latsize..., new_na)
        new_κs               = zeros(Float64, new_latsize..., new_na)
        new_gs               = zeros(Mat3, new_latsize..., new_na)
        new_extfield         = zeros(Vec3, new_latsize..., new_na)
        new_dipoles          = zeros(Vec3, new_latsize..., new_na)
        new_coherents        = zeros(CVec{N}, new_latsize..., new_na)
        new_dipole_buffers   = Array{Vec3, 4}[]
        new_coherent_buffers = Array{CVec{N}, 4}[]

        # Begin with empty interactions
        new_ints             = empty_interactions(origin.mode, new_na, N)
        new_ewald            = nothing

        new_sys = System(origin, origin.mode, new_cryst, new_latsize, new_Ns, new_κs, new_gs, new_ints, new_ewald,
            new_extfield, new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))

        # Fill in interactions. In the case of an inhomogeneous system, the
        # interactions in `sys` have detached from `orig`, so we use the latest
        # ones.
        set_interactions_from_origin!(new_sys, is_homogeneous(sys) ? origin : sys)
    end

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
    # an interaction matrix that depends on `new_latsize`.
    if !isnothing(sys.ewald)
        enable_dipole_dipole!(new_sys)
    end

    return new_sys
end


# Shape of a possibly reshaped unit cell, given in multiples of the original
# unit cell.
function cell_shape(sys)
    return orig_crystal(sys).latvecs \ sys.crystal.latvecs
end

"""
    resize_supercell(sys::System{N}, latsize::NTuple{3,Int}) where N

Creates a [`System`](@ref) with a given number of conventional unit cells in
each lattice vector direction. Interactions and other settings will be inherited
from `sys`.

Convenience function for:
```julia
reshape_supercell(sys, [latsize[1] 0 0; 0 latsize[2] 0; 0 0 latsize[3]])
```

See also [`reshape_supercell`](@ref).
"""
function resize_supercell(sys::System{N}, latsize::NTuple{3,Int}) where N
    return reshape_supercell(sys, diagm(collect(latsize)))
end

"""
    repeat_periodically(sys::System{N}, counts::NTuple{3,Int}) where N

Creates a [`System`](@ref) identical to `sys` but repeated a given number of
times in each dimension, specified by the tuple `counts`.

See also [`reshape_supercell`](@ref).
"""
function repeat_periodically(sys::System{N}, counts::NTuple{3,Int}) where N
    is_homogeneous(sys) || error("Cannot reshape system with inhomogeneous interactions.")
    all(>=(1), counts) || error("Require at least one count in each direction.")

    # Scale each column by `counts` and reshape
    return reshape_supercell_aux(sys, counts .* sys.latsize, cell_shape(sys))
end
