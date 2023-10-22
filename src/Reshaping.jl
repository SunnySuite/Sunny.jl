
"""
    reshape_supercell(sys::System, shape)

Maps an existing [`System`](@ref) to a new one that has the shape and
periodicity of a requested supercell. The columns of the ``3×3`` integer matrix
`shape` represent the supercell lattice vectors measured in units of the
original crystal lattice vectors.
"""
function reshape_supercell(sys::System{N}, shape) where N
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


# Transfer homogeneous interactions from `sys.origin` to reshaped `sys`
function set_interactions_from_origin!(sys::System{N}) where N
    origin = sys.origin
    new_ints  = interactions_homog(sys)
    orig_ints = interactions_homog(origin)

    for new_i in 1:natoms(sys.crystal)
        i = map_atom_to_crystal(sys.crystal, new_i, origin.crystal)
        new_ints[new_i].onsite = orig_ints[i].onsite

        new_pc = PairCoupling[]
        for pc in orig_ints[i].pair
            new_bond = transform_bond(sys.crystal, new_i, origin.crystal, pc.bond)
            isculled = bond_parity(new_bond)
            push!(new_pc, PairCoupling(isculled, new_bond, pc.scalar, pc.bilin, pc.biquad, pc.general))
        end
        new_pc = sort!(new_pc, by=c->c.isculled)
        new_ints[new_i].pair = new_pc
    end
end


function reshape_supercell_aux(sys::System{N}, new_latsize::NTuple{3, Int}, new_cell_shape::Mat3) where N
    is_homogeneous(sys) || error("Cannot reshape system with inhomogeneous interactions.")

    # `origin` describes the unit cell of the original system. For sequential
    # reshapings, `sys.origin` keeps its original meaning. Make a deep copy so
    # that the new system fully owns `origin`, and mutable updates to the
    # previous system won't affect this one.
    origin = clone_system(isnothing(sys.origin) ? sys : sys.origin)

    # If `new_cell_shape == I`, we can effectively restore the unit cell of
    # `origin`. Otherwise, we will need to reshape the crystal, map the
    # interactions, and keep a reference to the original system.
    if new_cell_shape ≈ I
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

        new_ints = interactions_homog(origin)

        new_sys = System(nothing, origin.mode, new_cryst, new_latsize, new_Ns, new_κs, new_gs, new_ints, nothing,
            new_extfield, new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))
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

        new_ints = empty_interactions(origin.mode, new_na, N)

        new_sys = System(origin, origin.mode, new_cryst, new_latsize, new_Ns, new_κs, new_gs, new_ints, nothing,
            new_extfield, new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))

        set_interactions_from_origin!(new_sys)
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

    # Restore dipole-dipole interactions if present
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
    @assert all(>=(1), counts)
    # Scale each column by `counts` and reshape
    return reshape_supercell_aux(sys, counts .* sys.latsize, cell_shape(sys))
end
