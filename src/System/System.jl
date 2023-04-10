"""
    System(crystal::Crystal, latsize, infos, mode; units=Units.meV, seed::Int)

Construct a `System` of spins for a given [`Crystal`](@ref) symmetry. The
`latsize` parameter determines the number of unit cells in each lattice vector
direction. The `infos` parameter is a list of [`SpinInfo`](@ref) objects, which
determine the magnitude ``S`` and ``g``-tensor of each spin.

The three possible options for `mode` are `:SUN`, `:dipole`, and `:large_S`. The
most variationally accurate choice is `:SUN`, in which each spin-``S`` degree of
freedom is described as an SU(_N_) coherent state, where ``N = 2S + 1``. Note
that an SU(_N_) coherent state fully describes any local spin state; this
description includes expected dipole components ``⟨Ŝᵅ⟩``, quadrupole components
``⟨ŜᵅŜᵝ+ŜᵝŜᵅ⟩``, etc.

The mode `:dipole` projects the SU(_N_) dynamics onto the space of pure dipoles.
In practice this means that Sunny will simulate Landau-Lifshitz dynamics, but
all single-ion anisotropy and biquadratic exchange interactions will be
automatically renormalized for maximum accuracy.

To disable such renormalization, e.g. to reproduce results using the historical
large-``S`` classical limit, use the experimental mode `:large_S`. Modes `:SUN`
or `:dipole` are strongly preferred for the development of new models.

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

An optional `seed` may be provided to achieve reproducible random number
generation.

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, infos::Vector{SpinInfo}, mode::Symbol;
                    units=Units.meV, seed=nothing)
    if mode ∉ [:SUN, :dipole, :large_S]
        error("Mode must be one of [:SUN, :dipole, :large_S].")
    end

    na = natoms(crystal)

    infos = propagate_site_info(crystal, infos)
    Ss = [si.S for si in infos]
    gs = [si.g for si in infos]
    Ns = @. Int(2Ss+1)

    if mode == :SUN
        allequal(Ns) || error("Currently all spins S must be equal in SU(N) mode.")
        N = first(Ns)
        κs = fill(1.0, na)
    else
        N = 0 # acts as marker for :dipole
        κs = Ss
    end

    # Repeat such that `A[:]` → `A[cell, :]` for every `cell`
    repeat_to_lattice(A) = permutedims(repeat(A, 1, latsize...), (2, 3, 4, 1))

    Ns = repeat_to_lattice(Ns)
    κs = repeat_to_lattice(κs)
    gs = repeat_to_lattice(gs)

    interactions = empty_interactions(na, N)
    ewald = nothing

    extfield = zeros(Vec3, latsize..., na)
    dipoles = fill(zero(Vec3), latsize..., na)
    coherents = fill(zero(CVec{N}), latsize..., na)
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = System(nothing, mode, crystal, latsize, Ns, κs, gs, interactions, ewald,
                 extfield, dipoles, coherents, dipole_buffers, coherent_buffers, units, rng)
    polarize_spins!(ret, (0,0,1))
    return ret
end

function Base.show(io::IO, ::MIME"text/plain", sys::System{N}) where N
    modename = if sys.mode==:SUN
        "SU($N)"
    elseif sys.mode==:dipole
        "Dipole mode"
    elseif sys.mode==:large_S
        "Large-S classical limit"
    else
        error("Unreachable")
    end
    printstyled(io, "System [$modename]\n"; bold=true, color=:underline)
    println(io, "Cell size $(natoms(sys.crystal)), Lattice size $(sys.latsize)")
    if !isnothing(sys.origin)
        println(io, "Reshaped cell geometry $(cell_dimensions(sys))")
    end
end

# Per Julia developers, `deepcopy` is memory unsafe, especially in conjunction
# with C libraries. We were observing very confusing crashes that surfaced in
# the FFTW library, https://github.com/JuliaLang/julia/issues/48722. To prevent
# this from happening again, avoid all uses of `deepcopy`, and create our own
# stack of `clone` functions instead.
Base.deepcopy(_::System) = error("Use `clone_system` instead of `deepcopy`.")

# Creates a clone of the system where all the mutable internal data is copied.
# It is intended to be thread-safe to use the original and the copied systems,
# without any restrictions, but see caveats in `clone_ewald()`.
function clone_system(sys::System{N}) where N
    (; origin, mode, crystal, latsize, Ns, gs, κs, extfield, interactions_union, ewald, dipoles, coherents, units, rng) = sys

    origin_clone = isnothing(origin) ? nothing : clone_system(origin)
    ewald_clone  = isnothing(ewald)  ? nothing : clone_ewald(ewald)

    # Dynamically dispatch to the correct `map` function for either homogeneous
    # (Vector) or inhomogeneous interactions (4D Array)
    interactions_clone = map(clone_interactions, interactions_union)
    
    # Empty buffers are required for thread safety.
    empty_dipole_buffers = Array{Vec3, 4}[]
    empty_coherent_buffers = Array{CVec{N}, 4}[]

    System(origin_clone, mode, crystal, latsize, Ns, copy(κs), copy(gs),
           interactions_clone, ewald_clone, copy(extfield), copy(dipoles), copy(coherents),
           empty_dipole_buffers, empty_coherent_buffers, units, copy(rng))
end


"""
    (cell1, cell2, cell3, i) :: Site

Four indices identifying a single site in a [`System`](@ref). The first three
indices select the lattice cell and the last selects the sublattice (i.e., the
atom within the unit cell).

This object can be used to index `dipoles` and `coherents` fields of a `System`.
A `Site` is also required to specify inhomogeneous interactions via functions
such as [`set_external_field_at!`](@ref) or [`set_exchange_at!`](@ref).

Note that the definition of a cell may change when a system is reshaped. In this
case, it is convenient to construct the `Site` using [`position_to_site`](@ref),
which always takes a position in fractional coordinates of the original lattice
vectors.
"""
const Site = NTuple{4, Int}

@inline to_cartesian(i::CartesianIndex{N}) where N = i
@inline to_cartesian(i::NTuple{N, Int})    where N = CartesianIndex(i)

# kbtodo: offsetcell ?
# Offset a `cell` by `ncells`
@inline offsetc(cell::CartesianIndex{3}, ncells, latsize) = CartesianIndex(mod1.(Tuple(cell) .+ Tuple(ncells), latsize))

# Split a site `site` into its cell and sublattice parts
@inline to_cell(site) = CartesianIndex((site[1],site[2],site[3]))
@inline to_atom(site) = site[4]

# An iterator over all unit cells using CartesianIndices
@inline all_cells(sys::System) = CartesianIndices(sys.latsize)

"""
    all_sites(sys::System)

An iterator over all [`Site`](@ref)s in the system. 
"""
@inline all_sites(sys::System) = CartesianIndices(sys.dipoles)

"""
    global_position(sys::System, site::Site)

Position of a [`Site`](@ref) in global coordinates.

To precompute a full list of positions, one can use [`all_sites`](@ref) as
below:

```julia
pos = [global_position(sys, site) for site in all_sites(sys)]
```
"""
function global_position(sys::System, site)
    r = sys.crystal.positions[site[4]] + Vec3(site[1]-1, site[2]-1, site[3]-1)
    return sys.crystal.latvecs * r
end

"""
    magnetic_moment(sys::System, site::Site)

Get the magnetic moment for a [`Site`](@ref). This is the spin dipole multiplied
by the Bohr magneton and the local g-tensor.
"""
magnetic_moment(sys::System, site) = sys.units.μB * sys.gs[site] * sys.dipoles[site]

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)

# The original crystal for a system, invariant under reshaping
orig_crystal(sys) = isnothing(sys.origin) ? sys.crystal : sys.origin.crystal

# Position of a site in fractional coordinates of the original crystal
function position(sys::System, site)
    return orig_crystal(sys).latvecs \ global_position(sys, site)
end

"""
    position_to_site(sys::System, r)

Converts a position `r` to four indices of a [`Site`](@ref). The coordinates of
`r` are given in units of the lattice vectors for the original crystal. This
function can be useful for working with systems that have been reshaped using
[`reshape_geometry`](@ref).

# Example

```julia
# Find the `site` at the center of a unit cell which is displaced by four
# multiples of the first lattice vector
site = position_to_site(sys, [4.5, 0.5, 0.5])

# Print the dipole at this site
println(sys.dipoles[site])
```
"""
function position_to_site(sys::System, r)
    # convert to fractional coordinates of possibly reshaped crystal
    r = Vec3(r)
    new_r = sys.crystal.latvecs \ orig_crystal(sys).latvecs * r
    i, offset = position_to_atom_and_offset(sys.crystal, new_r)
    cell = @. mod1(offset+1, sys.latsize) # 1-based indexing with periodicity
    return to_cartesian((cell..., i))
end


# Given two [`Site`](@ref)s for a possibly reshaped system, return the
# corresponding [`Bond`](@ref) for the original (unreshaped) crystal. This
# `bond` can be used for symmetry analysis, e.g., as input to
# [`print_bond`](@ref). See also [`site_to_atom`](@ref).
#
# Warning: This function will not be useful until site2 is wrapped properly.
function sites_to_bond(sys::System{N}, site1, site2) where N
    site1, site2 = to_cartesian.((site1, site2))

    # Position of sites in multiples of original latvecs.
    r1 = position(sys, site1)
    r2 = position(sys, site2)

    # Map to Bond for original crystal
    return Bond(orig_crystal(sys), BondPos(r1, r2))
end

# Given a [`Site`](@ref)s for a possibly reshaped system, return the
# corresponding atom index for the original (unreshaped) crystal. See also
# [`sites_to_bond`](@ref).
function site_to_atom(sys::System{N}, site) where N
    site = to_cartesian(site)
    r = position(sys, site)
    return position_to_atom(orig_crystal(sys), r)
end

# Maps atom `i` in `cryst` to the corresponding atom in `orig_cryst`
function map_atom_to_crystal(cryst, i, orig_cryst)
    r = cryst.positions[i]
    orig_r = orig_cryst.latvecs \ cryst.latvecs * r
    return position_to_atom(orig_cryst, orig_r)
end

# Given a `bond` for `cryst`, return a corresponding new bond for the reshaped
# `new_cryst`. The new bond will begin at atom `new_i`.
function transform_bond(new_cryst::Crystal, new_i::Int, cryst::Crystal, bond::Bond)
    # Positions in new fractional coordinates
    new_ri = new_cryst.positions[new_i]
    new_rj = new_ri + new_cryst.latvecs \ global_displacement(cryst, bond)

    # Verify that new_i (indexed into new_cryst) is consistent with bond.i
    # (indexed into original cryst).
    @assert bond.i == position_to_atom(cryst, cryst.latvecs \ new_cryst.latvecs * new_ri)

    # Construct bond using new indexing system
    new_j, new_n = position_to_atom_and_offset(new_cryst, new_rj)
    return Bond(new_i, new_j, new_n)
end


"""
    symmetry_equivalent_bonds(sys::System, bond::Bond)

Given a [`Bond`](@ref) for the original (unreshaped) crystal, return all
symmetry equivalent bonds in the [`System`](@ref). Each returned bond is
represented as a pair of [`Site`](@ref)s, which may be used as input to
[`set_exchange_at!`](@ref) or [`set_biquadratic_at!`](@ref). Reverse bonds are
not included (no double counting).

# Example
```julia
for (site1, site2) in symmetry_equivalent_bonds(sys, bond)
    @assert site1 < site2
    set_exchange_at!(sys, J, site1, site2)
end
```
"""
function symmetry_equivalent_bonds(sys::System, bond::Bond)
    ret = Tuple{Site, Site}[]

    for new_i in 1:natoms(sys.crystal)
        # atom index in original crystal
        i = map_atom_to_crystal(sys.crystal, new_i, orig_crystal(sys))

        # loop over symmetry equivalent bonds in original crystal
        for bond′ in all_symmetry_related_bonds_for_atom(orig_crystal(sys), i, bond)

            # map to a bond with indexing for new crystal
            new_bond = transform_bond(sys.crystal, new_i, orig_crystal(sys), bond′)

            # loop over all new crystal cells and push site pairs
            for new_cell_i in all_cells(sys)
                new_cell_j = offsetc(new_cell_i, new_bond.n, sys.latsize)
                site_i = (Tuple(new_cell_i)..., new_bond.i)
                site_j = (Tuple(new_cell_j)..., new_bond.j)
                site_i < site_j && push!(ret, (site_i, site_j))
            end
        end
    end

    return ret
end


struct SpinState{N}
    s::Vec3
    Z::CVec{N}
end

# Returns √κ * normalize(Z)
@inline function normalize_ket(Z::CVec{N}, κ) where N
    return iszero(κ) ? zero(CVec{N}) : Z/sqrt(dot(Z,Z)/κ)
end

# Returns κ * normalize(s)
@inline function normalize_dipole(s::Vec3, κ)
    return iszero(κ) ? zero(Vec3) : κ*normalize(s)
end

@inline function getspin(sys::System{N}, site) where N
    return SpinState(sys.dipoles[site], sys.coherents[site])
end

@inline function setspin!(sys::System{N}, spin::SpinState{N}, site) where N
    sys.dipoles[site] = spin.s
    sys.coherents[site] = spin.Z
    return
end

@inline function flip(spin::SpinState{N}) where N
    return SpinState(-spin.s, flip_ket(spin.Z))
end

@inline function randspin(sys::System{0}, site)
    s = normalize_dipole(randn(sys.rng, Vec3), sys.κs[site])
    return SpinState(s, CVec{0}())
end
@inline function randspin(sys::System{N}, site) where N
    Z = normalize_ket(randn(sys.rng, CVec{N}), sys.κs[site])
    s = expected_spin(Z)
    return SpinState(s, Z)
end

@inline function dipolarspin(sys::System{0}, site, dir)
    s = normalize_dipole(Vec3(dir), sys.κs[site])
    Z = CVec{0}()
    return SpinState(s, Z)
end
@inline function dipolarspin(sys::System{N}, site, dir) where N
    Z = normalize_ket(ket_from_dipole(Vec3(dir), Val(N)), sys.κs[site])
    s = expected_spin(Z)
    return SpinState(s, Z)
end


function randomize_spins!(sys::System{N}) where N
    for site in all_sites(sys)
        setspin!(sys, randspin(sys, site), site)
    end
end

"""
    set_coherent_state!(sys::System, Z, site::Site)

Set a coherent spin state at a [`Site`](@ref) using the ``N`` complex amplitudes
in `Z`, to be interpreted in the eigenbasis of ``𝒮̂ᶻ``. That is, `Z[1]`
represents the amplitude for the basis state fully polarized along the
``ẑ``-direction, and subsequent components represent states with decreasing
angular momentum along this axis (``m = S, S-1, …, -S``).
"""
function set_coherent_state!(sys::System{N}, Z, site) where N
    length(Z) != N && error("Length of coherent state does not match system.")
    iszero(N)      && error("Cannot set zero-length coherent state.")
    site = to_cartesian(site)
    Z = normalize_ket(CVec{N}(Z), sys.κs[site])
    sys.coherents[site] = Z
    sys.dipoles[site] = expected_spin(Z)
end


"""
    polarize_spin!(sys::System, dir, site::Site)

Polarize the spin at a [`Site`](@ref) along the direction `dir`.
"""
function polarize_spin!(sys::System{N}, dir, site) where N
    site = to_cartesian(site)
    setspin!(sys, dipolarspin(sys, site, dir), site)
end

"""
    polarize_spins!(sys::System, dir)

Polarize all spins in the system along the direction `dir`.
"""
function polarize_spins!(sys::System{N}, dir) where N
    for site in all_sites(sys)
        polarize_spin!(sys, dir, site)
    end
end


function get_dipole_buffers(sys::System{N}, numrequested) where N
    numexisting = length(sys.dipole_buffers)
    if numexisting < numrequested
        for _ in 1:(numrequested-numexisting)
            push!(sys.dipole_buffers, zero(sys.dipoles))
        end
    end
    return sys.dipole_buffers[1:numrequested]
end

function get_coherent_buffers(sys::System{N}, numrequested) where N
    numexisting = length(sys.coherent_buffers)
    if numexisting < numrequested
        for _ in 1:(numrequested-numexisting)
            push!(sys.coherent_buffers, zero(sys.coherents))
        end
    end
    return sys.coherent_buffers[1:numrequested]
end
