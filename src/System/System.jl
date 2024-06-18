"""
    System(crystal::Crystal, latsize, infos, mode; units=Units.meV, seed::Int)

Construct a `System` of spins for a given [`Crystal`](@ref) symmetry. The
`latsize` parameter determines the number of unit cells in each lattice vector
direction. The `infos` parameter is a list of [`SpinInfo`](@ref) objects, which
determine the magnitude ``S`` and ``g``-tensor of each spin.

The two primary options for `mode` are `:SUN` and `:dipole`. In the former, each
spin-``S`` degree of freedom is described as an SU(_N_) coherent state, i.e. a
quantum superposition of ``N = 2S + 1`` levels. This formalism can be useful to
capture multipolar spin fluctuations or local entanglement effects. 

Mode `:dipole` projects the SU(_N_) dynamics onto the restricted space of pure
dipoles. In practice this means that Sunny will simulate Landau-Lifshitz
dynamics, but single-ion anisotropy and biquadratic exchange interactions will
be renormalized to improve accuracy. To disable this renormalization, use the
mode `:dipole_large_S` which applies the ``S → ∞`` classical limit. For details,
see the documentation page: [Interaction Strength Renormalization](@ref).

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

An optional `seed` may be provided to achieve reproducible random number
generation.

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, infos::Vector{SpinInfo}, mode::Symbol;
                units=Units.meV, seed=nothing)
    if !in(mode, (:SUN, :dipole, :dipole_large_S))
        error("Mode must be `:SUN`, `:dipole`, or `:dipole_large_S`.")
    end

    # The lattice vectors of `crystal` must be conventional (`crystal` cannot be
    # reshaped).
    if !isnothing(crystal.root)
        @assert crystal.latvecs == crystal.root.latvecs
    end
    
    na = natoms(crystal)

    infos = propagate_site_info(crystal, infos)
    Ss = [si.S for si in infos]
    gs = [si.g for si in infos]

    # TODO: Label SU(2) rep instead
    Ns = @. Int(2Ss+1)

    if mode == :SUN
        allequal(Ns) || error("Currently all spins S must be equal in SU(N) mode.")
        N = first(Ns)
        κs = fill(1.0, na)
    elseif mode in (:dipole, :dipole_large_S)
        N = 0 # marker for :dipole mode
        κs = copy(Ss)
    end

    # Repeat such that `A[:]` → `A[cell, :]` for every `cell`
    repeat_to_lattice(A) = permutedims(repeat(A, 1, latsize...), (2, 3, 4, 1))

    Ns = repeat_to_lattice(Ns)
    κs = repeat_to_lattice(κs)
    gs = repeat_to_lattice(gs)

    interactions = empty_interactions(mode, na, N)
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

function mode_to_str(sys::System{N}) where N
    if sys.mode == :SUN
        return "[SU($N)]"
    elseif sys.mode == :dipole
        return "[Dipole mode]"
    elseif sys.mode == :dipole_large_S
        return "[Dipole mode, large-S]"
    else
        error()
    end
end

function lattice_to_str(sys::System)
    lat_str = isnothing(sys.origin) ? "Lattice" : "Reshaped lattice"
    return lat_str * " ($(join(sys.latsize, "×")))×$(natoms(sys.crystal))"
end

function energy_to_str(sys::System)
    if is_valid_normalization(sys)
        "Energy per site "*number_to_math_string(energy_per_site(sys))
    else
        "[Incorrectly normalized spin state!]"
    end
end

function Base.show(io::IO, sys::System{N}) where N
    print(io, "System($(mode_to_str(sys)), $(lattice_to_str(sys)), $(energy_to_str(sys)))")
end

function Base.show(io::IO, ::MIME"text/plain", sys::System{N}) where N
    printstyled(io, "System $(mode_to_str(sys))\n"; bold=true, color=:underline)
    println(io, lattice_to_str(sys))
    if !isnothing(sys.origin)
        shape = number_to_math_string.(cell_shape(sys))
        println(io, formatted_matrix(shape; prefix="Reshaped cell "))
    end
    println(io, energy_to_str(sys))
end

# Per Julia developers, `deepcopy` is memory unsafe, especially in conjunction
# with C libraries. We were observing very confusing crashes that surfaced in
# the FFTW library, https://github.com/JuliaLang/julia/issues/48722. To prevent
# this from happening again, avoid all uses of `deepcopy`, and create our own
# stack of `clone` functions instead.
Base.deepcopy(_::System) = error("Use `clone_system` instead of `deepcopy`.")

"""
    clone_system(sys::System)

Creates a full clone of the system, such that mutable updates to one copy will
not affect the other, and thread safety is guaranteed.
"""
function clone_system(sys::System{N}) where N
    (; origin, mode, crystal, latsize, Ns, gs, κs, extfield, interactions_union, ewald, dipoles, coherents, units, rng) = sys

    origin_clone = isnothing(origin) ? nothing : clone_system(origin)
    ewald_clone = nothing # TODO: use clone_ewald(ewald)

    # Dynamically dispatch to the correct `map` function for either homogeneous
    # (Vector) or inhomogeneous interactions (4D Array)
    interactions_clone = map(clone_interactions, interactions_union)
    
    # Empty buffers are required for thread safety.
    empty_dipole_buffers = Array{Vec3, 4}[]
    empty_coherent_buffers = Array{CVec{N}, 4}[]

    ret = System(origin_clone, mode, crystal, latsize, Ns, copy(κs), copy(gs),
                 interactions_clone, ewald_clone, copy(extfield), copy(dipoles), copy(coherents),
                 empty_dipole_buffers, empty_coherent_buffers, units, copy(rng))

    if !isnothing(ewald)
        # At the moment, clone_ewald is unavailable, so instead rebuild the
        # Ewald data structures from scratch. This might be fixed eventually.
        # See https://github.com/JuliaMath/FFTW.jl/issues/261.
        enable_dipole_dipole!(ret)
    end

    return ret
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
const Site = Union{NTuple{4, Int}, CartesianIndex{4}}

@inline to_cartesian(i::CartesianIndex{N}) where N = i
@inline to_cartesian(i::NTuple{N, Int})    where N = CartesianIndex(i)

# Like mod1(x, L), but short-circuits early in the common case. See
# https://github.com/SunnySuite/Sunny.jl/pull/184 for discussion.
@inline function altmod1(x::Int, L::Int)
    1 <= x <= L ? x : mod1(x, L)
end

# Offset a `cell` by `ncells`
@inline offsetc(cell::CartesianIndex{3}, ncells, latsize) = CartesianIndex(altmod1.(Tuple(cell) .+ Tuple(ncells), latsize))

# Split a site `site` into its cell and sublattice parts
@inline to_cell(site) = CartesianIndex((site[1],site[2],site[3]))
@inline to_atom(site) = site[4]

# An iterator over all unit cells using CartesianIndices
@inline eachcell(sys::System) = CartesianIndices(sys.latsize)

"""
    spin_label(sys::System, i::Int)

If atom `i` carries a single spin-``S`` moment, then returns the half-integer
label ``S``. Otherwise, throws an error.
"""
function spin_label(sys::System, i::Int)
    if sys.mode == :dipole_large_S
        return Inf
    else
        @assert sys.mode in (:dipole, :SUN)
        allequal(sys.Ns[:,:,:,i]) || error("Spin varies between chemical cells")
        return (sys.Ns[1,1,1,i]-1)/2
    end
end


"""
    eachsite(sys::System)

An iterator over all [`Site`](@ref)s in the system. 
"""
@inline eachsite(sys::System) = CartesianIndices(sys.dipoles)

"""
    global_position(sys::System, site::Site)

Position of a [`Site`](@ref) in global coordinates.

To precompute a full list of positions, one can use [`eachsite`](@ref) as
below:

```julia
pos = [global_position(sys, site) for site in eachsite(sys)]
```
"""
function global_position(sys::System, site)
    r = sys.crystal.positions[site[4]] + Vec3(site[1]-1, site[2]-1, site[3]-1)
    return sys.crystal.latvecs * r
end

"""
    magnetic_moment(sys::System, site::Site)

Get the magnetic moment for a [`Site`](@ref). This is the spin dipole ``𝐒``
multiplied by the local ``g``-tensor and the Bohr magneton ``μ_B``.

```math
μ = - μ_B g 𝐒.
```

The constant ``μ_B`` depends on the choice of [`Units`](@ref). By default,
energy is in meV and external field is in tesla, such that ``μ_B = 0.05788…``
(meV/T). If `Units.theory` is selected instead, the constant ``μ_B = -1``
effectively aligns magnetic and spin dipoles, which is a common convention when
specifying theoretical model systems.
"""
function magnetic_moment(sys::System, site)
    site = to_cartesian(site)
    return - sys.units.μB * sys.gs[site] * sys.dipoles[site]
end

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)

# The crystal originally used to construct a system. It is guaranteed to be
# un-reshaped, and its lattice vectors define the "conventional" unit cell. It
# may, however, be a subcrystal of `orig_crystal(sys).root`.
orig_crystal(sys) = something(sys.origin, sys).crystal

"""
    position(sys::System, site::Site)

Position of a [`Site`](@ref) in units of lattice vectors for the original
crystal.
"""
position(sys::System, site) = orig_crystal(sys).latvecs \ global_position(sys, site)

"""
    position_to_site(sys::System, r)

Converts a position `r` to four indices of a [`Site`](@ref). The coordinates of
`r` are given in units of the lattice vectors for the original crystal. This
function can be useful for working with systems that have been reshaped using
[`reshape_supercell`](@ref).

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


# Given a [`Site`](@ref)s for a possibly reshaped system, return the
# corresponding atom index for the original (unreshaped) crystal.
function site_to_atom(sys::System{N}, site) where N
    site = to_cartesian(site)
    r = position(sys, site)
    return position_to_atom(orig_crystal(sys), r)
end

# Maps atom `i` in `cryst` to the corresponding atom in `other_cryst`
function map_atom_to_other_crystal(cryst::Crystal, i, other_cryst::Crystal)
    global_r = cryst.latvecs * cryst.positions[i]
    orig_r = other_cryst.latvecs \ global_r
    return position_to_atom(other_cryst, orig_r)
end

# Maps atom `i` in `cryst` to the corresponding site in `other_sys`
function map_atom_to_other_system(cryst::Crystal, i, other_sys::System)
    global_r = cryst.latvecs * cryst.positions[i]
    other_r = orig_crystal(other_sys).latvecs \ global_r
    return position_to_site(other_sys, other_r)
end


# Given a `bond` for `cryst`, return a corresponding new bond for the reshaped
# `other_cryst`. The new bond will begin at atom `other_i`.
function map_bond_to_other_crystal(cryst::Crystal, bond::Bond, other_cryst::Crystal, other_i::Int)
    # Positions in new fractional coordinates
    other_ri = other_cryst.positions[other_i]
    other_rj = other_ri + other_cryst.latvecs \ global_displacement(cryst, bond)

    # Verify that new_i (indexed into new_cryst) is consistent with bond.i
    # (indexed into original cryst).
    @assert bond.i == position_to_atom(cryst, cryst.latvecs \ other_cryst.latvecs * other_ri)

    # Construct bond using new indexing system
    other_j, other_n = position_to_atom_and_offset(other_cryst, other_rj)
    return Bond(other_i, other_j, other_n)
end


"""
    symmetry_equivalent_bonds(sys::System, bond::Bond)

Given a [`Bond`](@ref) for the original (unreshaped) crystal, return all
symmetry equivalent bonds in the [`System`](@ref). Each returned bond is
represented as a pair of [`Site`](@ref)s, which may be used as input to
[`set_exchange_at!`](@ref) or [`set_pair_coupling_at!`](@ref). Reverse bonds are
not included in the iterator (no double counting).

# Example
```julia
for (site1, site2, offset) in symmetry_equivalent_bonds(sys, bond)
    @assert site1 < site2
    set_exchange_at!(sys, J, site1, site2; offset)
end
```
"""
function symmetry_equivalent_bonds(sys::System, bond::Bond)
    ret = Tuple{Site, Site, SVector{3, Int}}[]

    for new_i in 1:natoms(sys.crystal)
        # atom index in original crystal
        i = map_atom_to_other_crystal(sys.crystal, new_i, orig_crystal(sys))

        # loop over symmetry equivalent bonds in original crystal
        for bond′ in all_symmetry_related_bonds_for_atom(orig_crystal(sys), i, bond)

            # map to a bond with indexing for new crystal
            new_bond = map_bond_to_other_crystal(orig_crystal(sys), bond′, sys.crystal, new_i)

            # loop over all new crystal cells and push site pairs
            for new_cell_i in eachcell(sys)
                new_cell_j = offsetc(new_cell_i, new_bond.n, sys.latsize)
                site_i = (Tuple(new_cell_i)..., new_bond.i)
                site_j = (Tuple(new_cell_j)..., new_bond.j)
                site_i < site_j && push!(ret, (site_i, site_j, new_bond.n))
            end
        end
    end

    return ret
end

"""
    remove_ion_at!(sys::System, site::Site)

Remove all interactions associated with the magnetic ion at `site`. The system
must support inhomogeneous interactions via [`to_inhomogeneous`](@ref).
"""
function remove_ion_at!(sys::System{N}, site) where N
    is_homogeneous(sys) && error("Use `to_inhomogeneous` first.")

    site = to_cartesian(site)

    # Remove onsite coupling
    ints = interactions_inhomog(sys)
    ints[site].onsite = empty_anisotropy(sys.mode, N)

    # Remove this site from neighbors' pair lists
    for (; bond) in ints[site].pair
        cell′ = offsetc(to_cell(site), bond.n, sys.latsize)
        pair′ = ints[cell′, bond.j].pair
        deleteat!(pair′, only(findall(pc′ -> pc′.bond == reverse(bond), pair′)))
    end

    # Empty pair list from this site
    empty!(ints[site].pair)

    # Set g=0 to disable Zeeman and dipole-dipole interactions
    sys.gs[site] = Mat3(0I)

    # Set κ=0 to disable all expectation values
    sys.κs[site] = 0
    sys.dipoles[site] = zero(Vec3)
    sys.coherents[site] = zero(CVec{N})
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

@inline function coherent_state(sys::System{N}, site, Z) where N
    Z = normalize_ket(CVec{N}(Z), sys.κs[site])
    s = expected_spin(Z)
    return SpinState(s, Z)
end

@inline function dipolar_state(sys::System{0}, site, dir)
    s = normalize_dipole(Vec3(dir), sys.κs[site])
    Z = CVec{0}()
    return SpinState(s, Z)
end
@inline function dipolar_state(sys::System{N}, site, dir) where N
    return coherent_state(sys, site, ket_from_dipole(Vec3(dir), Val(N)))
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

@inline function getspin(sys::System{N}, site) where N
    return SpinState(sys.dipoles[site], sys.coherents[site])
end

@inline function setspin!(sys::System{N}, spin::SpinState{N}, site) where N
    sys.dipoles[site] = spin.s
    sys.coherents[site] = spin.Z
    return
end

function is_valid_normalization(sys::System{0})
    all(zip(sys.dipoles, sys.κs)) do (s, κ)
        norm2(s) ≈ κ^2
    end
end
function is_valid_normalization(sys::System{N}) where N
    all(zip(sys.coherents, sys.κs)) do (Z, κ)
        norm2(Z) ≈ κ
    end
end

function validate_normalization(sys::System)
    is_valid_normalization(sys) || error("Detected non-normalized spin state.")
end

"""
    randomize_spins!(sys::System)

Randomizes all spins under appropriate the uniform distribution.
"""
function randomize_spins!(sys::System{N}) where N
    for site in eachsite(sys)
        setspin!(sys, randspin(sys, site), site)
    end
end


"""
    set_coherent!(sys::System, Z, site::Site)

Set a coherent spin state at a [`Site`](@ref) using the ``N`` complex amplitudes
in `Z`.

For a standard [`SpinInfo`](@ref), these amplitudes will be interpreted in the
eigenbasis of ``𝒮̂ᶻ``. That is, `Z[1]` represents the amplitude for the basis
state fully polarized along the ``ẑ``-direction, and subsequent components
represent states with decreasing angular momentum along this axis (``m = S, S-1,
…, -S``).
"""
function set_coherent!(sys::System{N}, Z, site) where N
    site = to_cartesian(site)
    length(Z) != N && error("Length of coherent state does not match system.")
    iszero(N)      && error("Cannot set zero-length coherent state.")
    setspin!(sys, coherent_state(sys, site, Z), site)
end


"""
    set_dipole!(sys::System, dir, site::Site)

Polarize the spin at a [`Site`](@ref) along the direction `dir`.
"""
function set_dipole!(sys::System{N}, dir, site) where N
    site = to_cartesian(site)
    setspin!(sys, dipolar_state(sys, site, dir), site)
end

"""
    polarize_spins!(sys::System, dir)

Polarize all spins in the system along the direction `dir`.
"""
function polarize_spins!(sys::System{N}, dir) where N
    for site in eachsite(sys)
        set_dipole!(sys, dir, site)
    end
end

"""
    set_spin_rescaling!(sys, α)

In dipole mode, rescale all spin magnitudes ``S → α S``. In SU(N) mode, rescale
all SU(N) coherent states ``Z → √α Z`` such that every expectation value
rescales like ``⟨A⟩ → α ⟨A⟩``.
"""
function set_spin_rescaling!(sys::System{0}, α)
    for site in eachsite(sys)
        sys.κs[site] = α * (sys.Ns[site]-1)/2
        set_dipole!(sys, sys.dipoles[site], site)
    end
end

function set_spin_rescaling!(sys::System{N}, κ) where N
    sys.κs .= κ
    for site in eachsite(sys)
        set_coherent!(sys, sys.coherents[site], site)
    end
end


function get_dipole_buffers(sys::System{N}, numrequested) where N
    numexisting = length(sys.dipole_buffers)
    if numexisting < numrequested
        for _ in 1:(numrequested-numexisting)
            push!(sys.dipole_buffers, zero(sys.dipoles))
        end
    end
    return view(sys.dipole_buffers, 1:numrequested)
end

function get_coherent_buffers(sys::System{N}, numrequested) where N
    numexisting = length(sys.coherent_buffers)
    if numexisting < numrequested
        for _ in 1:(numrequested-numexisting)
            push!(sys.coherent_buffers, zero(sys.coherents))
        end
    end
    return view(sys.coherent_buffers, 1:numrequested)
end
