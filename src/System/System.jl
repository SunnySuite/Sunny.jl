"""
    System(crystal::Crystal, moments, mode; dims=(1, 1, 1), seed=nothing)

A spin system is constructed from the [`Crystal`](@ref) unit cell, a
specification of the spin `moments` symmetry-distinct sites, and a calculation
`mode`. Interactions can be added to the system using, e.g.,
[`set_exchange!`](@ref). The default supercell dimensions are 1×1×1 chemical
cells, but this can be changed with `dims`.

Spin `moments` comprise a list of pairs, `[i1 => Moment(...), i2 => ...]`, where
`i1, i2, ...` are a complete set of symmetry-distinct atoms. Each
[`Moment`](@ref) contains spin and ``g``-factor information.

The two primary options for `mode` are `:SUN` and `:dipole`. In the former, each
spin-``s`` degree of freedom is described as an SU(_N_) coherent state, i.e. a
quantum superposition of ``N = 2s + 1`` levels. This formalism can be useful to
capture multipolar spin fluctuations or local entanglement effects. 

Mode `:dipole` projects the SU(_N_) dynamics onto the restricted space of pure
dipoles. In practice this means that Sunny will simulate Landau-Lifshitz
dynamics, but single-ion anisotropy and biquadratic exchange interactions will
be renormalized to improve accuracy. To disable this renormalization, use the
mode `:dipole_uncorrected`, which corresponds to the formal ``s → ∞`` limit. For
details, see the documentation page: [Interaction Renormalization](@ref).

Stochastic operations on this system can be made reproducible by selecting a
`seed` for this system's pseudo-random number generator. The default system seed
will be generated from Julia's task-local RNG, which can itself be seeded using
`Random.seed!`.

All spins are initially polarized in the global ``z``-direction.
"""
function System(crystal::Crystal, moments::Vector{Pair{Int, Moment}}, mode::Symbol;
                dims::NTuple{3,Int}=(1, 1, 1), seed=nothing, units=nothing)
    if !isnothing(units)
        @warn "units argument to System is deprecated and will be ignored!"
    end
    if mode in (:dipole_large_S, :dipole_large_s)
        @warn "Deprecation warning! Use :dipole_uncorrected instead of $mode"
        mode = :dipole_uncorrected
    end

    if !in(mode, (:SUN, :dipole, :dipole_uncorrected))
        error("Mode must be `:SUN`, `:dipole`, or `:dipole_uncorrected`.")
    end

    # The lattice vectors of `crystal` must be conventional (`crystal` cannot be
    # reshaped).
    if !isnothing(crystal.root)
        @assert crystal.latvecs == crystal.root.latvecs
    end

    na = natoms(crystal)

    moments = propagate_moments(crystal, moments)
    Ss = [m.s for m in moments]
    gs = [m.g for m in moments]

    # TODO: Label SU(2) rep instead
    Ns = @. Int(2Ss+1)

    if mode == :SUN
        allequal(Ns) || error("Currently all spins S must be equal in SU(N) mode.")
        N = first(Ns)
        κs = fill(1.0, na)
    elseif mode in (:dipole, :dipole_uncorrected)
        N = 0 # marker for :dipole mode
        κs = copy(Ss)
    end

    Ns = reshape(Ns, 1, 1, 1, :)
    κs = reshape(κs, 1, 1, 1, :)
    gs = reshape(gs, 1, 1, 1, :)

    interactions = empty_interactions(mode, na, N)
    ewald = nothing

    extfield = zeros(Vec3, 1, 1, 1, na)
    dipoles = fill(zero(Vec3), 1, 1, 1, na)
    coherents = fill(zero(CVec{N}), 1, 1, 1, na)
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]

    rng = isnothing(seed) ? Random.Xoshiro(rand(UInt64, 4)...) : Random.Xoshiro(seed)

    ret = System(nothing, mode, crystal, (1, 1, 1), Ns, κs, gs, interactions, ewald,
                 extfield, dipoles, coherents, dipole_buffers, coherent_buffers, rng)
    polarize_spins!(ret, (0,0,1))
    return dims == (1, 1, 1) ? ret : repeat_periodically(ret, dims)
end

function mode_to_str(sys::System{N}) where N
    if sys.mode == :SUN
        return "[SU($N)]"
    elseif sys.mode == :dipole
        return "[Dipole mode]"
    elseif sys.mode == :dipole_uncorrected
        return "[Dipole mode, large-s]"
    else
        error()
    end
end

function supercell_to_str(dims, cryst)
    return "Supercell (" * join(dims, "×") * ")×" * string(natoms(cryst))
end

function energy_to_str(sys::System)
    if is_valid_normalization(sys)
        "Energy per site "*number_to_math_string(energy_per_site(sys))
    else
        "[Incorrectly normalized spin state!]"
    end
end

function Base.show(io::IO, sys::System{N}) where N
    print(io, "System($(mode_to_str(sys)), $(supercell_to_str(sys.dims, sys.crystal)), $(energy_to_str(sys)))")
end

function Base.show(io::IO, ::MIME"text/plain", sys::System{N}) where N
    printstyled(io, "System $(mode_to_str(sys))\n"; bold=true, color=:underline)
    println(io, supercell_to_str(sys.dims, sys.crystal))
    if !isnothing(sys.origin) && cell_shape(sys) != cell_shape(sys.origin)
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
    (; origin, mode, crystal, dims, Ns, gs, κs, extfield, interactions_union, ewald, dipoles, coherents, rng) = sys

    origin_clone = isnothing(origin) ? nothing : clone_system(origin)
    ewald_clone = nothing # TODO: use clone_ewald(ewald)

    # Dynamically dispatch to the correct `map` function for either homogeneous
    # (Vector) or inhomogeneous interactions (4D Array)
    interactions_clone = map(clone_interactions, interactions_union)
    
    # Empty buffers are required for thread safety.
    empty_dipole_buffers = Array{Vec3, 4}[]
    empty_coherent_buffers = Array{CVec{N}, 4}[]

    ret = System(origin_clone, mode, crystal, dims, Ns, copy(κs), copy(gs),
                 interactions_clone, ewald_clone, copy(extfield), copy(dipoles), copy(coherents),
                 empty_dipole_buffers, empty_coherent_buffers, copy(rng))

    if !isnothing(ewald)
        # At the moment, clone_ewald is unavailable, so instead rebuild the
        # Ewald data structures from scratch. This might be fixed eventually.
        # See https://github.com/JuliaMath/FFTW.jl/issues/261.
        enable_dipole_dipole!(ret, ewald.μ0_μB²)
    end

    return ret
end


"""
    (cell1, cell2, cell3, i) :: Site

Four indices identifying a single site in a [`System`](@ref). The first three
indices select the unit cell and the last index selects the sublattice, i.e.,
the ``i``th atom within the unit cell.

This object can be used to index `dipoles` and `coherents` fields of a `System`.
A `Site` is also required to specify inhomogeneous interactions via functions
such as [`set_field_at!`](@ref) or [`set_exchange_at!`](@ref).

Note that the definition of a cell may change when a system is reshaped. In this
case, it is convenient to construct the `Site` using [`position_to_site`](@ref),
which always takes a position in fractional coordinates of the original lattice
vectors.
"""
const Site = Union{NTuple{4, Int}, CartesianIndex{4}}

@inline to_cartesian(i::CartesianIndex{N}) where N = i
@inline to_cartesian(i::NTuple{N, Int})    where N = CartesianIndex(i)

# Split a site `site` into its cell and sublattice parts
@inline to_cell(site) = (site[1], site[2], site[3])
@inline to_atom(site) = site[4]

# Like mod1(x, L), but short-circuits early in the common case. See
# https://github.com/SunnySuite/Sunny.jl/pull/184 for discussion.
@inline function altmod1(x::Int, L::Int)
    1 <= x <= L ? x : mod1(x, L)
end

# Get the neighboring site associated with site and bond
@inline function bonded_site(site, bond, dims)
    # @assert to_atom(site) == bond.i
    CartesianIndex(altmod1.(to_cell(site) .+ bond.n, dims)..., bond.j)
end

"""
    spin_label(sys::System, i::Int)

If atom `i` carries a single spin-``s`` moment, then returns the half-integer
label ``s``. Otherwise, throws an error.
"""
function spin_label(sys::System, i::Int)
    if sys.mode == :dipole_uncorrected
        return Inf
    else
        @assert sys.mode in (:dipole, :SUN)
        allequal(sys.Ns[:,:,:,i]) || error("Spin varies between chemical cells")
        return (sys.Ns[1,1,1,i]-1)/2
    end
end


"""
    eachsite(sys::System)
    eachsite(sys::System, i)

An iterator over all [`Site`](@ref)s in the system. Restrict to one sublattice
`i` with an optional second argument.
"""
@inline eachsite(sys::System) = CartesianIndices(size(sys.dipoles))
@inline eachsite(sys::System, i) = CartesianIndices((sys.dims..., i:i))

"""
nsites(sys::System) = length(eachsite(sys))
"""
nsites(sys::System) = length(eachsite(sys))

# Number of (original) crystal cells in the system
ncells(sys::System) = nsites(sys) / natoms(orig_crystal(sys))


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

Returns ``- g 𝐒``, the local magnetic moment in units of the Bohr magneton. The
spin dipole ``𝐒`` and ``g``-tensor may both be [`Site`](@ref) dependent.
"""
function magnetic_moment(sys::System, site)
    site = to_cartesian(site)
    return - sys.gs[site] * sys.dipoles[site]
end

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.dims)

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
    cell = @. mod1(offset+1, sys.dims) # 1-based indexing with periodicity
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
represented as a pair of [`Site`](@ref)s and an `offset`, which may be used as
input to [`set_exchange_at!`](@ref) or [`set_pair_coupling_at!`](@ref). Reverse
bonds are not included in the iterator (no double counting).

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
            for site_i in eachsite(sys, new_bond.i)
                site_j = bonded_site(site_i, new_bond, sys.dims)
                site_i < site_j && push!(ret, (site_i, site_j, new_bond.n))
            end
        end
    end

    return ret
end


struct SpinState{N}
    S::Vec3
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
    S = expected_spin(Z)
    return SpinState(S, Z)
end

@inline function dipolar_state(sys::System{0}, site, dir)
    S = normalize_dipole(Vec3(dir), sys.κs[site])
    Z = CVec{0}()
    return SpinState(S, Z)
end
@inline function dipolar_state(sys::System{N}, site, dir) where N
    return coherent_state(sys, site, ket_from_dipole(Vec3(dir), Val(N)))
end

@inline function flip(spin::SpinState{N}) where N
    return SpinState(-spin.S, flip_ket(spin.Z))
end

@inline function randspin(sys::System{0}, site)
    S = normalize_dipole(randn(sys.rng, Vec3), sys.κs[site])
    return SpinState(S, CVec{0}())
end
@inline function randspin(sys::System{N}, site) where N
    Z = normalize_ket(randn(sys.rng, CVec{N}), sys.κs[site])
    S = expected_spin(Z)
    return SpinState(S, Z)
end

@inline function getspin(sys::System{N}, site) where N
    return SpinState(sys.dipoles[site], sys.coherents[site])
end

@inline function setspin!(sys::System{N}, spin::SpinState{N}, site) where N
    sys.dipoles[site] = spin.S
    sys.coherents[site] = spin.Z
    return
end

function is_valid_normalization(sys::System{0})
    all(zip(sys.dipoles, sys.κs)) do (S, κ)
        norm2(S) ≈ κ^2
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

For a single quantum spin-``s``, these amplitudes will be interpreted in the
eigenbasis of ``Ŝ^z``. That is, `Z[1]` represents the amplitude for the basis
state fully polarized along the ``ẑ``-direction, and subsequent components
represent states with decreasing angular momentum along this axis (``m = s, s-1,
…, -s``).
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
