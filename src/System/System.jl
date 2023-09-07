"""
    System(crystal::Crystal, latsize, infos, mode; units=Units.meV, seed::Int)

Construct a `System` of spins for a given [`Crystal`](@ref) symmetry. The
`latsize` parameter determines the number of unit cells in each lattice vector
direction. The `infos` parameter is a list of [`SpinInfo`](@ref) objects, which
determine the magnitude ``S`` and ``g``-tensor of each spin.

The two possible options for `mode` are `:SUN` and `:dipole`. In the former,
each spin-``S`` degree of freedom is described as an SU(_N_) coherent state,
i.e. a quantum superposition of ``N = 2S + 1`` levels. This approach can be
important, e.g., to capture multipolar spin fluctuations when there is a strong
single-ion anisotropy, or to explicitly resolve spin-orbit coupling. 

Mode `:dipole` projects the SU(_N_) dynamics onto the restricted space of pure
dipoles. In practice this means that Sunny will simulate Landau-Lifshitz
dynamics, but all single-ion anisotropy and biquadratic exchange interactions
will be renormalized for maximum accuracy. To disable this renormalization
(e.g., to model systems in the large-spin limit) construct anisotropy operators
using the special function [`large_S_spin_operators`](@ref) and define any
biquadratic exchange using the `large_S=true` in [`set_exchange!`](@ref).

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

An optional `seed` may be provided to achieve reproducible random number
generation.

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, infos::Vector{SpinInfo}, mode::Symbol;
                    units=Units.meV, seed=nothing)
    if !(mode == :SUN || mode == :dipole)
        error("Mode must be `:SUN` or `:dipole`.")
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
    elseif mode == :dipole
        N = 0 # acts as marker for :dipole
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

function Base.show(io::IO, sys::System{N}) where N
    modename = if sys.mode==:SUN
        "SU($N)"
    elseif sys.mode==:dipole
        "Dipole"
    else
        error()
    end
    print(io,"System{$modename}[$(sys.latsize)×$(natoms(sys.crystal))]")
    if !isnothing(sys.origin)
        print(io,"[Reshape = $(cell_dimensions(sys))]")
    end
end

function Base.show(io::IO, ::MIME"text/plain", sys::System{N}) where N
    modename = if sys.mode==:SUN
        "SU($N)"
    elseif sys.mode==:dipole
        "Dipole mode"
    else
        error()
    end
    printstyled(io, "System [$modename]\n"; bold=true, color=:underline)
    println(io, "Lattice: $(sys.latsize)×$(natoms(sys.crystal))")
    if !isnothing(sys.origin)
        #println(io, "Reshaped cell geometry $(cell_dimensions(sys))")
        is_enlarged = abs(det(sys.crystal.latvecs)) > abs(det(sys.origin.crystal.latvecs))
        println(io)
        println(io, "Unit cell has been $(is_enlarged ? "enlarged" : "reshaped") from original such that:")
        printstyled(io, "  [Original lattice vectors]";bold=true,color=:red)
        print(io," * $(cell_dimensions(sys)) = ")
        printstyled(io,"[Reshaped lattice vectors]\n";bold=true,color=:blue)
        println(io, "where")
        printstyled(io, "  [Original] ";bold=true,color=:red)
        show(io, sys.origin.crystal)
        printstyled(io, "\n  [Reshaped] ";bold=true,color=:blue)
        show(io, sys.crystal)
        println(io)
    else
        show(io, sys.crystal)
    end
    println(io)

    if is_homogeneous(sys)
        ints = interactions_homog(sys)
        if isempty(ints)
            println(io, "No interactions")
        else
            print(io, "Homogeneous interactions by atom:\n")
            for (i,int) in enumerate(ints)
                if sys.crystal.types[i] != ""
                    print(io,"  $i. '$(sys.crystal.types[i])' atom has ")
                else
                    print(io,"  Atom $i has ")
                end
                show(io,int)
                println(io)
            end
        end
    else
        print(io, "Inhomogeneous interactions (may differ at every site)")
    end

    if !isnothing(sys.ewald)
        println(io, "Long range dipole-dipole interactions enabled!")
    end

    if !iszero(sys.extfield)
        if allequal(sys.extfield)
            println(io, "Uniform magnetic field B = $(sys.extfield[1]) applied")
        else
            mean_field = sum(sys.extfield) ./ length(sys.extfield)
            rms_field = sqrt.(sum([(B .- mean_field) .^ 2 for B in sys.extfield]) ./ length(sys.extfield))

            mean_field_str = @sprintf "[%.4g %.4g %.4g]" mean_field[1] mean_field[2] mean_field[3]
            rms_field_str = @sprintf "[%.4g %.4g %.4g]" rms_field[1] rms_field[2] rms_field[3]

            println(io, "Spatially periodic magnetic field B = (mean $mean_field_str ± $rms_field_str RMS) applied")
        end
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
const Site = Union{NTuple{4, Int}, CartesianIndex{4}}

@inline to_cartesian(i::CartesianIndex{N}) where N = i
@inline to_cartesian(i::NTuple{N, Int})    where N = CartesianIndex(i)

# kbtodo: offsetcell ?
# Offset a `cell` by `ncells`
@inline offsetc(cell::CartesianIndex{3}, ncells, latsize) = CartesianIndex(mod1.(Tuple(cell) .+ Tuple(ncells), latsize))

# Split a site `site` into its cell and sublattice parts
@inline to_cell(site) = CartesianIndex((site[1],site[2],site[3]))
@inline to_atom(site) = site[4]

# An iterator over all unit cells using CartesianIndices
@inline eachcell(sys::System) = CartesianIndices(sys.latsize)

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

Get the magnetic moment for a [`Site`](@ref). This is the spin dipole multiplied
by the Bohr magneton and the local g-tensor.
"""
magnetic_moment(sys::System, site) = sys.units.μB * sys.gs[site] * sys.dipoles[site]

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)

# The original crystal for a system, invariant under reshaping
orig_crystal(sys) = something(sys.origin, sys).crystal

# Position of a site in fractional coordinates of the original crystal
function position(sys::System, site)
    return orig_crystal(sys).latvecs \ global_position(sys, site)
end

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
[`set_exchange_at!`](@ref). Reverse bonds are not included (no double counting).

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
        i = map_atom_to_crystal(sys.crystal, new_i, orig_crystal(sys))

        # loop over symmetry equivalent bonds in original crystal
        for bond′ in all_symmetry_related_bonds_for_atom(orig_crystal(sys), i, bond)

            # map to a bond with indexing for new crystal
            new_bond = transform_bond(sys.crystal, new_i, orig_crystal(sys), bond′)

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
    length(Z) != N && error("Length of coherent state does not match system.")
    iszero(N)      && error("Cannot set zero-length coherent state.")
    setspin!(sys, coherent_state(sys, site, Z), to_cartesian(site))
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


function spin_operators_pair(sys::System{N}, b::Bond) where N
    Si = spin_matrices(N=sys.Ns[b.i])
    Sj = spin_matrices(N=sys.Ns[b.j])
    return local_quantum_operators(Si, Sj)
end
