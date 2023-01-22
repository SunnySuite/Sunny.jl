"""
    System(crystal::Crystal, latsize, infos, mode; units=Units.meV, seed::Int)

Construct a `System` of spins for a given `crystal` symmetry. The `latsize`
parameter determines the number of unit cells in each lattice vector direction.
The `infos` parameter is a list of [`SpinInfo`](@ref) objects, which determine
the magnitude ``S`` and ``g``-tensor of each spin.

The three possible options for `mode` are `:dipole`, `:SUN`, and `:projected`.
The choice `:dipole` restricts each a spin to an angular momentum dipole. In
contrast, the choice `:SUN` describes each spin as a full SU(_N_) coherent
state, where ``N = 2S + 1``. The SU(_N_) approach captures more dynamical
degrees of freedom, e.g.,  quadrupolar moments ``⟨ŜᵅŜᵝ+ŜᵝŜᵅ⟩``. Finally, the
choice `:projected` projects the SU(_N_) dynamics onto the space of dipoles. In
practice, the distinction between `:projected` and `:dipoles` is that the former
will automatically apply appropriate renormalizations to the single-ion
anisotropy and biquadratic exchange interactions for maximum accuracy.

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, infos::Vector{SpinInfo}, mode::Symbol;
                    units=Units.meV, seed=nothing)
    mode==:projected && error("SU(N) projected mode not yet implemented.")

    if mode ∉ [:dipole, :SUN, :projected]
        error("Mode must be one of [:dipole, :SUN, :projected].")
    end

    nb = nbasis(crystal)

    infos = propagate_site_info(crystal, infos)
    Ss = [si.S for si in infos]
    gs = [si.g for si in infos]
    Ns = @. Int(2Ss+1)

    if mode == :SUN
        if !allequal(Ns)
            error("Currently all spins S must be equal in SU(N) mode.")
        end
        N = first(Ns)
        κs = fill(1.0, latsize..., nb)
    else
        N = 0
        # Repeat such that `κs[cell, :] == Ss` for every `cell`
        κs = permutedims(repeat(Ss, 1, latsize...), (2, 3, 4, 1))
    end

    interactions = Interactions(nb, latsize, N)

    dipoles = fill(zero(Vec3), latsize..., nb)
    coherents = fill(zero(CVec{N}), latsize..., nb)
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = System(mode, crystal, latsize, Ns, gs, κs, interactions, dipoles, coherents,
                    dipole_buffers, coherent_buffers, units, rng)
    polarize_spins!(ret, (0,0,1))
    return ret
end

function Base.show(io::IO, ::MIME"text/plain", sys::System{N}) where N
    modename = if sys.mode==:SUN
        "SU($N)"
    elseif sys.mode==:projected
        "Projected SU($N)"
    else
        "Dipole mode"
    end
    printstyled(io, "System [$modename]\n"; bold=true, color=:underline)
    println(io, "Cell size $(nbasis(sys.crystal)), Lattice size $(sys.latsize)")
end


function clone_spin_state(sys::System{N}) where N
    System(sys.mode, sys.crystal, sys.latsize, sys.Ns, sys.gs, sys.κs, sys.interactions,
        copy(sys.dipoles), copy(sys.coherents), sys.dipole_buffers, sys.coherent_buffers,
        sys.units, copy(sys.rng))
end

"""
    extend_periodically(sys::System{N}, mults::NTuple{3, Int64}) where N

Creates a new System identical to `sys` but with each dimension multiplied
by the corresponding factor given in the tuple `mults`. The original spin configuration
is simply repeated periodically.
"""
function extend_periodically(sys::System{N}, factors::NTuple{3, Int64}) where N
    @assert all(>=(1), factors)
    new_latsize = factors .* sys.latsize
    new_κs        = repeat(sys.κs, factors..., 1)
    new_dipoles   = repeat(sys.dipoles, factors..., 1)
    new_coherents = repeat(sys.coherents, factors..., 1)
    return System(sys.mode, sys.crystal, new_latsize, sys.Ns, sys.gs, new_κs, sys.interactions,
                    new_dipoles, new_coherents, Array{Vec3, 4}[], Array{CVec{N}, 4}[],
                    sys.units, copy(sys.rng))
end


"""
    Site(n1, n2, n3, i)

References a single site in a `System` via its unit cell `(n1,n2,n3)` and its
sublattice `i`. Can be used to index `dipoles` and `coherents` fields of a
`System`, or to set inhomogeneous interactions.

See also [`set_vacancy_at!`](@ref), [`set_external_field_at!`](@ref).
"""
@inline Site(idx::CartesianIndex{4})            = idx
@inline Site(idx::NTuple{4, Int})               = CartesianIndex(idx)
@inline Site(n1::Int, n2::Int, n3::Int, b::Int) = CartesianIndex(n1, n2, n3, b)

# Element-wise application of mod1(cell+off, latsize), returning CartesianIndex
@inline offsetc(cell::CartesianIndex{3}, off, latsize) = CartesianIndex(mod1.(Tuple(cell).+Tuple(off), latsize))

# Split a Cartesian index (cell,i) into its parts cell and i.
@inline splitidx(idx::CartesianIndex{4}) = (CartesianIndex((idx[1],idx[2],idx[3])), idx[4])

# An iterator over all sites using CartesianIndices
@inline all_sites(sys::System) = CartesianIndices(sys.dipoles)

# An iterator over all unit cells using CartesianIndices
@inline all_cells(sys::System) = CartesianIndices(sys.latsize)

# Position of a site in global coordinates
position(sys::System, idx) = sys.crystal.lat_vecs * (sys.crystal.positions[idx[4]] .+ (idx[1]-1, idx[2]-1, idx[3]-1))

# Magnetic moment at a site
magnetic_moment(sys::System, idx) = sys.units.μB * sys.gs[idx[4]] * sys.dipoles[idx]

# Positions of all sites in global coordinates
positions(sys::System) = [position(sys, idx) for idx in all_sites(sys)]

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)


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

@inline function getspin(sys::System{N}, idx::CartesianIndex{4}) where N
    return SpinState(sys.dipoles[idx], sys.coherents[idx])
end

@inline function setspin!(sys::System{N}, spin::SpinState{N}, idx::CartesianIndex{4}) where N
    sys.dipoles[idx] = spin.s
    sys.coherents[idx] = spin.Z
    nothing
end

@inline function flip(spin::SpinState{N}) where N
    return SpinState(-spin.s, flip_ket(spin.Z))
end

@inline function randspin(sys::System{0}, idx)
    s = normalize_dipole(randn(sys.rng, Vec3), sys.κs[idx])
    return SpinState(s, CVec{0}())
end
@inline function randspin(sys::System{N}, idx) where N
    Z = normalize_ket(randn(sys.rng, CVec{N}), sys.κs[idx])
    s = expected_spin(Z)
    return SpinState(s, Z)
end

@inline function dipolarspin(sys::System{0}, idx, dir)
    s = normalize_dipole(Vec3(dir), sys.κs[idx])
    Z = CVec{0}()
    return SpinState(s, Z)
end
@inline function dipolarspin(sys::System{N}, idx, dir) where N
    Z = normalize_ket(ket_from_dipole(Vec3(dir), Val(N)), sys.κs[idx])
    s = expected_spin(Z)
    return SpinState(s, Z)
end


function randomize_spins!(sys::System{N}) where N
    for idx in all_sites(sys)
        setspin!(sys, randspin(sys, idx), idx)
    end
end

function polarize_spins!(sys::System{N}, dir) where N
    for idx = all_sites(sys)
        spin = dipolarspin(sys, idx, dir)
        setspin!(sys, spin, idx)
    end
end

"""
    set_vacancy_at!(sys::System, idx::Site)

Make a single site nonmagnetic. [`Site`](@ref) includes a unit cell and a
sublattice index.
"""
function set_vacancy_at!(sys::System{N}, idx) where N
    idx = Site(idx)
    sys.κs[idx] = 0.0
    sys.dipoles[idx] = zero(Vec3)
    sys.coherents[idx] = zero(CVec{N})
end


function get_dipole_buffers(sys::System, numrequested)
    numexisting = length(sys.dipole_buffers)
    if numexisting < numrequested
        for _ in 1:(numrequested-numexisting)
            push!(sys.dipole_buffers, zero(sys.dipoles))
        end
    end
    return sys.dipole_buffers[1:numrequested]
end

function get_coherent_buffers(sys::System, numrequested)
    numexisting = length(sys.coherent_buffers)
    if numexisting < numrequested
        for _ in 1:(numrequested-numexisting)
            push!(sys.coherent_buffers, zero(sys.coherents))
        end
    end
    return sys.coherent_buffers[1:numrequested]
end
