"""
    System(crystal::Crystal, latsize, siteinfos; mode, units=Units.meV, seed::Int)

Construct a `System` of spins for a given `crystal` symmetry. The `latsize`
parameter determines the number of unit cells in each lattice vector direction.
The `siteinfos` parameter is a list of [`SiteInfo`](@ref) objects, which
determine the magnitude ``S`` and ``g``-tensor of each spin.

The three possible options for `mode` are `:dipole`, `:SUN`, and `:projected`.
The choice `:dipole` restricts the description of a spin to its angular momentum
dipole. The choice `:SUN` will expand the description of a spin to a full
SU(_N_) coherent state, where ``N = 2S + 1``. This approach captures more
quantum mechanical degrees of freedom, e.g., the dynamics of quadrupolar
fluctuations ``⟨ŜᵅŜᵝ+ŜᵝŜᵅ⟩``. Finally, the choice `:projected` produces the
SU(_N_) dynamics but projected to the space of dipoles. In practice, the
distinction between `:projected` and `:dipoles` is that the former will
automatically apply appropriate renormalizations to the single-ion anisotropy
and biquadratic exchange interactions for maximum accuracy.

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, siteinfos::Vector{SiteInfo}=SiteInfo[];
                    mode::Symbol, units=Units.meV, seed=nothing)
    mode==:projected && error("SU(N) projected mode not yet implemented.")

    if mode ∉ [:dipole, :SUN, :projected]
        error("Mode must be one of [:dipole, :SUN, :projected].")
    end

    siteinfos = if isempty(siteinfos)
        [SiteInfo(i; S=1.0) for i in 1:nbasis(crystal)]
    else
        propagate_site_info(crystal, siteinfos)
    end
    Ss = [si.S for si in siteinfos]
    gs = [si.g for si in siteinfos]

    # Determine dimension N of the local Hilbert space, or 0 if in dipole-only mode
    N, κs = if mode == :SUN
        if !allequal(Ss)
            error("Currently all spins S must be equal in SU(N) mode.")
        end
        N = Int(2first(Ss) + 1)
        (N, fill(1.0, nbasis(crystal)))
    else
        (0, Ss)
    end

    interactions = Interactions(crystal, N)

    # Initialize sites to all spins along +z
    dipoles = fill(zero(Vec3), latsize..., nbasis(crystal))
    coherents = fill(zero(CVec{N}), latsize..., nbasis(crystal))
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = System(mode, crystal, latsize, interactions, dipoles, coherents, κs, gs,
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


"""
    extend_periodically(sys::System{N}, mults::NTuple{3, Int64}) where N

Creates a new System identical to `sys` but with each dimension multiplied
by the corresponding factor given in the tuple `mults`. The original spin configuration
is simply repeated periodically.
"""
function extend_periodically(sys::System{N}, factors::NTuple{3, Int64}) where N
    @assert all(>=(1), factors)
    latsize = factors .* sys.latsize
    dipoles   = repeat(sys.dipoles, factors..., 1)
    coherents = repeat(sys.coherents, factors..., 1)
    #KBTODO: repeat κs
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    return System(sys.mode, sys.crystal, latsize, sys.interactions, dipoles, coherents, sys.κs, sys.gs,
                      dipole_buffers, coherent_buffers, sys.units, copy(sys.rng))
end



const Cell = CartesianIndex{3}
const Idx = CartesianIndex{4}

# Element-wise application of mod1(cell+off, latsize), returning CartesianIndex
@inline offsetc(cell::Cell, off, latsize) = CartesianIndex(mod1.(Tuple(cell).+Tuple(off), latsize))

# Split a Cartesian index (cell,i) into its parts cell and i.
@inline splitidx(idx::Idx) = (CartesianIndex((idx[1],idx[2],idx[3])), idx[4])

@inline convert_idx(idx::Idx) = idx
@inline convert_idx(idx::NTuple{4,Int}) = CartesianIndex(idx)

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

@inline function getspin(sys::System{N}, idx::Idx) where N
    return SpinState(sys.dipoles[idx], sys.coherents[idx])
end

@inline function setspin!(sys::System{N}, spin::SpinState{N}, idx::Idx) where N
    sys.dipoles[idx] = spin.s
    sys.coherents[idx] = spin.Z
    nothing
end

@inline function flip(spin::SpinState{N}) where N
    return SpinState(-spin.s, flip_ket(spin.Z))
end

@inline function randspin(sys::System{0}, idx)
    s = sys.κs[idx[4]] * normalize(randn(sys.rng, Vec3))
    return SpinState(s, CVec{0}())
end
@inline function randspin(sys::System{N}, idx) where N
    Z = normalize(randn(sys.rng, CVec{N}))
    s = sys.κs[idx[4]] * expected_spin(Z)
    return SpinState(s, Z)
end

@inline function dipolarspin(sys::System{0}, idx, dir)
    s = sys.κs[idx[4]] * normalize(Vec3(dir))
    Z = CVec{0}()
    return SpinState(s, Z)
end
@inline function dipolarspin(sys::System{N}, idx, dir) where N
    Z = ket_from_dipole(Vec3(dir), Val(N))
    s = sys.κs[idx[4]] * expected_spin(Z)
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


function get_dipole_buffers(sys::System, numrequested)
    numexisting = length(sys.dipole_buffers)
    if numexisting < numrequested
        for _ ∈ 1:(numrequested-numexisting)
            push!(sys.dipole_buffers, zero(sys.dipoles))
        end
    end
    return sys.dipole_buffers[1:numrequested]
end

function get_coherent_buffers(sys::System, numrequested)
    numexisting = length(sys.coherent_buffers)
    if numexisting < numrequested
        for _ ∈ 1:(numrequested-numexisting)
            push!(sys.coherent_buffers, zero(sys.coherents))
        end
    end
    return sys.coherent_buffers[1:numrequested]
end
