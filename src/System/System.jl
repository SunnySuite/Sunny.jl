"""
    System(crystal::Crystal, latsize, infos, mode; units=Units.meV, seed::Int)

Construct a `System` of spins for a given `crystal` symmetry. The `latsize`
parameter determines the number of unit cells in each lattice vector direction.
The `infos` parameter is a list of [`SpinInfo`](@ref) objects, which determine
the magnitude ``S`` and ``g``-tensor of each spin.

The three possible options for `mode` are `:SUN`, `:dipole`, and `:large_S`. The
most variationally accurate choice is `:SUN`, in which each spin-``S`` degree of
freedom is described as an SU(_N_) coherent state, where ``N = 2S + 1``. Note
that an SU(_N_) coherent state fully describes any local spin state; this
description includes expected dipole components ``⟨Ŝᵅ⟩``, quadrupole components
``⟨ŜᵅŜᵝ+ŜᵝŜᵅ⟩``, etc.

Classical spin models have historically used the expected dipoles alone.  The
choice `mode=:dipole` projects the SU(_N_) dynamics onto the space of pure
dipoles. In practice this means that Sunny will simulate Landau-Lifshitz
dynamics, but all single-ion anisotropy or biquadratic exchange interactions
will be automatically renormalized to achieve maximum accuracy.

To disable such renormalization, e.g. to reproduce results collected using the
historical large-``S`` classical limit, use `mode=:large_S`. We emphasize,
however, that `mode=:SUN` or `mode=:dipole` should be preferred for the
development of new models.

The default units system of (meV, Å, tesla) can be overridden by with the
`units` parameter; see [`Units`](@ref). 

All spins are initially polarized in the ``z``-direction.
"""
function System(crystal::Crystal, latsize::NTuple{3,Int}, infos::Vector{SpinInfo}, mode::Symbol;
                    units=Units.meV, seed=nothing)
    if mode ∉ [:SUN, :dipole, :large_S]
        error("Mode must be one of [:SUN, :dipole, :large_S].")
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
    elseif sys.mode==:dipole
        "Dipole mode"
    elseif sys.mode==:large_S
        "Large-S classical limit"
    else
        error("Unreachable")
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
    new_field     = repeat(sys.interactions.extfield, factors..., 1)
    (; anisos, pairexch, ewald) = sys.interactions
    new_interactions = Interactions(new_field, anisos, pairexch, nothing)
    new_sys =  System(sys.mode, sys.crystal, new_latsize, sys.Ns, sys.gs, new_κs, new_interactions,
                      new_dipoles, new_coherents, Array{Vec3, 4}[], Array{CVec{N}, 4}[],
                      sys.units, copy(sys.rng))
    if !isnothing(ewald)
        new_sys.interactions.ewald = Ewald(new_sys)
    end
    return new_sys
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


wavevec_str(q) = "["*join(number_to_math_string.(q), ", ")*"]"

"""
    print_dominant_wavevectors(sys::System; nmax=10)

Prints a list of wavevectors in the first Brillouin zone according to their
contributions to the structure factor weights, accumulated over individual
sublattices. Coordinates are measured in units of reciprocal lattice vectors.
These dominant wavevectors may be used as input to
[`suggest_magnetic_supercell`](@ref). Unlike in [`static_intensities`](@ref),
here the structure factor weights do not incorporate phase averaging between the
sublattices.
"""
function print_dominant_wavevectors(sys::System{N}; nmax=10) where N
    s = reinterpret(reshape, Float64, sys.dipoles)
    V = prod(sys.latsize) # number of spins in sublattice
    cell_dims = (2,3,4)
    sk = FFTW.fft(s, cell_dims) / √V
    Sk = real.(conj.(sk) .* sk)
    dims = (1,5) # sum over spin index and basis index
    # In Julia 1.9 this becomes: sum(eachslice(dat; dims)))
    Sk = dropdims(sum(Sk; dims); dims)
    @assert sum(Sk) ≈ norm(sys.dipoles)^2

    weight = Sk .* (100 / norm(sys.dipoles)^2)
    p = sortperm(-weight[:])

    println("Dominant wavevectors for spin sublattices:\n")
    for (i, m) in enumerate(CartesianIndices(sys.latsize)[p])
        q = (Tuple(m) .- 1) ./ sys.latsize
        q = [qi > 1/2 ? qi-1 : qi for qi in q]

        qstr = wavevec_str(q)
        
        if weight[m] < 0.01
            break
        end

        wstr = @sprintf "%6.2f" weight[m]
        nspaces = 20
        spacing = " " ^ max(0, nspaces - length(qstr))
        print("    $qstr $spacing $wstr%")
        println(i == 1 ? " weight" : "")

        if i > nmax
            spacing = " " ^ (nspaces-1)
            println("    ... $spacing ...")
            break
        end
    end
end

"""
    suggest_magnetic_supercell(qs, latsize)

Suggests a magnetic supercell, in units of the crystal lattice vectors, that is
consistent with periodicity of the wavevectors in `qs`. The wavevectors should
be commensurate with `latsize` in units of the lattice vectors.
"""
function suggest_magnetic_supercell(qs, latsize)
    foreach(qs) do q
        m = round.(Int, q .* latsize)
        if norm(m - q .* latsize) > 1e-12
            error("Wavevector $q incommensurate with lattice size $latsize")
        end
    end

    # All possible periodic offsets
    nmax = latsize .÷ 2
    ns = [Vec3(n1, n2, n3) for n1 in -nmax[1]:nmax[1], n2 in -nmax[2]:nmax[2], n3 in -nmax[3]:nmax[3]][:]

    # Filter out elements of ns that are not consistent with q ∈ qs
    for q in qs
        ns = filter(ns) do n            
            # Proposed offset is `r = nᵅ aᵅ` in lattice vectors `aᵅ`. Wave
            # vector is `k = qᵝ bᵝ` in reciprocal vectors `bᵝ`. Using
            # orthogonality property `aᵅ⋅bᵝ = 2πδᵅᵝ`, it follows `r⋅k = 2π n⋅q`.
            # The integer offsets `n` are allowed if `n⋅q` is exactly integer,
            # which confirms periodicity, `exp(-i r⋅k) = 1`.
            is_approx_integer(n⋅q; atol=1e-12)
        end
    end

    # Naive supercell is the entire lattice
    best_score = Inf
    best_A = diagm(collect(latsize))

    # Goodness of supervectors A1, A2, A3. Lower is better.
    function score(A)
        (A1, A2, A3) = eachcol(A)
        # Minimize supercell volume
        V = (A1 × A2) ⋅ A3
        # println(A)
        # println(V)
        # To split ties, maximize orthogonality
        c1 = (normalize(A1) × normalize(A2)) ⋅ normalize(A3)
        # To split ties, maximize alignment with Cartesian system
        c2 = normalize(A1)[1] + normalize(A2)[2] + normalize(A3)[3]
        V <= 0 ? Inf : V - 1e-3*c1 - 1e-6c2
    end

    any_changed = true
    while any_changed
        any_changed = false

        # try all possible replacements of individual columns
        for j in 1:3, n in ns
            A = copy(best_A)
            A[:, j] = n
            if best_score > score(A)
                best_score = score(A)
                best_A = A
                any_changed = true
            end
        end

        # try permuting columns with various sign changes
        (A1, A2, A3) = eachcol(best_A)
        for A in ([A1 -A3 A2], [A1 A3 -A2],
                  [-A3 A2 A1], [A3 A2 -A1],
                  [-A2 A1 A3], [A2 -A1 A3])
            if best_score > score(A)
                best_score = score(A)
                best_A = A
                any_changed = true
            end
        end
    end

    # # Alternative brute force implementation, for reference
    # push!(ns, [latsize[1],0,0])
    # push!(ns, [0,latsize[2],0])
    # push!(ns, [0,0,latsize[3]])
    # for A1 in ns, A2 in ns, A3 in ns
    #     A = [A1 A2 A3]
    #     if best_score > score(A)
    #         best_score = score(A)
    #         best_A = A
    #     end
    # end


    qstrs = join(map(wavevec_str, qs), ", ")

    println("Suggested magnetic supercell in multiples of lattice vectors:")
    println()
    println("    A = $best_A")
    println()
    println("for wavevectors [$qstrs].")
end
