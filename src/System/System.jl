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

The choice `:dipole` projects the SU(_N_) dynamics onto the space of pure
dipoles. In practice this means that Sunny will simulate Landau-Lifshitz
dynamics, but all single-ion anisotropy or biquadratic exchange interactions
will be automatically renormalized for maximum accuracy. [IN PROGRESS]

To disable such renormalization, e.g. to reproduce results collected using the
historical large-``S`` classical limit, use `mode=:large_S`. [IN PROGRESS] Modes
`:SUN` or `:dipole` should be preferred for the development of new models.

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
    extfield = zeros(Vec3, latsize..., nb)
    interactions = empty_interactions(nb, N)
    ewald = nothing
    dipoles = fill(zero(Vec3), latsize..., nb)
    coherents = fill(zero(CVec{N}), latsize..., nb)
    dipole_buffers = Array{Vec3, 4}[]
    coherent_buffers = Array{CVec{N}, 4}[]
    rng = isnothing(seed) ? Random.Xoshiro() : Random.Xoshiro(seed)

    ret = System(nothing, mode, crystal, latsize, Ns, gs, κs, extfield, interactions, ewald,
                 dipoles, coherents, dipole_buffers, coherent_buffers, units, rng)
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
# It should be thread-safe to use the original and the copied systems, without
# any restrictions.
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

    System(origin_clone, mode, crystal, latsize, Ns, copy(gs), copy(κs), copy(extfield),
           interactions_clone, ewald_clone, copy(dipoles), copy(coherents),
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

@inline convert_idx(site::CartesianIndex{4})           = site
@inline convert_idx(site::NTuple{4, Int})              = CartesianIndex(site)
@inline convert_idx(c::CartesianIndex{3}, i::Int)      = CartesianIndex(c[1], c[2], c[3], i)
@inline convert_idx(c1::Int, c2::Int, c3::Int, i::Int) = CartesianIndex(c1, c2, c3, i)

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
    return sys.crystal.lat_vecs * r
end

"""
    magnetic_moment(sys::System, site::Site)

Get the magnetic moment for a [`Site`](@ref). The result is `sys.dipoles[site]`
multiplied by the Bohr magneton and the ``g``-tensor for `site`.
"""
magnetic_moment(sys::System, site) = sys.units.μB * sys.gs[site[4]] * sys.dipoles[site]

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)

# The original crystal for a system, invariant under reshaping
orig_crystal(sys) = isnothing(sys.origin) ? sys.crystal : sys.origin.crystal

# Position of a site in fractional coordinates
function position(sys::System, site)
    return orig_crystal(sys).lat_vecs \ global_position(sys, site)
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
    new_r = sys.crystal.lat_vecs \ orig_crystal(sys).lat_vecs * r
    b, offset = position_to_index_and_offset(sys.crystal, new_r)
    cell = @. mod1(offset+1, sys.latsize) # 1-based indexing with periodicity
    return convert_idx(cell..., b)
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
    polarize_spin!(sys::System, dir, site::Site)

Polarize the spin at a [`Site`](@ref) along the direction `dir`.
"""
function polarize_spin!(sys::System{N}, dir, site) where N
    site = convert_idx(site)
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


"""
    reshape_geometry(sys::System, A)

Maps an existing [`System`](@ref) to a new one that has the shape and
periodicity of a requested supercell. The columns of the ``3×3`` integer matrix
`A` represent the supercell lattice vectors measured in units of the original
crystal lattice vectors.

The crystal unit cell may also need to be reshaped to achieve the desired
periodicity of the requested supercell. If this is the case, the returned
`System` object will be missing symmetry information. Consequently, certain
operations will be unavailable for this system, e.g., setting interactions by
symmetry propagation. In practice, one can set all interactions using the
original system, and then reshape as a final step.
"""
function reshape_geometry(sys::System{N}, A) where N
    # latsize for new system
    new_latsize = NTuple{3, Int}(gcd.(eachcol(A)))
    # Unit cell for new system, in units of original unit cell. Obtained by
    # dividing each column of A by corresponding new_latsize component.
    new_cell_size = Int.(A / diagm(collect(new_latsize)))
    return reshape_geometry_aux(sys, new_latsize, new_cell_size)
end


# Map all "unit cell" quantities from `origin` to `new_sys`. These are: Ns, gs,
# and homogeneous interactions. 
function transfer_unit_cell!(new_sys::System{N}, origin::System{N}) where N
    # Both systems must be homogeneous
    @assert is_homogeneous(new_sys) && is_homogeneous(origin)

    origin_ints = interactions_homog(origin)
    new_ints    = interactions_homog(new_sys)
    new_cryst   = new_sys.crystal

    for new_i in 1:nbasis(new_cryst)
        new_ri = new_cryst.positions[new_i]
        ri = origin.crystal.lat_vecs \ new_cryst.lat_vecs * new_ri
        i = position_to_index(origin.crystal, ri)

        # Spin descriptors
        new_sys.Ns[new_i] = origin.Ns[i]
        new_sys.gs[new_i] = origin.gs[i]

        # Interactions
        function map_couplings(couplings::Vector{Coupling{T}}) where T
            new_couplings = Coupling{T}[]
            for (; bond, J) in couplings
                new_bond = transform_bond(new_cryst, new_i, origin.crystal, bond)
                isculled = bond_parity(new_bond)
                push!(new_couplings, Coupling(isculled, new_bond, J))
            end
            return sort!(new_couplings, by=c->c.isculled)
        end
        new_ints[new_i].aniso    = origin_ints[i].aniso
        new_ints[new_i].heisen   = map_couplings(origin_ints[i].heisen)
        new_ints[new_i].exchange = map_couplings(origin_ints[i].exchange)
        new_ints[new_i].biquad   = map_couplings(origin_ints[i].biquad)
    end
end

function reshape_geometry_aux(sys::System{N}, new_latsize::NTuple{3, Int}, new_cell_size::Matrix{Int}) where N
    is_homogeneous(sys) || error("Cannot reshape system with inhomogeneous interactions.")

    # `origin` describes the unit cell of the original system. For sequential
    # reshapings, `sys.origin` keeps its original meaning. Make a deep copy so
    # that the new system fully owns `origin`, and mutable updates to the
    # previous system won't affect this one.
    origin = clone_system(isnothing(sys.origin) ? sys : sys.origin)

    # If `new_cell_size == I`, we can effectively restore the unit cell of
    # `origin`, but with `new_latsize`
    if new_cell_size == I
        new_cryst = origin.crystal
        new_nb = nbasis(new_cryst)

        new_κs               = zeros(Float64, new_latsize..., new_nb)
        new_extfield         = zeros(Vec3, new_latsize..., new_nb)
        new_ints             = interactions_homog(origin) # homogeneous only
        new_dipoles          = zeros(Vec3, new_latsize..., new_nb)
        new_coherents        = zeros(CVec{N}, new_latsize..., new_nb)
        new_dipole_buffers   = Array{Vec3, 4}[]
        new_coherent_buffers = Array{CVec{N}, 4}[]
        new_sys = System(nothing, origin.mode, origin.crystal, new_latsize, origin.Ns, origin.gs, new_κs, new_extfield, new_ints, nothing,
                         new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))
    
    # Else, rebuild the unit cell for the new crystal
    else
        new_cryst = reshape_crystal(origin.crystal, Mat3(new_cell_size))
        new_nb = nbasis(new_cryst)
        
        new_Ns               = zeros(Int, new_nb)
        new_gs               = zeros(Mat3, new_nb)
        new_κs               = zeros(Float64, new_latsize..., new_nb)
        new_extfield         = zeros(Vec3, new_latsize..., new_nb)
        new_ints             = empty_interactions(new_nb, N)
        new_dipoles          = zeros(Vec3, new_latsize..., new_nb)
        new_coherents        = zeros(CVec{N}, new_latsize..., new_nb)
        new_dipole_buffers   = Array{Vec3, 4}[]
        new_coherent_buffers = Array{CVec{N}, 4}[]
        
        new_sys = System(origin, origin.mode, new_cryst, new_latsize, new_Ns, new_gs, new_κs, new_extfield, new_ints, nothing,
                        new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))
        transfer_unit_cell!(new_sys, origin)
    end

    # Copy per-site quantities from `sys`
    for new_site in all_sites(new_sys)
        site = position_to_site(sys, position(new_sys, new_site))
        new_sys.κs[new_site] = sys.κs[site]
        new_sys.extfield[new_site] = sys.extfield[site]
        new_sys.dipoles[new_site] = sys.dipoles[site]
        new_sys.coherents[new_site] = sys.coherents[site]
    end

    # Enable dipole-dipole
    if !isnothing(sys.ewald)
        enable_dipole_dipole!(new_sys)
    end

    return new_sys
end

# Dimensions of a possibly reshaped unit cell, given in multiples of the
# original unit cell.
function cell_dimensions(sys)
    A = orig_crystal(sys).lat_vecs \ sys.crystal.lat_vecs
    @assert norm(A - round.(A)) < 1e-12
    return round.(Int, A)
end

"""
    resize_periodically(sys::System{N}, latsize) where N

Creates a [`System`](@ref) identical to `sys` but enlarged to a given number of
unit cells in each lattice vector direction.

An error will be thrown if `sys` is incommensurate with `latvecs`. Use
[`reshape_geometry`](@ref) instead to reduce the volume, or to perform an
incommensurate reshaping.
"""
function resize_periodically(sys::System{N}, latsize::NTuple{3,Int}) where N
    # Shape of the original system, in multiples of the original unit cell.
    sysdims = cell_dimensions(sys) * diagm(collect(sys.latsize))
    # Proposed system shape, given in fractional coordinates of original system
    # geometry
    A = sysdims \ diagm(collect(latsize))
    # All matrix elements must be integer
    if norm(A - round.(A)) > 1e-12
        error("Incommensurate system size.")
    end
    return reshape_geometry(sys, diagm(collect(latsize)))
end

"""
    repeat_periodically(sys::System{N}, counts) where N

Creates a [`System`](@ref) identical to `sys` but repeated a given number of
times in each dimension, specified by the tuple `counts`.
"""
function repeat_periodically(sys::System{N}, counts::NTuple{3,Int}) where N
    counts = NTuple{3,Int}(counts)
    @assert all(>=(1), counts)
    # Scale each column by `counts` and reshape
    return reshape_geometry_aux(sys, counts .* sys.latsize, Matrix(cell_dimensions(sys)))
end


wavevec_str(q) = "[" *join(number_to_math_string.(q), ", ")*"]"

"""
    print_wrapped_intensities(sys::System; nmax=10)

For Bravais lattices: Prints up to `nmax` wavevectors according to their
instantaneous (static) structure factor intensities, listed in descending order.
For non-Bravais lattices: Performs the same analysis for each spin sublattice
independently; the output weights are naïvely averaged over sublattices, without
incorporating phase shift information. Only wavevectors within the first
Brillouin zone are printed. Wavevector coordinates are given in reciprocal
lattice units, such that each coordinate is between ``-1/2`` and ``1/2``.  The
output from this function will typically be used as input to
[`suggest_magnetic_supercell`](@ref).

Because this function does not incorporate phase information in its averaging
over sublattices, the printed weights are not directly comparable with
experiment. For that purpose, use [`InstantStructureFactor`](@ref) instead.

The weights printed by `print_wrapped_intensities` may be given a physical
interpretation as follows: All possible ``q``-vectors are periodically wrapped
into the first Brillouin zone, and the average over their corresponding
instantaneous structure factor intensities produce the output weights.
"""
function print_wrapped_intensities(sys::System{N}; nmax=10) where N
    isnothing(sys.origin) || error("Cannot perform this analysis on reshaped system.")

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
consistent with periodicity of the wavevectors in `qs`. An upper bound for the
supercell is given by `latsize`, which is measured in units of lattice vectors,
and must be commensurate with the wavevectors.
"""
function suggest_magnetic_supercell(qs, latsize)
    foreach(qs) do q
        m = round.(Int, q .* latsize)
        if norm(m - q .* latsize) > 1e-12
            error("Wavevector $q incommensurate with lattice size $latsize")
        end
    end

    # All possible periodic offsets, sorted by length
    nmax = div.(latsize, 2)
    ns = [[n1, n2, n3] for n1 in -nmax[1]:nmax[1], n2 in -nmax[2]:nmax[2], n3 in -nmax[3]:nmax[3]][:]
    sort!(ns, by=n->n'*n)

    # Remove zero vector
    @assert iszero(ns[1])
    deleteat!(ns, 1)

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

    # Add vectors that wrap the entire lattice to ensure that a subset of the ns
    # span a nonzero volume.
    push!(ns, [latsize[1], 0, 0])
    push!(ns, [0, latsize[2], 0])
    push!(ns, [0, 0, latsize[3]])

    # Goodness of supervectors A1, A2, A3. Lower is better.
    function score(A)
        (A1, A2, A3) = eachcol(A)
        # Minimize supercell volume
        V = (A1 × A2) ⋅ A3
        # To split ties, maximize orthogonality
        c1 = (normalize(A1) × normalize(A2)) ⋅ normalize(A3)
        # To split ties, maximize alignment with Cartesian system
        c2 = normalize(A1)[1] + normalize(A2)[2] + normalize(A3)[3]
        V <= 0 ? Inf : V - 1e-3*c1 - 1e-6c2
    end

    # Find three vectors that span a nonzero volume which is hopefully small.
    # This will be our initial guess for the supercell.
    i1 = 1
    A1 = ns[i1]
    i2 = findfirst(n -> !iszero(n×A1), ns[i1+1:end])::Int
    A2 = ns[i1+i2]
    i3 = findfirst(n -> !iszero(n⋅(A1×A2)), ns[i1+i2+1:end])::Int
    A3 = ns[i1+i2+i3]
    best_A = [A1 A2 A3]
    best_score = score(best_A)

    # Iteratively search for an improved supercell. For efficiency, restrict the
    # search to some maximum number of displacement vectors
    ns = first(ns, 1000)

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

        # try all possible rotations
        (A1, A2, A3) = eachcol(best_A)
        for A in ([A1 -A3 A2], [A1 -A2 -A3], [A1 A3 -A2], # x-axis
                  [A3 A2 -A1], [-A1 A2 -A3], [-A3 A2 A1], # y-axis
                  [-A2 A1 A3], [-A1 -A2 A3], [A2 -A1 A3]) # z-axis
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
    println("    $best_A")
    println()
    println("for wavevectors [$qstrs].")
end
