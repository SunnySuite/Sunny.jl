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

Classical spin models have historically used the expected dipoles alone.  The
choice `mode=:dipole` projects the SU(_N_) dynamics onto the space of pure
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


function clone_spin_state(sys::System{N}) where N
    System(sys.origin, sys.mode, sys.crystal, sys.latsize, sys.Ns, sys.gs, sys.κs, sys.extfield,
        sys.interactions, sys.ewald, copy(sys.dipoles), copy(sys.coherents), sys.dipole_buffers, sys.coherent_buffers,
        sys.units, copy(sys.rng))
end


"""
    Site(n1, n2, n3, i)

Construct an index to a single site in a [`System`](@ref) via its unit cell
`(n1,n2,n3)` and its atom `i`. Can be used to index `dipoles` and `coherents`
fields of a `System`, or to set inhomogeneous interactions. This function is
effectively an alias for the four-component `CartesianIndex` constructor.
"""
@inline Site(idx::CartesianIndex{4})            = idx
@inline Site(idx::NTuple{4, Int})               = CartesianIndex(idx)
@inline Site(n1::Int, n2::Int, n3::Int, i::Int) = CartesianIndex(n1, n2, n3, i)

# Element-wise application of mod1(cell+off, latsize), returning CartesianIndex
@inline offsetc(cell::CartesianIndex{3}, off, latsize) = CartesianIndex(mod1.(Tuple(cell).+Tuple(off), latsize))

# Split a Cartesian index (cell,i) into its parts cell and i.
@inline splitidx(idx::CartesianIndex{4}) = (CartesianIndex((idx[1],idx[2],idx[3])), idx[4])

# An iterator over all sites using CartesianIndices
@inline all_sites(sys::System) = CartesianIndices(sys.dipoles)

# An iterator over all unit cells using CartesianIndices
@inline all_cells(sys::System) = CartesianIndices(sys.latsize)

# Position of a site in global coordinates
function global_position(sys::System, idx)
    r = sys.crystal.positions[idx[4]] + Vec3(idx[1]-1, idx[2]-1, idx[3]-1)
    return sys.crystal.lat_vecs * r
end

# Positions of all sites in global coordinates
global_positions(sys::System) = [global_position(sys, idx) for idx in all_sites(sys)]

# Magnetic moment at a site
magnetic_moment(sys::System, idx) = sys.units.μB * sys.gs[idx[4]] * sys.dipoles[idx]

# Total volume of system
volume(sys::System) = cell_volume(sys.crystal) * prod(sys.latsize)

# The original crystal for a system, invariant under reshaping
orig_crystal(sys) = isnothing(sys.origin) ? sys.crystal : sys.origin.crystal

# Position of a site in fractional coordinates
function position(sys::System, idx)
    return orig_crystal(sys).lat_vecs \ global_position(sys, idx)
end

# Convert a fractional position `r` to a Cartesian{4} site index
function position_to_site(sys::System, r::Vec3)
    # convert to fractional coordinates of possibly reshaped crystal
    new_r = sys.crystal.lat_vecs \ orig_crystal(sys).lat_vecs * r
    b, offset = position_to_index_and_offset(sys.crystal, new_r)
    cell = @. mod1(offset+1, sys.latsize) # 1-based indexing with periodicity
    return Site(cell..., b)
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

    origin_ints = interactions(origin)
    new_ints    = interactions(new_sys)
    new_cryst   = new_sys.crystal

    # Matrix that maps from fractional positions in `new_cryst` to fractional
    # positions in `cryst`
    to_orig_pos = origin.crystal.lat_vecs \ new_cryst.lat_vecs
    # Inverse mapping
    to_new_pos = inv(to_orig_pos)

    for new_i in 1:nbasis(new_cryst)
        new_ri = new_cryst.positions[new_i]
        i = position_to_index(origin.crystal, to_orig_pos * new_ri)

        # Spin descriptors
        new_sys.Ns[new_i] = origin.Ns[i]
        new_sys.gs[new_i] = origin.gs[i]

        # Interactions
        function map_couplings(couplings::Vector{Coupling{T}}) where T
            new_couplings = Coupling{T}[]
            for (; bond, J) in couplings
                displacement = origin.crystal.positions[bond.j] + bond.n - origin.crystal.positions[bond.i]
                new_rj = new_ri + to_new_pos * displacement
                new_j, new_n = position_to_index_and_offset(new_cryst, new_rj)
                new_bond = Bond(new_i, new_j, new_n)
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
    origin = deepcopy(isnothing(sys.origin) ? sys : sys.origin)

    # If `new_cell_size == I`, we can effectively restore the unit cell of
    # `origin`, but with `new_latsize`
    if new_cell_size == I
        new_cryst = origin.crystal
        new_nb = nbasis(new_cryst)

        new_κs               = zeros(Float64, new_latsize..., new_nb)
        new_extfield         = zeros(Vec3, new_latsize..., new_nb)
        new_ints             = interactions(origin) # homogeneous only
        new_dipoles          = zeros(Vec3, new_latsize..., new_nb)
        new_coherents        = zeros(CVec{N}, new_latsize..., new_nb)
        new_dipole_buffers   = Array{Vec3, 4}[]
        new_coherent_buffers = Array{CVec{N}, 4}[]
        new_sys = System(nothing, origin.mode, origin.crystal, new_latsize, origin.Ns, origin.gs, new_κs, new_extfield, new_ints, nothing,
                         new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))
    
    # Else we must rebuild the unit cell for the new crystal
    else
        new_cryst = resize_crystal(origin.crystal, Mat3(new_cell_size))
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
    for new_idx in all_sites(new_sys)
        idx = position_to_site(sys, position(new_sys, new_idx))
        new_sys.κs[new_idx] = sys.κs[idx]
        new_sys.extfield[new_idx] = sys.extfield[idx]
        new_sys.dipoles[new_idx] = sys.dipoles[idx]
        new_sys.coherents[new_idx] = sys.coherents[idx]
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
    if isnothing(sys.origin)
        return [1 0 0; 0 1 0; 0 0 1]
    else
        A = sys.origin.crystal.lat_vecs \ sys.crystal.lat_vecs
        @assert norm(A - round.(A)) < 1e-12
        return round.(Int, A)
    end
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
    return reshape_geometry(sys, cell_dimensions(sys) * diagm(collect(counts)))
end


wavevec_str(q) = "["*join(number_to_math_string.(q), ", ")*"]"

"""
    print_dominant_wavevectors(sys::System; nmax=10)

Prints a list of wavevectors according to their weights in the static structure
factor. Coordinates are given in units of reciprocal lattice vectors. These
dominant wavevectors may be used as input to
[`suggest_magnetic_supercell`](@ref).

Unlike in [`instant_intensities`](@ref), here the structure factor weights do not
incorporate phase averaging between sublattices. Instead, intensities are
calculated for each sublattice independently, and naïvely summed. This means
that wavevectors beyond the first Brillouin zone will be missing.
"""
function print_dominant_wavevectors(sys::System{N}; nmax=10) where N
    if !isnothing(sys.origin)
        error("Cannot perform this analysis on reshaped system.")
    end

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
