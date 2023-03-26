
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


# Transfer homogeneous interactions from `sys.origin` to reshaped `sys`
function set_interactions_from_origin!(sys::System{N}) where N
    origin = sys.origin
    new_ints  = interactions_homog(sys)
    orig_ints = interactions_homog(origin)

    for new_i in 1:natoms(sys.crystal)
        function map_couplings(couplings::Vector{Coupling{T}}) where T
            new_couplings = Coupling{T}[]
            for (; bond, J) in couplings
                new_bond = transform_bond(sys.crystal, new_i, origin.crystal, bond)
                isculled = bond_parity(new_bond)
                push!(new_couplings, Coupling(isculled, new_bond, J))
            end
            return sort!(new_couplings, by=c->c.isculled)
        end
        
        i = map_atom_to_crystal(sys.crystal, new_i, origin.crystal)
        new_ints[new_i].aniso    = orig_ints[i].aniso
        new_ints[new_i].heisen   = map_couplings(orig_ints[i].heisen)
        new_ints[new_i].exchange = map_couplings(orig_ints[i].exchange)
        new_ints[new_i].biquad   = map_couplings(orig_ints[i].biquad)
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
    # `origin`. Otherwise, we will need to reshape the crystal, map the
    # interactions, and keep a reference to the original system.
    if new_cell_size == I
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
        new_cryst = reshape_crystal(origin.crystal, Mat3(new_cell_size))

        new_na               = natoms(new_cryst)
        new_Ns               = zeros(Int, new_latsize..., new_na)
        new_κs               = zeros(Float64, new_latsize..., new_na)
        new_gs               = zeros(Mat3, new_latsize..., new_na)
        new_extfield         = zeros(Vec3, new_latsize..., new_na)
        new_dipoles          = zeros(Vec3, new_latsize..., new_na)
        new_coherents        = zeros(CVec{N}, new_latsize..., new_na)

        new_dipole_buffers   = Array{Vec3, 4}[]
        new_coherent_buffers = Array{CVec{N}, 4}[]

        new_ints = empty_interactions(new_na, N)

        new_sys = System(origin, origin.mode, new_cryst, new_latsize, new_Ns, new_κs, new_gs, new_ints, nothing,
                    new_extfield, new_dipoles, new_coherents, new_dipole_buffers, new_coherent_buffers, origin.units, copy(sys.rng))

        set_interactions_from_origin!(new_sys)
    end

    # Copy per-site quantities
    for new_site in all_sites(new_sys)
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

# Dimensions of a possibly reshaped unit cell, given in multiples of the
# original unit cell.
function cell_dimensions(sys)
    A = orig_crystal(sys).latvecs \ sys.crystal.latvecs
    @assert norm(A - round.(A)) < 1e-12
    return round.(Int, A)
end

"""
    resize_periodically(sys::System{N}, latsize) where N

Creates a [`System`](@ref) identical to `sys` but enlarged to a given number of
unit cells in each lattice vector direction.

An error will be thrown if `sys` is incommensurate with `latsize`. Use
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
    dims = (1,5) # sum over spin index and atom (sublattice) index
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
