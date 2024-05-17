"""
    print_wrapped_intensities(sys::System; nmax=10)

For Bravais lattices: Prints up to `nmax` wavevectors according to their
instantaneous (static) structure factor intensities, listed in descending order.
For non-Bravais lattices: Performs the same analysis for each spin sublattice
independently; the output weights are naïvely averaged over sublattices, without
incorporating phase shift information. This procedure therefore wraps all
wavevectors into the first Brillouin zone. Each wavevector coordinate is given
between ``-1/2`` and ``1/2`` in reciprocal lattice units (RLU).  The output from
this function will typically be used as input to
[`suggest_magnetic_supercell`](@ref).

Because this function does not incorporate phase information in its averaging
over sublattices, the printed weights are not directly comparable with
experiment. For that purpose, use [`instant_correlations`](@ref) instead.
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

        qstr = fractional_vec3_to_string(q)
        
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

function rationalize_simultaneously(xs; tol, maxsize)
    for denom = 1:maxsize
        numers = @. round(Int, xs * denom)
        errs = @. xs - numers / denom
        if all(e -> abs(e) < tol, errs)
            return numers, denom
        end
    end
    error("Wavevectors are incommensurate for lattice sizes less than $maxsize. Try increasing `tol` parameter.")
end


"""
    suggest_magnetic_supercell(qs; tol=1e-12, maxsize=100)

Suggests a magnetic supercell, in units of the crystal lattice vectors, that is
consistent with periodicity of the wavevectors `qs` in RLU. If the wavevectors
are incommensurate (with respect to the maximum supercell size `maxsize`), one
can select a larger error tolerance `tol` to find a supercell that is almost
commensurate.

Prints a ``3×3`` matrix of integers that is suitable for use in
[`reshape_supercell`](@ref).

# Examples

```julia
# A magnetic supercell for a single-Q structure. Will print
q1 = [0, -1/4, 1/4]
suggest_magnetic_supercell([q1])       # [1 0 0; 0 2 1; 0 -2 1]

# A larger magnetic supercell for a double-Q structure
q2 = [1/4, 0, 1/4]
suggest_magnetic_supercell([q1, q2])   # [1 2 2; -1 2 -2; -1 2 2]

# If given incommensurate wavevectors, find an approximate supercell that
# is exactly commensurate for nearby wavevectors.
suggest_magnetic_supercell([[0, 0, 1/√5], [0, 0, 1/√7]]; tol=1e-2)

# This prints [1 0 0; 0 1 0; 0 0 16], which becomes commensurate under the
# approximations `1/√5 ≈ 7/16` and `1/√7 ≈ 3/8`.
```
"""
function suggest_magnetic_supercell(qs; tol=1e-12, maxsize=100)
    new_qs = zeros(Float64, 3, length(qs))
    denoms = zeros(Int, 3)

    for i in 1:3
        xs = [q[i] for q in qs]
        numers, denom = rationalize_simultaneously(xs; tol, maxsize)
        new_qs[i, :] = numers ./ denom
        denoms[i] = denom
    end

    new_qs = eachcol(new_qs)

    if !isapprox(qs, new_qs; atol=1e-12)
        @info "Using adapted q values: " * join(fractional_vec3_to_string.(new_qs), ", ")
    end

    suggest_magnetic_supercell_aux(new_qs, denoms)
end

function suggest_magnetic_supercell_aux(qs, denoms)
    foreach(qs) do q
        numers = @. round(Int, q * denoms)
        @assert norm(numers - q .* denoms) < 1e-12
    end

    # All possible periodic offsets, sorted by length
    nmax = div.(denoms, 2)
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
    push!(ns, [denoms[1], 0, 0])
    push!(ns, [0, denoms[2], 0])
    push!(ns, [0, 0, denoms[3]])

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
    # push!(ns, [denoms[1],0,0])
    # push!(ns, [0,denoms[2],0])
    # push!(ns, [0,0,denoms[3]])
    # for A1 in ns, A2 in ns, A3 in ns
    #     A = [A1 A2 A3]
    #     if best_score > score(A)
    #         best_score = score(A)
    #         best_A = A
    #     end
    # end

    qstrs = join(map(fractional_vec3_to_string, qs), ", ")
    println("""Suggested magnetic supercell in multiples of lattice vectors:
               
                   $(repr(best_A))
               
               for wavevectors [$qstrs].""")
end


function check_commensurate(sys, q)
    q_reshaped = sys.crystal.recipvecs \ orig_crystal(sys).recipvecs * q
    commensurate = true
    for i = 1:3
        denom = denominator(rationalize(q_reshaped[i]))
        commensurate = commensurate && iszero(mod(sys.latsize[i], denom))
    end
    if !commensurate
        @warn "Wavevector $(fractional_vec3_to_string(q)) is incommensurate with system."
    end
end

"""
    set_spiral_order!(sys; q, axis, S0)

Initializes the system with a spiral order described by the wavevector `q`, an
axis of rotation `axis`, and an initial dipole direction `S0` at the real-space
origin. The wavevector is expected in repicrocal lattice units (RLU), while the
direction vectors `axis` and `S0` are expected in global Cartesian coordinates.

# Example

```julia
# Spiral order for a wavevector propagating in the direction of the first
# reciprocal lattice vector (i.e., orthogonal to the lattice vectors ``𝐚_2``
# and ``𝐚_3``), repeating with a period of 10 lattice constants, and spiraling
# about the ``ẑ``-axis. The spin at the origin will point in the direction
# ``𝐒_0 = ŷ + ẑ``.  Here, ``(x̂, ŷ, ẑ)`` are the axes of Cartesian coordinate
# system in the global frame.
set_spiral_order!(sys; q=[1/10, 0, 0], axis=[0, 0, 1], S0=[0, 1, 1])
```

See also [`set_spiral_order_on_sublattice!`](@ref).
"""
function set_spiral_order!(sys; q, axis, S0)
    check_commensurate(sys, q)
    q_absolute = orig_crystal(sys).recipvecs * q

    for site in eachsite(sys)
        r = global_position(sys, site)
        θ = q_absolute ⋅ r
        set_dipole!(sys, axis_angle_to_matrix(axis, θ) * S0, site)
    end
end


"""
    set_spiral_order_on_sublattice!(sys, i; q, axis, S0)

Initializes sublattice `i` with a spiral order described by the wavevector `q`,
an axis of rotation `axis`, and an initial dipole direction `S0`. The phase is
selected such that the spin at `sys.dipole[1,1,1,i]` will point in the direction
of `S0`. The wavevector is expected in repicrocal lattice units (RLU), while the
direction vectors `axis` and `S0` are expected in global Cartesian coordinates.

This function is not available for systems with reshaped unit cells.

See also [`set_spiral_order!`](@ref).
"""
function set_spiral_order_on_sublattice!(sys, i; q, axis, S0)
    if orig_crystal(sys) != sys.crystal
        error("Cannot operate on a reshaped crystal. Atom indices may have changed.")
    end

    check_commensurate(sys, q)
    q_absolute = orig_crystal(sys).recipvecs * q

    r0 = global_position(sys, (1, 1, 1, i))
    for cell in eachcell(sys)
        site = (Tuple(cell)..., i)
        r = global_position(sys, site)
        θ = q_absolute ⋅ (r - r0)
        set_dipole!(sys, axis_angle_to_matrix(axis, θ) * S0, site)
    end
end
