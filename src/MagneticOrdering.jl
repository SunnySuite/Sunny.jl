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
experiment. For that purpose, use [`SampledCorrelationsStatic`](@ref) instead.
"""
function print_wrapped_intensities(sys::System{N}; nmax=10) where N
    sys.crystal == orig_crystal(sys) || error("Cannot perform this analysis on reshaped system.")

    s = reinterpret(reshape, Float64, sys.dipoles)
    V = prod(sys.dims) # number of spins in sublattice
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
    for (i, m) in enumerate(CartesianIndices(sys.dims)[p])
        k = (Tuple(m) .- 1) ./ sys.dims
        k = [ki > 1/2 ? ki-1 : ki for ki in k]

        kstr = fractional_vec3_to_string(k)
        
        if weight[m] < 0.01
            break
        end

        wstr = @sprintf "%6.2f" weight[m]
        nspaces = 20
        spacing = " " ^ max(0, nspaces - length(kstr))
        print("    $kstr $spacing $wstr%")
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
            return numers .// denom
        end
    end
    error("Wavevectors are incommensurate for lattice sizes less than $maxsize. Try increasing `tol` parameter.")
end


"""
    suggest_magnetic_supercell(ks; tol=1e-12, maxsize=100)

Suggests a magnetic supercell, in units of the crystal lattice vectors, that is
consistent with periodicity of the wavevectors `ks` in RLU. If the wavevectors
are incommensurate (with respect to the maximum supercell size `maxsize`), one
can select a larger error tolerance `tol` to find a supercell that is almost
commensurate.

Prints a ``3×3`` matrix of integers that is suitable for use in
[`reshape_supercell`](@ref).

# Examples

```julia
# A magnetic supercell for a single-Q structure. Will print
k1 = [0, -1/4, 1/4]
suggest_magnetic_supercell([k1])       # [1 0 0; 0 2 1; 0 -2 1]

# A larger magnetic supercell for a double-Q structure
k2 = [1/4, 0, 1/4]
suggest_magnetic_supercell([k1, k2])   # [1 2 2; -1 2 -2; -1 2 2]

# If given incommensurate wavevectors, find an approximate supercell that
# is exactly commensurate for nearby wavevectors.
suggest_magnetic_supercell([[0, 0, 1/√5], [0, 0, 1/√7]]; tol=1e-2)

# This prints [1 0 0; 0 1 0; 0 0 16], which becomes commensurate under the
# approximations `1/√5 ≈ 7/16` and `1/√7 ≈ 3/8`.
```
"""
function suggest_magnetic_supercell(ks; tol=1e-12, maxsize=100)
    eltype(ks) <: AbstractVector{<: Number} || error("Pass a list of wavevectors")

    rational_ks = zeros(Rational{Int}, 3, length(ks))
    for i in 1:3
        xs = [k[i] for k in ks]
        rational_ks[i, :] = rationalize_simultaneously(xs; tol, maxsize)
    end

    suggest_magnetic_supercell_aux(eachcol(rational_ks))
end

function suggest_magnetic_supercell_aux(ks)
    # A supercell of shape A = diagm(nmax) is guaranteed to be commensurate
    nmax = [lcm(denominator.(row)...) for row in zip(ks...)]
    ns = Vec3.(eachcol(diagm(nmax)))

    # Candidate lattice vectors for a smaller supercell
    rg = div.(nmax, 2)
    append!(ns, Vec3.(Iterators.product(-rg[1]:rg[1], -rg[2]:rg[2], -rg[3]:rg[3])))

    # Sort by length
    sort!(ns, by=n->n'*n)

    # Remove zero vector
    @assert iszero(ns[1])
    deleteat!(ns, 1)

    # Filter out elements of ns that are not consistent with k ∈ ks
    for k in ks
        ns = filter(ns) do n
            # Wave vector `k` in RLU is commensurate if `n⋅k` is integer,
            # corresponding to the condition `exp(-i 2π n⋅k) = 1`.
            isinteger(n⋅k)
        end
    end

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

    # Initial guess for the supercell
    best_A = diagm(nmax)
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
    # all_As = ([A1 A2 A3] for (A1, A2, A3) in Iterators.product(ns, ns, ns))
    # @assert minimum(score, all_As) == best_score

    kstrs = join(map(fractional_vec3_to_string, ks), ", ")
    println("""Possible magnetic supercell in multiples of lattice vectors:
               
                   $(repr(best_A))
               
               for the rationalized wavevectors:
               
                   [$kstrs]""")
end


function check_commensurate(sys; k)
    q_reshaped = sys.crystal.recipvecs \ orig_crystal(sys).recipvecs * k
    commensurate = true
    for i = 1:3
        denom = denominator(rationalize(q_reshaped[i]))
        commensurate = commensurate && iszero(mod(sys.dims[i], denom))
    end
    if !commensurate
        @warn "Wavevector $(fractional_vec3_to_string(k)) is incommensurate with system."
    end
end


function deprecate_small_q(; q, k)
    if !isnothing(q)
        @warn "`q` argument is deprecated, use `k` to denote propagation wavevector instead!"
        k = q
    end
    isnothing(k) && error("Must provide wavevector `k` as named argument")
    return k
end

