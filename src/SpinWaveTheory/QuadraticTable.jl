# The quadratic Hamiltonian compiled for evaluation at many wavevectors. Its
# source is the K = 2 term of the boson expansion, a list of monomials, which is
# the natural form for Wick's theorem but a slow one to Fourier transform. Here
# the monomials are merged by cell offset 𝐧 and matrix entry, so that building
# the 2L×2L dynamical matrix at 𝐪 costs one `cis` per distinct offset and one
# multiply-add per nonzero entry,
#
#     H(𝐪)[r, c] = Σ_𝐧 exp(2πi 𝐪⋅𝐧) vals[r, c, 𝐧].
#
# Only the lower triangle is stored, already symmetrized, and the upper triangle
# is filled by conjugation. The result is then Hermitian exactly, not merely to
# within round-off, which a sum of the full matrix in a different order for each
# triangle would not guarantee.
struct QuadraticTable
    offsets :: Vector{Vec3}
    ranges  :: Vector{UnitRange{Int}}  # Entries carrying each offset
    rows    :: Vector{Int32}
    cols    :: Vector{Int32}
    vals    :: Vector{ComplexF64}
end

# Compiles monomials in the convention of `accum_quadratic!`, which this
# reproduces, Hermitian part included, for an `H` that is initially zero.
function QuadraticTable(terms, L)
    acc = Dict{Tuple{NTuple{3, Int}, Int32, Int32}, ComplexF64}()
    # Adds v exp(2πi 𝐪⋅𝐧) at (r, c) to the Hermitian part, of which only the
    # lower triangle is kept.
    function add!(r, c, v, n)
        if r >= c
            k = (n, Int32(r), Int32(c))
            acc[k] = get(acc, k, 0.0im) + v/2
        end
        if c >= r
            k = (.-n, Int32(c), Int32(r))
            acc[k] = get(acc, k, 0.0im) + conj(v)/2
        end
    end
    for (; c, as, ns) in terms
        n = round.(Int, Tuple(ns[1] - ns[2]))
        add!(nambu_conj(as[2], L), as[1], c, n)
        add!(nambu_conj(as[1], L), as[2], c, .-n)
    end

    entries = sort!([(k, v) for (k, v) in acc if !iszero(v)]; by=first)
    offsets = Vec3[]
    ranges = UnitRange{Int}[]
    for (e, ((n, _, _), _)) in enumerate(entries)
        if isempty(offsets) || offsets[end] != Vec3(n)
            push!(offsets, Vec3(n))
            push!(ranges, e:e)
        else
            ranges[end] = first(ranges[end]):e
        end
    end
    rows = [r for ((_, r, _), _) in entries]
    cols = [c for ((_, _, c), _) in entries]
    vals = [v for (_, v) in entries]
    return QuadraticTable(offsets, ranges, rows, cols, vals)
end

# Adds the lower triangle of the table at `q_reshaped` to `H`
function accum_lower!(H, t::QuadraticTable, q_reshaped)
    for (n, rg) in zip(t.offsets, t.ranges)
        φ = iszero(n) ? 1.0+0im : cis(2π * dot(q_reshaped, n))
        @inbounds for e in rg
            H[t.rows[e], t.cols[e]] += φ * t.vals[e]
        end
    end
end

# Fills the upper triangle of `H` from its lower, and makes its diagonal real
function hermitian_from_lower!(H)
    @inbounds for c in axes(H, 2)
        H[c, c] = real(H[c, c])
        for r in c+1:size(H, 1)
            H[c, r] = conj(H[r, c])
        end
    end
end

# Accumulates Y[q, :] += H(qs[q]) X[q, :], the matrix-free product, for each
# row q of X and Y.
function accum_mul!(Y, X, t::QuadraticTable, qs_reshaped, phases)
    for (n, rg) in zip(t.offsets, t.ranges)
        @. phases = cis(2π * dot(qs_reshaped, (n,)))
        @inbounds for e in rg
            (r, c, v) = (t.rows[e], t.cols[e], t.vals[e])
            for q in axes(Y, 1)
                Y[q, r] += phases[q] * v * X[q, c]
            end
            if r != c
                for q in axes(Y, 1)
                    Y[q, c] += conj(phases[q] * v) * X[q, r]
                end
            end
        end
    end
end
