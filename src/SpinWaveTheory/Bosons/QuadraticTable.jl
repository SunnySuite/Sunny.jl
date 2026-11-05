# The quadratic Hamiltonian compiled for evaluation at many wavevectors. Its
# source is the K = 2 term of the boson expansion, a list of monomials, which is
# the natural form for Wick's theorem but a slow one to Fourier transform. Here
# the monomials are merged by cell offset 𝐧 and matrix entry, and stored in
# compressed sparse row form,
#
#     H(𝐪)[r, cols[e]] = exp(2πi 𝐪⋅𝐧ₑ) vals[e],
#
# summed over the entries e of row r, where 𝐧ₑ = offsets[offset_index[e]]. Both
# triangles are stored. Each entry above the diagonal is the exact conjugate of
# its partner below, at the reversed offset, and the partners of each matrix
# element are summed in the same order. The diagonal entries pair up among
# themselves, and their sum is real. The dense matrix at 𝐪 is then Hermitian
# exactly, not merely to within round-off.
#
# The table is immutable. It is evaluated at one wavevector either densely, by
# `accum_table!`, or for matrix-free products, through a `HamiltonianOperator`.
struct QuadraticTable
    offsets      :: Vector{Vec3}
    rowptr       :: Vector{Int32}
    cols         :: Vector{Int32}
    offset_index :: Vector{Int32}
    vals         :: Vector{ComplexF64}
end

# Compiles monomials in the convention of `accum_quadratic!`, which this
# reproduces, Hermitian part included, for an `H` that is initially zero.
function QuadraticTable(terms, L)
    # The lower triangle of the Hermitian part. An entry v exp(2πi 𝐪⋅𝐧) at (r,
    # c) contributes v/2 there, and its adjoint conj(v)/2 at (c, r) and -𝐧.
    lower = Dict{Tuple{NTuple{3, Int}, Int32, Int32}, ComplexF64}()
    function add!(r, c, v, n)
        if r >= c
            k = (n, Int32(r), Int32(c))
            lower[k] = get(lower, k, 0.0im) + v/2
        end
        if c >= r
            k = (.-n, Int32(c), Int32(r))
            lower[k] = get(lower, k, 0.0im) + conj(v)/2
        end
    end
    for (; c, as, ns) in terms
        n = round.(Int, Tuple(ns[1] - ns[2]))
        add!(nambu_conj(as[2], L), as[1], c, n)
        add!(nambu_conj(as[1], L), as[2], c, .-n)
    end

    # Each entry of the lower triangle, and below the diagonal also its exact
    # conjugate in the upper, at the reversed offset. The diagonal already holds
    # both partners, at 𝐧 and -𝐧. Sorting the upper entries by -𝐧 sums each
    # element in the order of its partner.
    entries = Tuple{Int32, Int32, NTuple{3, Int}, ComplexF64}[]
    for ((n, r, c), v) in lower
        iszero(v) && continue
        push!(entries, (r, c, n, v))
        r == c || push!(entries, (c, r, .-n, conj(v)))
    end
    sort!(entries; by=((r, c, n, _),) -> (r, c, r >= c ? n : .-n))

    n = 2L
    rowptr = zeros(Int32, n+1)
    rowptr[1] = 1
    for (r, _, _, _) in entries
        rowptr[r+1] += 1
    end
    cumsum!(rowptr, rowptr)
    offsets = unique(Vec3(n) for (_, _, n, _) in entries)
    index = Dict(o => Int32(k) for (k, o) in enumerate(offsets))
    cols = [c for (_, c, _, _) in entries]
    offset_index = [index[Vec3(n)] for (_, _, n, _) in entries]
    vals = [v for (_, _, _, v) in entries]
    return QuadraticTable(offsets, rowptr, cols, offset_index, vals)
end

Base.size(t::QuadraticTable, _) = length(t.rowptr) - 1

# The phases exp(2πi 𝐪⋅𝐧) of the offsets of the table
offset_phases(t::QuadraticTable, q_supercell) = [cis(2π * dot(q_supercell, n)) for n in t.offsets]

# Adds the table at `q_supercell` to the dense matrix `H`. On the diagonal, the
# imaginary parts of partners cancel, and are dropped so that they cancel
# exactly.
function accum_table!(H, t::QuadraticTable, q_supercell)
    (; rowptr, cols, offset_index, vals) = t
    φ = offset_phases(t, q_supercell)
    @inbounds for r in 1:size(t, 1), e in rowptr[r]:rowptr[r+1]-1
        c = cols[e]
        v = φ[offset_index[e]] * vals[e]
        H[r, c] += r == c ? real(v) : v
    end
    return H
end

# A bound on the spectral radius of ĨH(𝐪) at every 𝐪, for the matrix H(𝐪) of
# the table: the largest Gershgorin row sum, whose phases have unit magnitude
function gershgorin_bound(t::QuadraticTable)
    (; rowptr, vals) = t
    rowsum(r) = sum(e -> abs(vals[e]), rowptr[r]:rowptr[r+1]-1; init=0.0)
    return maximum(rowsum, 1:size(t, 1); init=0.0)
end
