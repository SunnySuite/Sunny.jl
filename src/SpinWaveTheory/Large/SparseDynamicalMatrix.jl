# The dynamical matrix H(𝐪) of a `SpinWaveTheory` in compressed sparse row form,
# for the matrix-free products of `SpinWaveTheoryLarge`. Both triangles of the
# `QuadraticTable` are stored, so that row r of the product is a gather,
#
#     (H(𝐪) x)[r] = ε x[r] + Σₑ exp(2πi 𝐪⋅𝐧ₑ) vals[e] x[cols[e]],
#
# over the entries e of row r, where ε is the regularization. Each entry refers to
# its cell offset 𝐧ₑ by index. Per wavevector, `values_at!` folds the phases into
# a copy of the values, so that a product costs one multiply-add per entry.
#
# Vectors are the rows of a matrix X, so that the values X[:, c] of a batch of
# vectors are contiguous. Each entry is then loaded once for the whole batch.
struct SparseDynamicalMatrix
    offsets        :: Vector{Vec3}
    rowptr         :: Vector{Int32}
    cols           :: Vector{Int32}
    offset_index   :: Vector{Int32}
    vals           :: Vector{ComplexF64}
    regularization :: Float64
end

function SparseDynamicalMatrix(swt::SpinWaveTheory)
    isnothing(uncontracted_system(swt.sys).ewald) || error("Dipole-dipole interactions are not yet supported in matrix-free products")
    t = swt.quadratic_table
    n = 2nbands(swt)
    # Adding zero normalizes -0.0, so that 𝐧 and -𝐧 of a zero offset coincide
    offsets = unique([t.offsets; .-t.offsets] .+ (zero(Vec3),))
    index = Dict(o => Int32(k) for (k, o) in enumerate(offsets))

    # The lower triangle, and its conjugate transpose, which is H[c, r] =
    # conj(H[r, c]) at the reversed offset
    entries = Tuple{Int32, Int32, Int32, ComplexF64}[]
    for (o, rg) in zip(t.offsets, t.ranges), e in rg
        (r, c, v) = (t.rows[e], t.cols[e], t.vals[e])
        push!(entries, (r, c, index[o + zero(Vec3)], v))
        r != c && push!(entries, (c, r, index[-o + zero(Vec3)], conj(v)))
    end
    sort!(entries; by=e -> (e[1], e[2]))

    rowptr = zeros(Int32, n+1)
    rowptr[1] = 1
    for (r, _, _, _) in entries
        rowptr[r+1] += 1
    end
    cumsum!(rowptr, rowptr)
    cols = [c for (_, c, _, _) in entries]
    offset_index = [k for (_, _, k, _) in entries]
    vals = [v for (_, _, _, v) in entries]
    return SparseDynamicalMatrix(offsets, rowptr, cols, offset_index, vals, swt.regularization)
end

Base.size(H::SparseDynamicalMatrix, _) = length(H.rowptr) - 1

# The values of H(𝐪) itself, each entry multiplied by its phase exp(2πi 𝐪⋅𝐧).
# Computed once per wavevector, and reused by every product.
function values_at!(Hq, H::SparseDynamicalMatrix, q_supercell)
    φ = [cis(2π * dot(q_supercell, o)) for o in H.offsets]
    @inbounds for e in eachindex(Hq)
        Hq[e] = φ[H.offset_index[e]] * H.vals[e]
    end
    return Hq
end

# Y = H(𝐪) X for each row of X, given the values Hq of `values_at!`. A column of
# the batch is a static vector, so that each output accumulates in registers.
mul_rows!(Y::Matrix, H::SparseDynamicalMatrix, X::Matrix, Hq) = mul_rows!(Y, H, X, Hq, Val{size(X, 1)}())

function mul_rows!(Y, H::SparseDynamicalMatrix, X, Hq, ::Val{B}) where B
    (; rowptr, cols, regularization) = H
    @assert size(Y) == size(X) == (B, size(H, 1))
    Xs = reinterpret(reshape, SVector{B, ComplexF64}, X)
    Ys = reinterpret(reshape, SVector{B, ComplexF64}, Y)
    @inbounds for r in eachindex(Xs)
        acc = regularization * Xs[r]
        for e in rowptr[r]:rowptr[r+1]-1
            acc += Hq[e] * Xs[cols[e]]
        end
        Ys[r] = acc
    end
    return Y
end

# Y = Ĩ X for each row of X, where Ĩ = diag(I, -I) is the Bogoliubov identity
function mul_Ĩ_rows!(Y, X)
    L = size(X, 2) ÷ 2
    @views Y[:, 1:L] .= X[:, 1:L]
    @views Y[:, L+1:2L] .= .-X[:, L+1:2L]
    return Y
end
