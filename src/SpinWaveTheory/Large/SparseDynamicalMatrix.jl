# The dynamical matrix H(𝐪) of a `SpinWaveTheory` in compressed sparse row form,
# for the matrix-free products of `SpinWaveTheoryLarge`. Both triangles of the
# `QuadraticTable` are stored, so that row r of the product is a gather,
#
#     (H(𝐪) x)[r] = ε x[r] + Σₑ exp(2πi 𝐪⋅𝐧ₑ) vals[e] x[cols[e]],
#
# over the entries e of row r, where ε is the regularization. Each entry refers to
# its cell offset 𝐧ₑ by index. Per wavevector, `values_at!` folds the phases into
# a copy of the values, so that a product costs one multiply-add per entry.
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

# A bound on the spectral radius of ĨH(𝐪) at every 𝐪, by the Gershgorin theorem.
# The phases of `values_at!` have unit magnitude.
function spectral_bound(H::SparseDynamicalMatrix)
    (; rowptr, vals, regularization) = H
    rowsum(r) = sum(e -> abs(vals[e]), rowptr[r]:rowptr[r+1]-1; init=0.0)
    return regularization + maximum(rowsum, 1:size(H, 1))
end

# The values of H(𝐪) itself, each entry multiplied by its phase exp(2πi 𝐪⋅𝐧).
# Computed once per wavevector, and reused by every product.
function values_at!(Hq, H::SparseDynamicalMatrix, q_supercell)
    φ = [cis(2π * dot(q_supercell, o)) for o in H.offsets]
    @inbounds for e in eachindex(Hq)
        Hq[e] = φ[H.offset_index[e]] * H.vals[e]
    end
    return Hq
end

# y = H(𝐪) x, given the values Hq of `values_at!`
function mul_sparse!(y, H::SparseDynamicalMatrix, x, Hq)
    (; rowptr, cols, regularization) = H
    @assert length(y) == length(x) == size(H, 1)
    @inbounds for r in eachindex(y)
        acc = regularization * x[r]
        for e in rowptr[r]:rowptr[r+1]-1
            acc += Hq[e] * x[cols[e]]
        end
        y[r] = acc
    end
    return y
end

# Y = Ĩ X for each column of X, where Ĩ = diag(I, -I) is the Bogoliubov identity
function mul_Ĩ!(Y, X)
    L = size(X, 1) ÷ 2
    @views Y[1:L, :] .= X[1:L, :]
    @views Y[L+1:2L, :] .= .-X[L+1:2L, :]
    return Y
end
