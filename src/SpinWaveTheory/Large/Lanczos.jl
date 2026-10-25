# Block Lanczos for A = Ĩ, the Bogoliubov identity, in the inner product of the
# positive definite S applied by `mulS!` [1]. The p rows of `V` are the starting
# block, and the iteration builds an S-orthonormal basis of block Krylov space,
#
#     Q₁ R₀ = V,    Q_{j+1} B_j = Ĩ S Q_j - Q_j A_j - Q_{j-1} B_{j-1}†,
#
# with A_j = (S Q_j)† Ĩ (S Q_j), each block orthogonalized only against its two
# predecessors. Here blocks are the rows of matrices with n columns. A block may
# have fewer rows than its predecessor, when the new directions of the Krylov
# space are linearly dependent, and the iteration is exhausted when none remain.
# The returned Lanczos object holds the blocks A_j and B_j, and the factor R₀.
# Call `extend!` to iterate further.
#
# Every iteration is one product with S, one fused pass for the p×p Gram
# matrices, and one fused update pass. The loops are written out, rather than
# calling BLAS, because this runs in each of many threads and OpenBLAS scales
# poorly under Julia threads.
#
# 1. G. H. Golub and R. Underwood, The block Lanczos method for computing
#    eigenvalues, in Mathematical Software III, 361-377 (1977).
mutable struct BlockLanczos{F}
    const mulS! :: F
    const R₀    :: Matrix{ComplexF64}
    const As    :: Vector{Matrix{ComplexF64}}
    const Bs    :: Vector{Matrix{ComplexF64}}
    # Buffers of the starting block size, whose leading `rank` rows hold the
    # current and previous blocks, the current block times S, and the residual
    # and its image under S. Rows beyond the rank are zero, so that products
    # with S act on whole buffers.
    Q           :: Matrix{ComplexF64}
    Qprev       :: Matrix{ComplexF64}
    SQ          :: Matrix{ComplexF64}
    W           :: Matrix{ComplexF64}
    SW          :: Matrix{ComplexF64}
    rank        :: Int
    rank_prev   :: Int
    exhausted   :: Bool
end

function BlockLanczos(mulS!, V::Matrix{ComplexF64})
    W = copy(V)
    SW = zero(V)
    mulS!(SW, W)
    G = gram_rows(W, SW)
    (P, R₀) = orthonormalize_rows(G; scale=maximum(real, diag(G)))
    r = size(P, 1)
    (Q, SQ) = (zero(V), zero(V))
    apply_rows!(view(Q, 1:r, :), P, W)
    apply_rows!(view(SQ, 1:r, :), P, SW)
    return BlockLanczos(mulS!, R₀, Matrix{ComplexF64}[], Matrix{ComplexF64}[], Q, zero(V), SQ, W, SW, r, 0, false)
end

# Given the Gram matrix G = W S W† of rows W, the r×p transform P and r×p factor
# B such that the rows of Q = P W are S-orthonormal and span those of W, with W
# = Bᵀ Q, where r is the numerical rank. With G = U D U†, P = D^{-1/2} Uᵀ and B =
# D^{1/2} U†. Directions with weight below a relative threshold of `scale` are
# dropped. A direction of negative norm means that S, the dynamical matrix, is
# not positive definite, so that the magnetic structure is unstable.
function orthonormalize_rows(G; scale)
    (; values, vectors) = eigen(Hermitian(G))
    minimum(values) < -sqrt(eps()) * scale && throw(InstabilityError("S is not positive definite"))
    keep = values .> sqrt(eps()) * scale
    (D, U) = (values[keep], vectors[:, keep])
    return (Diagonal(1 ./ sqrt.(D)) * transpose(U), Diagonal(sqrt.(D)) * U')
end

# Y = P X for rows, where P is small and dense
function apply_rows!(Y, P, X)
    @inbounds for k in axes(X, 2), a in axes(P, 1)
        acc = zero(ComplexF64)
        for b in axes(P, 2)
            acc += P[a, b] * X[b, k]
        end
        Y[a, k] = acc
    end
    return Y
end

nsteps(lz::BlockLanczos) = length(lz.As)

# G[a, b] = X[a, :] ⋅ Y[b, :]
function gram_rows(X, Y)
    p = size(X, 1)
    G = zeros(ComplexF64, p, p)
    @inbounds for k in axes(X, 2), b in 1:p, a in 1:p
        G[a, b] += conj(X[a, k]) * Y[b, k]
    end
    return G
end

# Performs `m` more iterations, or fewer if the Krylov space is exhausted
function extend!(lz::BlockLanczos, m)
    (; mulS!, As, Bs) = lz
    n = size(lz.Q, 2)
    L = n ÷ 2
    for _ in 1:m
        lz.exhausted && break
        (p, p′) = (lz.rank, lz.rank_prev)
        Q = view(lz.Q, 1:p, :)
        Qprev = view(lz.Qprev, 1:p′, :)
        SQ = view(lz.SQ, 1:p, :)
        W = view(lz.W, 1:p, :)
        SW = view(lz.SW, 1:p, :)
        # W = Ĩ S Q, and A = (S Q)† Ĩ (S Q)
        A = zeros(ComplexF64, p, p)
        @inbounds for k in 1:n
            sgn = k <= L ? 1 : -1
            for b in 1:p
                W[b, k] = sgn * SQ[b, k]
                for a in 1:p
                    A[a, b] += sgn * conj(SQ[a, k]) * SQ[b, k]
                end
            end
        end
        A = (A + A') / 2
        push!(As, A)
        # W -= Aᵀ Q + B̄ Q_prev as rows, where the previous coupling B has
        # dimensions p × p′
        Bprev = isempty(Bs) ? zeros(ComplexF64, p, 0) : Bs[end]
        @inbounds for k in 1:n, b in 1:p
            acc = W[b, k]
            for a in 1:p
                acc -= A[a, b] * Q[a, k]
            end
            for a in 1:p′
                acc -= conj(Bprev[b, a]) * Qprev[a, k]
            end
            W[b, k] = acc
        end
        # Local reorthogonalization against the current block, which costs one
        # more Gram matrix and restores orthogonality lost to cancellation
        C = gram_rows(SQ, W)
        @inbounds for k in 1:n, b in 1:p
            acc = W[b, k]
            for a in 1:p
                acc -= C[a, b] * Q[a, k]
            end
            W[b, k] = acc
        end
        # The next block, with W = Bᵀ Q_{j+1}. The current block becomes the
        # previous one by swapping buffers.
        mulS!(lz.SW, lz.W)
        G = gram_rows(W, SW)
        (P, B) = orthonormalize_rows(G; scale=maximum(abs, A))
        r = size(P, 1)
        if iszero(r)
            lz.exhausted = true
            break
        end
        push!(Bs, B)
        (lz.Qprev, lz.Q) = (lz.Q, lz.Qprev)
        apply_rows!(view(lz.Q, 1:r, :), P, W)
        apply_rows!(view(lz.SQ, 1:r, :), P, SW)
        if r < p
            for M in (lz.Q, lz.SQ, lz.W, lz.SW)
                view(M, r+1:size(M, 1), :) .= 0
            end
        end
        (lz.rank_prev, lz.rank) = (p, r)
    end
    return lz
end

# The Hermitian block tridiagonal Lanczos matrix T, with diagonal blocks A_j
# and subdiagonal blocks B_j, such that Q† S Ĩ S Q = T in the S-orthonormal
# block Krylov basis Q
function lanczos_matrix(lz::BlockLanczos)
    (; As, Bs) = lz
    m = length(As)
    offsets = cumsum([0; size.(As, 1)])
    T = zeros(ComplexF64, offsets[end], offsets[end])
    for j in 1:m
        r = offsets[j]+1:offsets[j+1]
        T[r, r] = As[j]
        if j < m
            s = offsets[j+1]+1:offsets[j+2]
            T[s, r] = Bs[j]
            T[r, s] = Bs[j]'
        end
    end
    return Hermitian(T)
end

# The leading block X₁(z) = [(z - T)⁻¹]₁₁ of the resolvent of the block
# tridiagonal Lanczos matrix, by the block continued fraction
#
#     X_m = (z - A_m)⁻¹,    X_j = (z - A_j - B_j† X_{j+1} B_j)⁻¹,
#
# which eliminates one block of z - T at a time, from the last to the first.
# Each X_j is the leading block of the resolvent of a Hermitian matrix, so its
# norm is at most 1/|Im z|, and the recursion is stable for Im z ≠ 0. It costs
# a few p×p operations per step, and avoids diagonalizing T.
#
# For speed, blocks are static p×p matrices. A block that shrank to r < p rows
# is padded with zeros, i.e., decoupled rows, which leave the leading block of
# the continued fraction unchanged.
struct StaticBlocks{p, p2}
    As :: Vector{SMatrix{p, p, ComplexF64, p2}}
    Bs :: Vector{SMatrix{p, p, ComplexF64, p2}}
    # Diagonal mask of the padded rows of each A
    pads :: Vector{SVector{p, Bool}}
end

function StaticBlocks(lz::BlockLanczos)
    p = size(lz.As[1], 1)
    return StaticBlocks(lz, Val{p}())
end

function StaticBlocks(lz::BlockLanczos, ::Val{p}) where p
    function pad(X)
        Y = zeros(ComplexF64, p, p)
        Y[axes(X)...] = X
        return SMatrix{p, p, ComplexF64}(Y)
    end
    As = pad.(lz.As)
    Bs = pad.(lz.Bs)
    pads = [SVector{p, Bool}(ntuple(i -> i > size(A, 1), p)) for A in lz.As]
    return StaticBlocks(As, Bs, pads)
end

function resolvent_leading_block(sb::StaticBlocks{p}, z) where p
    (; As, Bs, pads) = sb
    # Padded rows of A are zero and decoupled, so their diagonal entry of z - A
    # may be replaced by 1, which keeps the inverse finite without affecting the
    # leading block
    zI(pad) = Diagonal(map(b -> b ? one(z) : z, pad))
    m = length(As)
    X = inv(zI(pads[m]) - As[m])
    @inbounds for j in m-1:-1:1
        X = inv(zI(pads[j]) - As[j] - Bs[j]' * X * Bs[j])
    end
    return X
end
