# The resolvent y = (z - T)⁻¹ e₁ at each of the nodes `zs`, for the tridiagonal
# Lanczos matrix T of A = Ĩ, the Bogoliubov identity, in the inner product of
# the positive definite S applied by `mulS!` [1]. From the starting vector v,
# Lanczos builds an S-orthonormal basis of Krylov space,
#
#     q₁ β₀ = v,    q_{j+1} β_j = Ĩ S q_j - α_j q_j - β_{j-1} q_{j-1},
#
# with α_j = (S q_j)† Ĩ (S q_j). The real coefficients α_j and β_j are the
# diagonal and subdiagonal of T = Q† S Ĩ S Q. Elimination from the first row, z
# - T = L D Lᵀ, with pivots D = diag(u) and g = L⁻¹ e₁, gives the leading entry
# x₁ = y₁, and the probe sums sₐ = dₐᵀ y for the overlaps dₐⱼ = vₐ† S q_j of the
# columns vₐ of `probes`,
#
#     x₁ = Σⱼ gⱼ² / uⱼ,    sₐ = Σⱼ gⱼ hₐⱼ,    u_j = z - α_j - β_{j-1}² / u_{j-1},
#     g_j = β_{j-1} g_{j-1} / u_{j-1},    hₐⱼ = (dₐⱼ + β_{j-1} hₐ,ⱼ₋₁) / u_j.
#
# Each term depends only on the coefficients up to j, so the sums are
# accumulated as the iteration proceeds, and need no history. For Im z > 0, Im
# u_j ≥ Im z, so that no pivot vanishes. The Krylov space is exhausted when the
# residual vanishes, which makes the sums exact.
#
# After m rows, x_m = g_m / u_m is the last entry of y. The further rows of T,
# coupled by β_m, change x₁ by δx₁ = β_m² x_m² τ / (1 - β_m² τ / u_m), by the
# Sherman–Morrison formula, where τ is the leading entry of the resolvent of the
# further rows. Taking τ ≈ 1/u_m, as if the further coefficients resembled those
# so far, gives the estimate |δx₁| ≈ β_m² |x_m|² |τ|. It is large near any mode
# that T does not yet resolve, i.e., where the Krylov basis has not yet
# separated eigenvalues closer than Im z.
#
# The probe sums converge more slowly, because a probe lies mostly outside the
# Krylov space. In full space, the error of Q y is e = β_m x_m (z - ĨS)⁻¹ q_{m+1},
# so that |δsₐ| = |vₐ† S e| ≤ ‖vₐ‖ ‖e‖ by Cauchy–Schwarz in the S inner product.
# For S-normalized q, ‖(z - ĨS)⁻¹ q‖² = -Im[q† S (z - ĨS)⁻¹ q] / Im z, by the
# spectral theorem for ĨS, which is self-adjoint in this inner product. With the
# same τ ≈ 1/u_m, the estimate is ‖e‖ ≈ β_m |x_m| (-Im τ / Im z)^{1/2}.
#
# Every step is one product with S, two passes over vectors of length n, one
# pass per probe, and a pass over the nodes. The loops are written out, rather
# than calling BLAS, because this runs in each of many threads and OpenBLAS
# scales poorly under Julia threads.
#
# 1. T. Chen, The Lanczos algorithm for matrix functions: a handbook for
#    scientists (2024), Sec. 3.
mutable struct LanczosResolvent{F}
    const mulS!  :: F
    const probes :: Matrix{ComplexF64}
    const zs     :: Vector{ComplexF64}
    const β₀     :: Float64
    const G      :: Vector{ComplexF64}  # v† Ĩ v, then vₐ† Ĩ v for each probe
    # The current and previous Lanczos vectors, the current vector times S, and
    # the residual and its image under S
    q            :: Vector{ComplexF64}
    qprev        :: Vector{ComplexF64}
    Sq           :: Vector{ComplexF64}
    w            :: Vector{ComplexF64}
    Sw           :: Vector{ComplexF64}
    # The rows eliminated, and the coupling β_m of the last to the next
    nrows        :: Int
    β            :: Float64
    # Per node, 1/u_m, g_m, and the sums x₁ and sₐ, with hₐₘ for each probe
    const invu   :: Vector{ComplexF64}
    const g      :: Vector{ComplexF64}
    const x₁     :: Vector{ComplexF64}
    const h      :: Matrix{ComplexF64}
    const s      :: Matrix{ComplexF64}
end

function LanczosResolvent(mulS!, v, probes, zs)
    Sv = mulS!(similar(v), v)
    β₀² = real(dot_noblas(v, Sv))
    β₀² > 0 || throw(InstabilityError("S is not positive definite"))
    β₀ = sqrt(β₀²)
    Ĩv = mul_Ĩ!(similar(v), v)
    G = [dot_noblas(v, Ĩv); [dot_noblas(view(probes, :, a), Ĩv) for a in axes(probes, 2)]]
    (nz, np) = (length(zs), size(probes, 2))
    z(dims...) = zeros(ComplexF64, dims...)
    return LanczosResolvent(mulS!, probes, zs, β₀, G, v / β₀, zero(v), Sv / β₀, zero(v), zero(v),
                            0, 0.0, z(nz), ones(ComplexF64, nz), z(nz), z(np, nz), z(np, nz))
end

exhausted(r::LanczosResolvent) = r.nrows > 0 && iszero(r.β)

# x† y, written out because BLAS scales poorly under Julia threads
function dot_noblas(x, y)
    acc = zero(ComplexF64)
    @inbounds for k in eachindex(x, y)
        acc += conj(x[k]) * y[k]
    end
    return acc
end

# Eliminates rows of T up to `m` in total, or fewer if the Krylov space is
# exhausted. Each row j is eliminated with α_j and the overlaps of q_j, and the
# Lanczos step that follows gives the coupling β_j to the next row.
function extend!(r::LanczosResolvent, m)
    (; mulS!, probes, zs, invu, g, x₁, h, s) = r
    n = length(r.q)
    L = n ÷ 2
    while r.nrows < m && !exhausted(r)
        (; q, qprev, Sq, w, Sw) = r
        β′ = r.β
        # α_j = (S q_j)† Ĩ (S q_j), and the overlaps of the probes
        α = 0.0
        @inbounds for k in 1:n
            α += (k <= L ? 1 : -1) * abs2(Sq[k])
        end
        d = [dot_noblas(view(probes, :, a), Sq) for a in axes(probes, 2)]
        j = r.nrows += 1
        @inbounds for i in eachindex(zs)
            if j > 1
                g[i] *= β′ * invu[i]
            end
            invu[i] = 1 / (zs[i] - α - β′^2 * invu[i])
            x₁[i] += g[i]^2 * invu[i]
            for a in eachindex(d)
                h[a, i] = (d[a] + β′ * h[a, i]) * invu[i]
                s[a, i] += g[i] * h[a, i]
            end
        end

        # The residual w = Ĩ S q - α q - β′ q_prev, followed by a local
        # reorthogonalization against q, which restores orthogonality lost to
        # cancellation
        c = zero(ComplexF64)
        @inbounds for k in 1:n
            wk = (k <= L ? 1 : -1) * Sq[k] - α * q[k] - β′ * qprev[k]
            w[k] = wk
            c += conj(Sq[k]) * wk
        end
        @inbounds for k in 1:n
            w[k] -= c * q[k]
        end
        mulS!(Sw, w)
        # A residual of negative norm means that S, the dynamical matrix, is
        # not positive definite, so that the magnetic structure is unstable. A
        # residual small relative to Ĩ S q exhausts the space.
        β² = real(dot_noblas(w, Sw))
        scale = α^2 + β′^2 + abs(β²)
        β² < -sqrt(eps()) * scale && throw(InstabilityError("S is not positive definite"))
        r.β = β² <= sqrt(eps()) * scale ? 0.0 : sqrt(β²)
        exhausted(r) && break
        # The next vector. The current one becomes the previous one by swapping
        # buffers.
        (r.q, r.qprev, r.w, r.Sq, r.Sw) = (w, q, qprev, Sw, Sq)
        r.q ./= r.β
        r.Sq ./= r.β
    end
    return r
end

# The estimates of |δx₁| and ‖e‖ at node i
function error_estimates(r::LanczosResolvent, i)
    (τ, xm) = (r.invu[i], r.g[i] * r.invu[i])
    return (r.β^2 * abs2(xm) * abs(τ), r.β * abs(xm) * sqrt(-imag(τ) / imag(r.zs[i])))
end
