# The Ewald pair term T† Ã(𝐪) T of `QuadraticEwald`, applied by FFT. With the
# Bloch phase φ(𝐜) = exp(-2πi 𝐪⋅𝐜/dims) of each cell, the winding phase is
# exp(-2πi 𝐪⋅w) = φ(Δ)/φ(Δ mod dims). The phase φ(Δ) = φ(𝐜ⱼ)/φ(𝐜ᵢ) is applied
# with T, and what remains is a circular convolution over the cells, with
# kernel K[𝐦] = A(𝐪)[𝐦] conj(φ(𝐦)), which an FFT diagonalizes.
#
# The operator is set to a wavevector by `set_wavevector!`, and applied by
# `mul!`. It holds the kernel and phases at that wavevector, scratch space, and
# FFT plans, and is owned by one thread. The field transforms run at every
# product and are planned with MEASURE. FFTW caches what it learns, so that
# planning again for each thread is fast. The coefficients of T, the kernel, and
# the fields are in precision FT. The Ewald tensor A(𝐪) is calculated in Float64
# and rounded.
struct EwaldOperator{FT}
    quadratic :: QuadraticEwald
    coefs     :: Vector{SVector{3, Complex{FT}}}  # The coefficients of T
    K         :: Array{Complex{FT}, 6}            # Transformed kernel [αβ, 𝐦, i, j]
    φ         :: Vector{Complex{FT}}              # Bloch phase of each bare site
    v         :: Array{Complex{FT}, 5}            # Fields [α, 𝐦, i]
    u         :: Array{Complex{FT}, 5}
    fft       :: FFTW.cFFTWPlan{Complex{FT}, -1, true, 5, UnitRange{Int}}
    bfft      :: FFTW.cFFTWPlan{Complex{FT}, 1, true, 5, UnitRange{Int}}
    kfft      :: FFTW.cFFTWPlan{Complex{FT}, 1, true, 6, UnitRange{Int}}
end

function EwaldOperator(quadratic::QuadraticEwald, ::Type{FT}) where FT
    (; dims, cryst) = quadratic.cache
    na = natoms(cryst)
    v = zeros(Complex{FT}, 3, dims..., na)
    K = zeros(Complex{FT}, 9, dims..., na, na)
    fft = FFTW.plan_fft!(v, 2:4; flags=FFTW.MEASURE, num_threads=1)
    bfft = FFTW.plan_bfft!(v, 2:4; flags=FFTW.MEASURE, num_threads=1)
    kfft = FFTW.plan_bfft!(K, 2:4; num_threads=1)
    coefs = SVector{3, Complex{FT}}.(quadratic.coefs)
    φ = zeros(Complex{FT}, length(quadratic.cells))
    return EwaldOperator(quadratic, coefs, K, φ, v, zero(v), fft, bfft, kfft)
end

function set_wavevector!(op::EwaldOperator{FT}, q_supercell) where FT
    (; cache, cells) = op.quadratic
    (; dims) = cache
    na = natoms(cache.cryst)
    q = Vec3(q_supercell)
    A = ewald_interaction_tensor(cache, q)
    K = reinterpret(reshape, SMatrix{3, 3, Complex{FT}, 9}, op.K)
    for m in CartesianIndices(dims), j in 1:na, i in 1:na
        # The normalization 1/Ncells of the inverse transform
        K[m, i, j] = cis(2π * dot(q, Vec3(m.I .- 1) ./ dims)) * A[m, i, j] / prod(dims)
    end
    op.kfft * op.K
    map!(c -> cis(-2π * dot(q, Vec3(c ./ dims))), op.φ, cells)
    return op
end

# y = α H_E(𝐪) x + β y, at the wavevector of `set_wavevector!`
function LinearAlgebra.mul!(y, op::EwaldOperator{FT}, x, α::Number, β::Number) where FT
    (; ptr, labels, cache) = op.quadratic
    (; dims) = cache
    (; coefs, φ, fft, bfft) = op
    na = natoms(cache.cryst)
    v = reinterpret(reshape, SVector{3, Complex{FT}}, op.v)
    u = reinterpret(reshape, SVector{3, Complex{FT}}, op.u)
    K = reinterpret(reshape, SMatrix{3, 3, Complex{FT}, 9}, op.K)
    @inbounds for i in eachindex(φ)
        acc = zero(SVector{3, Complex{FT}})
        for e in ptr[i]:ptr[i+1]-1
            acc += coefs[e] * x[labels[e]]
        end
        v[i] = φ[i] * acc
    end
    fft * op.v
    @inbounds for i in 1:na, m in CartesianIndices(dims)
        acc = zero(SVector{3, Complex{FT}})
        for j in 1:na
            acc += K[m, i, j] * v[m, j]
        end
        u[m, i] = acc
    end
    bfft * op.u
    if iszero(β)
        fill!(y, 0)
    elseif !isone(β)
        y .*= β
    end
    @inbounds for i in eachindex(φ)
        ui = FT(α) * conj(φ[i]) * u[i]
        for e in ptr[i]:ptr[i+1]-1
            y[labels[e]] += dot(coefs[e], ui)
        end
    end
    return y
end

# The Hamiltonian matrix H(𝐪) of `hamiltonian_matrix!`, as a matrix-free linear
# operator for `SpinWaveTheoryLarge`. It is set to a wavevector by
# `set_wavevector!`, and applied by `mul!`. It holds the values of the quadratic
# table with their phases exp(2πi 𝐪⋅𝐧ₑ) folded in, and the Ewald state at 𝐪, so
# it is owned by one thread. The products are in precision FT. The phased values
# are calculated in Float64 and rounded.
struct HamiltonianOperator{FT}
    swt    :: SpinWaveTheory
    vals_q :: Vector{Complex{FT}}
    ewald  :: Union{Nothing, EwaldOperator{FT}}
end

function HamiltonianOperator(swt::SpinWaveTheory, ::Type{FT}=Float64) where FT
    ewald = isnothing(swt.ewald) ? nothing : EwaldOperator(swt.ewald, FT)
    vals_q = zeros(Complex{FT}, length(swt.quadratic_table.vals))
    return HamiltonianOperator{FT}(swt, vals_q, ewald)
end

Base.size(H::HamiltonianOperator, _) = size(H.swt.quadratic_table, 1)

function set_wavevector!(H::HamiltonianOperator, q_supercell)
    (; swt, vals_q, ewald) = H
    t = swt.quadratic_table
    φ = offset_phases(t, q_supercell)
    @inbounds for e in eachindex(vals_q)
        vals_q[e] = φ[t.offset_index[e]] * t.vals[e]
    end
    isnothing(ewald) || set_wavevector!(ewald, q_supercell)
    return H
end

# y = H(𝐪) x, at the wavevector of `set_wavevector!`
function LinearAlgebra.mul!(y, H::HamiltonianOperator{FT}, x) where FT
    (; swt, vals_q, ewald) = H
    (; rowptr, cols) = swt.quadratic_table
    @assert length(y) == length(x) == size(H, 1)
    ε = FT(swt.regularization)
    @inbounds for r in eachindex(y)
        acc = ε * x[r]
        for e in rowptr[r]:rowptr[r+1]-1
            acc += vals_q[e] * x[cols[e]]
        end
        y[r] = acc
    end
    isnothing(ewald) || mul!(y, ewald, x, true, true)
    return y
end
