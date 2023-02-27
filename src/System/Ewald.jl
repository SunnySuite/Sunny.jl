
function Ewald(sys::System{N}) where N
    (; crystal, latsize, gs, units) = sys

    (; μ0, μB) = units
    nb = nbasis(crystal)
    A = (μ0/4π) * μB^2 .* precompute_dipole_ewald(crystal, latsize)
    # Scale g tensors into pair interactions A
    for b1 in 1:nb, b2 in 1:nb
        for cell in all_cells(sys)
            A[cell, b1, b2] = gs[b1]' * A[cell, b1, b2] * gs[b2]
        end
    end

    ϕ  = zeros(Vec3, latsize..., nb)

    Ar = reshape(reinterpret(Float64, A), 3, 3, size(A)...) # dims: [α,β,cell,b1,b2]
    FA = FFTW.rfft(Ar, 3:5) # FFT on cell indices
    sz_rft = size(FA)[3:5]  # First FT dimension (dimension 3) will be ~ halved
    Fs = zeros(ComplexF64, 3, sz_rft..., nb)
    Fϕ = zeros(ComplexF64, 3, sz_rft..., nb)

    mock_spins = zeros(3, latsize..., nb)
    plan     = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(Fs, latsize[1], 2:4; flags=FFTW.MEASURE)

    return Ewald(A, ϕ, FA, Fs, Fϕ, plan, ift_plan)
end

# Clone all mutable state within Ewald. Note that `A`, `FA`, and plans should be
# immutable data.
function clone_ewald(ewald::Ewald)
    (; A, ϕ, FA, Fs, Fϕ, plan, ift_plan) = ewald
    return Ewald(A, copy(ϕ), FA, copy(Fs), copy(Fϕ), plan, ift_plan)
end

# Tensor product of 3-vectors
(⊗)(a::Vec3,b::Vec3) = reshape(kron(a,b), 3, 3)

# Precompute the dipole-dipole interaction matrix A and its Fourier transform
# F[A]
function precompute_dipole_ewald(cryst::Crystal, latsize::NTuple{3,Int}) :: Array{Mat3, 5}
    nb = nbasis(cryst)
    A = zeros(Mat3, latsize..., nb, nb)

    # Superlattice vectors that describe periodicity of system and their inverse
    supervecs = cryst.lat_vecs .* Vec3(latsize)'
    recipvecs = 2π * inv(supervecs)
    # Split into individual vectors
    supervecs = collect(eachcol(supervecs))
    recipvecs = collect(eachrow(recipvecs))

    # Precalculate constants
    I₃ = Mat3(I)
    V = (supervecs[1] × supervecs[2]) ⋅ supervecs[3]
    L = cbrt(V)
    # Roughly balances the real and Fourier space costs
    σ = L/3
    σ² = σ*σ
    σ³ = σ^3
    # Corresponding to c0=6 in Ewalder.jl. Should give ~13 digits of accuracy.
    rmax = 6√2 * σ
    kmax = 6√2 / σ

    nmax = map(supervecs, recipvecs) do a, b
        round(Int, rmax / (a⋅normalize(b)) + 1e-6) + 1
    end
    mmax = map(supervecs, recipvecs) do a, b
        round(Int, kmax / (b⋅normalize(a)) + 1e-6)
    end

    # nmax and mmax should be balanced here
    # println("nmax $nmax mmax $mmax")

    for idx in CartesianIndices(latsize), b2 in 1:nb, b1 in 1:nb
        acc = zero(Mat3)
        cell_offset = Vec3(idx[1]-1, idx[2]-1, idx[3]-1)
        Δr = cryst.lat_vecs * (cell_offset + cryst.positions[b2] - cryst.positions[b1])
        
        #####################################################
        ## Real space part
        for n1 = -nmax[1]:nmax[1], n2 = -nmax[2]:nmax[2], n3 = -nmax[3]:nmax[3]
            n = Vec3(n1, n2, n3)
            rvec = Δr + n' * supervecs
            r² = rvec⋅rvec
            if 0 < r² <= rmax*rmax
                r = √r²
                r³ = r²*r
                rhat = rvec/r
                erfc0 = erfc(r/(√2*σ))
                gauss0 = √(2/π) * (r/σ) * exp(-r²/2σ²)    
                acc += (1/2) * ((I₃/r³) * (erfc0 + gauss0) - (3(rhat⊗rhat)/r³) * (erfc0 + (1+r²/3σ²) * gauss0))
            end
        end

        #####################################################
        ## Fourier space part
        for m1 = -mmax[1]:mmax[1], m2 = -mmax[2]:mmax[2], m3 = -mmax[3]:mmax[3]
            k = Vec3(m1, m2, m3)' * recipvecs
            k² = k⋅k
            if 0 < k² <= kmax*kmax
                # Replace exp(-ikr) -> cos(kr). It's valid to drop the imaginary
                # component because it will cancel in the interaction that exchanges
                # i ↔ j.
                acc += (4π/2V) * (exp(-σ²*k²/2) / k²) * (k⊗k) * cos(k⋅Δr)
            end
        end

        #####################################################
        ## Remove self energies
        if iszero(Δr)
            acc += - I₃/(3√(2π)*σ³)
        end

        # For sites idx1=(cell1, b1) and idx2=(cell2, b2) offset by an amount
        # (off = cell2-cell1), the pair-energy is (s1 ⋅ A[off, b1, b2] ⋅ s2).
        # Julia arrays start at one, so we index A using (idx = off .+ 1).
        A[idx, b1, b2] = acc
    end

    return A
end


function energy(dipoles::Array{Vec3, 4}, ewald::Ewald)
    (; FA, Fs, plan) = ewald
    latsize = size(dipoles)[1:3]
    even_rft_size = latsize[1] % 2 == 0

    E = 0.0
    mul!(Fs, plan, reinterpret(reshape, Float64, dipoles))

    # rfft() is missing half the elements of the first Fourier transformed
    # dimension (here, dimension 2). Account for these missing values by scaling
    # the output by 2.
    if even_rft_size
        @views Fs[:, 2:end-1, :, :, :] .*= √2
    else
        @views Fs[:, 2:end, :, :, :] .*= √2
    end
    # * The field is a cross correlation, h(r) = - 2 A(Δr) s(r+Δr) = - 2A⋆s, or
    #   in Fourier space, F[h] = - 2 conj(F[A]) F[s]
    # * The energy is an inner product, E = - (1/2)s⋅h, or using Parseval's
    #   theorem, E = - (1/2) conj(F[s]) F[h] / N
    # * Combined, the result is: E = conj(F[s]) conj(F[A]) F[s] / N
    (_, m1, m2, m3, nb) = size(Fs)
    ms = CartesianIndices((m1, m2, m3))
    @inbounds for b2 in 1:nb, b1 in 1:nb, m in ms, α in 1:3, β in 1:3
        E += real(conj(Fs[α, m, b1]) * conj(FA[α, β, m, b1, b2]) * Fs[β, m, b2])
    end
    return E / prod(latsize)
end

# Use FFT to accumulate the entire field -dE/ds for long-range dipole-dipole
# interactions
function accum_force!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ewald::Ewald)
    (; FA, Fs, Fϕ, ϕ, plan, ift_plan) = ewald

    fill!(Fϕ, 0.0)
    mul!(Fs, plan, reinterpret(reshape, Float64, dipoles))
    (_, m1, m2, m3, nb) = size(Fs)
    ms = CartesianIndices((m1, m2, m3))
    # Without @inbounds, performance degrades by ~50%
    @inbounds for b2 in 1:nb, b1 in 1:nb, m in ms, α in 1:3, β in 1:3
        Fϕ[α,m,b1] += conj(FA[α,β,m,b1,b2]) * Fs[β,m,b2]
    end
    ϕr = reinterpret(reshape, Float64, ϕ)
    mul!(ϕr, ift_plan, Fϕ)
    for i in eachindex(B)
        B[i] -= 2ϕ[i]
    end
end

# Calculate the field -dE/ds at idx1 generated by a dipole at idx2.
function pairwise_force_at(idx1, idx2, dipole::Vec3, ewald::Ewald)
    latsize = size(ewald.ϕ)[1:3]
    (cell1, b1) = splitidx(idx1)
    (cell2, b2) = splitidx(idx2)
    offset = mod.(Tuple(cell2-cell1), latsize)
    idx = CartesianIndex(offset .+ (1,1,1))

    # A prefactor of 2 is always appropriate here. If idx1 == idx2, it accounts
    # for the quadratic dependence on the dipole. If idx1 != idx2, it accounts
    # for energy contributions from both ordered pairs (idx1,idx2) and
    # (idx2,idx1).
    return - 2ewald.A[idx, b1, b2] * dipole
end

# Calculate the field -dE/ds at idx generated by all `dipoles`.
function force_at(dipoles::Array{Vec3, 4}, ewald::Ewald, idx)
    acc = zero(Vec3)
    for idx2 in CartesianIndices(dipoles)
        acc += pairwise_force_at(idx, idx2, dipoles[idx2], ewald)
    end
    return acc
end

# Calculate the change in dipole-dipole energy when the spin at `idx` is updated
# to `s`
function energy_delta(dipoles::Array{Vec3, 4}, ewald::Ewald, idx, s::Vec3)
    Δs = s - dipoles[idx]
    b = idx[4]
    h = Sunny.force_at(dipoles, ewald, idx)
    return - Δs⋅h + dot(Δs, ewald.A[1, 1, 1, b, b], Δs)
end
