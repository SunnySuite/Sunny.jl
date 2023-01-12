# Types of precomputed FFT plans
const rFTPlan = FFTW.rFFTWPlan{Float64, -1, false, 5, UnitRange{Int64}}
const rBFTPlan = FFTW.rFFTWPlan{ComplexF64, 1, false, 5, UnitRange{Int64}}
const rIFTPlan = FFTW.AbstractFFTs.ScaledPlan{ComplexF64, rBFTPlan, Float64}

struct EwaldCPU <: AbstractInteractionCPU
    A        :: Array{Mat3, 5}        # Interaction matrices in real-space, dims: [i,j,k,b1,b2]
    ϕ        :: Array{Vec3, 4}        # Space for holding real-space fields, dims: [i,j,k,b]
    FA       :: Array{ComplexF64, 7}  # Fourier transformed A, flattened dims: [α,β,m1,m2,m3,b1,b2]
    Fs       :: Array{ComplexF64, 5}  # Fourier transformed spins, dims: [α,m1,m2,m3,b]
    Fϕ       :: Array{ComplexF64, 5}  # Fourier-transformed fields, dims: [α,m1,m2,m3,b]
    plan     :: rFTPlan
    ift_plan :: rIFTPlan
end

function EwaldCPU(cryst::Crystal, sz::NTuple{3,Int}, gs::Vector{Mat3}, units::PhysicalConsts)
    (; μ0, μB) = units
    nb = nbasis(cryst)
    A = (μ0/4π) * μB^2 .* precompute_dipole_ewald(cryst, sz)
    # Scale interaction tensors by pair of g tensors for each site
    for b1 in 1:nb, b2 in 1:nb
        for ijk in CartesianIndices(sz)
            A[ijk, b1, b2] = gs[b1]' * A[ijk, b1, b2] * gs[b2]
        end
    end

    ϕ  = zeros(Vec3, sz..., nb)

    Ar = reshape(reinterpret(reshape, Float64, A), 3, 3, size(A)...) # dims: [α,β,i,j,k,b1,b2]
    FA = FFTW.rfft(Ar, 3:5) # FFT on cell indices (i,j,k)
    sz_rft = size(FA)[3:5]  # First FT dimension (dimension 3) will be ~ halved
    Fs = zeros(ComplexF64, 3, sz_rft..., nb)
    Fϕ = zeros(ComplexF64, 3, sz_rft..., nb)

    mock_spins = zeros(3, sz..., nb)
    plan     = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(Fs, sz[1], 2:4; flags=FFTW.MEASURE)

    return EwaldCPU(A, ϕ, FA, Fs, Fϕ, plan, ift_plan)
end


"Precompute the dipole interaction matrix, in ± compressed form."
function precompute_dipole_ewald(cryst::Crystal, sz::NTuple{3,Int}) :: Array{Mat3, 5}
    nb = nbasis(cryst)
    A = zeros(Mat3, sz..., nb, nb)

    # Superlattice vectors that describe periodicity of system and their inverse
    supervecs = cryst.lat_vecs .* Vec3(sz)'
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

    for idx in CartesianIndices(sz), b2 in 1:nb, b1 in 1:nb
        acc = zero(Mat3)
        # Note that `idx .- 1` represents a 3-vector offset
        Δr = position(cryst, b2, idx) - position(cryst, b1, (1,1,1))

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


function energy(dipoles::Array{Vec3, 4}, ewald::EwaldCPU)
    (; FA, Fs) = ewald
    sz = size(dipoles)[1:3]
    even_rft_size = sz[1] % 2 == 0

    E = 0.0
    mul!(Fs, ewald.plan, reinterpret(reshape, Float64, dipoles))

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
    @tullio E += real(
        conj(Fs[α, j, k, l, b1]) * conj(FA[α, β, j, k, l, b1, b2]) * Fs[β, j, k, l, b2]
    )
    return E / prod(sz)
end

"Use FFT to accumulate the entire field -dE/ds for long-range dipole-dipole interactions"
function accum_force!(B::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, ewald::EwaldCPU)
    (; FA, Fs, Fϕ, ϕ) = ewald

    fill!(Fϕ, 0.0)
    mul!(Fs, ewald.plan, reinterpret(reshape, Float64, dipoles))
    @tullio grad=false Fϕ[α,i,j,k,b1] += conj(FA[α,β,i,j,k,b1,b2]) * Fs[β,i,j,k,b2]
    ϕr = reinterpret(reshape, Float64, ϕ)
    mul!(ϕr, ewald.ift_plan, Fϕ)
    for i in eachindex(B)
        B[i] -= 2ϕ[i]
    end
end

"Calculate the field -dE/ds at idx1 generated by a dipole at idx2."
function pairwise_force_at(idx1, idx2, dipole::Vec3, ewald::EwaldCPU)
    sz = size(ewald.ϕ)[1:3]
    (cell1, b1) = splitidx(idx1)
    (cell2, b2) = splitidx(idx2)
    offset = mod.(Tuple(cell2-cell1), sz)
    idx = CartesianIndex(offset .+ (1,1,1))

    # A prefactor of 2 is always appropriate here. If idx1 == idx2, it accounts
    # for the quadratic dependence on the dipole. If idx1 != idx2, it accounts
    # for energy contributions from both ordered pairs (idx1,idx2) and
    # (idx2,idx1).
    return - 2ewald.A[idx, b1, b2] * dipole
end

"Calculate the field -dE/ds at idx generated by all `dipoles`."
function force_at(dipoles::Array{Vec3, 4}, ewald::EwaldCPU, idx)
    acc = zero(Vec3)
    for idx2 in CartesianIndices(dipoles)
        acc += pairwise_force_at(idx, idx2, dipoles[idx2], ewald)
    end
    return acc
end

# Calculate the change in dipole-dipole energy when the spin at `idx` is updated
# to `s`
function energy_delta(dipoles::Array{Vec3, 4}, ewald::EwaldCPU, idx, s::Vec3)
    Δs = s - dipoles[idx]
    b = idx[4]
    h = Sunny.force_at(dipoles, ewald, idx)
    return - Δs⋅h + dot(Δs, ewald.A[1, 1, 1, b, b], Δs)
end
