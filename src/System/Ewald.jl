
function Ewald(sys::System{N}) where N
    (; crystal, latsize, units) = sys

    A = (units.μ0/4π) * precompute_dipole_ewald(crystal, latsize)

    na = natoms(crystal)
    μ = zeros(Vec3, latsize..., na)
    ϕ = zeros(Vec3, latsize..., na)

    Ar = reshape(reinterpret(Float64, A), 3, 3, size(A)...) # dims: [α,β,cell,i,j]
    FA = FFTW.rfft(Ar, 3:5) # FFT on cell indices
    sz_rft = size(FA)[3:5]  # First FT dimension (dimension 3) will be ~ halved
    Fμ = zeros(ComplexF64, 3, sz_rft..., na)
    Fϕ = zeros(ComplexF64, 3, sz_rft..., na)

    mock_spins = zeros(3, latsize..., na)
    plan     = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(Fμ, latsize[1], 2:4; flags=FFTW.MEASURE)

    return Ewald(A, μ, ϕ, FA, Fμ, Fϕ, plan, ift_plan)
end

# Clone all mutable state within Ewald. Note that `A`, `FA` are immutable data.
# We're also treating FFTW plans as immutable. This is not 100% correct, because
# plans mutably cache inverse plans, which may possibly lead to data races in a
# multithreaded context. Unfortunately, FFTW does not yet provide a copy
# function. For context, see https://github.com/JuliaMath/FFTW.jl/issues/261.
function clone_ewald(ewald::Ewald)
    (; A, μ, ϕ, FA, Fμ, Fϕ, plan, ift_plan) = ewald
    return Ewald(A, copy(μ), copy(ϕ), FA, copy(Fμ), copy(Fϕ), plan, ift_plan)
end

# Tensor product of 3-vectors
(⊗)(a::Vec3,b::Vec3) = reshape(kron(a,b), 3, 3)

# Precompute the dipole-dipole interaction matrix A and its Fourier transform
# F[A]
function precompute_dipole_ewald(cryst::Crystal, latsize::NTuple{3,Int}) :: Array{Mat3, 5}
    na = natoms(cryst)
    A = zeros(Mat3, latsize..., na, na)

    # Superlattice vectors and reciprocals for the full system volume
    sys_size = diagm(Vec3(latsize))
    latvecs = cryst.latvecs * sys_size
    recipvecs = sys_size \ cryst.recipvecs

    # Precalculate constants
    I₃ = Mat3(I)
    V = det(latvecs)
    L = cbrt(V)
    # Roughly balances the real and Fourier space costs
    σ = L/3
    σ² = σ*σ
    σ³ = σ^3
    # Corresponding to c0=6 in Ewalder.jl. Should give ~13 digits of accuracy.
    rmax = 6√2 * σ
    kmax = 6√2 / σ

    nmax = map(eachcol(latvecs), eachcol(recipvecs)) do a, b
        round(Int, rmax / (a⋅normalize(b)) + 1e-6) + 1
    end
    mmax = map(eachcol(latvecs), eachcol(recipvecs)) do a, b
        round(Int, kmax / (b⋅normalize(a)) + 1e-6)
    end

    # nmax and mmax should be balanced here
    # println("nmax $nmax mmax $mmax")

    for cell in CartesianIndices(latsize), j in 1:na, i in 1:na
        acc = zero(Mat3)
        cell_offset = Vec3(cell[1]-1, cell[2]-1, cell[3]-1)
        Δr = cryst.latvecs * (cell_offset + cryst.positions[j] - cryst.positions[i])
        
        #####################################################
        ## Real space part
        for n1 = -nmax[1]:nmax[1], n2 = -nmax[2]:nmax[2], n3 = -nmax[3]:nmax[3]
            rvec = Δr + latvecs * Vec3(n1, n2, n3)
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
            k = recipvecs * Vec3(m1, m2, m3)
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

        # For sites site1=(cell1, i) and site2=(cell2, j) offset by an amount
        # (off = cell2-cell1), the pair-energy is (s1 ⋅ A[off, i, j] ⋅ s2).
        # Julia arrays start at one, so we index A using (cell = off .+ 1).
        A[cell, i, j] = acc
    end

    return A
end


function ewald_energy(sys::System{N}) where N
    (; μ, FA, Fμ, plan) = sys.ewald
    latsize = size(sys.dipoles)[1:3]
    even_rft_size = latsize[1] % 2 == 0

    for site in all_sites(sys)
        μ[site] = magnetic_moment(sys, site)
    end

    E = 0.0
    mul!(Fμ, plan, reinterpret(reshape, Float64, μ))

    # rfft() is missing half the elements of the first Fourier transformed
    # dimension (here, dimension 2). Account for these missing values by scaling
    # the output by 2.
    if even_rft_size
        @views Fμ[:, 2:end-1, :, :, :] .*= √2
    else
        @views Fμ[:, 2:end, :, :, :] .*= √2
    end
    # * The field is a cross correlation, h(r) = - 2 A(Δr) s(r+Δr) = - 2A⋆s, or
    #   in Fourier space, F[h] = - 2 conj(F[A]) F[s]
    # * The energy is an inner product, E = - (1/2)s⋅h, or using Parseval's
    #   theorem, E = - (1/2) conj(F[s]) F[h] / N
    # * Combined, the result is: E = conj(F[s]) conj(F[A]) F[s] / N
    (_, m1, m2, m3, na) = size(Fμ)
    ms = CartesianIndices((m1, m2, m3))
    @inbounds for j in 1:na, i in 1:na, m in ms, α in 1:3, β in 1:3
        E += real(conj(Fμ[α, m, i]) * conj(FA[α, β, m, i, j]) * Fμ[β, m, j])
    end
    return E / prod(latsize)
end

# Use FFT to accumulate the entire field dE/ds for long-range dipole-dipole
# interactions
function accum_ewald_grad!(∇E::Array{Vec3, 4}, dipoles::Array{Vec3, 4}, sys::System{N}) where N
    (; ewald, units, gs) = sys
    (; μ, FA, Fμ, Fϕ, ϕ, plan, ift_plan) = ewald

    # Fourier transformed magnetic moments
    for site in all_sites(sys)
        μ[site] = units.μB * gs[site] * dipoles[site]
    end
    mul!(Fμ, plan, reinterpret(reshape, Float64, μ))

    # Calculate magneto-potential ϕ in Fourier space. Without @inbounds,
    # performance degrades by ~50%
    fill!(Fϕ, 0.0)
    (_, m1, m2, m3, na) = size(Fμ)
    ms = CartesianIndices((m1, m2, m3))
    @inbounds for j in 1:na, i in 1:na, m in ms, α in 1:3, β in 1:3
        Fϕ[α,m,i] += conj(FA[α,β,m,i,j]) * Fμ[β,m,j]
    end

    # Inverse Fourier transform to get ϕ in real space
    ϕr = reinterpret(reshape, Float64, ϕ)
    mul!(ϕr, ift_plan, Fϕ)
    
    for site in all_sites(sys)
        ∇E[site] += 2 * units.μB * (gs[site]' * ϕ[site])
    end
end

# Calculate the field dE/ds at site1 generated by a dipole at site2.
function ewald_pairwise_grad_at(sys::System{N}, site1, site2) where N
    (; gs, ewald, units) = sys
    latsize = size(ewald.ϕ)[1:3]
    cell_offset = mod.(Tuple(to_cell(site2)-to_cell(site1)), latsize)
    cell = CartesianIndex(cell_offset .+ (1,1,1))

    # A prefactor of 2 is always appropriate here. If site1 == site2, it
    # accounts for the quadratic dependence on the dipole. If site1 != site2, it
    # accounts for energy contributions from both ordered pairs (site1, site2)
    # and (site2, site1).
    return 2 * units.μB^2 * gs[site1]' * ewald.A[cell, to_atom(site1), to_atom(site2)] * gs[site2] * sys.dipoles[site2]
end

# Calculate the field dE/ds at `site` generated by all `dipoles`.
function ewald_grad_at(sys::System{N}, site) where N
    acc = zero(Vec3)
    for site2 in all_sites(sys)
        acc += ewald_pairwise_grad_at(sys, site, site2)
    end
    return acc
end

# Calculate the change in dipole-dipole energy when the spin at `site` is
# updated to `s`
function ewald_energy_delta(sys::System{N}, site, s::Vec3) where N
    (; dipoles, ewald, units) = sys
    Δs = s - dipoles[site]
    Δμ = units.μB * (sys.gs[site] * Δs)
    i = to_atom(site)
    ∇E = ewald_grad_at(sys, site)
    return Δs⋅∇E + dot(Δμ, ewald.A[1, 1, 1, i, i], Δμ)
end
