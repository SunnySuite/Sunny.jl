
function Ewald(sys::System{N}, μ0_μB², demag) where N
    (; crystal, dims) = sys

    A = precompute_dipole_ewald(crystal, dims, demag) * μ0_μB²

    na = natoms(crystal)
    μ = zeros(Vec3, dims..., na)
    ϕ = zeros(Vec3, dims..., na)

    Ar = reshape(reinterpret(Float64, A), 3, 3, size(A)...) # dims: [α,β,cell,i,j]
    FA = FFTW.rfft(Ar, 3:5) # FFT on cell indices
    sz_rft = size(FA)[3:5]  # First FT dimension (dimension 3) will be ~ halved
    Fμ = zeros(ComplexF64, 3, sz_rft..., na)
    Fϕ = zeros(ComplexF64, 3, sz_rft..., na)

    mock_spins = zeros(3, dims..., na)
    plan     = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(Fμ, dims[1], 2:4; flags=FFTW.MEASURE)

    return Ewald(μ0_μB², demag, A, μ, ϕ, FA, Fμ, Fϕ, plan, ift_plan)
end

# Ideally, this would clone all mutable state within Ewald. Note that `A`, `FA`
# are immutable data. A blocker is that FFTW plans cannot currently be copied,
# and it is not 100% clear whether they can be treated as immutable. For
# example, they cache inverse plans, which may possibly lead to data races in a
# multithreaded context. See discussion at
# https://github.com/JuliaMath/FFTW.jl/issues/261.
function clone_ewald(ewald::Ewald)
    error("Not supported")
    (; μ0_μB², demag, A, μ, ϕ, FA, Fμ, Fϕ, plan, ift_plan) = ewald
    return Ewald(μ0_μB², demag, A, copy(μ), copy(ϕ), FA, copy(Fμ), copy(Fϕ), copy(plan), copy(ift_plan))
end

# Tensor product of 3-vectors
(⊗)(a::Vec3,b::Vec3) = reshape(kron(a,b), 3, 3)


function precompute_dipole_ewald(cryst::Crystal, dims::NTuple{3,Int}, demag::Mat3)
    precompute_dipole_ewald_aux(cryst, dims, demag, Vec3(0,0,0), cos, Val{Float64}())
end

function precompute_dipole_ewald_at_wavevector(cryst::Crystal, dims::NTuple{3,Int}, demag::Mat3, q_reshaped::Vec3)
    precompute_dipole_ewald_aux(cryst, dims, demag, q_reshaped, cis, Val{ComplexF64}())
end

# Precompute the pairwise interaction matrix A between magnetic moments μ. For
# q_reshaped = 0, this yields the usual Ewald energy, E = μᵢ Aᵢⱼ μⱼ / 2. Nonzero
# q_reshaped is useful in spin wave theory. Physically, this amounts to a
# modification of the periodic boundary conditions, such that μ(q) can be
# incommensurate with the magnetic cell. In all cases, the energy is E = μᵢ(-q)
# Aᵢⱼ(-q) μⱼ(q) / 2 in Fourier space, where q should be interpreted as a Fourier
# transform of the cell offset.
#
# As an optimization, this function returns real values when q_reshaped is zero.
# Effectively, one can replace `exp(i (q+k)⋅r) → cos(k⋅r)` because the imaginary
# part cancels in the symmetric sum over ±k. Specifically, replace `cis(x) ≡
# exp(i x) = cos(x) + i sin(x)` with just `cos(x)` for efficiency. The parameter
# `T ∈ {Float64, ComplexF64}` controls the return type in a type-stable way.
function precompute_dipole_ewald_aux(cryst::Crystal, dims::NTuple{3,Int}, demag, q_reshaped, cis, ::Val{T}) where T
    na = natoms(cryst)
    A = zeros(SMatrix{3, 3, T, 9}, dims..., na, na)

    # Superlattice vectors and reciprocals for the full system volume
    sys_size = diagm(Vec3(dims))
    latvecs = cryst.latvecs * sys_size
    recipvecs = cryst.recipvecs / sys_size

    # Precalculate constants
    I₃ = Mat3(I)
    V = det(latvecs)
    L = cbrt(V)
    # Roughly balances the real and Fourier space costs. Note that σ = 1/(√2 λ)
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

    for cell in CartesianIndices(dims), j in 1:na, i in 1:na
        acc = zero(eltype(A))
        cell_offset = Vec3(cell[1]-1, cell[2]-1, cell[3]-1)
        Δr = cryst.latvecs * (cell_offset + cryst.positions[j] - cryst.positions[i])

        #####################################################
        ## Real space part
        for n1 = -nmax[1]:nmax[1], n2 = -nmax[2]:nmax[2], n3 = -nmax[3]:nmax[3]
            n = Vec3(n1, n2, n3)
            rvec = Δr + latvecs * n
            r² = rvec⋅rvec
            if 0 < r² <= rmax*rmax
                r = √r²
                r³ = r²*r
                rhat = rvec/r
                erfc0 = erfc(r/(√2*σ))
                gauss0 = √(2/π) * (r/σ) * exp(-r²/2σ²)
                phase = cis(2π * dot(q_reshaped, n))
                acc += phase * (1/4π) * ((I₃/r³) * (erfc0 + gauss0) - (3(rhat⊗rhat)/r³) * (erfc0 + (1+r²/3σ²) * gauss0))
            end
        end

        #####################################################
        ## Fourier space part
        for m1 = -mmax[1]:mmax[1], m2 = -mmax[2]:mmax[2], m3 = -mmax[3]:mmax[3]
            m = Vec3(m1, m2, m3)
            k = recipvecs * (m + q_reshaped - round.(q_reshaped))
            k² = k⋅k

            ϵ² = 1e-16
            if k² <= ϵ²
                # Surface term Eₛ = μ₀ M⋅N M / 2V gives rise to demagnetization
                # effect. Net magnetization M is associated with mode k = 0.
                # Demagnetization factor tensor N, denoted `demag`, depends on
                # sample geometry and has trace 1 in vacuum background. This
                # Ewald correction was originally derived in S. W. DeLeeuw et
                # al., Proc. R. Soc. Lond. A 373, 27-56 (1980). See Ballenegger,
                # J. Chem. Phys. 140, 161102 (2014) for a pedagogical review.
                acc += demag / V
            elseif ϵ² < k² <= kmax*kmax
                phase = cis(-k⋅Δr)
                acc += phase * (1/V) * (exp(-σ²*k²/2) / k²) * (k⊗k)
            end
        end

        #####################################################
        ## Remove self energies
        if iszero(Δr)
            acc += - I₃/(3(2π)^(3/2)*σ³)
        end

        # For sites site1=(cell1, i) and site2=(cell2, j) offset by an amount
        # (off = cell2-cell1), the pair-energy is (s1 ⋅ A[off, i, j] ⋅ s2).
        # Julia arrays start at one, so we index A using (cell = off .+ 1).
        A[cell, i, j] = acc
    end

    # TODO: Verify that A[off, i, j] ≈ A[-off, j, i]'

    return A
end


function ewald_energy(sys::System{N}) where N
    (; μ, FA, Fμ, plan) = sys.ewald
    dims = size(sys.dipoles)[1:3]
    even_rft_size = dims[1] % 2 == 0

    E = 0.0
    μ .= magnetic_moments(sys)
    mul!(Fμ, plan, reinterpret(reshape, Float64, μ))

    # rfft() is missing half the elements of the first Fourier transformed
    # dimension (here, dimension 2). Account for these missing values by scaling
    # the output by 2.
    if even_rft_size
        @views Fμ[:, 2:end-1, :, :, :] .*= √2
    else
        @views Fμ[:, 2:end, :, :, :] .*= √2
    end

    # In real space, E = μ (A ⋆ μ) / 2. In Fourier space, the convolution
    # becomes an ordinary product using Parseval's theorem.
    (_, m1, m2, m3, na) = size(Fμ)
    ms = CartesianIndices((m1, m2, m3))
    @inbounds for j in 1:na, i in 1:na, m in ms, α in 1:3, β in 1:3
        E += (1/2) * real(conj(Fμ[α, m, i]) * conj(FA[α, β, m, i, j]) * Fμ[β, m, j])
    end
    return E / prod(dims)
end

# Use FFT to accumulate the entire field dE/dS for long-range dipole-dipole
# interactions
function accum_ewald_grad!(∇E, dipoles, sys::System{N}) where N
    (; ewald, gs) = sys
    (; μ, FA, Fμ, Fϕ, ϕ, plan, ift_plan) = ewald

    # Fourier transformed magnetic moments
    μ .= magnetic_moments(sys)
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

    for site in eachsite(sys)
        ∇E[site] -= gs[site]' * ϕ[site]
    end
end

# Calculate the field dE/dS at site1 generated by a dipole at site2.
function ewald_pairwise_grad_at(sys::System{N}, site1, site2) where N
    (; gs, ewald) = sys
    dims = size(ewald.ϕ)[1:3]
    cell_offset = mod.(to_cell(site2) .- to_cell(site1), dims)
    cell = CartesianIndex(cell_offset .+ (1,1,1))

    # The factor of 1/2 in the energy formula `E = μ (A ⋆ μ) / 2` disappears due
    # to quadratic appearance of μ = - g S.
    return gs[site1]' * ewald.A[cell, to_atom(site1), to_atom(site2)] * gs[site2] * sys.dipoles[site2]
end

# Calculate the field dE/dS at `site` generated by all `dipoles`.
function ewald_grad_at(sys::System{N}, site) where N
    acc = zero(Vec3)
    for site2 in eachsite(sys)
        acc += ewald_pairwise_grad_at(sys, site, site2)
    end
    return acc
end

# Calculate the change in dipole-dipole energy when the spin at `site` is
# updated to `S`
function ewald_energy_delta(sys::System{N}, site, S::Vec3) where N
    (; dipoles, ewald) = sys
    ΔS = S - dipoles[site]
    Δμ = sys.gs[site] * ΔS
    i = to_atom(site)
    ∇E = ewald_grad_at(sys, site)
    return ΔS⋅∇E + dot(Δμ, ewald.A[1, 1, 1, i, i], Δμ) / 2
end

"""
    modify_exchange_with_truncated_dipole_dipole!(sys::System, cutoff, μ0_μB²)

Like [`enable_dipole_dipole!`](@ref), the purpose of this function is to
introduce long-range dipole-dipole interactions between magnetic moments.
Whereas `enable_dipole_dipole!` employs Ewald summation, this function instead
employs real-space pair couplings with truncation at the specified `cutoff`
distance. The implicit demagnetization factor is 1/3, as appropriate for a
spherical sample in vacuum. If the cutoff is relatively small, then this
function may be faster than `enable_dipole_dipole!`.
"""
function modify_exchange_with_truncated_dipole_dipole!(sys::System{N}, cutoff, μ0_μB²=nothing) where N
    if isnothing(μ0_μB²)
        @warn "Deprecated syntax! Consider `modify_exchange_with_truncated_dipole_dipole!(sys, cutoff, units.vacuum_permeability)` where `units = Units(:meV, :angstrom)`."
        μ0_μB² = Units(:meV, :angstrom).vacuum_permeability
    end

    if !isnothing(sys.origin)
        modify_exchange_with_truncated_dipole_dipole!(sys.origin, cutoff, μ0_μB²)
        transfer_params_from_origin!(sys)
        return
    end

    # To support inhomogeneous systems, we would need a code path that modifies
    # the interactions on each site. See previous implementation in
    # https://github.com/SunnySuite/Sunny.jl/pull/416).
    is_homogeneous(sys) || error("System must be homogeneous")

    pairs = PairCoupling[]
    for bond in reference_bonds(sys.crystal, cutoff)
        for i in 1:natoms(sys.crystal)
            for bond′ in all_symmetry_related_bonds_for_atom(sys.crystal, i, bond)
                (; j) = bond′
                r = global_displacement(sys.crystal, bond′)
                iszero(r) && continue
                r̂ = normalize(r)
                bilin = (μ0_μB²/4π) * sys.gs[i]' * ((I - 3r̂⊗r̂) / norm(r)^3) * sys.gs[j]
                pc = PairCoupling(bond′, 0.0, Mat3(bilin), 0.0, zero(TensorDecomposition))
                push!(pairs, pc)
            end
        end
    end

    # Add to model params and repopulate couplings
    replace_model_param!(sys, :TruncatedDipoleDipole => 1.0; pairs)
    repopulate_couplings_from_params!(sys)
end
