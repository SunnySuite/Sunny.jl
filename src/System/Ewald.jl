
function Ewald(sys::System{N}, μ0_μB², demag) where N
    (; crystal, dims) = sys

    # Collect pieces to build the Ewald interaction tensor A(q) cheaply. This
    # step can be a bottleneck for spin wave theory.
    cache = EwaldTensorCache(crystal, dims, μ0_μB², demag)

    # Build and store the interaction tensor A(q=0). The imaginary part cancels
    # in the symmetric sum over ±k, so keep only the real part.
    A = real.(ewald_interaction_tensor(cache, zero(Vec3)))
    Ar = reshape(reinterpret(Float64, A), 3, 3, size(A)...) # dims: [α,β,cell,i,j]
    FA = FFTW.rfft(Ar, 3:5) # FFT on cell indices

    # Scratch space for calculating interactions in Fourier space
    na = natoms(crystal)
    μ = zeros(Vec3, dims..., na)
    ϕ = zeros(Vec3, dims..., na)
    sz_rft = size(FA)[3:5]  # First FT dimension (dimension 3) will be ~ halved
    Fμ = zeros(ComplexF64, 3, sz_rft..., na)
    Fϕ = zeros(ComplexF64, 3, sz_rft..., na)

    mock_spins = zeros(3, dims..., na)
    plan     = FFTW.plan_rfft(mock_spins, 2:4; flags=FFTW.MEASURE)
    ift_plan = FFTW.plan_irfft(Fμ, dims[1], 2:4; flags=FFTW.MEASURE)

    return Ewald(cache, A, μ, ϕ, FA, Fμ, Fϕ, plan, ift_plan)
end

# Ideally, this would clone all mutable state within Ewald. Note that `A`, `FA`
# are immutable data. A blocker is that FFTW plans cannot currently be copied,
# and it is not 100% clear whether they can be treated as immutable. For
# example, they cache inverse plans, which may possibly lead to data races in a
# multithreaded context. See discussion at
# https://github.com/JuliaMath/FFTW.jl/issues/261.
function clone_ewald(ewald::Ewald)
    error("Not supported")
    (; cache, A, μ, ϕ, FA, Fμ, Fϕ, plan, ift_plan) = ewald
    return Ewald(cache, A, copy(μ), copy(ϕ), FA, copy(Fμ), copy(Fϕ), copy(plan), copy(ift_plan))
end

# Outer product of 3-vectors as static Mat3
(⊗)(a::Vec3, b::Vec3) = a * b'

# Bare (single-image) point dipole-dipole coupling tensor for a displacement
# `r`, i.e. the r → 0 limit of the Ewald kernel, without the μ0_μB² prefactor.
function point_dipole_tensor(r::Vec3)
    r̂ = normalize(r)
    return (I - 3(r̂⊗r̂)) / (4π * norm(r)^3)
end

@inline cell_offset(cell) = Vec3(cell[1]-1, cell[2]-1, cell[3]-1)

function EwaldTensorCache(cryst::Crystal, dims::NTuple{3,Int}, μ0_μB², demag::Mat3)
    na = natoms(cryst)

    # Superlattice vectors and reciprocals for the full system volume
    sys_size = diagm(Vec3(dims))
    latvecs = cryst.latvecs * sys_size
    recipvecs = cryst.recipvecs / sys_size

    # Precalculate constants
    I₃ = Mat3(I)
    V = abs(det(latvecs))
    L = cbrt(V)

    # Selected to roughly balance the real and Fourier space costs
    σ = L/3
    σ² = σ*σ

    # Gives about 13 digits of accuracy (equivalent to c0=6 in Ewalder.jl)
    rmax = 6√2 * σ
    kmax = 6√2 / σ
    nmax = map(eachcol(latvecs), eachcol(recipvecs)) do a, b
        round(Int, rmax / (a⋅normalize(b)) + 1e-6) + 1
    end
    mmax = map(eachcol(latvecs), eachcol(recipvecs)) do a, b
        round(Int, kmax / (b⋅normalize(a)) + 1e-6)
    end

    # Precalculate q-independent part of real-space terms
    real_terms = Array{Vector{Tuple{NTuple{3, Int}, Mat3}}, 5}(undef, dims..., na, na)
    for cell in CartesianIndices(dims), j in 1:na, i in 1:na
        Δr = cryst.latvecs * (cell_offset(cell) + cryst.positions[j] - cryst.positions[i])

        terms = Tuple{NTuple{3, Int}, Mat3}[]
        for n1 = -nmax[1]:nmax[1], n2 = -nmax[2]:nmax[2], n3 = -nmax[3]:nmax[3]
            rvec = Δr + latvecs * Vec3(n1, n2, n3)
            r² = rvec⋅rvec
            0 < r² <= rmax*rmax || continue
            r = √r²
            r³ = r²*r
            rhat = rvec/r
            erfc0 = erfc(r/(√2*σ))
            gauss0 = √(2/π) * (r/σ) * exp(-r²/2σ²)
            A = (1/4π) * ((I₃/r³) * (erfc0 + gauss0) - (3(rhat⊗rhat)/r³) * (erfc0 + (1+r²/3σ²) * gauss0))
            push!(terms, ((n1, n2, n3), A))
        end
        real_terms[cell, i, j] = terms
    end

    return EwaldTensorCache(dims, cryst, μ0_μB², demag, σ², Tuple(mmax), kmax^2, Tuple(nmax), real_terms)
end

# Calculate the Ewald interaction tensor A[cell, i, j] at wavevector
# `q_reshaped`, from the q-independent pieces held in `cache`. For q_reshaped =
# 0, this yields the usual Ewald energy, E = μᵢ Aᵢⱼ μⱼ / 2. Nonzero q_reshaped
# is useful in spin wave theory. Physically, this amounts to a modification of
# the periodic boundary conditions, such that μ(q) can be incommensurate with
# the magnetic cell. In all cases, the energy is E = μᵢ(-q) Aᵢⱼ(-q) μⱼ(q) / 2 in
# Fourier space, where q should be interpreted as a Fourier transform of the
# cell offset.
function ewald_interaction_tensor(cache::EwaldTensorCache, q_reshaped::Vec3)
    (; dims, cryst, μ0_μB², demag, σ², mmax, kmax², nmax, real_terms) = cache
    (; positions) = cryst
    na = natoms(cryst)
    recipvecs = cryst.recipvecs / diagm(Vec3(dims))
    V = abs(det(cryst.latvecs)) * prod(dims)
    self_energy = -Mat3(I) / (3(2π)^(3/2) * σ²^(3/2))
    q0 = q_reshaped - round.(q_reshaped)
    A = zeros(CMat3, dims..., na, na)

    # Both sections below build the same set of (cell, i, j) entries: every
    # nonzero-offset block in full, but only the upper triangle i ≤ j of the
    # zero-offset block. That block is Hermitian at any q (its phase factorizes as
    # conj(siteφᵢ)⋅siteφⱼ, independent of the k-grid), so a final pass mirrors its
    # lower triangle. `imax` is the upper bound on i. Every other (nonzero-offset)
    # block is built in full, so the result is correct for any q, not just the
    # q = 0 case where dims ≠ 1 today.
    is_zero_offset(cell) = all(isone, cell.I)
    imax(cell, j) = is_zero_offset(cell) ? j : na

    # The k = 0 Fourier mode is singular and omitted from the reciprocal sum; in
    # its place it contributes the demag surface term Eₛ = μ₀ M⋅N M / 2V, added to
    # every pair below. This term is present only when q0 = 0 (the sole mode with
    # k = recipvecs(m + q0) = 0). The factor tensor N (`demag`) has trace 1 in
    # vacuum. See S. W. DeLeeuw et al., Proc. R. Soc. Lond. A 373, 27-56 (1980)
    # and Ballenegger, J. Chem. Phys. 140, 161102 (2014). The same q0_is_zero flag
    # drives both the demag term and the mode skip below, so they stay consistent.
    q0_is_zero = iszero(q0)
    demag_term = q0_is_zero ? demag / V : zero(Mat3)

    #####################################################
    ## Fourier space part
    # With k = recipvecs (m + q0), the phase cis(-k⋅Δr) of a displacement
    # Δr = latvecs (off + rⱼ - rᵢ) factorizes over the three axes, since
    # recipvecsᵀ⋅latvecs = 2π I. So tabulate the 1D per-axis phases e^{-2πi(m+q0)x}
    # (for site fractions x = rᵢ/dims and cell fractions x = off/dims) and take
    # products; each surviving mode then costs only `na` products, no
    # transcendentals. The mode loop runs outermost so each mode's site phases
    # `siteφ` and tensor `Aₖ` are computed once and reused across cells and pairs.
    siteφ_1d = ntuple(3) do a
        [cis(-2π*(m+q0[a]) * positions[i][a]/dims[a]) for i in 1:na, m in centered(-mmax[a]:mmax[a])]
    end
    cellφ_1d = ntuple(3) do a
        [cis(-2π*(m+q0[a]) * (c-1)/dims[a]) for c in 1:dims[a], m in centered(-mmax[a]:mmax[a])]
    end
    # On the zero-offset diagonal (i = j) the phase is identically 1 (cellφ = 1 and
    # conj(siteφᵢ)·siteφᵢ = 1), so its Fourier value is the same mode sum ∑ₖ Aₖ for
    # every site — independent of cell and pair. Accumulate it once as `diagF` (a
    # single register add per mode, vs `na` scattered writes into A), skip those
    # entries in the inner loop, and write `diagF` to the diagonal afterward.
    diagF = zero(CMat3)
    siteφ = zeros(ComplexF64, na)
    for m1 = -mmax[1]:mmax[1], m2 = -mmax[2]:mmax[2], m3 = -mmax[3]:mmax[3]
        # The m = 0 mode at q0 = 0 is the singular k = 0 term (demag stands in for it).
        q0_is_zero && all(iszero, (m1, m2, m3)) && continue
        k = recipvecs * (Vec3(m1, m2, m3) + q0)
        k² = k⋅k
        k² <= kmax² || continue  # reciprocal-space cutoff
        Aₖ = ((1/V) * exp(-σ²*k²/2) / k²) * (k⊗k)  # real, symmetric, pair-independent
        diagF += Aₖ
        @inbounds for i in 1:na
            siteφ[i] = siteφ_1d[1][i, m1] * siteφ_1d[2][i, m2] * siteφ_1d[3][i, m3]
        end
        @inbounds for cell in CartesianIndices(dims)
            cellφ = cellφ_1d[1][cell[1], m1] * cellφ_1d[2][cell[2], m2] * cellφ_1d[3][cell[3], m3]
            for j in 1:na, i in 1:imax(cell, j)
                is_zero_offset(cell) && i == j && continue  # constant; written from diagF below
                A[cell, i, j] += (cellφ * conj(siteφ[i]) * siteφ[j]) * Aₖ
            end
        end
    end
    @inbounds for i in 1:na
        A[1, 1, 1, i, i] = diagF
    end

    #####################################################
    ## Real-space part (q-independent tensors cached in `real_terms`), added onto
    # the Fourier tensor already in `A`, then the demag, self-energy, and μ0_μB²
    # prefactor. For sites site1=(cell1, i) and site2=(cell2, j) offset by
    # (off = cell2-cell1), the pair-energy is (s1 ⋅ A[off, i, j] ⋅ s2); Julia arrays
    # start at one, so we index A using (cell = off .+ 1). The self-energy corrects
    # only the on-site diagonal. The phase cis(2π q⋅n) depends only on the shift n
    # and factorizes over axes, so tabulate the three per-axis phases e^{2πi qₐ nₐ}
    # (indexed directly by nₐ) and take their product for each term's stored shift n.
    real_phases = ntuple(3) do a
        [cis(2π * q_reshaped[a] * n) for n in centered(-nmax[a]:nmax[a])]
    end
    @inbounds for cell in CartesianIndices(dims), j in 1:na, i in 1:imax(cell, j)
        acc = A[cell, i, j] + demag_term
        for ((n1, n2, n3), Aⁿ) in real_terms[cell, i, j]
            acc += (real_phases[1][n1] * real_phases[2][n2] * real_phases[3][n3]) * Aⁿ
        end
        (is_zero_offset(cell) && i == j) && (acc += self_energy)
        A[cell, i, j] = μ0_μB² * acc
    end

    # Fill the lower triangle of the zero-offset block by Hermitian symmetry.
    @inbounds for j in 1:na, i in 1:j-1
        A[1, 1, 1, j, i] = A[1, 1, 1, i, j]'
    end
    return A
end

# The @nospecialize(sys) hint satisfies JET when Hilbert size N is not known
# statically.
function ewald_energy(@nospecialize(sys::System))
    (; μ, FA, Fμ, plan) = sys.ewald
    dims = size(sys.dipoles)[1:3]
    even_rft_size = dims[1] % 2 == 0

    E = 0.0
    @. μ = - sys.gs * sys.dipoles # i.e., magnetic_moments(sys)
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
# interactions. The @nospecialize(sys) hint satisfies JET when Hilbert size N is
# not known statically.
function accum_ewald_grad!(∇E, dipoles, @nospecialize(sys::System))
    (; gs, ewald) = sys
    (; μ, FA, Fμ, Fϕ, ϕ, plan, ift_plan) = ewald

    # Fourier transformed magnetic moments for the provided trial dipoles
    @. μ = - gs * dipoles
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
    off = mod.(to_cell(site2) .- to_cell(site1), dims)
    cell = CartesianIndex(off .+ (1,1,1))

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
    Δμ = - sys.gs[site] * ΔS
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
                bilin = μ0_μB² * sys.gs[i]' * point_dipole_tensor(r) * sys.gs[j]
                pc = PairCoupling(bond′, 0.0, Mat3(bilin), 0.0, zero(TensorDecomposition))
                push!(pairs, pc)
            end
        end
    end

    # Add to model params and repopulate couplings
    replace_model_param!(sys, :TruncatedDipoleDipole => 1.0; pairs)
    repopulate_couplings_from_params!(sys)
end
