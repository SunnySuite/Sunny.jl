# Bogoliubov transformation that diagonalizes a quadratic bosonic Hamiltonian,
# allowing for anomalous terms. The general procedure derives from Colpa,
# Physica A, 93A, 327-353 (1978).
function bogoliubov!(V::Matrix{ComplexF64}, H::Matrix{ComplexF64})
    L = div(size(H, 1), 2)
    @assert size(V) == size(H) == (2L, 2L)

    # Initialize V to the para-unitary identity Ĩ = diagm([ones(L), -ones(L)])
    V .= 0
    for i in 1:L
        V[i, i] = 1
        V[i+L, i+L] = -1
    end

    # Solve generalized eigenvalue problem, Ĩ t = λ H t, for columns t of V.
    # Eigenvalues are sorted such that positive values appear first, and are
    # otherwise ascending in absolute value.
    sortby(x) = (-sign(x), abs(x))
    λ, V0 = eigen!(Hermitian(V), Hermitian(H); sortby)

    # Note that V0 and V refer to the same data.
    @assert V0 === V

    # Normalize columns of V so that para-unitarity holds, V† Ĩ V = Ĩ.
    for j in axes(V, 2)
        c = 1 / sqrt(abs(λ[j]))
        view(V, :, j) .*= c
    end

    # Inverse of λ are eigenvalues of Ĩ H. A factor of 2 yields the physical
    # quasiparticle energies.
    energies = λ        # reuse storage
    @. energies = 2 / λ

    # The first L elements are positive, while the next L energies are negative.
    # Their absolute values are excitation energies for the wavevectors q and
    # -q, respectively.
    @assert all(>(0), view(energies, 1:L)) && all(<(0), view(energies, L+1:2L))
    
    # Disable tests below for speed. Note that the data in H has been
    # overwritten by eigen!, so H0 should refer to an original copy of H.
    #=
    Ĩ = Diagonal([ones(L); -ones(L)])
    @assert V' * Ĩ * V ≈ Ĩ
    @assert diag(V' * H0 * V) ≈ Ĩ * energies / 2
    # If H(q) = H(-q) (reflection symmetry), eigenvalues come in pairs
    if H0[1:L, 1:L] ≈ H0[L+1:2L, L+1:2L]
        @assert energies[1:L] ≈ -energies[L+1:2L]
    end
    =#

    return energies
end


# Returns |1 + nB(ω)| where nB(ω) = 1 / (exp(βω) - 1) is the Bose function. See
# also `classical_to_quantum` which additionally "undoes" the classical
# Boltzmann distribution.
function thermal_prefactor(kT, ω)
    if iszero(kT)
        return ω >= 0 ? 1 : 0
    else
        @assert kT > 0
        return abs(1 / (1 - exp(-ω/kT)))
    end
end


function calculate_excitations!(V, H, swt::SpinWaveTheory, q)
    (; sys) = swt
    q_global = orig_crystal(sys).recipvecs * q
    q_reshaped = sys.crystal.recipvecs \ q_global

    if sys.mode == :SUN
        swt_hamiltonian_SUN!(H, swt, q_reshaped)
    else
        @assert sys.mode in (:dipole, :dipole_large_S)
        swt_hamiltonian_dipole!(H, swt, q_reshaped)
    end

    try
        return bogoliubov!(V, H)
    catch _
        error("Instability at wavevector q = $q")
    end
end

"""
    excitations(swt::SpinWaveTheory, q)

Given a wavevector `q`, returns a pair `(energies, V)`. The first half of
`energies` are the quasi-particle excitation energies for `q`. The second half
are the negated quasiparticle energies for `-q`. The columns of `V`, to be
contracted with the Holstein-Primakoff bosons ``[𝐛_𝐪, 𝐛_{-𝐪}^†]``, are the
corresponding eigenvectors of the quadratic spin wave Hamiltonian.
"""
function excitations(swt::SpinWaveTheory, q)
    L = nbands(swt)
    V = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    energies = calculate_excitations!(V, H, swt, q)
    return (energies, V)
end

"""
    dispersion(swt::SpinWaveTheory, qpts)

Given a list of wavevectors `qpts` in reciprocal lattice units (RLU), returns
excitation energies for each band. The return value `ret` is 2D array, and
should be indexed as `ret[band_index, q_index]`.
"""
function dispersion(swt::SpinWaveTheory, qpts)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    disp = [view(excitations(swt, q)[1], 1:L) for q in qpts.qs]
    return reduce(hcat, disp)
end

"""
    intensities_bands(swt::SpinWaveTheory, qpts; formfactors=nothing)

Calculate spin wave excitation bands for a set of q-points in reciprocal space.
"""
function intensities_bands(swt::SpinWaveTheory, qpts; formfactors=nothing)
    (; sys, measure) = swt
    isempty(measure.observables) && error("No observables! Construct SpinWaveTheory with a `measure` argument.")

    qpts = convert(AbstractQPoints, qpts)
    cryst = orig_crystal(sys)

    # Number of atoms in magnetic cell
    @assert sys.latsize == (1,1,1)
    Na = length(eachsite(sys))
    # Number of chemical cells in magnetic cell
    Ncells = Na / natoms(cryst)
    # Number of quasiparticle modes
    L = nbands(swt)
    # Number of wavevectors
    Nq = length(qpts.qs)

    # Preallocation
    V = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    Avec_pref = zeros(ComplexF64, Na)
    disp = zeros(Float64, L, Nq)
    intensity = zeros(eltype(measure), L, Nq)

    # Temporary storage for pair correlations
    Ncorr = length(measure.corr_pairs)
    corrbuf = zeros(ComplexF64, Ncorr)

    Nobs = size(measure.observables, 1)

    # Expand formfactors for symmetry classes to formfactors for all atoms in
    # crystal
    ff_atoms = propagate_form_factors_to_atoms(formfactors, sys.crystal)

    for (iq, q) in enumerate(qpts.qs)
        q_global = cryst.recipvecs * q
        view(disp, :, iq) .= view(calculate_excitations!(V, H, swt, q), 1:L)

        for i in 1:Na
            r_global = global_position(sys, (1,1,1,i))
            Avec_pref[i] = exp(- im * dot(q_global, r_global))
            Avec_pref[i] *= compute_form_factor(ff_atoms[i], norm2(q_global))
        end

        Avec = zeros(ComplexF64, Nobs)

        # Fill `intensity` array
        for band = 1:L
            fill!(Avec, 0)
            if sys.mode == :SUN
                data = swt.data::SWTDataSUN
                N = sys.Ns[1]
                v = reshape(view(V, :, band), N-1, Na, 2)
                for i in 1:Na, μ in 1:Nobs
                    O = data.observables_localized[μ, i]
                    for α in 1:N-1
                        Avec[μ] += Avec_pref[i] * (O[α, N] * v[α, i, 2] + O[N, α] * v[α, i, 1])
                    end
                end
            else
                @assert sys.mode in (:dipole, :dipole_large_S)
                data = swt.data::SWTDataDipole
                v = reshape(view(V, :, band), Na, 2)
                for i in 1:Na, μ in 1:Nobs
                    O = data.observables_localized[μ, i]
                    # This is the Avec of the two transverse and one
                    # longitudinal directions in the local frame. (In the
                    # local frame, z is longitudinal, and we are computing
                    # the transverse part only, so the last entry is zero)
                    displacement_local_frame = SA[v[i, 2] + v[i, 1], im * (v[i, 2] - v[i, 1]), 0.0]
                    Avec[μ] += Avec_pref[i] * (data.sqrtS[i]/sqrt(2)) * (O' * displacement_local_frame)[1]
                end
            end

            map!(corrbuf, measure.corr_pairs) do (α, β)
                Avec[α] * conj(Avec[β]) / Ncells
            end
            intensity[band, iq] = measure.combiner(q_global, corrbuf)
        end
    end

    return BandIntensities(cryst, qpts, disp, intensity)
end

"""
    intensities!(data, swt::SpinWaveTheory, qpts; energies, kernel, formfactors=nothing)
    intensities!(data, swt::SpiralSpinWaveTheory, qpts; energies, kernel, formfactors=nothing)
    intensities!(data, swt::SampledCorrelations, qpts; energies=nothing, kernel=nothing, formfactors=nothing)

Like [`intensities`](@ref), but makes use of storage space `data` to avoid
allocation costs.
"""
function intensities!(data, swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, formfactors=nothing)
    @assert size(data) == (length(energies), size(bands.data, 2))
    bands = intensities_bands(swt, qpts; formfactors)
    @assert eltype(bands) == eltype(data)
    broaden!(data, bands; energies, kernel)
    return BroadenedIntensities(bands.crystal, bands.qpts, collect(energies), data)
end

"""
    intensities(swt::SpinWaveTheory, qpts; energies, kernel, formfactors=nothing)
    intensities(swt::SpiralSpinWaveTheory, qpts; energies, kernel, formfactors=nothing)
    intensities(swt::SampledCorrelations, qpts; energies=nothing, kernel=nothing, formfactors=nothing)

Calculate spin wave intensities for a set of q-points in reciprocal space. A
broadening `kernel` is required for spin wave theory, but optional for `SampledCorrelations`.
"""
function intensities(swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, formfactors=nothing)
    return broaden(intensities_bands(swt, qpts; formfactors); energies, kernel)
end
