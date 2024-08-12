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
    # Eigenvalues are sorted such that positive values appear first, and are in
    # ascending order.
    λ, V0 = eigen!(Hermitian(V), Hermitian(H); sortby = x -> -1/real(x))

    # Note that V0 and V refer to the same data.
    @assert V0 === V

    # Normalize columns of V so that para-unitarity holds, V† Ĩ V = Ĩ.
    for j in axes(V, 2)
        c = 1 / sqrt(abs(λ[j]))
        view(V, :, j) .*= c
    end

    # Disable test for speed
    #=
    Ĩ = Diagonal([ones(L); -ones(L)])
    @assert V' * Ĩ * V ≈ Ĩ
    =#

    # Verify that half the eigenvalues are positive and the other half are
    # negative. The positive eigenvalues are quasiparticle energies for the
    # wavevector q that defines the dynamical matrix H(q). The absolute value of
    # the negative eigenvalues would be quasiparticle energies for H(-q), which
    # we are not considering in the present context.
    @assert all(>(0), view(λ, 1:L)) && all(<(0), view(λ, L+1:2L))
    
    # Inverse of λ are eigenvalues of Ĩ H. We only care about the first L
    # eigenvalues, which are positive. A factor of 2 is needed to get the
    # physical quasiparticle energies. These will be in descending order.
    disp = resize!(λ, L)
    @. disp = 2 / disp

    # In the special case that H(q) = H(-q) (i.e., a magnetic ordering with
    # reflection symmetry), the eigenvalues come in pairs. Note that the data in
    # H has been overwritten by eigen!, so H0 should refer to an original copy
    # of H.
    #=
    @assert diag(V' * H0 * V) ≈ [disp/2; reverse(disp)/2]
    =#

    return disp
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

Given a single wavevector `q`, returns a pair `(energies, V)`. Elements of
`energies` are the quasi-particle excitation energies. The columns of `V`, to be
contracted with the Holstein-Primakoff bosons ``[𝐛^†, b]``, are the
corresponding eigenvectors of the quadratic spin wave Hamiltonian.
"""
function excitations(swt::SpinWaveTheory, q)
    L = nbands(swt)
    V = zeros(ComplexF64, 2L, 2L)
    H = zeros(ComplexF64, 2L, 2L)
    disp = calculate_excitations!(V, H, swt, q)
    return (disp, view(V, :, 1:L))
end

"""
    dispersion(swt::SpinWaveTheory, qpts)

Given a list of wavevectors `qpts` in reciprocal lattice units (RLU), returns
excitation energies for each band. The return value `ret` is 2D array, and
should be indexed as `ret[band_index, q_index]`.
"""
function dispersion(swt::SpinWaveTheory, qpts)
    qpts = convert(AbstractQPoints, qpts)
    disp = [excitations(swt, q)[1] for q in qpts.qs]
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
        disp[:, iq] .= calculate_excitations!(V, H, swt, q)

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
    intensities(swt::SpinWaveTheory, qpts; energies, kernel, formfactors=nothing)

Calculate spin wave intensities for a set of q-points in reciprocal space. TODO.
"""
function intensities(swt::SpinWaveTheory, qpts; energies, kernel::AbstractBroadening, formfactors=nothing)
    return broaden(intensities_bands(swt, qpts; formfactors), energies; kernel)
end

