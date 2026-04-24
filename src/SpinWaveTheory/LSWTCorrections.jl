"""
    energy_per_site_lswt_correction(swt::SpinWaveTheory; opts...)

Computes a perturbative correction to the classical energy per site assuming
large ``s`` (spin magnitude in dipole mode) or large ``λ`` (representation label
in SU(N) mode). The correction is a sum of two terms. The first is a uniform
(``𝐪 = 0``) correction. The second integrates over the zero-point energy for
all spin-wave modes, i.e., 1/2 ∑ₙ ∫d³q ω(q, n), where q belongs to the first
magnetic Brillouin zone and n is the band index. The correction appears at
sub-leading order in ``s`` or ``λ``. For instance, if the classical energy is
``J s^2``, the correction appears at order ``J s``.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration. See the
[HCubature](https://github.com/JuliaMath/HCubature.jl) documentation for
details.
"""
function energy_per_site_lswt_correction(swt::SpinWaveTheory; opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    Natoms = natoms(sys.crystal)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    # The uniform correction to the classical energy (trace of the (1,1)-block
    # of the spin-wave Hamiltonian)
    dynamical_matrix!(H, swt, zero(Vec3))
    δE₁ = -real(tr(view(H, 1:L, 1:L))) / 2Natoms

    # Integrate zero-point energy over the first Brillouin zone 𝐪 ∈ [0, 1]³ for
    # magnetic cell in reshaped RLU
    δE₂ = hcubature((0,0,0), (1,1,1); opts...) do q_reshaped
        dynamical_matrix!(H, swt, q_reshaped)
        ωs = bogoliubov!(V, H)
        return sum(view(ωs, 1:L)) / 2Natoms
    end

    # Error bars in δE₂[2] are discarded
    return δE₁ + δE₂[1]
end

# Calculates the magnetization reduction for :SUN mode for all atoms
function magnetization_lswt_correction_sun(swt::SpinWaveTheory; opts...)
    (; sys, data) = swt

    N = sys.Ns[1]
    Natoms = natoms(sys.crystal)
    L = (N - 1) * Natoms

    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    # Construct angular momentum operators O = n⋅S aligned with quantization
    # axis.
    S = spin_matrices_of_dim(; N)
    O = zeros(ComplexF64, N, N, Natoms)
    for i in 1:Natoms
        n = normalize(swt.sys.dipoles[i])
        U = data.local_unitaries[i]
        O[:, :, i] += U' * (n' * S) * U
        @assert O[N, N, i] ≈ norm(swt.sys.dipoles[i])
    end

    δS = hcubature((0,0,0), (1,1,1); opts...) do q
        swt_hamiltonian_SUN!(H, swt, q)
        bogoliubov!(V, H)
        ret = zeros(Natoms)
        for band in L+1:2L
            v = reshape(view(V, :, band), N-1, Natoms, 2)
            for i in 1:Natoms, α in 1:N-1, β in 1:N-1
                ret[i] -= real((O[N, N, i]*δ(α, β) - O[α, β, i]) * conj(v[α, i, 1]) * v[β, i, 1])
            end
        end
        return SVector{Natoms}(ret)
    end

    # Error bars in δS[2] are discarded
    return δS[1]
end

# Calculates the magnetization reduction for :dipole mode for every site
function magnetization_lswt_correction_dipole(swt::SpinWaveTheory; opts...)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    δS = hcubature((0,0,0), (1,1,1); opts...) do q
        swt_hamiltonian_dipole!(H, swt, Vec3(q))
        bogoliubov!(V, H)
        return SVector{L}(-norm2(view(V, L+i, 1:L)) for i in 1:L)
    end

    # Error bars in δS[2] are discarded
    return δS[1]
end

"""
    magnetization_lswt_correction(swt::SpinWaveTheory; opts...)

Calculates the reduction in the classical dipole magnitude for all atoms in the
magnetic cell. In the case of `:dipole` and `:dipole_uncorrected` mode, the
classical dipole magnitude is constrained to spin-`s`. While in `:SUN` mode, the
classical dipole magnitude can be smaller than `s` due to anisotropic
interactions.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration. See the HCubature package documentation
for details.
"""
function magnetization_lswt_correction(swt::SpinWaveTheory; opts...)
    any(in(keys(opts)), (:rtol, :atol, :maxevals)) || error("Must specify one of `rtol`, `atol`, or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    if sys.mode == :SUN
        δS = magnetization_lswt_correction_sun(swt; opts...)
    else
        @assert sys.mode in (:dipole, :dipole_uncorrected)
        δS = magnetization_lswt_correction_dipole(swt; opts...)
    end
    return δS
end


"""
    intensities_free_two_magnon(swt::SpinWaveTheory, q, energies, η::Float64; opts...)

Calculates dynamical pair correlation intensities for a provided ``𝐪``-point in
reciprocal space.
"""
function intensities_free_two_magnon(swt::SpinWaveTheory, q, energies, η::Float64; opts...)
    kernel = lorentzian(; fwhm=2η)

    (; sys, data, measure) = swt
    cryst = orig_crystal(sys)
    q_global = cryst.recipvecs * q
    q_reshaped = to_reshaped_rlu(sys, q)

    # Number of atoms in magnetic cell
    Nm = length(sys.dipoles)
    # Number of chemical cells in magnetic cell
    Ncells = Nm / natoms(cryst)
    # Dimension of Hilbert space
    N = sys.Ns[1]
    # Number of quasiparticle modes
    L = nbands(swt)

    # Get some integer numbers for later use
    num_energies = length(energies)
    num_obs = num_observables(measure)
    num_corrs = num_correlations(measure)

    # Preallocation
    # H1, V1 for q+k
    H1 = zeros(ComplexF64, 2L, 2L)
    V1 = zeros(ComplexF64, 2L, 2L)
    # H2, V2 for -k
    H2 = zeros(ComplexF64, 2L, 2L)
    V2 = zeros(ComplexF64, 2L, 2L)

    Avec_pref = zeros(ComplexF64, num_obs, Nm)
    Avec = zeros(ComplexF64, num_obs, L, L)
    corrbuf = zeros(ComplexF64, num_corrs, num_energies)
    resbuf = zeros(num_energies)


    for i = 1:Nm, μ in 1:num_obs
        r_global = global_position(sys, (1,1,1,i))
        ff = get_swt_formfactor(measure, μ, i)
        Avec_pref[μ, i] = exp(-1im * dot(q_global, r_global))
        Avec_pref[μ, i] *= compute_form_factor(ff, norm2(q_global))
    end

    ints = hcubature((0,0,0), (1,1,1); opts...) do k_reshaped
        qpk_reshaped = q_reshaped + k_reshaped
        if sys.mode == :SUN
            swt_hamiltonian_SUN!(H1, swt, qpk_reshaped)
            swt_hamiltonian_SUN!(H2, swt, -k_reshaped)
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            swt_hamiltonian_dipole!(H1, swt, qpk_reshaped)
            swt_hamiltonian_dipole!(H2, swt, -k_reshaped)
        end

        disp1 = bogoliubov!(V1, H1)
        disp2 = bogoliubov!(V2, H2)

        # Fill the buffers with zeros
        Avec .= 0.0
        resbuf .= 0.0
        corrbuf .= 0.0

        if sys.mode == :SUN
            for band1 = 1:L
                v1 = reshape(view(V1, :, band1), N-1, Nm, 2)
                for band2 = 1:L
                    v2 = reshape(view(V2, :, band2), N-1, Nm, 2)
                    for i = 1:Nm
                        for μ = 1:num_obs
                            O = data.observables_localized[μ, i]
                            for α = 1:N-1
                                for β = 1:N-1
                                    Avec[μ, band1, band2] += Avec_pref[μ, i] * (O[α, β] - δ(α, β) * O[N, N]) * (v1[α, i, 2]*v2[β, i, 1] + v1[β, i, 1]*v2[α, i, 2])
                                end
                            end
                        end
                    end
                end
            end
        else
            @assert sys.mode in (:dipole, :dipole_large_S)
            for band1 = 1:L
                v1 = reshape(view(V1, :, band1), Nm, 2)
                for band2 = 1:L
                    v2 = reshape(view(V2, :, band2), Nm, 2)
                    for i = 1:Nm
                        for μ = 1:num_obs
                            O = data.observables_localized[μ, i]
                            Avec[μ, band1, band2] += Avec_pref[μ, i] * O[3] * (v1[i, 2]*v2[i, 1] + v1[i, 1]*v2[i, 2])
                        end
                    end
                end
            end
        end

        for (ie, energy) in enumerate(energies)
            for (i, (α, β)) in enumerate(measure.corr_pairs)
                for band1 in 1:L, band2 in 1:L
                    corrbuf[i, ie] += Avec[α, band1, band2] * conj(Avec[β, band1, band2]) * kernel(disp1[band1]+disp2[band2], energy) / Ncells
                end
            end

            resbuf[ie] += measure.combiner(q_global, corrbuf[:, ie])
        end

        return SVector{num_energies}(resbuf)
    end

    return Vector(ints[1])
end
