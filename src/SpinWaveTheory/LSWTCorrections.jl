"""
    energy_per_site_lswt_correction(swt::SpinWaveTheory; opts...)

Computes the [𝒪(1/λ) or 𝒪(1/S)] correction to the classical energy **per
site** [𝒪(λ²) or 𝒪(S²)] given a [`SpinWaveTheory`](@ref). The correction
[𝒪(λ) or 𝒪(S)] includes a uniform term (For instance, if the classical energy
is αJS², the LSWT gives a correction like αJS) and the summation over the
zero-point energy for all spin-wave modes, i.e., 1/2 ∑ₙ ∫d³q ω(q, n), where q
belongs to the first magnetic Brillouin zone and n is the band index.

A keyword argument `rtol`, `atol`, or `maxevals` is required to control the
accuracy of momentum-space integration. See the HCubature package documentation
for details.
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
