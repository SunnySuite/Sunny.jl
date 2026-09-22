# Static, frequency-independent 1/s corrections: the ground-state energy and the
# ordered moments. Both are integrals of the harmonic zero-point fluctuations over
# the Brillouin zone, which is what makes them cheap compared to the dynamical
# corrections of CorrectedIntensities.jl. Conventions are collected in
# Corrections.jl.

"""
    energy_per_site_lswt_correction(swt::SpinWaveTheory; rtol=nothing, maxevals=nothing)

Correction to the classical energy per site at relative order ``1/s``, where
``s`` is the spin magnitude in dipole mode, or ``1/λ`` for the representation
label ``λ`` in SU(N) mode. If the classical energy is ``J s²``, this correction
appears at order ``J s``.

It is the zero-point energy of the harmonic magnons,
``(1/2) Σ_n ∫d³q ω(𝐪, n)`` over the first magnetic Brillouin zone, less the
uniform ``𝐪 = 0`` term that the Holstein-Primakoff normal ordering leaves behind,
together with the constant that [`anisotropy_correction`](@ref) generates from an
onsite coupling. The last of these vanishes identically in `:dipole` mode, where
`rcs_factors` makes the classical energy exact.

Not included are the corrections of [`tadpole_correction`](@ref) and
[`hartree_fock_correction`](@ref), which are smaller by a further power of ``1/s``
and so belong to the next order. To include them, add their `δE` fields to this
result.

The Brillouin-zone integral is performed by adaptive cubature, controlled by at
least one of `rtol` (a relative accuracy target) or `maxevals` (a budget of
integrand evaluations).
"""
function energy_per_site_lswt_correction(swt::SpinWaveTheory; rtol=nothing, maxevals=nothing)
    isnothing(rtol) && isnothing(maxevals) && error("Must specify `rtol` or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    # Normalize per physical site
    Nsites = nsites(uncontracted_system(sys))
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    # The uniform correction to the classical energy (trace of the (1,1)-block
    # of the spin-wave Hamiltonian)
    dynamical_matrix!(H, swt, zero(Vec3))
    δE₁ = -real(tr(view(H, 1:L, 1:L))) / 2Nsites

    # Integrate zero-point energy over the first Brillouin zone 𝐪 ∈ [0, 1]³ for
    # magnetic cell in reshaped RLU. Error bars are discarded.
    (δE₂, _) = hcubature((0, 0, 0), (1, 1, 1); rtol=@something(rtol, 0),
                         maxevals=@something(maxevals, typemax(Int))) do q_reshaped
        dynamical_matrix!(H, swt, q_reshaped)
        ωs = bogoliubov!(V, H)
        return sum(view(ωs, 1:L)) / 2Nsites
    end

    # The Stevens machinery behind this term is specific to dipole mode, where the
    # 1/s expansion of an onsite coupling leaves a constant; in :SUN mode the
    # coupling enters the boson Hamiltonian exactly and there is nothing to add.
    δE₃ = sys.mode == :SUN ? 0.0 : anisotropy_correction(swt).δE

    return δE₁ + δE₂ + δE₃
end

# Reduction in the magnitude of each classical dipole, for :SUN mode
function magnetization_lswt_correction_sun(swt::SpinWaveTheory; rtol, maxevals)
    (; sys, data) = swt

    # This correction measures the reduction of the classical dipole moment along
    # its own direction. It is undefined for an entangled system, whose ordered
    # object is a unit's product-space state (e.g. a dimer singlet has zero net
    # dipole), not a spin-(N-1)/2 dipole.
    is_entangled(sys) && error("Magnetization correction is not supported for entangled units.")

    N = sys.Ns[1]
    Natoms = natoms(sys.crystal)
    L = (N - 1) * Natoms

    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    # Construct angular momentum operators O = n⋅S aligned with quantization
    # axis, where S are the bare spin matrices and n = normalize(dipoles[i]).
    S = SVector{3}(spin_matrices_of_dim(; N))
    O = zeros(ComplexF64, N, N, Natoms)
    for i in 1:Natoms
        n = normalize(swt.sys.dipoles[i])
        U = data.local_unitaries[i]
        O[:, :, i] += U' * (n' * S) * U
        @assert O[N, N, i] ≈ norm(swt.sys.dipoles[i])
    end

    (δS, _) = hcubature((0, 0, 0), (1, 1, 1); rtol, maxevals) do q
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

    return δS
end

# Reduction in the magnitude of each classical dipole, for :dipole mode
function magnetization_lswt_correction_dipole(swt::SpinWaveTheory; rtol, maxevals)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    V = zeros(ComplexF64, 2L, 2L)

    (δS, _) = hcubature((0, 0, 0), (1, 1, 1); rtol, maxevals) do q
        swt_hamiltonian_dipole!(H, swt, Vec3(q))
        bogoliubov!(V, H)
        return SVector{L}(-norm2(view(V, L+i, 1:L)) for i in 1:L)
    end

    return δS
end

"""
    magnetization_lswt_correction(swt::SpinWaveTheory; rtol=nothing, maxevals=nothing)

Reduction in the magnitude of each classical dipole of the magnetic cell, caused
by zero-point fluctuations and appearing at relative order ``1/s``. Returns one
negative number per site. In `:dipole` and `:dipole_uncorrected` mode the
classical magnitude is constrained to be spin-`s`, whereas in `:SUN` mode it may
already be smaller than `s` because of anisotropic interactions.

This shortens each dipole without reorienting it. At the same order the ordered
structure also tilts, for which see [`corrected_dipoles`](@ref).

The Brillouin-zone integral is performed by adaptive cubature, controlled by at
least one of `rtol` (a relative accuracy target) or `maxevals` (a budget of
integrand evaluations).
"""
function magnetization_lswt_correction(swt::SpinWaveTheory; rtol=nothing, maxevals=nothing)
    isnothing(rtol) && isnothing(maxevals) && error("Must specify `rtol` or `maxevals` to control momentum-space integration.")
    opts = (; rtol=@something(rtol, 0), maxevals=@something(maxevals, typemax(Int)))

    (; sys) = swt
    if sys.mode == :SUN
        return magnetization_lswt_correction_sun(swt; opts...)
    else
        @assert sys.mode in (:dipole, :dipole_uncorrected)
        return magnetization_lswt_correction_dipole(swt; opts...)
    end
end

"""
    corrected_dipoles(swt::SpinWaveTheory; rtol=nothing, maxevals=nothing)

Classical dipoles of the magnetic cell, corrected at relative order ``1/s``.
Zero-point fluctuations act on the ordered structure in two independent ways at
this order, and both are applied here: each dipole is shortened by
[`magnetization_lswt_correction`](@ref), and the structure is tilted by the
zero-point pressure of [`tadpole_correction`](@ref), whose canonical example is
the change in canting angle of an antiferromagnet in an applied field.

The tilt vanishes for a collinear structure, and requires the cubic vertex, which
is available for a narrower class of models than the shortening is; see
[`tadpole_correction`](@ref). Where it is unavailable the dipoles are returned
shortened but untilted, which is still correct at this order for any structure the
tilt would not move.

The Brillouin-zone integrals are performed by adaptive cubature, controlled by at
least one of `rtol` (a relative accuracy target) or `maxevals` (a budget of
integrand evaluations).
"""
function corrected_dipoles(swt::SpinWaveTheory; rtol=nothing, maxevals=nothing)
    (; sys) = swt
    δS = magnetization_lswt_correction(swt; rtol, maxevals)
    # `tadpole_correction` already returns dipoles of the classical magnitude,
    # rotated but not shortened, so the two corrections compose by scaling.
    dipoles = if isnothing(corrections_unsupported_reason(swt))
        tadpole_correction(swt; rtol, maxevals).dipoles
    else
        # One dipole per element of δS, which is a site of the magnetic cell; note
        # that in :SUN mode this is not `nbands(swt)`.
        [sys.dipoles[1, 1, 1, i] for i in eachindex(δS)]
    end
    return map((d, δ) -> (1 + δ/norm(d)) * d, dipoles, δS)
end
