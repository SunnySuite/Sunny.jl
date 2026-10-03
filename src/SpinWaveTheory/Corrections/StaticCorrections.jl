# Static, frequency-independent 1/s corrections, each a contraction in the
# vacuum of Vacuum.jl: the mean field of H₄, the tadpole, the corrections to the
# observables, and the energy and ordered moments that they and the zero-point
# fluctuations determine. Conventions are collected in Corrections.jl.

# Mean-field (Hartree-Fock) treatment of the four-boson term H₄, which shifts
# the magnon dispersion at O(1/s) relative to LSWT. Iterating it to
# self-consistency is the vacuum `MagnonVacuum(swt, terms2)` that reproduces its
# own mean field, which resums one class of higher-order diagrams and drops
# others of the same order. That is uncontrolled, and in particular it violates
# the Ward identity of a broken continuous symmetry, gapping a mode that should
# be gapless.

"""
    hartree_fock_correction(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)

Decouples the four-boson term of the Holstein-Primakoff expansion into a
mean-field correction to the quadratic (LSWT) Hamiltonian, with mean fields
taken in the LSWT ground state. The correction is smaller than the LSWT
Hamiltonian by a factor of order ``1/s``. Returns `(; terms2, δE)`, where
`terms2` can be passed to [`corrected_dispersion`](@ref) and `δE` is a
correction to the energy per site.

Brillouin-zone integrals are controlled by `tol`, a relative accuracy target
for adaptive cubature, and/or `maxevals`, its budget of integrand evaluations.
Alternatively `grid` replaces the cubature by an average over a uniform grid of
the given dimensions.
"""
function hartree_fock_correction(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)
    check_corrections_supported(swt)
    terms4 = quartic_monomials(swt)
    (; terms2, δE) = mean_field(terms4, contractions(MagnonVacuum(swt), terms4, BZQuadrature(; tol, maxevals, grid)))
    return (; terms2, δE = δE / nsites(uncontracted_system(swt.sys)))
end

# Mean-field decoupling of the quartic monomials with contractions `g`. One
# contraction leaves the quadratic form, two the constant. Since ⟨terms2⟩ =
# 2⟨H₄⟩, the constant is subtracted once, which is what makes the mean-field
# form reproduce ⟨H₄⟩ itself.
function mean_field(terms4, g::Contractions)
    terms2 = wick_reduce(terms4, g, nothing, Val{2}())
    δE = -wick_expectation(terms4, g, nothing; nfactors=2)
    # δE is real only once the mean fields are exact, so its imaginary part
    # measures the error of their momentum integrals and shrinks with `tol`. A
    # mistake in the Wick decoupling, which is what this checks for, would
    # instead appear at O(1).
    @assert abs(imag(δE)) < g.noise * max(abs(δE), 1)
    return (; terms2, δE = real(δE))
end

# Tadpole correction. Zero-point fluctuations shift the ordered structure away
# from the classical energy minimum, e.g. they change the canting angle of an
# antiferromagnet in a field.
#
# The shift is realized as a uniform displacement of the Holstein-Primakoff
# bosons, b_{i,𝐫} → b_{i,𝐫} + v_i, rather than as a rotation of the classical
# dipoles. The two are equivalent, since ⟨S⁺_i⟩ = σ_i v_i tilts the moment away
# from the local ẑ axis, but the displacement is far better behaved numerically:
# away from a classical energy minimum the linear term H₁ no longer vanishes,
# and the quadratic form H₂ ceases to be positive definite near 𝐪 = 0, so
# `bogoliubov!` cannot even be applied there. Working with v keeps every
# calculation at the classical minimum, where LSWT is well defined.
#
# Substituting b → b + v produces a term linear in the bosons from H₂, and
# contracting two of the three legs of H₃ produces another (the tadpole itself).
# Requiring the total to vanish fixes v, which is of order s^(-1/2). The
# accompanying correction to H₂ comes from H₃ with a single leg replaced by v,
# and is of order s⁰, i.e. smaller than H₂ by 1/s, matching the mean-field
# correction from H₄.

"""
    tadpole_correction(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)

Computes the shift of the ordered magnetic structure caused by zero-point
fluctuations, which appears at relative order ``1/s``. An example is the change
in canting angle of an antiferromagnet in an applied field. Returns `(; terms2,
δE, v)`, where `terms2` is a correction to the quadratic Hamiltonian that can be
passed to [`corrected_dispersion`](@ref), `δE` is a correction to the energy per
site, and `v` is the boson displacement that realizes the shift, which
[`observable_corrections`](@ref) and [`corrected_magnetic_moments`](@ref) need.
The latter reads the shifted structure off `v`, together with the zero-point
depletion, and converts to ``μ = -g 𝐒``.

The displaced structure is the one that minimizes the energy reported by
[`corrected_energy_per_site`](@ref), i.e. the classical energy together with the
zero-point energy of the magnons, and `δE` is the resulting gain. The correction
is of the same size as the mean-field correction of
[`hartree_fock_correction`](@ref), so for a noncollinear structure the two
should be applied together. Restoring an exact Goldstone mode requires, in
addition, the self-energy generated by the cubic vertex, which enters at the
same order.

Brillouin-zone integrals are controlled by `tol`, a relative accuracy target for
adaptive cubature, and/or `maxevals`, its budget of integrand evaluations.
Alternatively `grid` replaces the cubature by an average over a uniform grid of
the given dimensions.
"""
function tadpole_correction(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)
    quad = BZQuadrature(; tol, maxevals, grid)
    check_corrections_supported(swt)
    terms3 = cubic_monomials(swt)
    (; terms2, δE, w) = tadpole(MagnonVacuum(swt), terms3, contractions(MagnonVacuum(swt), terms3, quad))
    return (; terms2, δE = δE / nsites(uncontracted_system(swt.sys)), v = w[1:nbands(swt)])
end

# The tadpole about the vacuum with contractions `g`: the Nambu-packed
# displacement `w`, the correction `terms2` it makes to the quadratic
# Hamiltonian, and the energy gain `δE` of the magnetic cell.
function tadpole(vac::MagnonVacuum, terms3, g::Contractions)
    (; swt) = vac
    L = nbands(swt)
    (; noise) = g

    # The linear term at the classical minimum: one contraction of H₃, the
    # tadpole itself, plus the one-boson word of an onsite anisotropy. The
    # latter is proportional to the classical energy gradient, so like the
    # contractions it vanishes at a classical minimum but not once tadpole
    # relaxation has moved the structure off one. In mode :dipole it vanishes
    # identically, because `rcs_factors` leaves the classical energy exact. See
    # `anisotropy_coefficients`.
    ℓ = nambu_vector([wick_reduce(terms3, g, nothing, Val{1}()); anisotropy_monomials(swt, Val{1}())], L, noise)

    # Packing v into the Nambu vector w, with w[i] = v_i and w[L+i] = conj(v_i),
    # the displacement of H₂ = Σ_𝐪 (1/2) x†_𝐪 H_𝐪 x_𝐪, that of LSWT rather
    # than of the vacuum, contributes
    #
    #     (√N/2) (w† H₀ x₀ + x₀† H₀ w) = √N Σ_a (H₀ w)[ā] x₀[a],
    #
    # where the two halves coincide by the para-symmetry H₀[ā,b̄] =
    # conj(H₀[a,b]). Since Σ_𝐫 O_a(𝐫) = √N x₀[a], cancellation of the total
    # linear term requires (H₀ w)[ā] = -ℓ[a] for every a.
    #
    # H₀ is singular whenever the structure has a continuous degeneracy: a
    # Goldstone mode is a direction along which the structure may be rotated
    # freely, and by symmetry ℓ has no component along it. The pseudo-inverse
    # selects the minimum-norm solution, which is the one that leaves that
    # freedom unused. Its cutoff must also discard `swt.regularization`, which
    # would otherwise turn the Goldstone direction into a huge displacement.
    H = zeros(ComplexF64, 2L, 2L)
    dynamical_matrix!(H, swt, zero(Vec3))
    w = -pinv(H; rtol=1e-6) * [ℓ[nambu_conj(a, L)] for a in 1:2L]
    @assert norm(w[L+1:2L] - conj(w[1:L])) < noise * max(norm(w), 1)

    # Correction to H₂, from each cubic monomial with one leg displaced
    terms2 = wick_reduce(terms3, nothing, w, Val{2}())

    # Half the linear response, the usual energy gain of a displaced harmonic
    # system, since at the stationary point w† H₀ w = -Σ_a ℓ[a] w[a].
    δE = sum(ℓ .* w) / 2
    @assert abs(imag(δE)) < noise * max(abs(δE), 1)

    return (; terms2, δE = real(δE), w)
end

"""
    observable_corrections(swt::SpinWaveTheory; v=nothing, tol=nothing, maxevals=nothing, grid=nothing)

Correction of relative order ``1/s`` to the amplitude for a magnon to be created
by each observable. Two effects contribute at this order: the cubic term of the
Holstein-Primakoff expansion of the transverse spin components, and the tilt of
the ordered structure by zero-point fluctuations. The latter requires the boson
displacement `v` of [`tadpole_correction`](@ref), and is omitted if `v` is
`nothing`. Returns the coefficients `δc[a, μ]` of the one-boson operators,
labeled as in [`accum_observable_corrections!`](@ref), which is what applies
them.

Brillouin-zone integrals are controlled by `tol`, a relative accuracy target for
adaptive cubature, and/or `maxevals`, its budget of integrand evaluations.
Alternatively `grid` replaces the cubature by an average over a uniform grid of
the given dimensions.
"""
function observable_corrections(swt::SpinWaveTheory; v=nothing, tol=nothing, maxevals=nothing, grid=nothing)
    check_corrections_supported(swt)
    w = isnothing(v) ? nothing : [v; conj(v)]
    g = contractions(MagnonVacuum(swt), observable_cubic_monomials(swt), BZQuadrature(; tol, maxevals, grid))
    return observable_corrections(swt, w, g)
end

# The same with contractions `g` and the Nambu-packed displacement `w`. One
# contraction of the cubic word, and the tilt as one displaced leg of the
# longitudinal word. Every contracted pair acts on the same site as the
# surviving operator, so only onsite correlations are needed and no wavevector
# dependence survives.
function observable_corrections(swt::SpinWaveTheory, w, g::Contractions)
    L = nbands(swt)
    return stack(1:num_observables(swt.measure)) do μ
        nambu_vector([wick_reduce(observable_monomials(swt, μ, Val{3}()), g, w, Val{1}());
                      wick_reduce(observable_monomials(swt, μ, Val{2}()), g, w, Val{1}())], L, g.noise)
    end
end

"""
    corrected_energy_per_site(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)

Energy per site of the magnetic structure, corrected at relative order ``1/s``,
where ``s`` is the spin magnitude in dipole mode, or ``1/λ`` for the
representation label ``λ`` in SU(N) mode. If the classical energy is ``J s²``,
the correction appears at order ``J s``.

The correction is the zero-point energy of the harmonic magnons, ``(1/2) Σ_n
∫d³q ω(𝐪, n)`` over the first magnetic Brillouin zone, less the uniform ``𝐪 =
0`` term that the Holstein-Primakoff normal ordering leaves behind, together
with the constant that [`anisotropy_correction`](@ref) generates from an onsite
coupling. The last of these vanishes identically in `:dipole` mode, where
`rcs_factors` makes the classical energy exact.

Not included are the corrections of [`tadpole_correction`](@ref) and
[`hartree_fock_correction`](@ref), which are smaller by a further power of
``1/s`` and so belong to the next order. To include them, add their `δE` fields
to this result.

Brillouin-zone integrals are controlled by `tol`, a relative accuracy target for
adaptive cubature, and/or `maxevals`, its budget of integrand evaluations.
Alternatively `grid` replaces the cubature by an average over a uniform grid of
the given dimensions.
"""
function corrected_energy_per_site(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)
    quad = BZQuadrature(; tol, maxevals, grid)

    (; sys) = swt
    # Normalize per physical site
    Nsites = nsites(uncontracted_system(sys))
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)

    # The uniform correction to the classical energy (trace of the (1,1)-block
    # of the spin-wave Hamiltonian)
    dynamical_matrix!(H, swt, zero(Vec3))
    δE₁ = -real(tr(view(H, 1:L, 1:L))) / 2Nsites

    # Zero-point energy, averaged over the magnetic Brillouin zone
    δE₂ = bz_average(quad) do q_reshaped
        dynamical_matrix!(H, swt, q_reshaped)
        ωs = bogoliubov!(ws, H)
        return sum(view(ωs, 1:L)) / 2Nsites
    end

    # Vanishes in :SUN mode, where an onsite coupling enters the boson
    # Hamiltonian exactly and there is nothing to add.
    δE₃ = anisotropy_correction(swt).δE

    return swt.classical_energy + δE₁ + δE₂ + δE₃
end

"""
    boson_density(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)

Zero-point density of Holstein-Primakoff bosons, ``n_i = Σ_α ⟨b^†_{iα}
b_{iα}⟩``, summed over the flavors ``α`` carried by each site of the magnetic
cell. This is the small parameter that controls the ``1/s`` expansion: every
correction in this module is a power of ``n``, so ``n ≪ 1`` is the statement
that the expansion is converging, and ``n`` of order one means it is not.

Use this as the health check for a spin-wave calculation. In `:dipole` mode
there is one boson per site, and ``n_i`` is also the shortening of the classical
dipole: ``⟨S^z_i⟩ = s - n_i``. In `:SUN` mode there are ``N-1`` bosons per site
and the two readings part company, because the depletion of a general
``N``-level state is not a reduction of a dipole length; ``n`` remains the
controlled parameter, and is well defined even where the ordered state carries
no dipole at all, as for a quadrupolar state. Use
[`corrected_magnetic_moments`](@ref) for the ordered moment itself, which is
available in either mode.

Brillouin-zone integrals are controlled by `tol`, a relative accuracy target
for adaptive cubature, and/or `maxevals`, its budget of integrand evaluations.
Alternatively `grid` replaces the cubature by an average over a uniform grid of
the given dimensions.
"""
function boson_density(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)
    quad = BZQuadrature(; tol, maxevals, grid)
    L = nbands(swt)
    nf = nflavors(swt)
    # Bosons are laid out as (flavor, atom) with flavor fastest, so each site owns a
    # contiguous run of `nf` flavors.
    o = zero(Vec3)
    ns = [BosonMonomial(1.0+0im, (L+a, a), (o, o)) for a in 1:L]
    g = contractions(MagnonVacuum(swt), ns, quad)
    return [real(wick_expectation(ns[(i-1)*nf .+ (1:nf)], g, nothing)) for i in 1:div(L, nf)]
end

"""
    corrected_magnetic_moments(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)

Magnetic moments ``μ = -g 𝐒`` in units of the Bohr magneton, corrected at
relative order ``1/s``, for each site of the magnetic cell. Compare to the
classical [`magnetic_moments`](@ref), which this reduces to when the correction
is switched off.

Zero-point fluctuations act on the ordered structure in two independent ways at
this order, and both are applied here: the moment is depleted by the bosons of
[`boson_density`](@ref), and the structure is tilted by the zero-point pressure
of [`tadpole_correction`](@ref), whose canonical example is the change in canting
angle of an antiferromagnet in an applied field. Summing these moments over the
magnetic cell gives the uniform magnetization, so a sweep over
[`set_field!`](@ref) yields a corrected ``M`` vs. ``H`` curve. Because an
anisotropic ``g`` need not commute with the tilt, ``μ`` and ``𝐒`` are corrected
by different amounts; use `boson_density` for a statement about the boson count
itself.

Both effects are read off the same expansion of ``𝐒`` in the local frame, the
tilt from its one-boson word and the depletion from its two-boson word, and the
two are summed. The result is thermodynamically consistent in either mode, ``Σ_i
μ_i = -∂E/∂𝐁`` against the corrected energy of
[`corrected_energy_per_site`](@ref).

The tilt vanishes for a collinear structure in `:dipole` mode, and requires the
cubic vertex, which is available for a narrower class of models than the
depletion is; see [`tadpole_correction`](@ref). Where it is unavailable the
moments are returned depleted but untilted, which is still correct at this order
for any structure the tilt would not move.

Brillouin-zone integrals are controlled by `tol`, a relative accuracy target
for adaptive cubature, and/or `maxevals`, its budget of integrand evaluations.
Alternatively `grid` replaces the cubature by an average over a uniform grid of
the given dimensions.
"""
function corrected_magnetic_moments(swt::SpinWaveTheory; tol=nothing, maxevals=nothing, grid=nothing)
    quad = BZQuadrature(; tol, maxevals, grid)
    (; sys) = swt
    L = nbands(swt)
    Na = nsites(sys)

    # The tilt needs the cubic vertex; without it the depletion below is still
    # the whole correction for any structure the tilt would not move.
    v = isnothing(corrections_unsupported_reason(swt)) ? tadpole_correction(swt; tol, maxevals, grid).v : zeros(ComplexF64, L)
    w = [v; conj(v)]

    # ⟨𝐒⟩ to O(1/s): the classical dipole, its tilt by one displaced leg, and
    # its depletion by one contraction. Every word acts on a single site, so
    # there is no wavevector dependence to keep.
    g = contractions(MagnonVacuum(swt), (t for α in 1:3, i in 1:Na for t in spin_monomials(swt, α, i, Val{2}())), quad)
    # Shaped like `magnetic_moments`, i.e. indexed by `Site`. The leading dims
    # are (1, 1, 1) because `SpinWaveTheory` flattens any supercell into its
    # cell. Each component is real, the spin operators being Hermitian.
    μs = map(1:Na) do i
        S = Vec3(ntuple(3) do α
            sum(K -> real(wick_expectation(spin_monomials(swt, α, i, K), g, w)), (Val{0}(), Val{1}(), Val{2}()))
        end)
        return -sys.gs[1, 1, 1, i] * S
    end
    return reshape(μs, 1, 1, 1, Na)
end

"""
    static_self_energy(swt::SpinWaveTheory, qpts, terms2; vacuum=MagnonVacuum(swt))

Band-resolved energy shift caused by a correction `terms2` to the quadratic
Hamiltonian, evaluated in the quasi-particles of `vacuum` and including its
counterterm, as returned by [`hartree_fock_correction`](@ref) or
[`tadpole_correction`](@ref). Corrections from multiple sources should be
concatenated. The result has the same shape as [`dispersion`](@ref), to which it
is to be added.

Whereas [`corrected_dispersion`](@ref) rediagonalizes the corrected Hamiltonian,
this function evaluates the correction to first order only. The two differ at
relative order ``1/s²``, but only the latter can be combined with
[`cubic_self_energy`](@ref), which is likewise a first-order energy shift. It is
also the appropriate choice when the structure supports a Goldstone mode: an
``O(1/s)`` correction to a Hamiltonian with a protected zero mode produces a gap
of order ``\\sqrt{1/s}`` upon rediagonalization, obscuring the cancellation
between the terms that keeps the mode gapless. That cancellation is exact when
`terms2` is integrated with `grid` equal to the loop grid of
[`cubic_self_energy`](@ref); otherwise each quadrature leaves its own error,
and the mode is gapped by the larger of the two.
"""
function static_self_energy(swt::SpinWaveTheory, qpts, terms2; vacuum=MagnonVacuum(swt))
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)
    ret = stack(qpts.qs) do q
        q_reshaped = to_reshaped_rlu(swt.sys, q)
        vacuum_bogoliubov!(ws, H, vacuum, q_reshaped)
        δH = zeros(ComplexF64, 2L, 2L)
        accum_quadratic!(δH, terms2, q_reshaped)
        accum_counterterm!(δH, vacuum, q_reshaped)
        # Writing the quadratic form as (1/2) y† (T† δH T) y in the quasi-particle
        # basis, the coefficient of α†_n α_n is the n-th diagonal element.
        real.(diag(ws.T' * δH * ws.T))[1:L]
    end
    return reshape(ret, L, size(qpts.qs)...)
end

"""
    corrected_dispersion(swt::SpinWaveTheory, qpts, terms2)

Excitation energies including a mean-field correction to the quadratic
Hamiltonian, as returned by [`hartree_fock_correction`](@ref) or
[`tadpole_correction`](@ref). Corrections from multiple sources should be
concatenated. Otherwise like [`dispersion`](@ref). See
[`static_self_energy`](@ref) for the alternative of applying the correction to
first order only.
"""
function corrected_dispersion(swt::SpinWaveTheory, qpts, terms2)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    vac = MagnonVacuum(swt, terms2)
    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)
    disp = stack(qpts.qs) do q
        vacuum_bogoliubov!(ws, H, vac, to_reshaped_rlu(swt.sys, q))[1:L]
    end
    return reshape(disp, L, size(qpts.qs)...)
end
