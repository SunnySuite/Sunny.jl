# Static, frequency-independent 1/s corrections, each a contraction in the
# vacuum `swt.state`: the mean field of H₄, the tadpole, the corrections to the
# observables, and the energy and ordered moments that they and the zero-point
# fluctuations determine. Conventions are collected in Corrections.jl.

# Mean-field (Hartree-Fock) treatment of the four-boson term H₄, which shifts
# the magnon dispersion at O(1/s) relative to LSWT. Iterating it to
# self-consistency gives the vacuum that reproduces its own mean field, which
# resums one class of higher-order diagrams and drops others of the same order.
# That is uncontrolled, and in particular it violates the Ward identity of a
# broken continuous symmetry, gapping a mode that should be gapless.

"""
    hartree_fock_correction(swt::SpinWaveTheory; bz=default_bz(swt))

Decouples the four-boson term of the Holstein-Primakoff expansion into a
mean-field correction to the quadratic (LSWT) Hamiltonian, with mean fields
taken in the Gaussian state of `swt`, by default the LSWT ground state. The
correction is smaller than the LSWT Hamiltonian by a factor of order ``1/s``.
Returns `(; terms2, δE)`, where `terms2` is a correction to the quadratic
Hamiltonian, e.g. for `Sunny.with_correction`, and `δE` is a correction to the
energy per site.

The LSWT Hamiltonian need not be positive definite. If its frequencies are real,
as for a structure that quantum fluctuations select from a classically
degenerate family (order by disorder), the mean fields are taken in the Gaussian
state of its positive-norm modes, and the corrected Hamiltonian may then be
stable where LSWT is not.

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function hartree_fock_correction(swt::SpinWaveTheory; bz::BZIntegration=default_bz(swt))
    check_corrections_supported(swt)
    check_grid(swt, bz)
    terms4 = quartic_monomials(swt)
    (; terms2, δE) = mean_field(terms4, contractions(swt, terms4, bz))
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

"""
    GoldstonePin(q)

A constraint for [`self_consistent`](@ref) that keeps the Goldstone mode at
wavevector `q` (RLU) exactly gapless. At each iterate a uniform onsite term
``μ Σ_i b^†_i b_i`` is added, with ``μ`` chosen so that the quadratic Hamiltonian
has a zero mode at `q`, like the Lagrange multiplier of modified spin-wave
theory. The resulting state is not stationary in
[`gaussian_energy_per_site`](@ref), and since [`corrected_intensities`](@ref)
subtracts every correction of the state as a counterterm, it removes the pin
again. The pinned spectrum itself is that of [`intensities`](@ref).
"""
struct GoldstonePin
    q :: Vec3
end

GoldstonePin(q::AbstractVector) = GoldstonePin(Vec3(q))

# Uniform onsite term c Σ_a b†_a b_a
uniform_shift(c, L) = [BosonMonomial(complex(c), (L+a, a), (zero(Vec3), zero(Vec3))) for a in 1:L]

# Terms that the constraint `pin` adds to the correction `terms2` of the
# harmonic `swt`, and the multiplier that fixes them. A uniform shift moves
# every eigenvalue of the dynamical matrix equally, so one evaluation at `q`
# solves for it exactly. The zero mode is left at `swt.regularization`, like any
# Goldstone mode of LSWT.
function constraint_terms(pin::GoldstonePin, swt, terms2)
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    dynamical_matrix!(H, with_correction(swt, terms2), to_supercell_rlu(swt.sys, pin.q))
    μ = swt.regularization - eigmin(Hermitian(H))
    return (uniform_shift(μ, L), (; μ))
end

"""
    self_consistent(swt::SpinWaveTheory; grid, constraints=(), guess=nothing, scf_tol=1e-8, maxiters=200)

The [`SpinWaveTheory`](@ref) about the Gaussian state that reproduces its own
mean field, i.e. self-consistent Hartree-Fock-Bogoliubov. Its `dispersion` and
`intensities` are those of the self-consistent quadratic Hamiltonian, and the
``1/s`` corrections take it as the vacuum about which the bosons are expanded.
Among Gaussian states it is stationary in the energy of
[`gaussian_energy_per_site`](@ref). It exists also for a structure that is
classically unstable but stabilized by quantum fluctuations, in which case LSWT
has no vacuum. The iteration starts from the state of `guess`, a
`SpinWaveTheory`, if given, otherwise from LSWT, and throws an
`InstabilityError` if the fixed point needs a stabilizing shift.
Self-consistency resums one class of higher-order terms and drops others of the
same order. It is uncontrolled, and it gaps the Goldstone mode of a broken
continuous symmetry. Each of `constraints`, e.g. [`GoldstonePin`](@ref), adds
terms to the correction at every iterate, and the multipliers that fix them are
reported in the `info` of the state.

[`self_consistent_tadpole`](@ref) additionally relaxes the structure.

Momentum integrals are discrete sums over the points of `grid`, a
[`BZGrid`](@ref). The state records it, and every ``1/s`` correction about
this state uses the same grid by default, so that the counterterm cancels the
mean field exactly in [`corrected_intensities`](@ref). A grid that resolves a
regulator `η` there is, e.g., `auto_bzgrid(swt; η, tol=0.01)`.
"""
function self_consistent(swt::SpinWaveTheory; grid::BZGrid, constraints=(), guess=nothing, scf_tol=1e-8, maxiters=200,
                         memory=5)
    # The fixed point of g ↦ contractions in the vacuum of H₂ plus the mean
    # field of g. The iteration is accelerated by Anderson mixing over the last
    # `memory` iterates, with real coefficients, so each iterate is an affine
    # combination of correlations of Gaussian states and its mean field stays
    # Hermitian. Where an iterate has no vacuum, a uniform shift s b†b restores
    # one, starting from the most negative eigenvalue on the grid. The shift
    # must vanish at the fixed point.
    check_corrections_supported(swt)
    swt = harmonic(swt)
    L = nbands(swt)
    terms4 = quartic_monomials(swt)
    probe = grid_points(grid)
    shifted(terms2, s) = with_correction(swt, [terms2; uniform_shift(s, L)])

    # The mean field of `g` with the terms of each constraint, and their
    # multipliers
    function correction(g)
        terms2 = mean_field(terms4, g).terms2
        info = (;)
        for c in constraints
            (extra, i) = constraint_terms(c, swt, terms2)
            terms2 = [terms2; extra]
            info = Base.merge(info, i)
        end
        return (terms2, info)
    end

    H = zeros(ComplexF64, 2L, 2L)
    function stabilized(terms2)
        s = 0.0
        while true
            try
                return (contractions(shifted(terms2, s), terms4, grid), s)
            catch err
                err isa InstabilityError || rethrow()
                sw = shifted(terms2, s)
                λ = minimum(probe) do q
                    dynamical_matrix!(H, sw, q)
                    eigmin(Hermitian(H))
                end
                s += max(-1.1λ, 1e-3)
            end
        end
    end

    (g, _) = stabilized(isnothing(guess) ? BosonMonomial{2}[] : guess.state.correction)
    (xs, rs) = (Vector{ComplexF64}[], Vector{ComplexF64}[])
    realify(v) = [real(v); imag(v)]
    for iter in 1:maxiters
        (terms2, info) = correction(g)
        (g′, s) = stabilized(terms2)
        r = g′.values - g.values
        norm(r) < scf_tol && return iszero(s) ? SpinWaveTheory(swt, GaussianState(terms2, grid, (; method=:hartree_fock, iter, info...))) :
            throw(InstabilityError("Quantum fluctuations do not stabilize the structure: the self-consistent \
                                    mean field leaves the quadratic Hamiltonian indefinite."))
        push!(xs, g.values)
        push!(rs, r)
        length(xs) > memory + 1 && (popfirst!(xs); popfirst!(rs))
        x = g.values + r
        if length(xs) > 1
            ΔX = stack(diff(xs))
            ΔR = stack(diff(rs))
            γ = stack(realify, eachcol(ΔR)) \ realify(r)
            x -= (ΔX + ΔR) * γ
        end
        g = Contractions(g, x)
    end
    error("Self-consistent mean field not converged in $maxiters iterations")
end

"""
    self_consistent_tadpole(swt::SpinWaveTheory; grid, constraints=(), guess=nothing, scf_tol=1e-8, maxiters=200, memory=5)

Like [`self_consistent`](@ref), but additionally relaxes the structure: the
dipoles are rotated until the energy of the state is stationary in their
directions too, so that zero-point fluctuations shift e.g. the canting angle in
a field, as [`tadpole_correction`](@ref) does at first order. The result is the
stationary point of [`gaussian_energy_per_site`](@ref) jointly in the Gaussian
state and the structure, so its derivatives obey the Hellmann-Feynman theorem,
e.g. the magnetization is the expectation value of the moment in the state.
Only mode `:dipole` is supported.

The energy of a Gaussian state about a structure,

    F(frames, g) = E_cl + (⟨H₂⟩ + ⟨H₄⟩) / N,

with contractions g, is the functional that `self_consistent` makes stationary
in g at fixed frames. Making it stationary in the frames as well, at fixed g,
completes the variational principle: the joint stationary point is the one
whose derivatives in any coupling, e.g. the field, are partial derivatives of F
at fixed state and structure. Relaxing to the stationary point of E_cl + ⟨H₂⟩
alone, which is the condition that H₁ cancel one contraction of H₃, would be
the first-order tadpole of `tadpole_correction` iterated, and breaks that by
the frame dependence of ⟨H₄⟩, about 1% in the magnetization of a canted
antiferromagnet at s = 1/2. That dependence is the two-fold contraction of a
five-boson term, which the truncated expansion has no other place for.

The frames rather than the displacement of the bosons are what relax. A
displacement is the first-order description of a rotation; at finite tilt it
would leave the structure off a classical extremum, where the expansion
carries H₁ and the vacuum of H₂ need not exist. Each relaxation step instead
rotates every frame by a small angle about its own transverse axes, which
continues the frames smoothly, so that the correction, written in local
labels, remains a good guess for the next solve.

A rotation of every frame by a symmetry of the Hamiltonian, e.g. about the
field, leaves F invariant only together with the state. At fixed contractions
the anomalous ones carry a phase in each local frame, which pins the azimuth,
so the relaxation cannot drift along the symmetry, and a converged state
rotated by it is again converged, with the same spectrum.
"""
function self_consistent_tadpole(swt::SpinWaveTheory; grid::BZGrid, constraints=(), guess=nothing, scf_tol=1e-8,
                                 maxiters=200, memory=5)
    swt.sys.mode == :dipole || error("Self-consistent tadpole requires mode :dipole")
    frames = swt.data.local_rotations
    for iter in 1:maxiters
        sc = self_consistent(with_frames(swt, frames); grid, constraints, guess, scf_tol, maxiters, memory)
        frames′ = relax_frames(sc, grid)
        δ = maximum(norm(R[:, 3] - R′[:, 3]) for (R, R′) in zip(frames, frames′))
        δ < 10scf_tol && return SpinWaveTheory(sc, GaussianState(sc.state.correction, grid,
                                                       Base.merge(sc.state.info, (; tadpole_iter=iter))))
        (frames, guess) = (frames′, sc)
    end
    error("Self-consistent tadpole not converged in $maxiters iterations")
end

# F of the comment above at `frames`, with contractions `g` of the labels of
# `frame_labels`
function frame_energy(swt::SpinWaveTheory, frames, g)
    sw = with_frames(swt, frames)
    E = wick_expectation(quadratic_monomials(sw.sys, sw.data), g, nothing) +
        wick_expectation(quartic_monomials(sw), g, nothing; nfactors=2)
    return sw.classical_energy + anisotropy_correction(sw).δE + real(E) / nsites(uncontracted_system(sw.sys))
end

# Every Nambu combination of the slot pairs of H₂ and H₄, as unit monomials. A
# coefficient that vanishes in one frame need not vanish in another, e.g. the
# anomalous terms of a collinear structure, so the labels must not depend on
# which happen to be nonzero.
function frame_labels(swt::SpinWaveTheory)
    L = nbands(swt)
    ret = BosonMonomial{2}[]
    for (; as, ns) in Iterators.flatten((quadratic_monomials(swt.sys, swt.data), quartic_monomials(swt)))
        for p in eachindex(as), q in p+1:lastindex(as)
            (i, j) = (mod1(as[p], L), mod1(as[q], L))
            for a in (i, L+i), b in (j, L+j)
                push!(ret, BosonMonomial(1.0+0im, (a, b), (ns[p], ns[q])))
            end
        end
    end
    return ret
end

# The theory of `swt` about the structure whose dipoles point along column 3 of
# the local `frames`, in the harmonic vacuum. The frames are taken as given,
# rather than rebuilt from the dipoles, so that the labels of the bosons, and
# with them the phases of the anomalous contractions, continue smoothly as
# `relax_frames` tilts each frame. The structure need not be a classical energy
# minimum.
function with_frames(swt::SpinWaveTheory, frames::Vector{Mat3})
    sys = clone_system(swt.sys)
    for (i, site) in enumerate(eachsite(sys))
        iszero(frames[i]) || set_dipole!(sys, frames[i][:, 3], site)
    end
    data = swt_data(sys, swt.measure; frames)
    return SpinWaveTheory(sys, data, swt.measure, swt.regularization, energy_per_site(sys))
end

# The local frames of `swt` rotated by angles θ[2i-1], θ[2i] about the two
# transverse axes of frame i
function rotate_frames(frames, θ)
    return map(eachindex(frames)) do i
        ω = frames[i] * Vec3(θ[2i-1], θ[2i], 0)
        iszero(ω) ? frames[i] : axis_angle_to_matrix(ω / norm(ω), norm(ω)) * frames[i]
    end
end

# Frames at which F is stationary for the contractions of the state of `sc`,
# found by Newton's method from the frames of `sc`, with derivatives by finite
# differences of step `h`. The F of a frozen state is smooth in the frames, and
# its Hessian is computed once. The pseudo-inverse only guards against a
# direction in which F is flat.
function relax_frames(sc::SpinWaveTheory, grid; h=1e-4, maxiters=20, tol=1e-12)
    g = contractions(sc, frame_labels(sc), grid)
    frames = sc.data.local_rotations
    n = 2length(frames)
    e(k) = h * (1:n .== k)
    F(θ) = frame_energy(sc, rotate_frames(frames, θ), g)
    ∇F(θ) = [(F(θ + e(k)) - F(θ - e(k))) / 2h for k in 1:n]
    ∇²F = Symmetric([(F(e(k) + e(l)) - F(e(k) - e(l)) - F(e(l) - e(k)) + F(-e(k) - e(l))) / 4h^2 for k in 1:n, l in 1:n])
    H⁺ = pinv(∇²F; rtol=1e-6)
    θ = zeros(n)
    for _ in 1:maxiters
        r = ∇F(θ)
        norm(r) < tol && break
        θ -= H⁺ * r
    end
    return rotate_frames(frames, θ)
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
    tadpole_correction(swt::SpinWaveTheory; bz=default_bz(swt))

Computes the shift of the ordered magnetic structure caused by zero-point
fluctuations, which appears at relative order ``1/s``. An example is the change
in canting angle of an antiferromagnet in an applied field. Returns `(; terms2,
δE, v)`, where `terms2` is a correction to the quadratic Hamiltonian that can be
passed to `Sunny.with_correction`, `δE` is a correction to the energy per
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

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function tadpole_correction(swt::SpinWaveTheory; bz::BZIntegration=default_bz(swt))
    check_corrections_supported(swt)
    check_grid(swt, bz)
    terms3 = cubic_monomials(swt)
    (; terms2, δE, w) = tadpole(swt, terms3, contractions(swt, terms3, bz))
    return (; terms2, δE = δE / nsites(uncontracted_system(swt.sys)), v = w[1:nbands(swt)])
end

# The tadpole about the vacuum with contractions `g`: the Nambu-packed
# displacement `w`, the correction `terms2` it makes to the quadratic
# Hamiltonian, and the energy gain `δE` of the magnetic cell.
function tadpole(swt::SpinWaveTheory, terms3, g::Contractions)
    L = nbands(swt)
    (; noise) = g

    # The linear term at the classical minimum: one contraction of H₃, the
    # tadpole itself, plus the sub-leading one-boson words. These are
    # proportional to the classical energy gradient, so like the contractions
    # they vanish at a classical minimum but not once tadpole relaxation has
    # moved the structure off one. In mode :dipole those of an onsite
    # anisotropy vanish identically, because `rcs_factors` leaves the classical
    # energy exact. See `monomials_dipole`.
    ℓ = nambu_vector([wick_reduce(terms3, g, nothing, Val{1}()); lswt_correction_monomials(swt, Val{1}())], L, noise)

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
    dynamical_matrix!(H, harmonic(swt), zero(Vec3))
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
    observable_corrections(swt::SpinWaveTheory; v=nothing, bz=default_bz(swt))

Correction of relative order ``1/s`` to the amplitude for a magnon to be created
by each observable. Two effects contribute at this order: the cubic term of the
Holstein-Primakoff expansion of the transverse spin components, and the tilt of
the ordered structure by zero-point fluctuations. The latter requires the boson
displacement `v` of [`tadpole_correction`](@ref), and is omitted if `v` is
`nothing`. Returns the coefficients `δc[a, μ]` of the one-boson operators,
labeled as in [`accum_observable_corrections!`](@ref), which is what applies
them.

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function observable_corrections(swt::SpinWaveTheory; v=nothing, bz::BZIntegration=default_bz(swt))
    check_corrections_supported(swt)
    check_grid(swt, bz)
    w = isnothing(v) ? nothing : [v; conj(v)]
    g = contractions(swt, observable_cubic_monomials(swt), bz)
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
    corrected_energy_per_site(swt::SpinWaveTheory; bz=BZAdaptive(tol=0.01))

Energy per site of the magnetic structure, corrected at relative order ``1/s``,
where ``s`` is the spin magnitude in dipole mode, or ``1/λ`` for the
representation label ``λ`` in SU(N) mode. If the classical energy is ``J s²``,
the correction appears at order ``J s``.

The correction is the zero-point energy of the harmonic magnons, ``(1/2) Σ_n
∫d³q ω(𝐪, n)`` over the first magnetic Brillouin zone, less the constant that
the Holstein-Primakoff normal ordering leaves behind, together with the constant
that [`anisotropy_correction`](@ref) generates from an onsite coupling. The last
of these vanishes identically in `:dipole` mode, where `rcs_factors` makes the
classical energy exact.

Not included are the corrections of [`tadpole_correction`](@ref) and
[`hartree_fock_correction`](@ref), which are smaller by a further power of
``1/s`` and so belong to the next order. To include them, add their `δE` fields
to this result. See [`gaussian_energy_per_site`](@ref) for the energy of a
self-consistent state.

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function corrected_energy_per_site(swt::SpinWaveTheory; bz::BZIntegration=BZAdaptive(tol=0.01))
    # Vanishes in :SUN mode, where an onsite coupling enters the boson
    # Hamiltonian exactly and there is nothing to add.
    return swt.classical_energy + zero_point_energy(harmonic(swt), bz) + anisotropy_correction(swt).δE
end

"""
    gaussian_energy_per_site(swt::SpinWaveTheory; bz=default_bz(swt))

Expectation value of the Hamiltonian per site, through its four-boson term, in
the Gaussian state of `swt`. This is a variational upper bound on the ground
state energy, which is stationary in the state of [`self_consistent`](@ref), and
so suited to comparing candidate structures under
self-consistent Hartree-Fock. In the harmonic vacuum it is the result of
[`corrected_energy_per_site`](@ref) plus ``⟨H₄⟩``, i.e. less the `δE` of
[`hartree_fock_correction`](@ref).

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function gaussian_energy_per_site(swt::SpinWaveTheory; bz::BZIntegration=default_bz(swt))
    check_corrections_supported(swt)
    check_grid(swt, bz)
    terms4 = quartic_monomials(swt)
    (; correction) = swt.state
    g = contractions(swt, Iterators.flatten((correction, terms4)), bz)
    # ⟨H₂⟩ is the zero-point energy of H₂ plus the correction, less
    # ⟨correction⟩, and ⟨H₄⟩ is minus the constant of its mean field
    E₄ = -mean_field(terms4, g).δE - real(wick_expectation(correction, g, nothing))
    return swt.classical_energy + zero_point_energy(swt, bz) + anisotropy_correction(swt).δE +
           E₄ / nsites(uncontracted_system(swt.sys))
end

# Zero-point energy per site of the vacuum, averaged over the magnetic Brillouin
# zone, less the normal-ordering constant: the zone average of the trace of the
# (1,1)-block of its quadratic Hamiltonian. This trace is q-dependent wherever a
# bond joins a site to its own periodic image, e.g. a ferromagnet, so it cannot
# be taken at q = 0.
function zero_point_energy(swt::SpinWaveTheory, bz)
    # Normalize per physical site
    Nsites = nsites(uncontracted_system(swt.sys))
    L = nbands(swt)
    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)
    return bz_average(bz) do q_supercell
        dynamical_matrix!(H, swt, q_supercell)
        trA = real(tr(view(H, 1:L, 1:L)))
        ωs = bogoliubov!(ws, H)
        return (sum(view(ωs, 1:L)) - trA) / 2Nsites
    end
end

"""
    boson_density(swt::SpinWaveTheory; bz=default_bz(swt))

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

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function boson_density(swt::SpinWaveTheory; bz::BZIntegration=default_bz(swt))
    is_entangled(swt.sys) && error("Corrections are not supported for entangled units")
    check_grid(swt, bz)
    L = nbands(swt)
    nf = nflavors(swt)
    # Bosons are laid out as (flavor, atom) with flavor fastest, so each site owns a
    # contiguous run of `nf` flavors.
    o = zero(Vec3)
    ns = [BosonMonomial(1.0+0im, (L+a, a), (o, o)) for a in 1:L]
    g = contractions(swt, ns, bz)
    return [real(wick_expectation(ns[(i-1)*nf .+ (1:nf)], g, nothing)) for i in 1:div(L, nf)]
end

"""
    corrected_magnetic_moments(swt::SpinWaveTheory; bz=default_bz(swt))

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

Brillouin zone integrals are evaluated according to `bz`, either
[`BZAdaptive`](@ref) parameters or an explicit [`BZGrid`](@ref).
"""
function corrected_magnetic_moments(swt::SpinWaveTheory; bz::BZIntegration=default_bz(swt))
    is_entangled(swt.sys) && error("Corrections are not supported for entangled units")
    check_grid(swt, bz)
    (; sys) = swt
    L = nbands(swt)
    Na = nsites(sys)

    # The tilt needs the cubic vertex; without it the depletion below is still
    # the whole correction for any structure the tilt would not move.
    v = isnothing(corrections_unsupported_reason(swt)) ? tadpole_correction(swt; bz).v : zeros(ComplexF64, L)
    w = [v; conj(v)]

    # ⟨𝐒⟩ to O(1/s): the classical dipole, its tilt by one displaced leg, and
    # its depletion by one contraction. Every word acts on a single site, so
    # there is no wavevector dependence to keep.
    g = contractions(swt, (t for α in 1:3, i in 1:Na for t in spin_monomials(swt, α, i, Val{2}())), bz)
    # Shaped like `magnetic_moments`, i.e. indexed by `Site`. Each component is
    # real, the spin operators being Hermitian.
    μs = map(1:Na) do i
        S = Vec3(ntuple(3) do α
            sum(K -> real(wick_expectation(spin_monomials(swt, α, i, K), g, w)), (Val{0}(), Val{1}(), Val{2}()))
        end)
        return -sys.gs[i] * S
    end
    return reshape(μs, size(sys.dipoles))
end

"""
    static_self_energy(swt::SpinWaveTheory, qpts, terms2)

Band-resolved energy shift caused by a correction `terms2` to the quadratic
Hamiltonian, evaluated in the quasi-particles of the Gaussian state of `swt` and
including its counterterm, as returned by [`hartree_fock_correction`](@ref) or
[`tadpole_correction`](@ref). Corrections from multiple sources should be
concatenated. The result has the same shape as [`dispersion`](@ref), to which it
is to be added.

Whereas `dispersion(Sunny.with_correction(swt, terms2), qpts)` rediagonalizes
the corrected Hamiltonian, this function evaluates the correction to first
order only. The two differ at relative order ``1/s²``, but only the latter can
be combined with
[`cubic_self_energy`](@ref), which is likewise a first-order energy shift. It is
also the appropriate choice when the structure supports a Goldstone mode: an
``O(1/s)`` correction to a Hamiltonian with a protected zero mode produces a gap
of order ``\\sqrt{1/s}`` upon rediagonalization, obscuring the cancellation
between the terms that keeps the mode gapless. That cancellation is exact when
`terms2` and [`cubic_self_energy`](@ref) are integrated with the same
[`BZGrid`](@ref); otherwise each quadrature leaves its own error, and the mode
is gapped by the larger of the two.
"""
function static_self_energy(swt::SpinWaveTheory, qpts, terms2)
    L = nbands(swt)
    qpts = convert(AbstractQPoints, qpts)
    H = zeros(ComplexF64, 2L, 2L)
    ws = BogoliubovWorkspace(L)
    ret = stack(qpts.qs) do q
        q_supercell = to_supercell_rlu(swt.sys, q)
        excitations_supercell!(ws, H, swt, q_supercell)
        # Writing the quadratic form as (1/2) y† Σ y in the quasi-particle
        # basis, the coefficient of α†_n α_n is the n-th diagonal element.
        real.(diag(static_self_energy_matrix(swt, terms2, q_supercell, ws.T)))[1:L]
    end
    return reshape(ret, L, size(qpts.qs)...)
end

# Nambu matrix of the static correction `terms2`, counterterm of the vacuum
# included, in the quasi-particle basis `T` of the vacuum at `q_supercell`
function static_self_energy_matrix(swt::SpinWaveTheory, terms2, q_supercell, T)
    δH = zeros(ComplexF64, size(T))
    accum_quadratic!(δH, terms2, q_supercell)
    accum_counterterm!(δH, swt, q_supercell)
    return T' * δH * T
end
