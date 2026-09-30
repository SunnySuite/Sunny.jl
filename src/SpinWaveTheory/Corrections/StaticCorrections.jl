# Static, frequency-independent 1/s corrections: the ground-state energy and the
# ordered moments. Both are integrals of the harmonic zero-point fluctuations over
# the Brillouin zone, which is what makes them cheap compared to the dynamical
# corrections of CorrectedIntensities.jl. Conventions are collected in
# Corrections.jl.

"""
    corrected_energy_per_site(swt::SpinWaveTheory; tol=nothing, maxevals=nothing)

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

The Brillouin-zone integral is performed by adaptive cubature, controlled by at
least one of `tol` (a relative accuracy target) or `maxevals` (a budget of
integrand evaluations).
"""
function corrected_energy_per_site(swt::SpinWaveTheory; tol=nothing, maxevals=nothing)
    isnothing(tol) && isnothing(maxevals) && error("Must specify `tol` or `maxevals` to control momentum-space integration.")

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

    # Integrate zero-point energy over the first Brillouin zone 𝐪 ∈ [0, 1]³ for
    # magnetic cell in reshaped RLU. Error bars are discarded.
    (δE₂, _) = hcubature((0, 0, 0), (1, 1, 1); rtol=@something(tol, 0),
                         maxevals=@something(maxevals, typemax(Int))) do q_reshaped
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
    boson_density(swt::SpinWaveTheory; tol=nothing, maxevals=nothing)

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

The Brillouin-zone integral is performed by adaptive cubature, controlled by at
least one of `tol` (a relative accuracy target) or `maxevals` (a budget of
integrand evaluations).
"""
function boson_density(swt::SpinWaveTheory; tol=nothing, maxevals=nothing)
    isnothing(tol) && isnothing(maxevals) && error("Must specify `tol` or `maxevals` to control momentum-space integration.")

    L = nbands(swt)
    nf = nflavors(swt)
    # ⟨b†_a b_a⟩ for every boson of the magnetic cell, in the Nambu labeling that
    # `nambu_correlations` expects
    gs = nambu_correlations(swt, [(L+a, a, (0, 0, 0)) for a in 1:L],
                            BosonMonomial{2}[]; tol, maxevals)
    # Bosons are laid out as (flavor, atom) with flavor fastest, so each site owns a
    # contiguous run of `nf` flavors.
    return [sum(α -> real(gs[(i-1)*nf + α]), 1:nf) for i in 1:div(L, nf)]
end

"""
    corrected_magnetic_moments(swt::SpinWaveTheory; tol=nothing, maxevals=nothing)

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

The Brillouin-zone integrals are performed by adaptive cubature, controlled by
at least one of `tol` (a relative accuracy target) or `maxevals` (a budget of
integrand evaluations).
"""
function corrected_magnetic_moments(swt::SpinWaveTheory; tol=nothing, maxevals=nothing)
    isnothing(tol) && isnothing(maxevals) && error("Must specify `tol` or `maxevals` to control momentum-space integration.")

    (; sys) = swt
    L = nbands(swt)
    Na = nsites(sys)

    # The tilt needs the cubic vertex; without it the depletion below is still
    # the whole correction for any structure the tilt would not move.
    v = if isnothing(corrections_unsupported_reason(swt))
        tadpole_correction(swt; tol, maxevals).v
    else
        zeros(ComplexF64, L)
    end
    # Nambu packing of the displacement, as `tadpole_correction` forms it
    w = [v; conj(v)]

    # The classical dipole, the tilt, and the depletion. Only onsite
    # correlations appear, every word acting on a single site, so there is no
    # wavevector dependence to keep.
    (words0, words1, words2) = map(K -> [spin_words(swt, α, i, K) for α in 1:3, i in 1:Na],
                                   (Val{0}(), Val{1}(), Val{2}()))
    ckeys = correlation_keys(L, reduce(vcat, words2))
    gs = nambu_correlations(swt, ckeys, BosonMonomial{2}[]; tol, maxevals)
    g = correlation_lookup(ckeys, gs, L)

    # Shaped like `magnetic_moments`, i.e. indexed by `Site`. The leading dims
    # are (1, 1, 1) because `SpinWaveTheory` flattens any supercell into its
    # cell.
    μs = map(1:Na) do i
        # Each component is real, the spin operators being Hermitian
        S = Vec3(ntuple(3) do α
            total(f, ws) = sum(f, ws[α, i]; init=zero(ComplexF64))
            real(total(t -> t.c, words0) +
                 total(t -> t.c * w[t.as[1]], words1) +
                 total(t -> t.c * g(t.as[1], t.as[2], t.ns[2] - t.ns[1]), words2))
        end)
        return -sys.gs[1, 1, 1, i] * S
    end
    return reshape(μs, 1, 1, 1, Na)
end
