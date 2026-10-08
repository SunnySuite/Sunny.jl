# Corrections to linear spin wave theory (LSWT) at sub-leading order in 1/s.
#
# Expanding the spin Hamiltonian in Holstein-Primakoff bosons b, b† yields
#
#     H = E_cl + H₁ + H₂ + H₃ + H₄ + …,
#
# where Hₙ collects the n-boson terms and scales as s^(2-n/2). At a classical
# energy minimum H₁ vanishes. LSWT retains H₂ only, which `bogoliubov!`
# diagonalizes into quasi-particle modes α, α† with energies ε. Relative to
# LSWT, both H₄ (at first order) and H₃ (at second order) contribute at O(1/s).
#
# Writing quasi-particle momenta as subscripts and letting N denote the number
# of magnetic cells, the cubic term splits into a decay and a source channel,
#
#     H₃ = (1/2!√N) Σ_{1+2=3} [Γ₁(3; 1, 2) α†₃ α₁ α₂ + h.c.]
#        + (1/3!√N) Σ_{1+2+3=0} [Γ₂(1, 2, 3) α†₁ α†₂ α†₃ + h.c.],
#
# whose self-energies SelfEnergy.jl writes down. At T = 0 only the decay channel
# resonates at ω > 0, giving magnons a finite lifetime.
#
# Sunny's Nambu conventions are reused throughout. In particular the columns
# `L+1:2L` of a Bogoliubov matrix `T` obtained at wavevector `q` are the
# eigenvectors at `-q` (see `excitations!`), so a band index ranging over the
# full Nambu space `1:2L` reaches both channels at once.
#
# Resummation. Collect the quasi-particles into y_𝐪 = [α_𝐪; α†_{-𝐪}], with
# metric Ĩ = diagm([ones(L); -ones(L)]), and write w = T†u for the observable
# amplitudes (u from `set_swt_observable_vectors!`, corrected by
# Observables.jl). LSWT gives the retarded response
#
#     χ^{μν}(z) = w_μ' G₀(z) w_ν,   G₀(z) = (zĨ - |ε|)⁻¹,
#
# whose anti-Hermitian part (χ' - χ)/2πi is the broadened `intensities` at z = ω
# + iη, including the mirror poles at ω < 0.
#
# At O(1/s) the static mean fields of StaticCorrections.jl, the mean field,
# tadpole and anisotropy corrections, add Σstat = T†δH T, for a perturbation
# (1/2)x†δH x of the quadratic Hamiltonian. The cubic vertex couples each magnon
# to pairs of magnons, and the observable creates pairs directly too (Sᶻ = s -
# b†b). Both are captured by an auxiliary quadratic model: magnons coupled to a
# bath of free two-magnon states. For each pair of internal lines (𝐩 a, 𝐪-𝐩
# b) at pair energy x, define
#
#     y = [√18 U[a, b, :]; β],
#
# the vertex to each of the 2L external Nambu legs and the amplitudes β for each
# observable to create the pair (`pair_amplitude`). Forward lines (x > 0) are
# bath particles, backward lines (x < 0) bath holes of the opposite metric sign.
# Integrating out the bath gives the Cauchy transform
#
#     K(z) = Σ_pairs ± y y† / (z - x),
#
# with blocks K_mm (the cubic self-energy), K_md, K_dm and K_dd (the direct
# two-magnon continuum). The exact response of the auxiliary model is then
#
#     χ = w'Gw + K_dm G w + w'G K_md + K_dm G K_md + K_dd,
#     G = (zĨ - |ε| - Σstat - K_mm)⁻¹.
#
# This is exact at O(1/s) and, being the resolvent of a quadratic model,
# inherits its structure: Nambu symmetry, Goldstone protection (the static and
# dynamic 1/ε divergences cancel in the full 2L inverse), η as a pure
# Lorentzian, and the commutator sum rule ∫dω S = w'Ĩw. It is positive at ω > 0
# whenever the auxiliary model is stable, i.e. |ε| + Σstat + K_mm(0) ⪰ 0. Where
# it is not, the 1/s correction to some mode is comparable to its harmonic
# energy, and a resummed pole moves onto the imaginary axis. This is a breakdown
# of the expansion for that mode rather than a physical instability, and that
# mode may carry little weight in the observable; `corrected_channels` reports
# it. The interference terms linear in K_md are absent from the standard 1/s
# treatment, which takes the magnon spectral function as the major component of
# S and adds the two-magnon continuum separately (PRB 79, 144416, Sec. VI). For
# a trace measure they cancel in a zone sum, but not pointwise.
#
# K depends on frequency only through the scalar x, so the masses y y†, summed
# over the wavevectors of `loop_wavevectors`, are accumulated into bins of x and
# the Cauchy transform is applied afterwards. Linear splitting between
# neighbouring bins keeps each channel's measure semidefinite and preserves its
# zeroth and first moments; the shape error is O((Δ/η)²).
#
# The test suite certifies every term by comparing to exact diagonalization of a
# cluster with anisotropic interactions and readouts.
#
# Entry points. Every function takes `swt`, a `tol` that is either a relative
# accuracy or a `BZGrid`, and a `vacuum`, harmonic by default. A vacuum that
# `self_consistent_vacuum` solves records its `tol`, which every function then
# shares.
#
#   Vacuum       MagnonVacuum, self_consistent_vacuum
#   Spectra      corrected_intensities, corrected_intensities_bands
#   Statics      corrected_energy_per_site, gaussian_energy_per_site,
#                boson_density, corrected_magnetic_moments
#   Pieces       hartree_fock_correction, tadpole_correction,
#                observable_corrections, anisotropy_correction,
#                static_self_energy, cubic_self_energy, corrected_dispersion
#   Internals    OneLoop, SelfEnergy, corrected_channels

# Why the 1/s corrections of this directory are unavailable for `swt`, or `nothing`
# if they are available. Returned rather than thrown so that a caller offering a
# correct but weaker result in the unsupported cases, such as
# `corrected_magnetic_moments`, can ask without catching.
function corrections_unsupported_reason(swt::SpinWaveTheory)
    (; sys) = swt
    @assert sys.mode in (:dipole, :dipole_uncorrected, :SUN)

    is_entangled(sys) && return "are not supported for entangled units"
    isnothing(sys.ewald) || return "do not yet support long-range dipole-dipole interactions"

    # In :SUN mode every coupling, biquadratic included, has been decomposed into
    # the tensor pairs that `sun_monomials` expands, so there is nothing to reject.
    if sys.mode != :SUN
        for int in sys.interactions_union
            for pc in int.pair
                pc.isculled && break
                iszero(pc.biquad) || return "do not yet support biquadratic exchange"
            end
        end
    end
    return nothing
end

# Errors unless `swt` describes a model for which the 1/s corrections in this
# directory are implemented.
function check_corrections_supported(swt::SpinWaveTheory)
    reason = corrections_unsupported_reason(swt)
    isnothing(reason) || error("1/s corrections $reason.")
end
